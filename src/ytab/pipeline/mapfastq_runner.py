"""Restartable, one-sample-at-a-time local MapFastq orchestration."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml


STATUS_FIELDS = [
    "sample", "status", "layout", "fastq_1", "fastq_2", "output_dir",
    "log_file", "elapsed_seconds", "message",
]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _resolve(value: str | Path | None, root: Path) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _index_complete(prefix: Path) -> bool:
    return any(
        all(Path(f"{prefix}{suffix}.{extension}").is_file() for suffix in (
            ".1", ".2", ".3", ".4", ".rev.1", ".rev.2"
        ))
        for extension in ("bt2", "bt2l")
    )


def load_project_for_mapping(project_config: Path) -> dict:
    config_path = Path(project_config).expanduser().resolve()
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"Project config is not a YAML mapping: {config_path}")

    repo_root = _resolve(config.get("repo_root"), config_path.parent) or Path(__file__).resolve().parents[3]
    config["_project_config"] = str(config_path)
    config["_repo_root"] = str(repo_root)

    sample_sheet = _resolve(config.get("sample_sheet"), repo_root)
    if sample_sheet:
        if not sample_sheet.is_file():
            raise FileNotFoundError(f"Sample sheet is missing: {sample_sheet}")
        with sample_sheet.open(newline="", encoding="utf-8") as handle:
            config["samples"] = list(csv.DictReader(handle))

    reference_path = config_path.with_name("reference_resolved.json")
    if reference_path.is_file():
        with reference_path.open(encoding="utf-8") as handle:
            reference = json.load(handle)
        if not isinstance(reference, dict):
            raise ValueError(f"Reference file is not a JSON object: {reference_path}")
        config["reference"] = reference
    config["_reference_resolved"] = str(reference_path)
    return config


def get_included_samples(config: dict) -> list[dict]:
    return [dict(sample) for sample in config.get("samples", []) if _as_bool(sample.get("include", True))]


def validate_mapping_inputs(config: dict) -> list[str]:
    errors: list[str] = []
    root = Path(config["_repo_root"])
    try:
        threads = int(config.get("threads", 0))
    except (TypeError, ValueError):
        threads = 0
    if threads < 2:
        errors.append("Threads must be at least 2.")

    reference = config.get("reference") or {}
    prefix = _resolve(reference.get("bowtie2_index_prefix"), root)
    if prefix is None or not _index_complete(prefix):
        errors.append("No complete Bowtie2 index found. Run reference preparation before alignment.")

    samples = get_included_samples(config)
    if not samples:
        errors.append("No included samples were found in the sample sheet.")
    for sample in samples:
        name = str(sample.get("sample") or "")
        if not name:
            errors.append("An included sample has no sample name.")
        elif Path(name).name != name:
            errors.append(f"Sample name contains path separators: {name}")
        layout = str(sample.get("layout") or "single").lower()
        fastq_1 = _resolve(sample.get("fastq_1"), root)
        fastq_2 = _resolve(sample.get("fastq_2"), root)
        if fastq_1 is None or not fastq_1.is_file():
            errors.append(f"FASTQ file missing for sample {name}: {fastq_1 or 'fastq_1 not set'}")
        if layout == "paired" and (fastq_2 is None or not fastq_2.is_file()):
            errors.append(f"FASTQ R2 missing for paired sample {name}: {fastq_2 or 'fastq_2 not set'}")
        if layout not in {"single", "paired"}:
            errors.append(f"Unsupported layout for sample {name}: {layout}")
    return errors


def _paths(config: dict, sample_name: str) -> dict[str, Path]:
    root = Path(config["_repo_root"])
    project_dir = _resolve(config.get("output_project_dir"), root)
    export_dir = _resolve(config.get("output_export_dir"), root)
    if project_dir is None or export_dir is None:
        raise ValueError("Project and export output directories must be configured.")
    return {
        "output": project_dir / "mapfastq" / sample_name,
        "log": project_dir / "logs" / "mapfastq" / f"{sample_name}.log",
        "manifest": project_dir / "manifests" / "mapfastq" / f"{sample_name}.mapfastq_manifest.json",
        "status": project_dir / "manifests" / "mapfastq" / "mapfastq_status.csv",
        "qc": export_dir / "qc" / "mapping",
    }


def build_mapfastq_command(sample: dict, config: dict, threads: int | None = None) -> list[str]:
    root = Path(config["_repo_root"])
    selected_threads = int(config.get("threads") if threads is None else threads)
    if selected_threads < 2:
        raise ValueError("Threads must be at least 2.")
    name = str(sample["sample"])
    layout = str(sample.get("layout") or "single").lower()
    prefix = _resolve((config.get("reference") or {}).get("bowtie2_index_prefix"), root)
    output = _paths(config, name)["output"]
    script = root / "src" / "ytab" / "mapping" / "MapFastq.py"
    command = [
        sys.executable, str(script), "--out-dir", str(output), "--bt2-index",
        str(prefix), "--threads", str(selected_threads), "--sample-name", name,
    ]
    fastq_1 = _resolve(sample.get("fastq_1"), root)
    fastq_2 = _resolve(sample.get("fastq_2"), root)
    if layout == "paired":
        command.extend(["--r1", str(fastq_1), "--r2", str(fastq_2)])
    else:
        command.extend(["--input-file-name", str(fastq_1)])
    return command


def _detected_outputs(output_dir: Path) -> list[str]:
    return sorted(str(path) for path in output_dir.iterdir() if path.is_file()) if output_dir.is_dir() else []


def _expected_outputs(output_dir: Path, sample: str) -> bool:
    bam = output_dir / f"{sample}.sorted.bam"
    return bam.is_file() and (Path(f"{bam}.bai").is_file() or bam.with_suffix(".bai").is_file())


def _write_manifest(path: Path, manifest: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def run_mapfastq_sample(
    sample: dict, config: dict, threads: int | None = None,
    force: bool = False, dry_run: bool = False,
) -> dict:
    root = Path(config["_repo_root"])
    name = str(sample["sample"])
    paths = _paths(config, name)
    selected_threads = int(config.get("threads") if threads is None else threads)
    command = build_mapfastq_command(sample, config, selected_threads)
    fastq_1 = _resolve(sample.get("fastq_1"), root)
    fastq_2 = _resolve(sample.get("fastq_2"), root)
    prefix = _resolve((config.get("reference") or {}).get("bowtie2_index_prefix"), root)
    start_time = _now()
    started = time.monotonic()
    warnings = [str(sample.get("warnings"))] if sample.get("warnings") else []
    manifest = {
        "project_id": config.get("project_id"), "sample": name,
        "fastq_1": str(fastq_1) if fastq_1 else None,
        "fastq_2": str(fastq_2) if fastq_2 else None,
        "layout": str(sample.get("layout") or "single").lower(),
        "species": config.get("species"), "bowtie2_index_prefix": str(prefix),
        "threads": selected_threads, "output_dir": str(paths["output"]),
        "log_file": str(paths["log"]), "status": "failed",
        "start_time": start_time, "end_time": None, "elapsed_seconds": 0.0,
        "command_run": command, "detected_outputs": [], "warnings": warnings,
        "error_message": None,
    }

    if not force and paths["manifest"].is_file() and _expected_outputs(paths["output"], name):
        try:
            previous = json.loads(paths["manifest"].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = {}
        if previous.get("status") in {"success", "skipped"}:
            manifest.update(status="skipped", end_time=_now())
            manifest["warnings"].append("Successful cached outputs found; mapping skipped.")
            manifest["detected_outputs"] = _detected_outputs(paths["output"])
            _write_manifest(paths["manifest"], manifest)
            return manifest

    if dry_run:
        manifest.update(status="skipped", end_time=_now())
        manifest["warnings"].append("Dry run; command was not executed.")
        _write_manifest(paths["manifest"], manifest)
        return manifest

    paths["output"].mkdir(parents=True, exist_ok=True)
    paths["log"].parent.mkdir(parents=True, exist_ok=True)
    try:
        with paths["log"].open("w", encoding="utf-8") as log:
            log.write("Command: " + " ".join(command) + "\n\n")
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
        if completed.returncode != 0:
            raise RuntimeError(f"MapFastq exited with status {completed.returncode}; see {paths['log']}")
        if not _expected_outputs(paths["output"], name):
            raise RuntimeError("MapFastq completed without the expected sorted BAM and index.")
        stats = paths["output"] / f"{name}.mapping_stats.csv"
        if stats.is_file():
            paths["qc"].mkdir(parents=True, exist_ok=True)
            shutil.copy2(stats, paths["qc"] / stats.name)
        manifest["status"] = "success"
    except Exception as exc:
        manifest["error_message"] = str(exc)
    manifest["end_time"] = _now()
    manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
    manifest["detected_outputs"] = _detected_outputs(paths["output"])
    _write_manifest(paths["manifest"], manifest)
    return manifest


def _write_status(path: Path, results: list[dict]) -> None:
    existing: dict[str, dict] = {}
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            existing = {row["sample"]: row for row in csv.DictReader(handle)}
    for result in results:
        existing[result["sample"]] = {
            "sample": result["sample"], "status": result["status"],
            "layout": result["layout"], "fastq_1": result["fastq_1"] or "",
            "fastq_2": result["fastq_2"] or "", "output_dir": result["output_dir"],
            "log_file": result["log_file"], "elapsed_seconds": result["elapsed_seconds"],
            "message": result.get("error_message") or "; ".join(result.get("warnings") or []),
        }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerows(existing[name] for name in sorted(existing))


def run_mapfastq_project(
    project_config: Path, samples: list[str] | None = None,
    threads: int | None = None, force: bool = False, dry_run: bool = False,
    keep_going: bool = False,
) -> dict:
    config = load_project_for_mapping(project_config)
    if threads is not None:
        config["threads"] = threads
    errors = validate_mapping_inputs(config)
    if errors:
        raise ValueError("\n".join(errors))
    included = get_included_samples(config)
    by_name = {str(sample["sample"]): sample for sample in included}
    selected_names = list(by_name) if samples is None else samples
    unknown = [name for name in selected_names if name not in by_name]
    if unknown:
        raise ValueError(f"Selected samples are not included or do not exist: {', '.join(unknown)}")
    results: list[dict] = []
    for name in selected_names:
        result = run_mapfastq_sample(by_name[name], config, threads, force, dry_run)
        results.append(result)
        _write_status(_paths(config, name)["status"], results)
        if result["status"] == "failed" and not keep_going:
            break
    return {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "bowtie2_index_prefix": str(_resolve((config.get("reference") or {}).get("bowtie2_index_prefix"), Path(config["_repo_root"]))),
        "threads": int(config["threads"]), "included_samples": len(included),
        "selected_samples": selected_names, "results": results,
        "success": sum(row["status"] == "success" for row in results),
        "failed": sum(row["status"] == "failed" for row in results),
        "skipped": sum(row["status"] == "skipped" for row in results),
        "dry_run": dry_run,
    }
