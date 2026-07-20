"""Restartable local orchestration for per-sample CreateHitFile runs."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from .mapfastq_runner import get_included_samples as _get_included_samples
from .mapfastq_runner import load_project_for_mapping


STATUS_FIELDS = [
    "sample", "status", "input_bam", "hits_file", "output_dir",
    "log_file", "elapsed_seconds", "message",
]
TRACK_SUFFIXES = (".bed", ".bedgraph", ".wig", ".bw", ".bigwig")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _resolve(value: str | Path | None, root: Path) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def load_project_for_create_hit_file(project_config: Path) -> dict:
    return load_project_for_mapping(project_config)


def get_included_samples(config: dict) -> list[dict]:
    return _get_included_samples(config)


def _reference_files(config: dict) -> dict[str, Path | str]:
    root = Path(config["_repo_root"])
    reference = config.get("reference") or {}
    fasta = _resolve(reference.get("fasta"), root)
    gff = _resolve(reference.get("gff"), root)
    gtf = _resolve(reference.get("gtf"), root)
    feature_table = _resolve(reference.get("feature_table"), root)
    if gff and gff.is_file():
        features, feature_format = gff, "gff"
    elif gtf and gtf.is_file():
        features, feature_format = gtf, "gtf"
    elif feature_table and feature_table.is_file():
        features = feature_table
        feature_format = "cgd_tab" if "chromosomal_feature" in feature_table.name else "ncbi_feature_table"
    else:
        features, feature_format = None, ""
    return {"fasta": fasta, "features": features, "feature_format": feature_format}


def _paths(config: dict, sample_name: str) -> dict[str, Path]:
    root = Path(config["_repo_root"])
    project_dir = _resolve(config.get("output_project_dir"), root)
    export_dir = _resolve(config.get("output_export_dir"), root)
    if project_dir is None or export_dir is None:
        raise ValueError("Project and export output directories must be configured.")
    return {
        "project": project_dir,
        "output": project_dir / "create_hit_file" / sample_name,
        "hits": project_dir / "create_hit_file" / sample_name / f"{sample_name}_hits.txt",
        "log": project_dir / "logs" / "create_hit_file" / f"{sample_name}.log",
        "manifest": project_dir / "manifests" / "create_hit_file" / f"{sample_name}.create_hit_file_manifest.json",
        "status": project_dir / "manifests" / "create_hit_file" / "create_hit_file_status.csv",
        "tracks": export_dir / "browser_tracks" / "create_hit_file" / sample_name,
    }


def find_sample_bam(sample: dict, config: dict) -> Path:
    name = str(sample["sample"])
    return _paths(config, name)["project"] / "mapfastq" / name / f"{name}.sorted.bam"


def _bam_index_exists(bam: Path) -> bool:
    return Path(f"{bam}.bai").is_file() or bam.with_suffix(".bai").is_file()


def validate_create_hit_file_inputs(config: dict) -> list[str]:
    errors: list[str] = []
    try:
        threads = int(config.get("threads", 0))
    except (TypeError, ValueError):
        threads = 0
    if threads < 1:
        errors.append("Threads must be at least 1.")
    reference = _reference_files(config)
    fasta = reference["fasta"]
    features = reference["features"]
    if not isinstance(fasta, Path) or not fasta.is_file():
        errors.append("Reference FASTA is missing. Fix or prepare the reference resources first.")
    if not isinstance(features, Path) or not features.is_file():
        errors.append("No feature table, GFF, or GTF is available. Fix the reference resources first.")
    samples = get_included_samples(config)
    if not samples:
        errors.append("No included samples were found in the sample sheet.")
    for sample in samples:
        name = str(sample.get("sample") or "")
        if not name:
            errors.append("An included sample has no sample name.")
            continue
        if Path(name).name != name:
            errors.append(f"Sample name contains path separators: {name}")
        bam = find_sample_bam(sample, config)
        if not bam.is_file():
            errors.append(f"Mapped BAM missing for sample {name}: {bam}")
        elif not _bam_index_exists(bam):
            errors.append(f"BAM index missing for sample {name}: {bam}.bai")
    return errors


def build_create_hit_file_command(
    sample: dict, config: dict, threads: int | None = None
) -> list[str]:
    selected_threads = int(config.get("threads") if threads is None else threads)
    if selected_threads < 1:
        raise ValueError("Threads must be at least 1.")
    name = str(sample["sample"])
    root = Path(config["_repo_root"])
    reference = _reference_files(config)
    command = [
        sys.executable,
        str(root / "src" / "ytab" / "insertions" / "CreateHitFile.py"),
        "--bam", str(find_sample_bam(sample, config)),
        "--out-dir", str(_paths(config, name)["output"]),
        "--out-prefix", name,
        "--flat-output",
        "--fasta", str(reference["fasta"]),
        "--features", str(reference["features"]),
        "--feature-format", str(reference["feature_format"]),
        "--threads", str(selected_threads),
        "--write-browser-tracks",
    ]
    return command


def _detected_outputs(output_dir: Path) -> list[str]:
    return sorted(str(path) for path in output_dir.iterdir() if path.is_file()) if output_dir.is_dir() else []


def _write_manifest(path: Path, manifest: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def run_create_hit_file_sample(
    sample: dict, config: dict, threads: int | None = None,
    force: bool = False, dry_run: bool = False,
) -> dict:
    name = str(sample["sample"])
    paths = _paths(config, name)
    reference = _reference_files(config)
    selected_threads = int(config.get("threads") if threads is None else threads)
    command = build_create_hit_file_command(sample, config, selected_threads)
    bam = find_sample_bam(sample, config)
    started = time.monotonic()
    manifest = {
        "project_id": config.get("project_id"), "sample": name,
        "input_bam": str(bam), "species": config.get("species"),
        "reference_files_used": {
            "fasta": str(reference["fasta"]), "features": str(reference["features"]),
            "feature_format": reference["feature_format"],
        },
        "output_dir": str(paths["output"]), "hits_file": str(paths["hits"]),
        "log_file": str(paths["log"]), "status": "failed",
        "start_time": _now(), "end_time": None, "elapsed_seconds": 0.0,
        "command_run": command, "detected_outputs": [], "warnings": [],
        "error_message": None,
    }
    if not force and paths["manifest"].is_file() and paths["hits"].is_file():
        try:
            previous = json.loads(paths["manifest"].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = {}
        if previous.get("status") in {"success", "skipped"}:
            manifest.update(status="skipped", end_time=_now())
            manifest["warnings"].append("Successful cached hit file found; CreateHitFile skipped.")
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
            raise RuntimeError(f"CreateHitFile exited with status {completed.returncode}; see {paths['log']}")
        if not paths["hits"].is_file():
            raise RuntimeError("CreateHitFile completed without the expected hits file.")
        tracks = [path for path in paths["output"].iterdir() if path.is_file() and path.suffix.lower() in TRACK_SUFFIXES]
        if tracks:
            paths["tracks"].mkdir(parents=True, exist_ok=True)
            for track in tracks:
                shutil.copy2(track, paths["tracks"] / track.name)
            if not any(track.suffix.lower() in {".bw", ".bigwig"} for track in tracks):
                manifest["warnings"].append(
                    "BigWig was not produced; bedGraph and WIG tracks remain available."
                )
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
            "input_bam": result["input_bam"], "hits_file": result["hits_file"],
            "output_dir": result["output_dir"], "log_file": result["log_file"],
            "elapsed_seconds": result["elapsed_seconds"],
            "message": result.get("error_message") or "; ".join(result.get("warnings") or []),
        }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerows(existing[name] for name in sorted(existing))


def run_create_hit_file_project(
    project_config: Path, samples: list[str] | None = None,
    threads: int | None = None, force: bool = False, dry_run: bool = False,
    keep_going: bool = False,
) -> dict:
    config = load_project_for_create_hit_file(project_config)
    if threads is not None:
        config["threads"] = threads
    errors = validate_create_hit_file_inputs(config)
    if errors:
        raise ValueError("\n".join(errors))
    included = get_included_samples(config)
    by_name = {str(sample["sample"]): sample for sample in included}
    selected_names = list(by_name) if samples is None else samples
    unknown = [name for name in selected_names if name not in by_name]
    if unknown:
        raise ValueError(f"Selected samples are not included or do not exist: {', '.join(unknown)}")
    results = []
    for name in selected_names:
        result = run_create_hit_file_sample(by_name[name], config, threads, force, dry_run)
        results.append(result)
        _write_status(_paths(config, name)["status"], results)
        if result["status"] == "failed" and not keep_going:
            break
    return {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "threads": int(config["threads"]), "included_samples": len(included),
        "selected_samples": selected_names, "results": results,
        "input_bam_dir": str(_paths(config, selected_names[0])["project"] / "mapfastq") if selected_names else "",
        "output_dir": str(_paths(config, selected_names[0])["project"] / "create_hit_file") if selected_names else "",
        "success": sum(row["status"] == "success" for row in results),
        "failed": sum(row["status"] == "failed" for row in results),
        "skipped": sum(row["status"] == "skipped" for row in results),
        "dry_run": dry_run,
    }
