"""Restartable local orchestration for per-sample SummaryTable runs."""

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
    "sample", "status", "input_hits_file", "output_dir", "feature_tables",
    "log_file", "elapsed_seconds", "message",
]
STATS_MAP = {
    "File name": "sample",
    "Total Reads": "total_reads",
    "Total Hits": "total_hits",
    "% of hits in features": "percent_hits_in_features",
    "% of intergenic hits": "percent_intergenic_hits",
    "% of features hit": "percent_features_hit",
    "Mean hits per feature": "mean_hits_per_feature",
    "Mean reads per feature": "mean_reads_per_feature",
    "Mean Reads Per Hit": "mean_reads_per_hit",
    "Mean reads per hit in feature": "mean_reads_per_hit_in_feature",
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _resolve(value: str | Path | None, root: Path) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def load_project_for_summary_table(project_config: Path) -> dict:
    return load_project_for_mapping(project_config)


def get_included_samples(config: dict) -> list[dict]:
    return _get_included_samples(config)


def _reference_files(config: dict) -> dict[str, Path | None]:
    root = Path(config["_repo_root"])
    reference = config.get("reference") or {}
    return {
        "fasta": _resolve(reference.get("fasta"), root),
        "gff": _resolve(reference.get("gff"), root),
        "gtf": _resolve(reference.get("gtf"), root),
        "feature_table": _resolve(reference.get("feature_table"), root),
        "orthology_file": _resolve(reference.get("orthology_file"), root),
    }


def _paths(config: dict, sample_name: str) -> dict[str, Path]:
    root = Path(config["_repo_root"])
    project_dir = _resolve(config.get("output_project_dir"), root)
    export_dir = _resolve(config.get("output_export_dir"), root)
    if project_dir is None or export_dir is None:
        raise ValueError("Project and export output directories must be configured.")
    return {
        "project": project_dir,
        "output": project_dir / "summary" / sample_name,
        "log": project_dir / "logs" / "summary" / f"{sample_name}.log",
        "manifest": project_dir / "manifests" / "summary" / f"{sample_name}.summary_manifest.json",
        "status": project_dir / "manifests" / "summary" / "summary_status.csv",
        "combined_stats": project_dir / "summary" / "summary_stats.all_samples.csv",
        "qc": export_dir / "qc" / "summary",
    }


def find_sample_hits_file(sample: dict, config: dict) -> Path:
    name = str(sample["sample"])
    return _paths(config, name)["project"] / "create_hit_file" / name / f"{name}_hits.txt"


def validate_summary_table_inputs(config: dict) -> list[str]:
    errors: list[str] = []
    species = str(config.get("species") or "")
    if species not in {"albicans", "pombe", "cerevisiae", "glabrata"}:
        errors.append(f"SummaryTable does not yet provide a feature DB for species: {species}")
    try:
        threads = int(config.get("threads", 0))
    except (TypeError, ValueError):
        threads = 0
    if threads < 1:
        errors.append("Threads must be at least 1.")
    reference = _reference_files(config)
    annotation = reference["gff"] or reference["gtf"] or reference["feature_table"]
    if annotation is None or not annotation.is_file():
        errors.append("No feature table, GFF, or GTF is available. Fix the reference resources first.")
    if species == "glabrata":
        if reference["fasta"] is None or not reference["fasta"].is_file():
            errors.append("Reference FASTA is missing. Fix the reference resources first.")
        if reference["gff"] is None or not reference["gff"].is_file():
            errors.append("Glabrata GFF is missing. Fix the reference resources first.")
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
        hits = find_sample_hits_file(sample, config)
        if not hits.is_file():
            errors.append(f"Hit file missing for sample {name}: {hits}")
    return errors


def build_summary_table_command(
    sample: dict, config: dict, threads: int | None = None
) -> list[str]:
    selected_threads = int(config.get("threads") if threads is None else threads)
    if selected_threads < 1:
        raise ValueError("Threads must be at least 1.")
    name = str(sample["sample"])
    root = Path(config["_repo_root"])
    hits = find_sample_hits_file(sample, config)
    reference = _reference_files(config)
    command = [
        sys.executable, str(root / "src" / "ytab" / "summary" / "SummaryTable.py"),
        "--input-dir", str(hits.parent), "--hit-glob", hits.name,
        "--output-dir", str(_paths(config, name)["output"]),
        "--feature-db", str(config["species"]), "--threads", str(selected_threads),
        "--overwrite",
    ]
    if config.get("species") == "glabrata":
        command.extend([
            "--feature-gff", str(reference["gff"]),
            "--reference-fasta", str(reference["fasta"]),
        ])
        if reference["orthology_file"] and reference["orthology_file"].is_file():
            command.extend(["--scer-orthologs", str(reference["orthology_file"])])
    return command


def _feature_tables(output_dir: Path) -> list[Path]:
    patterns = ("*.feature_table*.csv", "*.feature_table*.tsv", "*.feature_table*.txt")
    return sorted({path for pattern in patterns for path in output_dir.glob(pattern) if path.is_file()}) if output_dir.is_dir() else []


def _write_manifest(path: Path, manifest: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def run_summary_table_sample(
    sample: dict, config: dict, threads: int | None = None,
    force: bool = False, dry_run: bool = False,
) -> dict:
    name = str(sample["sample"])
    paths = _paths(config, name)
    hits = find_sample_hits_file(sample, config)
    reference = _reference_files(config)
    selected_threads = int(config.get("threads") if threads is None else threads)
    command = build_summary_table_command(sample, config, selected_threads)
    started = time.monotonic()
    manifest = {
        "project_id": config.get("project_id"), "sample": name,
        "input_hits_file": str(hits), "species": config.get("species"),
        "reference_files_used": {key: str(value) if value else None for key, value in reference.items()},
        "output_dir": str(paths["output"]), "detected_feature_tables": [],
        "log_file": str(paths["log"]), "status": "failed",
        "start_time": _now(), "end_time": None, "elapsed_seconds": 0.0,
        "command_run": command, "warnings": [], "error_message": None,
    }
    existing_tables = _feature_tables(paths["output"])
    if not force and paths["manifest"].is_file() and existing_tables:
        try:
            previous = json.loads(paths["manifest"].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = {}
        if previous.get("status") in {"success", "skipped"}:
            manifest.update(status="skipped", end_time=_now())
            manifest["warnings"].append("Successful cached feature table found; SummaryTable skipped.")
            manifest["detected_feature_tables"] = [str(path) for path in existing_tables]
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
            log.flush()
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
        if completed.returncode != 0:
            raise RuntimeError(f"SummaryTable exited with status {completed.returncode}; see {paths['log']}")
        tables = _feature_tables(paths["output"])
        if not tables:
            raise RuntimeError("SummaryTable completed without an expected feature table.")
        sample_qc = paths["qc"] / name
        sample_qc.mkdir(parents=True, exist_ok=True)
        for source in tables + ([paths["output"] / "stats.csv"] if (paths["output"] / "stats.csv").is_file() else []):
            shutil.copy2(source, sample_qc / source.name)
        manifest["detected_feature_tables"] = [str(path) for path in tables]
        manifest["status"] = "success"
    except Exception as exc:
        manifest["error_message"] = str(exc)
    manifest["end_time"] = _now()
    manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
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
            "input_hits_file": result["input_hits_file"], "output_dir": result["output_dir"],
            "feature_tables": ";".join(result["detected_feature_tables"]),
            "log_file": result["log_file"], "elapsed_seconds": result["elapsed_seconds"],
            "message": result.get("error_message") or "; ".join(result.get("warnings") or []),
        }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerows(existing[name] for name in sorted(existing))


def collect_summary_stats(config: dict) -> Path | None:
    included = get_included_samples(config)
    rows = []
    for sample in included:
        name = str(sample["sample"])
        stats_path = _paths(config, name)["output"] / "stats.csv"
        if not stats_path.is_file():
            continue
        with stats_path.open(newline="", encoding="utf-8") as handle:
            for raw in csv.DictReader(handle):
                row = {target: raw.get(source, "") for source, target in STATS_MAP.items()}
                row["sample"] = name
                rows.append(row)
    if not rows:
        return None
    first = _paths(config, str(included[0]["sample"]))
    fields = list(STATS_MAP.values())
    path = first["combined_stats"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    first["qc"].mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, first["qc"] / path.name)
    return path


def run_summary_table_project(
    project_config: Path, samples: list[str] | None = None,
    threads: int | None = None, force: bool = False, dry_run: bool = False,
    keep_going: bool = False,
) -> dict:
    config = load_project_for_summary_table(project_config)
    if threads is not None:
        config["threads"] = threads
    errors = validate_summary_table_inputs(config)
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
        result = run_summary_table_sample(by_name[name], config, threads, force, dry_run)
        results.append(result)
        _write_status(_paths(config, name)["status"], results)
        if result["status"] == "failed" and not keep_going:
            break
    stats_path = None if dry_run else collect_summary_stats(config)
    project_dir = _paths(config, selected_names[0])["project"] if selected_names else Path()
    return {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "threads": int(config["threads"]), "included_samples": len(included),
        "selected_samples": selected_names, "results": results,
        "input_hits_dir": str(project_dir / "create_hit_file"),
        "output_dir": str(project_dir / "summary"),
        "combined_stats": str(stats_path) if stats_path else None,
        "success": sum(row["status"] == "success" for row in results),
        "failed": sum(row["status"] == "failed" for row in results),
        "skipped": sum(row["status"] == "skipped" for row in results),
        "dry_run": dry_run,
    }
