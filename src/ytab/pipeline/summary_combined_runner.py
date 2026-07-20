"""Run SummaryTable on one combined normalized parent hit file."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from .combine_hits_runner import load_project_for_combine_hits, resolve_combine_target


STATUS_FIELDS = [
    "target", "target_tag", "status", "input_combined_hits_file",
    "stable_combined_feature_table", "summary_stats_file", "output_dir",
    "log_file", "elapsed_seconds", "message",
]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _resolve(value: str | Path | None, root: Path) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def _paths(config: dict, target_tag: str | None = None) -> dict[str, Path]:
    root = Path(config["_repo_root"])
    project = _resolve(config.get("output_project_dir"), root)
    export_root = _resolve(config.get("output_export_dir"), root)
    if project is None or export_root is None:
        raise ValueError("Project and export output directories must be configured.")
    paths = {
        "project": project,
        "status": project / "manifests" / "summary_combined" / "summary_combined_status.csv",
    }
    if target_tag:
        output = project / "summary_combined" / target_tag
        paths.update({
            "output": output,
            "stable_feature": output / f"combined_feature_table.{target_tag}.txt",
            "stable_stats": output / f"combined_summary_stats.{target_tag}.csv",
            "log": project / "logs" / "summary_combined" / f"{target_tag}.log",
            "manifest": project / "manifests" / "summary_combined" / f"{target_tag}.summary_combined_manifest.json",
            "export": export_root / "qc" / "summary_combined" / target_tag,
        })
    return paths


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


def load_project_for_summary_combined(project_config: Path) -> dict:
    return load_project_for_combine_hits(project_config)


def resolve_summary_combined_target(config: dict, target: str = "recommended") -> dict:
    return resolve_combine_target(config, target)


def find_combined_hits_file(config: dict, target_tag: str) -> Path:
    return _paths(config)["project"] / "combined_hits" / target_tag / f"combined_parent_hits.{target_tag}.txt"


def validate_summary_combined_inputs(config: dict, target_info: dict) -> list[str]:
    errors: list[str] = []
    species = str(config.get("species") or "")
    if species not in {"albicans", "pombe", "cerevisiae", "glabrata"}:
        errors.append(f"SummaryTable does not yet provide a feature DB for species: {species}")
    hits = find_combined_hits_file(config, target_info["target_tag"])
    if not hits.is_file() or hits.stat().st_size == 0:
        errors.append(f"Combined parent hit file missing or empty: {hits}. Run Step 8 first.")
    reference = _reference_files(config)
    annotation = reference["gff"] or reference["gtf"] or reference["feature_table"]
    if annotation is None or not annotation.is_file():
        errors.append("No feature table, GFF, or GTF is available. Fix the reference resources first.")
    if species == "glabrata":
        if reference["fasta"] is None or not reference["fasta"].is_file():
            errors.append("Reference FASTA is missing. Fix the reference resources first.")
        if reference["gff"] is None or not reference["gff"].is_file():
            errors.append("Glabrata GFF is missing. Fix the reference resources first.")
    return errors


def build_summary_combined_command(config: dict, target_info: dict, threads: int | None = None) -> list[str]:
    selected_threads = int(config.get("threads") if threads is None else threads)
    if selected_threads < 1:
        raise ValueError("Threads must be at least 1.")
    root = Path(config["_repo_root"])
    tag = target_info["target_tag"]
    hits = find_combined_hits_file(config, tag)
    reference = _reference_files(config)
    command = [
        sys.executable, str(root / "src" / "ytab" / "summary" / "SummaryTable.py"),
        "--input-dir", str(hits.parent), "--hit-glob", hits.name,
        "--output-dir", str(_paths(config, tag)["output"]),
        "--feature-db", str(config["species"]), "--threads", str(selected_threads),
        "--overwrite",
    ]
    if config.get("species") == "glabrata":
        command.extend(["--feature-gff", str(reference["gff"]), "--reference-fasta", str(reference["fasta"])])
        if reference["orthology_file"] and reference["orthology_file"].is_file():
            command.extend(["--scer-orthologs", str(reference["orthology_file"])])
    return command


def detect_combined_feature_tables(output_dir: Path) -> list[Path]:
    patterns = ("*.feature_table*.csv", "*.feature_table*.tsv", "*.feature_table*.txt")
    found = {path for pattern in patterns for path in output_dir.glob(pattern) if path.is_file()}
    return sorted(path for path in found if not path.name.startswith("combined_feature_table.")) if output_dir.is_dir() else []


def choose_main_combined_feature_table(feature_tables: list[Path]) -> Path:
    if not feature_tables:
        raise ValueError("SummaryTable produced no combined feature table.")
    return sorted(
        feature_tables,
        key=lambda path: (
            0 if ".feature_table.RDF_1" in path.name else 1,
            0 if path.suffix.lower() == ".csv" else 1,
            path.name,
        ),
    )[0]


def collect_combined_summary_stats(output_dir: Path, destination: Path) -> Path | None:
    source = output_dir / "stats.csv"
    if not source.is_file():
        return None
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def _write_status(path: Path, manifest: dict) -> None:
    existing: dict[str, dict] = {}
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            existing = {row["target_tag"]: row for row in csv.DictReader(handle)}
    existing[manifest["target_tag"]] = {
        "target": manifest["target"], "target_tag": manifest["target_tag"],
        "status": manifest["status"], "input_combined_hits_file": manifest["input_combined_hits_file"],
        "stable_combined_feature_table": manifest.get("stable_combined_feature_table") or "",
        "summary_stats_file": manifest.get("summary_stats_file") or "",
        "output_dir": manifest["output_dir"], "log_file": manifest["log_file"],
        "elapsed_seconds": manifest["elapsed_seconds"],
        "message": manifest.get("error_message") or "; ".join(manifest.get("warnings") or []),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerows(existing[tag] for tag in sorted(existing))


def run_summary_combined_project(
    project_config: Path, target: str = "recommended", threads: int | None = None,
    force: bool = False, dry_run: bool = False,
) -> dict:
    config = load_project_for_summary_combined(project_config)
    target_info = resolve_summary_combined_target(config, target)
    selected_threads = int(config.get("threads") if threads is None else threads)
    if selected_threads < 1:
        raise ValueError("Threads must be at least 1.")
    errors = validate_summary_combined_inputs(config, target_info)
    if errors:
        raise ValueError("\n".join(errors))
    tag = target_info["target_tag"]
    paths = _paths(config, tag)
    hits = find_combined_hits_file(config, tag)
    references = _reference_files(config)
    command = build_summary_combined_command(config, target_info, selected_threads)
    started = time.monotonic()
    manifest = {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "target": target_info["target"], "target_tag": tag,
        "target_resolution_source": target_info["source"],
        "input_combined_hits_file": str(hits), "output_dir": str(paths["output"]),
        "stable_combined_feature_table": None, "summary_stats_file": None,
        "reference_files_used": {key: str(value) if value else None for key, value in references.items()},
        "log_file": str(paths["log"]), "status": "failed", "start_time": _now(),
        "end_time": None, "elapsed_seconds": 0.0, "command_run": command,
        "detected_feature_tables": [], "warnings": list(target_info["warnings"]), "error_message": None,
    }
    cached = False
    if not force and paths["manifest"].is_file() and paths["stable_feature"].is_file():
        try:
            previous = json.loads(paths["manifest"].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = {}
        cached = previous.get("status") in {"success", "skipped"} and previous.get("input_combined_hits_file") == str(hits)
    if cached:
        manifest["status"] = "skipped"
        manifest["warnings"].append("Successful cached combined feature table found; SummaryTable skipped.")
    elif dry_run:
        manifest["status"] = "skipped"
        manifest["warnings"].append("Dry run; command was not executed.")
    else:
        paths["output"].mkdir(parents=True, exist_ok=True)
        paths["log"].parent.mkdir(parents=True, exist_ok=True)
        try:
            with paths["log"].open("w", encoding="utf-8") as log:
                log.write("Command: " + " ".join(command) + "\n\n")
                log.flush()
                completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
            if completed.returncode != 0:
                raise RuntimeError(f"SummaryTable exited with status {completed.returncode}; see {paths['log']}")
            tables = detect_combined_feature_tables(paths["output"])
            main_table = choose_main_combined_feature_table(tables)
            if len(tables) > 1:
                manifest["warnings"].append(f"Multiple feature tables detected; selected {main_table.name} as the stable table.")
            shutil.copy2(main_table, paths["stable_feature"])
            stats = collect_combined_summary_stats(paths["output"], paths["stable_stats"])
            paths["export"].mkdir(parents=True, exist_ok=True)
            for source in tables + [paths["stable_feature"]] + ([stats] if stats else []):
                shutil.copy2(source, paths["export"] / source.name)
            manifest["status"] = "success"
        except Exception as exc:
            manifest["error_message"] = str(exc)
    tables = detect_combined_feature_tables(paths["output"])
    manifest["detected_feature_tables"] = [str(path) for path in tables]
    if paths["stable_feature"].is_file():
        manifest["stable_combined_feature_table"] = str(paths["stable_feature"])
    if paths["stable_stats"].is_file():
        manifest["summary_stats_file"] = str(paths["stable_stats"])
    manifest["end_time"] = _now()
    manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
    _write_json(paths["manifest"], manifest)
    _write_status(_paths(config)["status"], manifest)
    return {**manifest, "manifest_path": str(paths["manifest"]), "dry_run": dry_run}
