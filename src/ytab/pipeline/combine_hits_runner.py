"""Combine normalized parent hit files for one selected normalization target."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from .normalization_runner import (
    count_hit_sites_and_reads,
    get_normalization_samples,
    load_project_for_normalization,
    normalize_target_tag,
)


STATUS_FIELDS = [
    "target", "target_tag", "status", "selected_samples", "input_hit_files_count",
    "combined_hits_file", "output_dir", "log_file", "elapsed_seconds", "message",
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
        "normalization": project / "sample_normalization",
        "combined_base": project / "combined_hits",
        "manifest_base": project / "manifests" / "combined_hits",
        "status": project / "manifests" / "combined_hits" / "combined_hits_status.csv",
        "export_base": export_root / "qc" / "combined_hits",
    }
    if target_tag:
        output = project / "combined_hits" / target_tag
        stem = f"combined_parent_hits.{target_tag}"
        paths.update({
            "output": output,
            "combined": output / f"{stem}.txt",
            "metadata": output / f"{stem}.metadata.json",
            "utility_summary": output / f"{stem}.combine_summary.csv",
            "utility_manifest": output / f"{stem}.combine_manifest.json",
            "log": project / "logs" / "combined_hits" / f"{target_tag}.log",
            "manifest": project / "manifests" / "combined_hits" / f"{target_tag}.combined_hits_manifest.json",
            "export": export_root / "qc" / "combined_hits" / target_tag,
        })
    return paths


def load_project_for_combine_hits(project_config: Path) -> dict:
    return load_project_for_normalization(project_config)


def get_parent_samples(config: dict, samples: list[str] | None = None) -> list[dict]:
    selected = get_normalization_samples(config, "parents", samples)
    if samples is not None:
        treated = [
            str(row["sample"]) for row in selected
            if str(row.get("condition") or row.get("guessed_condition") or row.get("treatment") or "").strip().lower()
            in {"treated", "treatment"}
        ]
        if treated:
            raise ValueError("Explicit combine samples cannot include treated libraries: " + ", ".join(treated))
    return selected


def _read_recommendation(path: Path, label: str) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"{label} recommendation is missing: {path}")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid recommendation JSON: {path}: {exc}") from exc
    tag = str(data.get("recommended_target_tag") or "")
    if not tag:
        raise ValueError(f"Recommendation lacks recommended_target_tag: {path}")
    return {"target": float(data["recommended_target"]), "target_tag": tag}


def resolve_combine_target(config: dict, target: str = "recommended") -> dict:
    value = str(target).strip()
    if not value:
        raise ValueError("Combine target cannot be empty.")
    normalization = _paths(config)["normalization"]
    mode = value.lower()
    warnings: list[str] = []
    feature_path = normalization / "normalization_feature_recommendation.json"
    site_path = normalization / "normalization_recommendation.json"
    if mode in {"recommended", "auto", "feature-recommended"}:
        if feature_path.is_file():
            info = _read_recommendation(feature_path, "Feature-level")
            source = str(feature_path)
        else:
            info = _read_recommendation(site_path, "Site-retention")
            source = str(site_path)
            warnings.append("Feature-level recommendation is missing; using the Step 6 site-retention recommendation.")
    elif mode == "site-recommended":
        info = _read_recommendation(site_path, "Site-retention")
        source = str(site_path)
    else:
        tag = value if value.startswith("T") else normalize_target_tag(value)
        try:
            numeric = float(tag[1:].replace("p", "."))
        except ValueError:
            raise ValueError(f"Invalid normalization target: {value}") from None
        if normalize_target_tag(numeric) != tag:
            raise ValueError(f"Invalid normalization target tag: {tag}")
        info, source = {"target": numeric, "target_tag": tag}, "manual"
    target_dir = normalization / info["target_tag"]
    if not target_dir.is_dir():
        raise FileNotFoundError(f"Normalization target directory is missing: {target_dir}. Run Step 6 first.")
    return {**info, "target_dir": str(target_dir), "source": source, "warnings": warnings}


def find_normalized_parent_hits(
    config: dict, target_tag: str, selected_samples: list[dict]
) -> list[Path]:
    base = _paths(config)["normalization"] / target_tag
    return [base / str(row["sample"]) / f"{row['sample']}_normalized_hits.txt" for row in selected_samples]


def validate_combine_hits_inputs(
    config: dict, target_info: dict, selected_samples: list[dict]
) -> list[str]:
    errors: list[str] = []
    if len(selected_samples) < 2:
        errors.append("At least two parent samples are required to build a combined parent library.")
    for sample, path in zip(selected_samples, find_normalized_parent_hits(config, target_info["target_tag"], selected_samples)):
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"Normalized parent hit file missing or empty for {sample['sample']}: {path}. Run Step 6 first.")
    return errors


def build_combine_hits_command(
    config: dict, target_info: dict, selected_samples: list[dict]
) -> list[str]:
    root = Path(config["_repo_root"])
    paths = _paths(config, target_info["target_tag"])
    sample_name = f"combined_parent_hits.{target_info['target_tag']}"
    return [
        sys.executable, str(root / "src" / "ytab" / "insertions" / "CombineHitFiles.py"),
        "--hits-txt", *[str(path) for path in find_normalized_parent_hits(config, target_info["target_tag"], selected_samples)],
        "--outdir", str(paths["output"]), "--sample-name", sample_name,
        "--output-file", paths["combined"].name, "--sort",
    ]


def collect_combined_hits_metadata(combined_hits_file: Path) -> dict:
    counts = count_hit_sites_and_reads(Path(combined_hits_file))
    return {"total_combined_sites": counts["hit_sites"], "total_combined_reads": counts["total_reads"]}


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
        "status": manifest["status"], "selected_samples": ";".join(manifest["selected_samples"]),
        "input_hit_files_count": len(manifest["input_normalized_hit_files"]),
        "combined_hits_file": manifest["combined_hits_file"], "output_dir": manifest["output_dir"],
        "log_file": manifest["log_file"], "elapsed_seconds": manifest["elapsed_seconds"],
        "message": manifest.get("error_message") or "; ".join(manifest.get("warnings") or []),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerows(existing[tag] for tag in sorted(existing))


def run_combine_hits_project(
    project_config: Path, target: str = "recommended", samples: list[str] | None = None,
    force: bool = False, dry_run: bool = False,
) -> dict:
    config = load_project_for_combine_hits(project_config)
    selected = get_parent_samples(config, samples)
    target_info = resolve_combine_target(config, target)
    errors = validate_combine_hits_inputs(config, target_info, selected)
    if errors:
        raise ValueError("\n".join(errors))
    paths = _paths(config, target_info["target_tag"])
    input_files = find_normalized_parent_hits(config, target_info["target_tag"], selected)
    command = build_combine_hits_command(config, target_info, selected)
    started = time.monotonic()
    manifest = {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "target": target_info["target"], "target_tag": target_info["target_tag"],
        "selected_samples": [str(row["sample"]) for row in selected],
        "input_normalized_hit_files": [str(path) for path in input_files],
        "output_dir": str(paths["output"]), "combined_hits_file": str(paths["combined"]),
        "log_file": str(paths["log"]), "status": "failed", "start_time": _now(),
        "end_time": None, "elapsed_seconds": 0.0, "command_run": command,
        "detected_outputs": [], "warnings": list(target_info["warnings"]), "error_message": None,
    }
    cached = False
    if not force and paths["manifest"].is_file() and paths["combined"].is_file():
        try:
            previous = json.loads(paths["manifest"].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = {}
        cached = previous.get("status") in {"success", "skipped"} and previous.get("selected_samples") == manifest["selected_samples"]
    if cached:
        manifest["status"] = "skipped"
        manifest["warnings"].append("Successful cached combined parent hit file found; combine skipped.")
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
                raise RuntimeError(f"CombineHitFiles exited with status {completed.returncode}; see {paths['log']}")
            if not paths["combined"].is_file() or paths["combined"].stat().st_size == 0:
                raise RuntimeError("CombineHitFiles completed without the expected combined hit file.")
            counts = collect_combined_hits_metadata(paths["combined"])
            metadata = {
                "project_id": config.get("project_id"), "species": config.get("species"),
                "target": target_info["target"], "target_tag": target_info["target_tag"],
                "combined_sample_name": f"combined_parent_hits.{target_info['target_tag']}",
                "selected_parent_samples": manifest["selected_samples"],
                "input_normalized_hit_files": manifest["input_normalized_hit_files"],
                "output_combined_hits_file": str(paths["combined"]),
                "total_input_files": len(input_files), **counts, "timestamp": _now(),
                "warnings": manifest["warnings"],
            }
            _write_json(paths["metadata"], metadata)
            paths["export"].mkdir(parents=True, exist_ok=True)
            for source in (paths["combined"], paths["metadata"], paths["utility_summary"], paths["utility_manifest"]):
                if source.is_file():
                    shutil.copy2(source, paths["export"] / source.name)
            manifest["status"] = "success"
        except Exception as exc:
            manifest["error_message"] = str(exc)
    manifest["end_time"] = _now()
    manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
    manifest["detected_outputs"] = sorted(str(path) for path in paths["output"].glob("*") if path.is_file()) if paths["output"].is_dir() else []
    if paths["combined"].is_file():
        manifest.update(collect_combined_hits_metadata(paths["combined"]))
    if not dry_run:
        _write_json(paths["manifest"], manifest)
        _write_status(_paths(config)["status"], manifest)
    return {**manifest, "manifest_path": str(paths["manifest"]), "metadata_path": str(paths["metadata"]), "dry_run": dry_run}
