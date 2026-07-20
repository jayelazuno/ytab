"""Restartable local runner for essentiality classification of a parent library."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from .summary_combined_runner import load_project_for_summary_combined, resolve_summary_combined_target


REQUIRED_COLUMNS = {
    "standard_name", "coding_length", "hits", "reads", "neighborhood_index",
    "freedom_index", "insertion_index", "upstream_hits_100",
}
STATUS_FIELDS = [
    "target", "target_tag", "status", "input_combined_feature_table",
    "input_feature_count", "stable_prediction_table", "output_prediction_count",
    "output_dir", "log_file", "elapsed_seconds", "message",
]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _resolve(value: str | Path | None, root: Path) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _paths(config: dict, tag: str | None = None) -> dict[str, Path]:
    root = Path(config["_repo_root"])
    project = _resolve(config.get("output_project_dir"), root)
    export = _resolve(config.get("output_export_dir"), root)
    if project is None or export is None:
        raise ValueError("Project and export output directories must be configured.")
    paths = {
        "project": project,
        "status": project / "manifests" / "classifier" / "classifier_status.csv",
        "final_txt": project / "config" / "final_classifier_target.txt",
        "final_json": project / "config" / "final_classifier_target.json",
    }
    if tag:
        output = project / "classifier" / tag
        paths.update({
            "output": output,
            "stable": output / f"essentiality_predictions.{tag}.csv",
            "summary": output / f"classifier_summary.{tag}.csv",
            "metadata": output / f"classifier_run_metadata.{tag}.json",
            "log": project / "logs" / "classifier" / f"{tag}.log",
            "manifest": project / "manifests" / "classifier" / f"{tag}.classifier_manifest.json",
            "export": export / "classifier" / tag,
        })
    return paths


def load_project_for_classifier(project_config: Path) -> dict:
    return load_project_for_summary_combined(project_config)


def resolve_classifier_target(config: dict, target: str = "recommended") -> dict:
    return resolve_summary_combined_target(config, target)


def find_combined_feature_table(config: dict, target_tag: str) -> Path:
    project = _paths(config)["project"]
    return project / "summary_combined" / target_tag / f"combined_feature_table.{target_tag}.txt"


def resolve_classifier_resources(config: dict) -> dict:
    root = Path(config["_repo_root"])
    reference = config.get("reference") or {}
    project = _paths(config)["project"]
    resolved_path = project / "config" / "reference_resolved.json"
    resolved = {}
    if resolved_path.is_file():
        try:
            resolved = json.loads(resolved_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            resolved = {}
    species_root = root / "resources" / "species"
    orthology = _resolve(resolved.get("orthology_file") or reference.get("orthology_file"), root)
    resources = {
        "orthology_file": orthology,
        "kornmann_wig_files": sorted((species_root / "Kornmann").glob("*WildType*.wig")),
        "cerevisiae_feature_table": species_root / "cerevisiae" / "SGD_features.tab",
        "cerevisiae_viable_annotations": species_root / "cerevisiae" / "cerevisiae_viable_annotations.txt",
        "cerevisiae_inviable_annotations": species_root / "cerevisiae" / "cerevisiae_inviable_annotations.txt",
        "cerevisiae_paralogs": species_root / "cerevisiae" / "hasParalogs_sc.txt",
    }
    return resources


def _read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        first = handle.readline()
        skip_rdf = first.strip().startswith("RDF")
        if not skip_rdf:
            handle.seek(0)
        sample = handle.readline()
        delimiter = "\t" if "\t" in sample and "," not in sample else ","
        handle.seek(0)
        if skip_rdf:
            handle.readline()
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = list(reader)
        return [str(name).strip() for name in (reader.fieldnames or [])], rows


def inspect_classifier_requirements(config: dict, feature_table: Path, resources: dict) -> dict:
    columns, rows = _read_table(feature_table)
    missing = sorted(REQUIRED_COLUMNS - set(columns))
    invalid_rows = 0
    for row in rows:
        try:
            for name in REQUIRED_COLUMNS - {"standard_name"}:
                float(row.get(name, ""))
            if not str(row.get("standard_name") or "").strip():
                raise ValueError
        except (TypeError, ValueError):
            invalid_rows += 1
    return {
        "required_columns": sorted(REQUIRED_COLUMNS), "columns": columns,
        "missing_columns": missing, "input_feature_count": len(rows),
        "invalid_feature_rows": invalid_rows,
    }


def validate_classifier_inputs(config: dict, target_info: dict, feature_table: Path, resources: dict) -> list[str]:
    errors: list[str] = []
    if str(config.get("species") or "") != "glabrata":
        errors.append("Classifier table mode currently supports species 'glabrata' only.")
    if not feature_table.is_file() or feature_table.stat().st_size == 0:
        errors.append(f"Combined feature table missing or empty: {feature_table}. Run Steps 8 and 9 for {target_info['target_tag']} first.")
        return errors
    inspection = inspect_classifier_requirements(config, feature_table, resources)
    if inspection["missing_columns"]:
        errors.append("Combined feature table is missing required classifier columns: " + ", ".join(inspection["missing_columns"]))
    if inspection["invalid_feature_rows"]:
        errors.append(f"Combined feature table has {inspection['invalid_feature_rows']} rows with missing or non-numeric required classifier values.")
    for name in ("orthology_file", "cerevisiae_feature_table", "cerevisiae_viable_annotations", "cerevisiae_inviable_annotations", "cerevisiae_paralogs"):
        path = resources.get(name)
        if not isinstance(path, Path) or not path.is_file():
            errors.append(f"Required classifier resource is missing ({name}): {path}")
    wigs = resources.get("kornmann_wig_files") or []
    if not wigs:
        errors.append("Required Kornmann WildType training WIG files were not found.")
    elif any(not path.is_file() for path in wigs):
        errors.append("One or more Kornmann WildType training WIG files are missing.")
    return errors


def build_classifier_command(config: dict, target_info: dict, feature_table: Path, resources: dict, seed: int | None = None) -> list[str]:
    root = Path(config["_repo_root"])
    command = [
        sys.executable, str(root / "src" / "ytab" / "essentiality" / "Classifier.py"),
        "--mode", "scer-train-gla-classify", "--gla-feature-table", str(feature_table),
        "--out-dir", str(_paths(config, target_info["target_tag"])["output"]),
        "--threads", str(int(config.get("threads", 1))), "--target-fpr", "0.10",
        "--orthology-file", str(resources["orthology_file"]), "--combine", "--overwrite",
    ]
    if seed is not None:
        command.extend(["--seed", str(int(seed))])
    return command


def detect_classifier_outputs(output_dir: Path) -> list[Path]:
    if not output_dir.is_dir():
        return []
    stable_prefixes = ("essentiality_predictions.", "classifier_summary.", "classifier_run_metadata.")
    return sorted(path for path in output_dir.rglob("*") if path.is_file() and not path.name.startswith(stable_prefixes))


def choose_primary_prediction_table(outputs: list[Path]) -> Path:
    candidates = [path for path in outputs if path.name.endswith(".predictions.csv")]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        candidates = [path for path in outputs if path.name == "combined_glabrata_RF-G4.tsv"]
    if len(candidates) != 1:
        names = ", ".join(str(path) for path in candidates) or "none"
        raise ValueError(f"Could not identify one unambiguous classifier prediction table; candidates: {names}")
    return candidates[0]


def _prediction_info(path: Path) -> dict:
    columns, rows = _read_table(path)
    label_columns = [name for name in columns if "ess. for FPR" in name or name.lower().endswith("verdict")]
    probability_columns = [name for name in columns if name == "RF - G4" or "probab" in name.lower() or "score" in name.lower()]
    label = label_columns[0] if len(label_columns) == 1 else None
    probability = probability_columns[0] if len(probability_columns) == 1 else None
    counts: dict[str, int] = {}
    if label:
        for row in rows:
            value = str(row.get(label) or "")
            counts[value] = counts.get(value, 0) + 1
    return {"columns": columns, "rows": rows, "label_column": label, "probability_column": probability, "class_counts": counts}


def summarize_classifier_predictions(prediction_table: Path, target_info: dict) -> Path:
    info = _prediction_info(prediction_table)
    values = []
    if info["probability_column"]:
        for row in info["rows"]:
            try:
                values.append(float(row[info["probability_column"]]))
            except (TypeError, ValueError):
                pass
    values.sort()
    mean = sum(values) / len(values) if values else ""
    median = values[len(values) // 2] if len(values) % 2 else ((values[len(values) // 2 - 1] + values[len(values) // 2]) / 2 if values else "")
    labels = info["class_counts"] or {"": len(info["rows"])}
    fields = ["target", "target_tag", "input_feature_count", "output_prediction_count", "prediction_label", "prediction_count", "prediction_fraction", "probability_column", "mean_probability", "median_probability", "minimum_probability", "maximum_probability"]
    destination = prediction_table.parent / f"classifier_summary.{target_info['target_tag']}.csv"
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        total = len(info["rows"])
        for label, count in sorted(labels.items()):
            writer.writerow({
                "target": target_info["target"], "target_tag": target_info["target_tag"],
                "input_feature_count": total, "output_prediction_count": total,
                "prediction_label": label, "prediction_count": count,
                "prediction_fraction": count / total if total else "",
                "probability_column": info["probability_column"] or "", "mean_probability": mean,
                "median_probability": median, "minimum_probability": min(values) if values else "",
                "maximum_probability": max(values) if values else "",
            })
    return destination


def _serial_resources(resources: dict, hashes: bool = False) -> dict:
    result = {}
    for key, value in resources.items():
        paths = value if isinstance(value, list) else [value]
        result[key] = [{"path": str(path), "sha256": _sha256(path) if hashes and path.is_file() else None} for path in paths]
    return result


def _package_versions() -> dict:
    versions = {}
    for module in ("numpy", "pandas", "sklearn"):
        try:
            imported = __import__(module)
            versions[module] = getattr(imported, "__version__", "unknown")
        except ImportError:
            versions[module] = None
    return versions


def write_classifier_run_metadata(config: dict, target_info: dict, feature_table: Path, resources: dict, outputs: list[Path], seed: int | None = None) -> Path:
    paths = _paths(config, target_info["target_tag"])
    versions = _package_versions()
    data = {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "target": target_info["target"], "target_tag": target_info["target_tag"],
        "target_resolution_source": target_info["source"], "input_combined_feature_table": str(feature_table),
        "input_sha256": _sha256(feature_table), "classifier_resources": _serial_resources(resources, hashes=True),
        "detected_outputs": [str(path) for path in outputs], "seed": 0 if seed is None else seed,
        "python_version": platform.python_version(), "package_versions": versions, "timestamp": _now(),
    }
    paths["metadata"].write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return paths["metadata"]


def persist_final_classifier_target(config: dict, target_info: dict, feature_table: Path, warning: str | None = None) -> dict:
    paths = _paths(config)
    smoke_warning = "This target was selected during workflow testing with a shallow FASTQ subset and should not be interpreted as a biological normalization recommendation."
    selected_warning = warning or (smoke_warning if "smoke" in str(config.get("project_id", "")).lower() else "")
    data = {
        "project_id": config.get("project_id"), "target": target_info["target"], "target_tag": target_info["target_tag"],
        "selection_source": target_info["source"], "feature_table": str(feature_table),
        "selected_by": "ytab_run_classifier --save-final-target", "timestamp": _now(), "warning": selected_warning,
    }
    paths["final_txt"].write_text(target_info["target_tag"] + "\n", encoding="utf-8")
    paths["final_json"].write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return {"text": str(paths["final_txt"]), "json": str(paths["final_json"])}


def _write_status(path: Path, manifest: dict) -> None:
    existing = {}
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            existing = {row["target_tag"]: row for row in csv.DictReader(handle)}
    existing[manifest["target_tag"]] = {
        "target": manifest["target"], "target_tag": manifest["target_tag"], "status": manifest["status"],
        "input_combined_feature_table": manifest["input_combined_feature_table"], "input_feature_count": manifest["input_feature_count"],
        "stable_prediction_table": manifest.get("stable_prediction_table") or "", "output_prediction_count": manifest.get("output_prediction_count") or "",
        "output_dir": manifest["output_dir"], "log_file": manifest["log_file"], "elapsed_seconds": manifest["elapsed_seconds"],
        "message": manifest.get("error_message") or "; ".join(manifest.get("warnings") or []),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS); writer.writeheader()
        writer.writerows(existing[tag] for tag in sorted(existing))


def run_classifier_project(project_config: Path, target: str = "recommended", seed: int | None = None, force: bool = False, dry_run: bool = False, save_final_target: bool = False) -> dict:
    config = load_project_for_classifier(project_config)
    target_info = resolve_classifier_target(config, target)
    feature_table = find_combined_feature_table(config, target_info["target_tag"])
    resources = resolve_classifier_resources(config)
    errors = validate_classifier_inputs(config, target_info, feature_table, resources)
    if errors:
        raise ValueError("\n".join(errors))
    inspection = inspect_classifier_requirements(config, feature_table, resources)
    input_hash = _sha256(feature_table)
    paths = _paths(config, target_info["target_tag"])
    command = build_classifier_command(config, target_info, feature_table, resources, seed)
    resource_snapshot = _serial_resources(resources, hashes=True)
    started = time.monotonic()
    manifest = {
        "project_id": config.get("project_id"), "species": config.get("species"), "target": target_info["target"],
        "target_tag": target_info["target_tag"], "target_resolution_source": target_info["source"],
        "input_combined_feature_table": str(feature_table), "input_sha256": input_hash,
        "input_feature_count": inspection["input_feature_count"], "classifier_script": command[1],
        "classifier_resources": resource_snapshot, "output_dir": str(paths["output"]), "stable_prediction_table": None,
        "detected_outputs": [], "classifier_summary_file": None, "run_metadata_file": None,
        "log_file": str(paths["log"]), "seed": 0 if seed is None else seed,
        "python_version": platform.python_version(), "numpy_version": _package_versions()["numpy"],
        "pandas_version": _package_versions()["pandas"], "scikit_learn_version": _package_versions()["sklearn"],
        "status": "failed", "start_time": _now(), "end_time": None,
        "elapsed_seconds": 0.0, "command_run": command, "warnings": list(target_info["warnings"]), "error_message": None,
        "output_prediction_count": None, "class_counts": {}, "excluded_feature_count": 0,
    }
    cached = False
    if not force and paths["manifest"].is_file() and paths["stable"].is_file():
        try:
            previous = json.loads(paths["manifest"].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = {}
        cached = previous.get("status") == "success" and previous.get("target_tag") == target_info["target_tag"] and previous.get("input_sha256") == input_hash and previous.get("classifier_resources") == resource_snapshot
    if cached:
        manifest["status"] = "skipped"
        manifest["warnings"].append("Successful cached classifier result matches the current input and resources; classifier skipped.")
    elif dry_run:
        manifest["status"] = "skipped"
        manifest["warnings"].append("Dry run; classifier was not executed.")
    else:
        paths["log"].parent.mkdir(parents=True, exist_ok=True)
        try:
            with paths["log"].open("w", encoding="utf-8") as log:
                log.write("Command: " + " ".join(command) + "\n\n"); log.flush()
                completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
            if completed.returncode != 0:
                raise RuntimeError(f"Classifier exited with status {completed.returncode}; see {paths['log']}")
            outputs = detect_classifier_outputs(paths["output"])
            primary = choose_primary_prediction_table(outputs)
            shutil.copy2(primary, paths["stable"])
            summary = summarize_classifier_predictions(paths["stable"], target_info)
            metadata = write_classifier_run_metadata(config, target_info, feature_table, resources, outputs, seed)
            paths["export"].mkdir(parents=True, exist_ok=True)
            for source in outputs + [paths["stable"], summary, metadata]:
                destination = paths["export"] / source.relative_to(paths["output"])
                destination.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, destination)
            manifest["status"] = "success"
        except Exception as exc:
            manifest["error_message"] = str(exc)
    outputs = detect_classifier_outputs(paths["output"])
    manifest["detected_outputs"] = [str(path) for path in outputs]
    if paths["stable"].is_file():
        info = _prediction_info(paths["stable"])
        manifest["stable_prediction_table"] = str(paths["stable"])
        manifest["output_prediction_count"] = len(info["rows"])
        manifest["class_counts"] = info["class_counts"]
    if paths["summary"].is_file(): manifest["classifier_summary_file"] = str(paths["summary"])
    if paths["metadata"].is_file(): manifest["run_metadata_file"] = str(paths["metadata"])
    if save_final_target and not dry_run and manifest["status"] in {"success", "skipped"}:
        manifest["final_target_files"] = persist_final_classifier_target(config, target_info, feature_table)
    manifest["end_time"] = _now(); manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
    paths["manifest"].parent.mkdir(parents=True, exist_ok=True)
    if not cached:
        paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    _write_status(_paths(config)["status"], manifest)
    return {**manifest, "manifest_path": str(paths["manifest"]), "classifier_resources_resolved": _serial_resources(resources), "save_final_target": save_final_target, "dry_run": dry_run}
