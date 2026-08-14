"""Run SummaryTable on existing normalized hit files and evaluate targets."""

from __future__ import annotations

import csv
import json
import math
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
from .summary_table_runner import STATS_MAP


STATUS_FIELDS = [
    "target", "target_tag", "sample", "condition", "status",
    "input_normalized_hits_file", "output_dir", "feature_tables", "log_file",
    "elapsed_seconds", "message",
]
EVALUATION_FIELDS = [
    "target", "target_tag", "sample", "condition", "raw_hit_sites",
    "normalized_hit_sites", "hit_site_retention_fraction", "raw_total_reads",
    "normalized_total_reads", "raw_percent_features_hit",
    "normalized_percent_features_hit", "feature_retention_fraction",
    "feature_loss_fraction", "raw_total_hits", "normalized_total_hits",
    "raw_mean_hits_per_feature", "normalized_mean_hits_per_feature",
    "normalized_hits_file", "normalized_feature_table", "status", "warning",
]
SUMMARY_FIELDS = [
    "target", "target_tag", "sample_count", "min_hit_site_retention_fraction",
    "mean_hit_site_retention_fraction", "min_feature_retention_fraction",
    "mean_feature_retention_fraction", "max_feature_loss_fraction",
    "mean_normalized_percent_features_hit", "status",
    "recommended_by_site_retention", "recommended_by_feature_retention",
]
RAW_WARNING = "Raw feature-level metrics unavailable; run raw SummaryTable first."
SMOKE_WARNING = "The 10k smoke-subset recommendation is a workflow test, not a biological recommendation."


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _resolve(value: str | Path | None, root: Path) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def _condition(sample: dict) -> str:
    value = str(sample.get("condition") or sample.get("guessed_condition") or "").strip().lower()
    if value:
        return value
    return "parent" if "parent" in str(sample.get("sample") or "").lower() else "unknown"


def _paths(config: dict, target_tag: str | None = None, sample: str | None = None) -> dict[str, Path]:
    root = Path(config["_repo_root"])
    project = _resolve(config.get("output_project_dir"), root)
    export_root = _resolve(config.get("output_export_dir"), root)
    if project is None or export_root is None:
        raise ValueError("Project and export output directories must be configured.")
    normalization = project / "sample_normalization"
    norm_export = export_root / "qc" / "sample_normalization"
    paths = {
        "project": project, "normalization": normalization,
        "summary_base": project / "summary_normalized",
        "manifest_base": project / "manifests" / "summary_normalized",
        "status": project / "manifests" / "summary_normalized" / "summary_normalized_status.csv",
        "summary_export": export_root / "qc" / "summary_normalized",
        "evaluation": normalization / "normalization_target_evaluation.csv",
        "target_summary": normalization / "normalization_target_summary.csv",
        "feature_recommendation": normalization / "normalization_feature_recommendation.json",
        "evaluation_export": norm_export / "normalization_target_evaluation.csv",
        "target_summary_export": norm_export / "normalization_target_summary.csv",
        "feature_recommendation_export": norm_export / "normalization_feature_recommendation.csv",
    }
    if target_tag:
        paths.update({
            "target_summary": project / "summary_normalized" / target_tag,
            "target_export": export_root / "qc" / "summary_normalized" / target_tag,
        })
    if target_tag and sample:
        paths.update({
            "output": project / "summary_normalized" / target_tag / sample,
            "log": project / "logs" / "summary_normalized" / target_tag / f"{sample}.log",
            "manifest": project / "manifests" / "summary_normalized" / target_tag / f"{sample}.summary_normalized_manifest.json",
        })
    return paths


def load_project_for_summary_normalized(project_config: Path) -> dict:
    return load_project_for_normalization(project_config)


def get_summary_normalized_samples(
    config: dict, sample_mode: str = "parents", samples: list[str] | None = None
) -> list[dict]:
    return get_normalization_samples(config, sample_mode, samples)


def _tag_value(tag: str) -> float:
    if not tag.startswith("T"):
        raise ValueError(f"Invalid normalization target tag: {tag}")
    try:
        value = float(tag[1:].replace("p", "."))
    except ValueError:
        raise ValueError(f"Invalid normalization target tag: {tag}") from None
    if not math.isfinite(value) or value <= 0 or normalize_target_tag(value) != tag:
        raise ValueError(f"Invalid normalization target tag: {tag}")
    return value


def resolve_normalization_targets(
    config: dict, targets: str | list[str] = "recommended"
) -> list[dict]:
    base = _paths(config)["normalization"]
    values = [part.strip() for part in targets.split(",")] if isinstance(targets, str) else [str(part).strip() for part in targets]
    if not values or any(not value for value in values):
        raise ValueError("Normalization target selection cannot be empty.")
    mode = values[0].lower() if len(values) == 1 else "manual"
    if mode in {"recommended", "auto"}:
        recommendation = base / "normalization_recommendation.json"
        if not recommendation.is_file():
            raise FileNotFoundError(
                f"Normalization recommendation is missing: {recommendation}. Run Step 6 first."
            )
        data = json.loads(recommendation.read_text(encoding="utf-8"))
        tag = str(data.get("recommended_target_tag") or "")
        if not tag:
            raise ValueError(f"Recommendation does not contain recommended_target_tag: {recommendation}")
        values = [tag]
    elif mode == "all":
        values = [path.name for path in sorted(base.glob("T*")) if path.is_dir()]
        if not values:
            raise FileNotFoundError(f"No normalization target directories found under: {base}")
    infos: list[dict] = []
    seen: set[str] = set()
    for value in values:
        tag = value if value.startswith("T") else normalize_target_tag(value)
        numeric = _tag_value(tag)
        if tag in seen:
            raise ValueError(f"Duplicate normalization target: {tag}")
        seen.add(tag)
        infos.append({"target": numeric, "target_tag": tag, "normalization_dir": str(base / tag)})
    return infos


def find_normalized_hits_file(config: dict, target_tag: str, sample: dict) -> Path:
    name = str(sample["sample"])
    return _paths(config)["normalization"] / target_tag / name / f"{name}_normalized_hits.txt"


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


def validate_summary_normalized_inputs(
    config: dict, target_infos: list[dict], selected_samples: list[dict]
) -> list[str]:
    errors: list[str] = []
    reference = _reference_files(config)
    annotation = reference["gff"] or reference["gtf"] or reference["feature_table"]
    if annotation is None or not annotation.is_file():
        errors.append("No feature table, GFF, or GTF is available. Fix reference resources first.")
    if config.get("species") == "glabrata":
        if not reference["fasta"] or not reference["fasta"].is_file():
            errors.append("Reference FASTA is missing. Fix reference resources first.")
        if not reference["gff"] or not reference["gff"].is_file():
            errors.append("Glabrata GFF is missing. Fix reference resources first.")
    for info in target_infos:
        for sample in selected_samples:
            path = find_normalized_hits_file(config, info["target_tag"], sample)
            if not path.is_file():
                errors.append(f"Normalized hit file missing for {info['target_tag']}/{sample['sample']}: {path}. Run Step 6 first.")
    return errors


def build_summary_normalized_command(
    config: dict, target_info: dict, sample: dict, threads: int | None = None
) -> list[str]:
    selected_threads = int(config.get("threads") if threads is None else threads)
    if selected_threads < 1:
        raise ValueError("Threads must be at least 1.")
    root = Path(config["_repo_root"])
    name = str(sample["sample"])
    hits = find_normalized_hits_file(config, target_info["target_tag"], sample)
    output = _paths(config, target_info["target_tag"], name)["output"]
    reference = _reference_files(config)
    command = [
        sys.executable, str(root / "src" / "ytab" / "summary" / "SummaryTable.py"),
        "--input-dir", str(hits.parent), "--hit-glob", hits.name,
        "--output-dir", str(output), "--feature-db", str(config["species"]),
        "--threads", str(selected_threads), "--overwrite",
    ]
    if config.get("species") == "glabrata":
        command.extend(["--feature-gff", str(reference["gff"]), "--reference-fasta", str(reference["fasta"])])
        if reference["orthology_file"] and reference["orthology_file"].is_file():
            command.extend(["--scer-orthologs", str(reference["orthology_file"])])
    return command


def _feature_tables(output: Path) -> list[Path]:
    patterns = ("*.feature_table*.csv", "*.feature_table*.tsv", "*.feature_table*.txt")
    return sorted({path for pattern in patterns for path in output.glob(pattern) if path.is_file()}) if output.is_dir() else []


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def run_summary_normalized_sample(
    config: dict, target_info: dict, sample: dict, threads: int | None = None,
    force: bool = False, dry_run: bool = False,
) -> dict:
    name = str(sample["sample"])
    tag = target_info["target_tag"]
    paths = _paths(config, tag, name)
    hits = find_normalized_hits_file(config, tag, sample)
    command = build_summary_normalized_command(config, target_info, sample, threads)
    started = time.monotonic()
    manifest = {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "target": target_info["target"], "target_tag": tag, "sample": name,
        "condition": _condition(sample), "input_normalized_hits_file": str(hits),
        "output_dir": str(paths["output"]), "log_file": str(paths["log"]),
        "status": "failed", "start_time": _now(), "end_time": None,
        "elapsed_seconds": 0.0, "command_run": command, "feature_tables": [],
        "warnings": [], "error_message": None,
    }
    tables = _feature_tables(paths["output"])
    cached = False
    if not force and paths["manifest"].is_file() and tables:
        try:
            previous = json.loads(paths["manifest"].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = {}
        cached = previous.get("status") in {"success", "skipped"} and previous.get("input_normalized_hits_file") == str(hits)
    if cached:
        manifest["status"] = "skipped"
        manifest["warnings"].append("Successful cached normalized feature table found; SummaryTable skipped.")
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
            tables = _feature_tables(paths["output"])
            if not tables:
                raise RuntimeError("SummaryTable completed without an expected normalized feature table.")
            export = _paths(config, tag)["target_export"] / name
            export.mkdir(parents=True, exist_ok=True)
            for source in tables + ([paths["output"] / "stats.csv"] if (paths["output"] / "stats.csv").is_file() else []):
                shutil.copy2(source, export / source.name)
            manifest["status"] = "success"
        except Exception as exc:
            manifest["error_message"] = str(exc)
    manifest["end_time"] = _now()
    manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
    manifest["feature_tables"] = [str(path) for path in _feature_tables(paths["output"])]
    if not dry_run:
        _write_json(paths["manifest"], manifest)
    manifest["manifest_path"] = str(paths["manifest"])
    return manifest


def _write_status(path: Path, results: list[dict]) -> None:
    existing: dict[tuple[str, str], dict] = {}
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                existing[(row["target_tag"], row["sample"])] = row
    for result in results:
        existing[(result["target_tag"], result["sample"])] = {
            "target": result["target"], "target_tag": result["target_tag"],
            "sample": result["sample"], "condition": result["condition"],
            "status": result["status"],
            "input_normalized_hits_file": result["input_normalized_hits_file"],
            "output_dir": result["output_dir"], "feature_tables": ";".join(result["feature_tables"]),
            "log_file": result["log_file"], "elapsed_seconds": result["elapsed_seconds"],
            "message": result.get("error_message") or "; ".join(result.get("warnings") or []),
        }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerows(existing[key] for key in sorted(existing))


def _read_stats(path: Path) -> dict[str, str] | None:
    if not path.is_file():
        return None
    with path.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle), None)
    if not row:
        return None
    return {target: str(row.get(source, "")) for source, target in STATS_MAP.items()}


def _as_number(value: object) -> float | None:
    text = str(value or "").strip().rstrip("%")
    if not text or text.upper() == "NA":
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def _raw_stats(config: dict) -> dict[str, dict[str, str]]:
    path = _paths(config)["project"] / "summary" / "summary_stats.all_samples.csv"
    rows: dict[str, dict[str, str]] = {}
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            rows = {str(row["sample"]): row for row in csv.DictReader(handle)}
    return rows


def _write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def collect_normalization_target_evaluation(
    config: dict, min_feature_retention: float = 0.95
) -> dict:
    if not math.isfinite(min_feature_retention) or not 0 < min_feature_retention <= 1:
        raise ValueError("Minimum feature retention must be greater than 0 and at most 1.")
    paths = _paths(config)
    raw = _raw_stats(config)
    step6_comparison = paths["normalization"] / "normalization_comparison.csv"
    site_rows: dict[tuple[str, str], dict] = {}
    if step6_comparison.is_file():
        with step6_comparison.open(newline="", encoding="utf-8") as handle:
            site_rows = {(row["target_tag"], row["sample"]): row for row in csv.DictReader(handle)}
    evaluation: list[dict] = []
    manifest_paths = sorted(paths["manifest_base"].glob("T*/*.summary_normalized_manifest.json")) if paths["manifest_base"].is_dir() else []
    for manifest_path in manifest_paths:
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if manifest.get("status") not in {"success", "skipped"} or not manifest.get("feature_tables"):
            continue
        tag, sample = str(manifest["target_tag"]), str(manifest["sample"])
        site = site_rows.get((tag, sample), {})
        raw_row = raw.get(sample)
        norm_stats = _read_stats(Path(manifest["output_dir"]) / "stats.csv")
        raw_percent = _as_number(raw_row.get("percent_features_hit")) if raw_row else None
        norm_percent = _as_number(norm_stats.get("percent_features_hit")) if norm_stats else None
        feature_retention = norm_percent / raw_percent if raw_percent and norm_percent is not None else None
        warnings: list[str] = []
        if raw_row is None:
            warnings.append(RAW_WARNING)
        if feature_retention is None:
            warnings.append("Feature retention could not be calculated from SummaryTable statistics.")
        evaluation.append({
            "target": manifest["target"], "target_tag": tag, "sample": sample,
            "condition": manifest.get("condition", "unknown"),
            "raw_hit_sites": site.get("raw_hit_sites", "NA"),
            "normalized_hit_sites": site.get("normalized_hit_sites", "NA"),
            "hit_site_retention_fraction": site.get("hit_site_retention_fraction", "NA"),
            "raw_total_reads": site.get("raw_total_reads", raw_row.get("total_reads", "NA") if raw_row else "NA"),
            "normalized_total_reads": site.get("normalized_total_reads", norm_stats.get("total_reads", "NA") if norm_stats else "NA"),
            "raw_percent_features_hit": raw_row.get("percent_features_hit", "NA") if raw_row else "NA",
            "normalized_percent_features_hit": norm_stats.get("percent_features_hit", "NA") if norm_stats else "NA",
            "feature_retention_fraction": round(feature_retention, 8) if feature_retention is not None else "NA",
            "feature_loss_fraction": round(1 - feature_retention, 8) if feature_retention is not None else "NA",
            "raw_total_hits": raw_row.get("total_hits", "NA") if raw_row else "NA",
            "normalized_total_hits": norm_stats.get("total_hits", "NA") if norm_stats else "NA",
            "raw_mean_hits_per_feature": raw_row.get("mean_hits_per_feature", "NA") if raw_row else "NA",
            "normalized_mean_hits_per_feature": norm_stats.get("mean_hits_per_feature", "NA") if norm_stats else "NA",
            "normalized_hits_file": manifest["input_normalized_hits_file"],
            "normalized_feature_table": manifest["feature_tables"][0],
            "status": manifest["status"], "warning": "; ".join(warnings),
        })
    evaluation.sort(key=lambda row: (float(row["target"]), row["sample"]))
    _write_csv(paths["evaluation"], EVALUATION_FIELDS, evaluation)
    paths["evaluation_export"].parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(paths["evaluation"], paths["evaluation_export"])

    site_recommended = ""
    recommendation_path = paths["normalization"] / "normalization_recommendation.json"
    if recommendation_path.is_file():
        try:
            site_recommended = str(json.loads(recommendation_path.read_text(encoding="utf-8")).get("recommended_target_tag") or "")
        except (OSError, json.JSONDecodeError):
            pass
    grouped: dict[str, list[dict]] = {}
    for row in evaluation:
        grouped.setdefault(str(row["target_tag"]), []).append(row)
    safe: list[tuple[float, str]] = []
    summaries: list[dict] = []
    for tag, rows in grouped.items():
        sites = [_as_number(row["hit_site_retention_fraction"]) for row in rows]
        features = [_as_number(row["feature_retention_fraction"]) for row in rows]
        percents = [_as_number(row["normalized_percent_features_hit"]) for row in rows]
        sites_ok = [value for value in sites if value is not None]
        features_ok = [value for value in features if value is not None]
        percents_ok = [value for value in percents if value is not None]
        if len(features_ok) == len(rows) and min(features_ok) >= min_feature_retention:
            safe.append((float(rows[0]["target"]), tag))
        summaries.append({
            "target": rows[0]["target"], "target_tag": tag, "sample_count": len(rows),
            "min_hit_site_retention_fraction": min(sites_ok) if sites_ok else "NA",
            "mean_hit_site_retention_fraction": sum(sites_ok) / len(sites_ok) if sites_ok else "NA",
            "min_feature_retention_fraction": min(features_ok) if features_ok else "NA",
            "mean_feature_retention_fraction": sum(features_ok) / len(features_ok) if features_ok else "NA",
            "max_feature_loss_fraction": max(1 - value for value in features_ok) if features_ok else "NA",
            "mean_normalized_percent_features_hit": sum(percents_ok) / len(percents_ok) if percents_ok else "NA",
            "status": "success" if len(features_ok) == len(rows) else "incomplete",
            "recommended_by_site_retention": tag == site_recommended,
            "recommended_by_feature_retention": False,
        })
    warnings = [SMOKE_WARNING]
    if safe:
        recommended_target, recommended_tag = max(safe)
        reason = "Highest evaluated target satisfying the minimum feature-retention threshold across selected samples."
    elif summaries:
        least_aggressive = max(summaries, key=lambda row: float(row["target"]))
        recommended_target, recommended_tag = float(least_aggressive["target"]), str(least_aggressive["target_tag"])
        reason = "No evaluated target satisfied the feature-retention threshold; selected the least aggressive target tested."
        warnings.append(reason)
    else:
        recommended_target, recommended_tag = None, None
        reason = "No successful normalized SummaryTable results were available for feature-level recommendation."
        warnings.append(reason)
    for row in summaries:
        row["recommended_by_feature_retention"] = row["target_tag"] == recommended_tag
    summaries.sort(key=lambda row: float(row["target"]))
    _write_csv(paths["target_summary"], SUMMARY_FIELDS, summaries)
    shutil.copy2(paths["target_summary"], paths["target_summary_export"])
    recommendation = {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "recommended_target": recommended_target, "recommended_target_tag": recommended_tag,
        "min_feature_retention": min_feature_retention, "evaluated_target_tags": list(grouped),
        "reason": reason, "warnings": warnings, "timestamp": _now(),
    }
    _write_json(paths["feature_recommendation"], recommendation)
    fields = ["recommended_target", "recommended_target_tag", "min_feature_retention", "reason"]
    _write_csv(paths["feature_recommendation_export"], fields, [{key: recommendation[key] for key in fields}])
    return {
        "evaluation_path": str(paths["evaluation"]), "target_summary_path": str(paths["target_summary"]),
        "feature_recommendation_path": str(paths["feature_recommendation"]),
        "feature_recommendation": recommendation, "evaluation_rows": evaluation,
        "target_summary_rows": summaries,
    }


def run_summary_normalized_project(
    project_config: Path, targets: str | list[str] = "recommended",
    sample_mode: str = "parents", samples: list[str] | None = None,
    threads: int | None = None, force: bool = False, dry_run: bool = False,
    keep_going: bool = False, min_feature_retention: float = 0.95,
) -> dict:
    config = load_project_for_summary_normalized(project_config)
    selected = get_summary_normalized_samples(config, sample_mode, samples)
    target_infos = resolve_normalization_targets(config, targets)
    errors = validate_summary_normalized_inputs(config, target_infos, selected)
    if errors:
        raise ValueError("\n".join(errors))
    results: list[dict] = []
    stop = False
    for target_info in target_infos:
        for sample in selected:
            result = run_summary_normalized_sample(config, target_info, sample, threads, force, dry_run)
            results.append(result)
            if not dry_run:
                _write_status(_paths(config)["status"], results)
            if result["status"] == "failed" and not keep_going:
                stop = True
                break
        if stop:
            break
    evaluation = None if dry_run else collect_normalization_target_evaluation(config, min_feature_retention)
    return {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "sample_mode": "explicit" if samples is not None else sample_mode,
        "selected_samples": [str(row["sample"]) for row in selected],
        "target_infos": target_infos, "results": results,
        "output_dir": str(_paths(config)["summary_base"]), "evaluation": evaluation,
        "failed": sum(row["status"] == "failed" for row in results),
        "success": sum(row["status"] == "success" for row in results),
        "skipped": sum(row["status"] == "skipped" for row in results), "dry_run": dry_run,
    }
