"""Local MidLC normalization sweeps and insertion-site retention summaries."""

from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path

from .library_diagnostics_runner import load_project_for_library_diagnostics
from .mapfastq_runner import get_included_samples


STATUS_FIELDS = [
    "target", "target_tag", "sample", "condition", "status", "input_hits_file",
    "normalized_hits_file", "output_dir", "log_file", "elapsed_seconds", "message",
]
COMPARISON_FIELDS = [
    "target", "target_tag", "sample", "condition", "raw_hit_sites",
    "normalized_hit_sites", "hit_site_retention_fraction", "hit_site_loss_fraction",
    "raw_total_reads", "normalized_total_reads", "raw_percent_features_hit",
    "normalized_percent_features_hit", "normalized_hits_file", "status", "warning",
]
FEATURE_WARNING = (
    "Normalized feature-level metrics require SummaryTable on normalized hits in a later step."
)
RETENTION_WARNING = (
    "Feature-level retention requires SummaryTable on normalized hits and will be confirmed in a later step."
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _resolve(value: str | Path | None, root: Path) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def _paths(config: dict, target: float | None = None) -> dict[str, Path]:
    root = Path(config["_repo_root"])
    project = _resolve(config.get("output_project_dir"), root)
    export_root = _resolve(config.get("output_export_dir"), root)
    if project is None or export_root is None:
        raise ValueError("Project and export output directories must be configured.")
    base = project / "sample_normalization"
    export = export_root / "qc" / "sample_normalization"
    manifests = project / "manifests" / "sample_normalization"
    paths = {
        "project": project, "base": base, "export": export, "manifests": manifests,
        "status": manifests / "sample_normalization_status.csv",
        "comparison": base / "normalization_comparison.csv",
        "export_comparison": export / "normalization_comparison.csv",
        "recommendation": base / "normalization_recommendation.json",
        "export_recommendation": export / "normalization_recommendation.csv",
    }
    if target is not None:
        tag = normalize_target_tag(target)
        paths.update({
            "output": base / tag,
            "log": project / "logs" / "sample_normalization" / f"{tag}.log",
            "manifest": manifests / f"{tag}.sample_normalization_manifest.json",
            "target_export": export / tag,
        })
    return paths


def load_project_for_normalization(project_config: Path) -> dict:
    return load_project_for_library_diagnostics(project_config)


def _condition(sample: dict) -> str:
    role_values = {
        str(sample.get("classifier_role") or "").strip().lower(),
        str(sample.get("library_role") or "").strip().lower(),
        str(sample.get("control_or_treated") or "").strip().lower(),
        str(sample.get("fitness_role") or "").strip().lower(),
    }
    if role_values.intersection({"classifier_control", "parent", "control"}):
        return "parent"
    if role_values.intersection({"exclude", "treated", "treatment"}):
        return "treated"
    value = str(sample.get("condition") or sample.get("guessed_condition") or "").strip().lower()
    if value in {"parent", "mock", "control", "untreated"}:
        return "parent"
    if value in {"treated", "treatment"} or "treated" in value:
        return "treated"
    name = str(sample.get("sample") or "").lower()
    if "treated" in name:
        return "treated"
    if "mock" in name or "parent" in name:
        return "parent"
    return "unknown"


def get_normalization_samples(
    config: dict, sample_mode: str = "parents", samples: list[str] | None = None
) -> list[dict]:
    included = get_included_samples(config)
    by_name = {str(row.get("sample") or ""): row for row in included}
    if samples is not None:
        if not samples or len(samples) != len(set(samples)):
            raise ValueError("Explicit sample names must be a non-empty unique list.")
        unknown = [name for name in samples if name not in by_name]
        if unknown:
            raise ValueError("Selected samples are not included or do not exist: " + ", ".join(unknown))
        return [by_name[name] for name in samples]
    if sample_mode not in {"parents", "all", "treated"}:
        raise ValueError(f"Unsupported sample mode: {sample_mode}")
    if sample_mode == "all":
        selected = included
    elif sample_mode == "treated":
        selected = [row for row in included if _condition(row) == "treated"]
    else:
        selected = [row for row in included if _condition(row) == "parent"]
    if not selected:
        raise ValueError(f"No included samples matched sample mode: {sample_mode}")
    return [dict(row) for row in selected]


def find_sample_hits_file(sample: dict, config: dict) -> Path:
    name = str(sample["sample"])
    return _paths(config)["project"] / "create_hit_file" / name / f"{name}_hits.txt"


def validate_normalization_inputs(config: dict, selected_samples: list[dict]) -> list[str]:
    errors: list[str] = []
    if not selected_samples:
        return ["No samples were selected for normalization."]
    for sample in selected_samples:
        name = str(sample.get("sample") or "")
        if not name:
            errors.append("An included sample has no sample name.")
        elif Path(name).name != name:
            errors.append(f"Sample name contains path separators: {name}")
        hits = find_sample_hits_file(sample, config)
        if not hits.is_file():
            errors.append(f"Hit file missing for sample {name}: {hits}")
    return errors


def _target_decimal(value: int | float | str | Decimal) -> Decimal:
    try:
        target = Decimal(str(value).strip())
    except (InvalidOperation, ValueError):
        raise ValueError(f"Invalid normalization target: {value!r}") from None
    if not target.is_finite() or target <= 0:
        raise ValueError(f"Normalization target must be a finite positive number: {value!r}")
    return target.normalize()


def parse_targets(targets: str | list[str] | None) -> str | list[float]:
    if targets is None:
        return "auto"
    if isinstance(targets, str):
        text = targets.strip()
        if text.lower() == "auto":
            return "auto"
        values = [part.strip() for part in text.split(",")]
    else:
        values = [str(part).strip() for part in targets]
    if not values or any(not value for value in values):
        raise ValueError("Normalization target list cannot be empty.")
    parsed = [float(_target_decimal(value)) for value in values]
    if len(parsed) != len(set(parsed)):
        raise ValueError("Normalization targets must be unique.")
    return parsed


def normalize_target_tag(target: int | float | str) -> str:
    value = _target_decimal(target)
    text = format(value, "f").rstrip("0").rstrip(".") if "." in format(value, "f") else format(value, "f")
    whole, _, fraction = text.partition(".")
    return "T" + whole.zfill(3) + ("p" + fraction if fraction else "")


def count_hit_sites_and_reads(hits_file: Path) -> dict:
    sites = reads = 0
    with Path(hits_file).open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Chromosome", "Hit position", "Hit count"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Hit file lacks required columns {sorted(required)}: {hits_file}")
        for row in reader:
            try:
                count = int(float(row["Hit count"]))
            except (TypeError, ValueError):
                continue
            if count > 0:
                sites += 1
                reads += count
    return {"hit_sites": sites, "total_reads": reads}


def build_normalization_command(
    config: dict, target: float, selected_samples: list[dict], threads: int | None = None
) -> list[str]:
    if threads is not None and int(threads) < 1:
        raise ValueError("Threads must be at least 1.")
    tag_paths = _paths(config, target)
    root = Path(config["_repo_root"])
    return [
        sys.executable, str(root / "src" / "ytab" / "qc" / "LibraryDiagnostics.py"),
        "--hits-txt", *[str(find_sample_hits_file(row, config)) for row in selected_samples],
        "--outdir", str(tag_paths["output"]),
        "--normalize-to-midlc", format(_target_decimal(target), "f"),
        "--normalized-hits-suffix", "_normalized_hits.txt",
    ]


def _normalized_file(output: Path, sample: dict) -> Path:
    name = str(sample["sample"])
    return output / name / f"{name}_normalized_hits.txt"


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def _write_status(path: Path, results: list[dict]) -> None:
    existing: dict[tuple[str, str], dict] = {}
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                existing[(row["target_tag"], row["sample"])] = row
    for result in results:
        for row in result["sample_results"]:
            existing[(row["target_tag"], row["sample"])] = {key: row.get(key, "") for key in STATUS_FIELDS}
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerows(existing[key] for key in sorted(existing))


def _comparison_row(config: dict, target: float, sample: dict, status: str) -> dict:
    paths = _paths(config, target)
    raw_file = find_sample_hits_file(sample, config)
    normalized = _normalized_file(paths["output"], sample)
    raw = count_hit_sites_and_reads(raw_file)
    norm = count_hit_sites_and_reads(normalized) if normalized.is_file() else {"hit_sites": 0, "total_reads": 0}
    retention = norm["hit_sites"] / raw["hit_sites"] if raw["hit_sites"] else 0.0
    return {
        "target": format(_target_decimal(target), "f"), "target_tag": normalize_target_tag(target),
        "sample": str(sample["sample"]), "condition": _condition(sample),
        "raw_hit_sites": raw["hit_sites"], "normalized_hit_sites": norm["hit_sites"],
        "hit_site_retention_fraction": round(retention, 8),
        "hit_site_loss_fraction": round(1.0 - retention, 8),
        "raw_total_reads": raw["total_reads"], "normalized_total_reads": norm["total_reads"],
        "raw_percent_features_hit": "NA", "normalized_percent_features_hit": "NA",
        "normalized_hits_file": str(normalized), "status": status, "warning": FEATURE_WARNING,
    }


def _copy_target_exports(paths: dict[str, Path]) -> None:
    for source in paths["output"].rglob("*"):
        if source.is_file():
            target = paths["target_export"] / source.relative_to(paths["output"])
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)


def run_normalization_target(
    config: dict, target: float, selected_samples: list[dict], threads: int | None = None,
    force: bool = False, dry_run: bool = False,
) -> dict:
    target = float(_target_decimal(target))
    paths = _paths(config, target)
    tag = normalize_target_tag(target)
    command = build_normalization_command(config, target, selected_samples, threads)
    started = time.monotonic()
    expected = [_normalized_file(paths["output"], row) for row in selected_samples]
    selected_names = [str(row["sample"]) for row in selected_samples]
    warnings = [RETENTION_WARNING]
    manifest = {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "target": format(_target_decimal(target), "f"), "target_tag": tag,
        "selected_samples": selected_names,
        "input_hit_files": [str(find_sample_hits_file(row, config)) for row in selected_samples],
        "output_dir": str(paths["output"]), "export_dir": str(paths["target_export"]),
        "log_file": str(paths["log"]), "status": "failed", "start_time": _now(),
        "end_time": None, "elapsed_seconds": 0.0, "command_run": command,
        "normalized_hit_files": [], "comparison_rows": [], "warnings": warnings,
        "error_message": None,
    }
    cached = False
    if not force and paths["manifest"].is_file() and all(path.is_file() for path in expected):
        try:
            previous = json.loads(paths["manifest"].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = {}
        cached = previous.get("status") in {"success", "skipped"} and previous.get("selected_samples") == selected_names
    if cached:
        manifest["status"] = "skipped"
        warnings.append("Successful cached normalized hit files found; target skipped.")
    elif dry_run:
        manifest["status"] = "skipped"
        warnings.append("Dry run; command was not executed.")
    else:
        paths["output"].mkdir(parents=True, exist_ok=True)
        paths["log"].parent.mkdir(parents=True, exist_ok=True)
        try:
            with paths["log"].open("w", encoding="utf-8") as log:
                log.write("Command: " + " ".join(command) + "\n\n")
                log.flush()
                completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
            if completed.returncode != 0:
                raise RuntimeError(f"LibraryDiagnostics exited with status {completed.returncode}; see {paths['log']}")
            missing = [str(path) for path in expected if not path.is_file()]
            if missing:
                raise RuntimeError("Expected normalized hit files were not created: " + ", ".join(missing))
            _copy_target_exports(paths)
            manifest["status"] = "success"
        except Exception as exc:
            manifest["error_message"] = str(exc)
    manifest["end_time"] = _now()
    manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
    manifest["normalized_hit_files"] = [str(path) for path in expected if path.is_file()]
    if manifest["status"] in {"success", "skipped"} and all(path.is_file() for path in expected):
        manifest["comparison_rows"] = [
            _comparison_row(config, target, row, manifest["status"]) for row in selected_samples
        ]
    manifest["sample_results"] = []
    for sample, normalized in zip(selected_samples, expected):
        ok = normalized.is_file()
        status = manifest["status"] if ok or dry_run else "failed"
        manifest["sample_results"].append({
            "target": manifest["target"], "target_tag": tag, "sample": str(sample["sample"]),
            "condition": _condition(sample), "status": status,
            "input_hits_file": str(find_sample_hits_file(sample, config)),
            "normalized_hits_file": str(normalized), "output_dir": str(normalized.parent),
            "log_file": str(paths["log"]), "elapsed_seconds": manifest["elapsed_seconds"],
            "message": manifest.get("error_message") or "; ".join(warnings),
        })
    if not dry_run:
        _write_json(paths["manifest"], {key: value for key, value in manifest.items() if key != "sample_results"})
    manifest["manifest_path"] = str(paths["manifest"])
    return manifest


def collect_normalization_comparison(config: dict) -> Path | None:
    paths = _paths(config)
    rows: list[dict] = []
    if paths["manifests"].is_dir():
        for manifest_path in sorted(paths["manifests"].glob("T*.sample_normalization_manifest.json")):
            try:
                data = json.loads(manifest_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            rows.extend(data.get("comparison_rows") or [])
    if not rows:
        return None
    paths["comparison"].parent.mkdir(parents=True, exist_ok=True)
    with paths["comparison"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COMPARISON_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    paths["export"].mkdir(parents=True, exist_ok=True)
    shutil.copy2(paths["comparison"], paths["export_comparison"])
    return paths["comparison"]


def _observed_midlc(config: dict, selected_samples: list[dict]) -> tuple[dict[str, float], list[str]]:
    summary = _paths(config)["project"] / "library_diagnostics" / "library_diagnostics.summary.csv"
    warnings: list[str] = []
    values: dict[str, float] = {}
    if summary.is_file():
        with summary.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row.get("sample") in {str(sample["sample"]) for sample in selected_samples}:
                    try:
                        value = float(row["midlc_est"])
                    except (KeyError, TypeError, ValueError):
                        continue
                    if math.isfinite(value) and value > 0:
                        values[str(row["sample"])] = value
    if len(values) != len(selected_samples):
        warnings.append("MidLC estimates were not available; using fallback candidate targets.")
        return {}, warnings
    return values, warnings


def _candidate_targets(
    observed: dict[str, float], minimum: float | None, maximum: float | None, step: float
) -> list[float]:
    step = float(_target_decimal(step))
    if not observed and minimum is None and maximum is None:
        return [20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0]
    min_target = float(_target_decimal(minimum)) if minimum is not None else max(1.0, math.floor(min(observed.values()) * 0.25))
    max_target = float(_target_decimal(maximum)) if maximum is not None else float(math.floor(min(observed.values())))
    if max_target < min_target:
        raise ValueError("Auto maximum target must be greater than or equal to auto minimum target.")
    values: list[float] = []
    current = Decimal(str(min_target))
    end = Decimal(str(max_target))
    increment = Decimal(str(step))
    while current <= end:
        values.append(float(current))
        current += increment
    if values[-1] != max_target:
        values.append(max_target)
    return values


def _simulate_retention(raw: dict, midlc: float | None, target: float) -> float:
    if not midlc or raw["total_reads"] <= 0 or raw["hit_sites"] <= 0:
        return 1.0
    p = min(1.0, (target * midlc) / raw["total_reads"])
    # Conservative singletons bound: every retained site survives with at least p.
    return p


def recommend_normalization_target(
    config: dict, selected_samples: list[dict], min_site_retention: float = 0.95,
    auto_min_target: float | None = None, auto_max_target: float | None = None,
    auto_step: float = 5.0, threads: int | None = None, force: bool = False,
    dry_run: bool = False,
) -> dict:
    if not math.isfinite(min_site_retention) or not 0 < min_site_retention <= 1:
        raise ValueError("Minimum site retention must be greater than 0 and at most 1.")
    observed, warnings = _observed_midlc(config, selected_samples)
    candidates = _candidate_targets(observed, auto_min_target, auto_max_target, auto_step)
    raw_counts = {
        str(sample["sample"]): count_hit_sites_and_reads(find_sample_hits_file(sample, config))
        for sample in selected_samples
    }
    scored: list[tuple[float, float, float]] = []
    for target in candidates:
        retained = [
            _simulate_retention(
                raw_counts[str(sample["sample"])],
                observed.get(str(sample["sample"])), target,
            )
            for sample in selected_samples
        ]
        scored.append((target, min(retained), sum(retained) / len(retained)))
    feasible = [row for row in scored if row[1] >= min_site_retention]
    if feasible:
        recommended, min_observed, mean_observed = feasible[-1]
        reason = "Highest candidate target satisfying the insertion-site retention threshold across selected samples."
    else:
        recommended, min_observed, mean_observed = scored[0]
        reason = "No candidate satisfied the insertion-site retention threshold; selected the least aggressive target tested."
        warnings.append(reason)
    target_result = run_normalization_target(
        config, recommended, selected_samples, threads, force=force, dry_run=dry_run
    )
    if target_result["comparison_rows"]:
        actual = [float(row["hit_site_retention_fraction"]) for row in target_result["comparison_rows"]]
        min_observed, mean_observed = min(actual), sum(actual) / len(actual)
    recommendation = {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "selected_samples": [str(row["sample"]) for row in selected_samples],
        "sample_mode": config.get("_sample_mode", "parents"),
        "observed_midlc_values": observed, "candidate_targets": candidates,
        "candidate_target_tags": [normalize_target_tag(value) for value in candidates],
        "recommended_target": recommended, "recommended_target_tag": normalize_target_tag(recommended),
        "retention_threshold": min_site_retention,
        "min_site_retention_observed": round(min_observed, 8),
        "mean_site_retention_observed": round(mean_observed, 8),
        "reason": reason, "warnings": warnings + [RETENTION_WARNING], "timestamp": _now(),
    }
    paths = _paths(config)
    if not dry_run:
        _write_json(paths["recommendation"], recommendation)
        paths["export"].mkdir(parents=True, exist_ok=True)
        with paths["export_recommendation"].open("w", newline="", encoding="utf-8") as handle:
            fields = ["recommended_target", "recommended_target_tag", "min_site_retention_observed", "mean_site_retention_observed", "retention_threshold", "reason"]
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerow({key: recommendation[key] for key in fields})
    return {**recommendation, "target_result": target_result, "recommendation_path": str(paths["recommendation"]), "recommendation_csv": str(paths["export_recommendation"])}


def run_normalization_project(
    project_config: Path, targets: str | list[str] | None = "auto",
    sample_mode: str = "parents", samples: list[str] | None = None,
    threads: int | None = None, force: bool = False, dry_run: bool = False,
    keep_going: bool = False, min_site_retention: float = 0.95,
    auto_min_target: float | None = None, auto_max_target: float | None = None,
    auto_step: float = 5.0,
) -> dict:
    config = load_project_for_normalization(project_config)
    selected = get_normalization_samples(config, sample_mode, samples)
    errors = validate_normalization_inputs(config, selected)
    if errors:
        raise ValueError("\n".join(errors))
    config["_sample_mode"] = "explicit" if samples is not None else sample_mode
    parsed = parse_targets(targets)
    results: list[dict] = []
    recommendation = None
    if parsed == "auto":
        recommendation = recommend_normalization_target(
            config, selected, min_site_retention, auto_min_target, auto_max_target,
            auto_step, threads, force, dry_run,
        )
        results.append(recommendation["target_result"])
    else:
        for target in parsed:
            result = run_normalization_target(config, target, selected, threads, force, dry_run)
            results.append(result)
            if result["status"] == "failed" and not keep_going:
                break
    if not dry_run:
        _write_status(_paths(config)["status"], results)
    comparison = None if dry_run else collect_normalization_comparison(config)
    return {
        "project_id": config.get("project_id"), "species": config.get("species"),
        "sample_mode": config["_sample_mode"], "selected_samples": [str(row["sample"]) for row in selected],
        "input_hit_files": [str(find_sample_hits_file(row, config)) for row in selected],
        "target_mode": "auto" if parsed == "auto" else "manual",
        "targets": [result["target"] for result in results],
        "target_tags": [result["target_tag"] for result in results], "results": results,
        "recommendation": recommendation, "output_dir": str(_paths(config)["base"]),
        "comparison_path": str(comparison) if comparison else None,
        "failed": sum(result["status"] == "failed" for result in results),
        "dry_run": dry_run,
    }
