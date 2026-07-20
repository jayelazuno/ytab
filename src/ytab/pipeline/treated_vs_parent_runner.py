"""Local raw-SummaryTable treated-versus-parent fitness orchestration.

The scientific calculation remains in ``scripts/ytab_treated_vs_parent_screen.R``.
This module deliberately refuses every normalized input mode.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

from .mapfastq_runner import load_project_for_mapping

DESIGN_COLUMNS = ["comparison_id", "pair_id", "parent_sample", "treated_sample",
                  "background", "pool", "parent_condition", "treated_condition",
                  "include", "notes"]
STATUS_COLUMNS = ["analysis_id", "status", "input_mode", "comparison_count",
                  "parent_samples", "treated_samples", "stable_result_table",
                  "output_dir", "log_file", "elapsed_seconds", "message"]
INPUT_MODE_ERROR = ("Treated-versus-parent analysis supports raw per-sample SummaryTable "
                    "inputs only. MidLC-normalized inputs are reserved for the essentiality classifier.")
CPM_ERROR = ("CPM normalization could not be confirmed in "
             "scripts/ytab_treated_vs_parent_screen.R. Step 11 requires the established "
             "raw-summary to CPM workflow. Review the R script before continuing.")


def _now():
    return datetime.now(timezone.utc).isoformat()


def _bool(value):
    return value if isinstance(value, bool) else str(value).strip().lower() in {"1", "true", "yes", "y"}


def _sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve(value, root: Path):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def _safe_id(value: str):
    result = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    if not result:
        raise ValueError("Analysis ID is empty after sanitization.")
    return result


def _paths(config: dict, analysis_id: str | None = None):
    root = Path(config["_repo_root"])
    project = _resolve(config["output_project_dir"], root)
    export = _resolve(config["output_export_dir"], root)
    paths = {"project": project, "design": project / "config" / "comparison_design.csv",
             "issues": project / "config" / "comparison_design_issues.csv",
             "status": project / "manifests" / "treated_vs_parent" / "treated_vs_parent_status.csv"}
    if analysis_id:
        aid = _safe_id(analysis_id)
        output = project / "treated_vs_parent" / aid
        paths.update(output=output, stable=output / "treated_vs_parent_results.csv",
                     summary=output / "treated_vs_parent_comparison_summary.csv",
                     metadata=output / "treated_vs_parent_run_metadata.json",
                     log=project / "logs" / "treated_vs_parent" / f"{aid}.log",
                     manifest=project / "manifests" / "treated_vs_parent" / f"{aid}.treated_vs_parent_manifest.json",
                     export=export / "treated_vs_parent" / aid)
    return paths


def load_project_for_treated_vs_parent(project_config: Path) -> dict:
    return load_project_for_mapping(project_config)


def _condition(row):
    explicit = str(row.get("condition") or "").lower()
    guessed = str(row.get("guessed_condition") or "").lower()
    treatment = str(row.get("treatment") or "").lower()
    name = str(row.get("sample") or "").lower()
    if explicit == "parent" or guessed == "parent" or treatment == "untreated" or "parent" in name:
        return "parent"
    if explicit == "treated" or guessed == "treated" or (treatment and treatment != "untreated") or "treated" in name:
        return "treated"
    return ""


def _metadata(row, name):
    value = row.get(name) or row.get(f"guessed_{name}") or ""
    if value:
        return str(value)
    sample = str(row.get("sample") or "")
    if name == "pool":
        hit = re.search(r"(?:^|[-_])pool[-_]?([A-Za-z0-9.]+)", sample, re.I)
        return hit.group(1) if hit else ""
    if name == "background":
        hit = re.search(r"(?:^|[-_])(y?H\d+)(?:[-_]|$)", sample, re.I)
        return hit.group(1) if hit else ""
    return ""


def write_comparison_design(path: Path, rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=DESIGN_COLUMNS)
        writer.writeheader(); writer.writerows({key: row.get(key, "") for key in DESIGN_COLUMNS} for row in rows)
    return path


def _write_design_issues(path: Path, rows: list[dict]) -> Path:
    fields = ["comparison_id", "pair_id", "issue", "message"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows({key: row.get(key, "") for key in fields} for row in rows)
    return path


def infer_comparison_design(config: dict, overwrite: bool = False) -> Path:
    path = _paths(config)["design"]
    if path.is_file() and not overwrite:
        return path
    samples = [row for row in config.get("samples", []) if _bool(row.get("include", True))]
    parents = [row for row in samples if _condition(row) == "parent"]
    treated = [row for row in samples if _condition(row) == "treated"]
    rows = []
    for child in treated:
        background, pool = _metadata(child, "background"), _metadata(child, "pool")
        matches = [parent for parent in parents if _metadata(parent, "background") == background and _metadata(parent, "pool") == pool]
        parent = matches[0] if len(matches) == 1 else {}
        treated_name = str(child.get("sample") or "")
        parent_name = str(parent.get("sample") or "")
        pair = f"{background}_pool{pool}" if background and pool else f"pair_{len(rows)+1}"
        rows.append({"comparison_id": f"{pair}_treated_vs_parent", "pair_id": pair,
                     "parent_sample": parent_name, "treated_sample": treated_name,
                     "background": background, "pool": pool, "parent_condition": "parent",
                     "treated_condition": "treated", "include": "true",
                     "notes": "" if len(matches) == 1 else f"Expected one matching parent; found {len(matches)}."})
    write_comparison_design(path, rows)
    validation = validate_comparison_design(config, rows)
    _write_design_issues(_paths(config)["issues"], validation["issues"])
    return path


def load_or_create_comparison_design(config: dict, overwrite: bool = False) -> Path:
    return infer_comparison_design(config, overwrite)


def load_comparison_design(design_file: Path) -> list[dict]:
    with Path(design_file).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def validate_comparison_design(config: dict, design_rows: list[dict]) -> dict:
    samples = {str(row.get("sample") or ""): row for row in config.get("samples", [])}
    issues = []
    seen_comparisons, seen_pairs, parent_use, treated_use = set(), set(), {}, {}
    def issue(row, code, message):
        issues.append({"comparison_id": row.get("comparison_id", ""), "pair_id": row.get("pair_id", ""),
                       "issue": code, "message": message})
    for row in design_rows:
        if not _bool(row.get("include", True)): continue
        cid, pid = row.get("comparison_id", ""), row.get("pair_id", "")
        parent, treated = row.get("parent_sample", ""), row.get("treated_sample", "")
        if not cid or cid in seen_comparisons: issue(row, "duplicate_comparison_id", "Comparison ID is empty or duplicated.")
        if not pid or pid in seen_pairs: issue(row, "duplicate_pair_id", "Pair ID is empty or duplicated.")
        seen_comparisons.add(cid); seen_pairs.add(pid)
        for role, sample in (("parent", parent), ("treated", treated)):
            if sample not in samples: issue(row, f"unknown_{role}_sample", f"{sample or '(empty)'} is not in sample_sheet.csv.")
            elif not _bool(samples[sample].get("include", True)): issue(row, f"excluded_{role}_sample", f"{sample} is excluded in sample_sheet.csv.")
        if parent and parent == treated: issue(row, "same_sample", "Parent and treated samples are identical.")
        if not row.get("background"): issue(row, "missing_background", "Background is missing.")
        if not row.get("pool"): issue(row, "missing_pool", "Pool is missing.")
        parent_use[parent] = parent_use.get(parent, 0) + 1; treated_use[treated] = treated_use.get(treated, 0) + 1
    for sample, count in parent_use.items():
        if sample and count > 1: issue({}, "parent_reused", f"Parent sample {sample} occurs in {count} included pairs.")
    for sample, count in treated_use.items():
        if sample and count > 1: issue({}, "treated_reused", f"Treated sample {sample} occurs in {count} included pairs.")
    return {"valid": not issues, "issues": issues, "included_count": sum(_bool(r.get("include", True)) for r in design_rows)}


def get_included_comparisons(design_rows: list[dict], comparisons: list[str] | None = None) -> list[dict]:
    rows = [dict(row) for row in design_rows if _bool(row.get("include", True))]
    if comparisons:
        requested = set(comparisons)
        rows = [row for row in rows if row.get("comparison_id") in requested or row.get("pair_id") in requested]
        found = {row.get("comparison_id") for row in rows} | {row.get("pair_id") for row in rows}
        missing = requested - found
        if missing: raise ValueError("Unknown comparison or pair IDs: " + ", ".join(sorted(missing)))
    if not rows: raise ValueError("No included comparisons were selected.")
    return rows


def _table_columns(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        first = handle.readline(); second = handle.readline()
    header = second if first.strip().startswith("RDF") else first
    delimiter = "\t" if "\t" in header and "," not in header else ","
    return [value.strip() for value in next(csv.reader([header], delimiter=delimiter))]


def detect_raw_summary_feature_table(config: dict, sample: str) -> Path:
    directory = _paths(config)["project"] / "summary" / sample
    candidates = sorted({p for pattern in ("*.feature_table*.csv", "*.feature_table*.tsv", "*.feature_table*.txt") for p in directory.glob(pattern) if p.is_file() and p.stat().st_size})
    preferred = [p for p in candidates if "feature_table.RDF_1" in p.name]
    selected = preferred if preferred else candidates
    if not selected:
        raise FileNotFoundError(f"Raw SummaryTable output is missing for sample {sample}. Run Step 4 before treated-versus-parent analysis.")
    if len(selected) != 1:
        raise ValueError(f"Ambiguous raw SummaryTable feature tables for sample {sample}: " + ", ".join(map(str, selected)))
    columns = {re.sub(r"[^a-z0-9]+", "_", x.lower()).strip("_") for x in _table_columns(selected[0])}
    if not columns.intersection({"standard_name", "feature_name", "feature", "gene", "gene_name", "locus_tag", "orf", "qng_id"}):
        raise ValueError(f"Raw SummaryTable for {sample} has no recognized feature identifier column.")
    if not columns.intersection({"reads", "read_count", "feature_reads", "reads_per_feature", "total_reads", "sum_reads"}):
        raise ValueError(f"Raw SummaryTable for {sample} has no recognized read-count column.")
    return selected[0]


def inspect_treated_vs_parent_script(script_path: Path) -> dict:
    text = Path(script_path).read_text(encoding="utf-8")
    evidence = []
    for pattern, label in ((r"reads\s*/\s*total_feature_reads\s*\*\s*1e6", "reads / total_feature_reads * 1e6"),
                           (r"counts?\s*/\s*[^\n]+\*\s*(?:1e6|1000000)", "count/library-size scaling to one million")):
        if re.search(pattern, text, re.I): evidence.append(label)
    def value(pattern):
        hit = re.search(pattern, text, re.I); return hit.group(1) if hit else None
    return {"script": str(Path(script_path).resolve()), "cpm_normalization_detected": bool(evidence),
            "cpm_detection_evidence": evidence, "normalization_method": "counts per million (CPM)" if evidence else None,
            "pseudocount": value(r"PSEUDOCOUNT[^\n]*[\"']([0-9.]+)"),
            "significance_thresholds": {"z_quantile": value(r"Z_QUANTILE[^\n]*[\"']([0-9.]+)")},
            "fold_change_thresholds": {}, "multiple_testing_method": None,
            "analysis_method": "background-specific parent-parent local noise model and z-score calls"}


def verify_cpm_normalization(script_inspection: dict) -> dict:
    if not script_inspection.get("cpm_normalization_detected"): raise ValueError(CPM_ERROR)
    return script_inspection


def resolve_treated_vs_parent_inputs(config: dict, comparisons: list[dict], input_mode: str = "auto") -> dict:
    if input_mode not in {"auto", "raw-summary"}: raise ValueError(INPUT_MODE_ERROR)
    all_rows = load_comparison_design(_paths(config)["design"])
    samples = sorted({r[k] for r in all_rows if _bool(r.get("include", True)) for k in ("parent_sample", "treated_sample") if r.get(k)})
    files = {sample: detect_raw_summary_feature_table(config, sample) for sample in samples}
    return {"requested_input_mode": input_mode, "input_mode": "raw-summary", "sample_files": files,
            "selected_comparisons": comparisons}


def resolve_optional_classifier_annotation(config: dict, target: str | None = None) -> dict | None:
    from .classifier_runner import find_combined_feature_table, resolve_classifier_target
    info = resolve_classifier_target(config, target or "recommended")
    table = _paths(config)["project"] / "classifier" / info["target_tag"] / f"essentiality_predictions.{info['target_tag']}.csv"
    if not table.is_file(): raise FileNotFoundError(f"Classifier prediction table is missing: {table}")
    return {**info, "prediction_table": table}


def validate_treated_vs_parent_inputs(config, comparisons, resolved_inputs, script_inspection):
    errors = []
    if not script_inspection.get("cpm_normalization_detected"): errors.append(CPM_ERROR)
    for sample, path in resolved_inputs.get("sample_files", {}).items():
        if not Path(path).is_file(): errors.append(f"Raw SummaryTable input is missing for {sample}: {path}")
    return errors


def build_treated_vs_parent_command(config, analysis_id, comparisons, resolved_inputs, output_dir,
                                    classifier_annotation=None, extra_args=None):
    root = Path(config["_repo_root"]); script = root / "scripts" / "ytab_treated_vs_parent_screen.R"
    rscript = shutil.which("Rscript")
    if not rscript: raise FileNotFoundError("Rscript was not found in the active environment or PATH.")
    ids = ",".join(row["comparison_id"] for row in comparisons)
    command = [rscript, str(script), "--project-id", str(config["project_id"]),
               "--comparison-design", str(_paths(config)["design"]), "--comparisons", ids,
               "--output-dir", str(output_dir), "--analysis-id", _safe_id(analysis_id)]
    for sample, path in sorted(resolved_inputs["sample_files"].items()):
        command.extend(["--sample-table", f"{sample}={path}"])
    command.extend(extra_args or [])
    return command


def detect_treated_vs_parent_outputs(output_dir: Path) -> list[Path]:
    stable = {"treated_vs_parent_results.csv", "treated_vs_parent_comparison_summary.csv", "treated_vs_parent_run_metadata.json"}
    return sorted(p for p in Path(output_dir).rglob("*") if p.is_file() and p.name not in stable)


def choose_primary_treated_vs_parent_table(outputs: list[Path]) -> Path:
    candidates = [p for p in outputs if p.name == "treated_vs_parent.summary_by_feature.csv"]
    if len(candidates) != 1:
        raise ValueError("Could not identify one treated_vs_parent.summary_by_feature.csv output; candidates: " + ", ".join(map(str, candidates)))
    return candidates[0]


def _read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle: return list(csv.DictReader(handle))


def create_stable_treated_vs_parent_outputs(config, analysis_id, comparisons, outputs):
    paths = _paths(config, analysis_id); primary = choose_primary_treated_vs_parent_table(outputs)
    shutil.copy2(primary, paths["stable"])
    pair_table = next((p for p in outputs if p.name == "treated_vs_parent.by_pool.log2fc_z.csv"), None)
    summary_rows = []
    if pair_table:
        rows = _read_csv(pair_table)
        for comp in comparisons:
            matches = [r for r in rows if r.get("contrast") == comp["comparison_id"]]
            counts = {}
            for row in matches: counts[row.get("call", "")] = counts.get(row.get("call", ""), 0) + 1
            summary_rows.append({**{k: comp.get(k, "") for k in DESIGN_COLUMNS[:6]}, "result_rows": len(matches),
                                 "call_counts": json.dumps(counts, sort_keys=True)})
    fields = DESIGN_COLUMNS[:6] + ["result_rows", "call_counts"]
    with paths["summary"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(summary_rows)
    return {"stable_result_table": paths["stable"], "comparison_summary_file": paths["summary"],
            "result_row_count": len(_read_csv(paths["stable"])), "scientific_columns": _table_columns(paths["stable"])}


def annotate_treated_vs_parent_results(result_table: Path, classifier_annotation: dict) -> dict:
    results, annotations = _read_csv(result_table), _read_csv(classifier_annotation["prediction_table"])
    result_cols = list(results[0]) if results else []; annotation_cols = list(annotations[0]) if annotations else []
    keys = [key for key in ("feature_id", "standard_name", "gene", "gene_id") if key in result_cols and key in annotation_cols]
    if len(keys) != 1: raise ValueError("Could not identify one shared result/classifier feature identifier.")
    key = keys[0]; lookup = {row[key]: row for row in annotations if row.get(key)}; matched = 0
    added = [f"classifier_{name}" for name in annotation_cols if name != key]
    for row in results:
        hit = lookup.get(row.get(key)); matched += int(hit is not None)
        for source, dest in zip([x for x in annotation_cols if x != key], added): row[dest] = hit.get(source, "") if hit else ""
    with result_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=result_cols + added); writer.writeheader(); writer.writerows(results)
    return {"classifier_annotation_match_count": matched, "classifier_annotation_unmatched_count": len(results)-matched, "join_key": key}


def write_treated_vs_parent_metadata(config, analysis_id, comparisons, resolved_inputs, script_inspection,
                                     outputs, classifier_annotation=None):
    paths = _paths(config, analysis_id); design = _paths(config)["design"]
    r_environment = _r_environment(outputs)
    metadata = {"project_id": config.get("project_id"), "species": config.get("species"), "analysis_id": _safe_id(analysis_id),
                "input_mode": "raw-summary", "comparison_design": str(design), "comparison_design_sha256": _sha256(design),
                "included_comparisons": comparisons, "parent_samples": sorted({r["parent_sample"] for r in comparisons}),
                "treated_samples": sorted({r["treated_sample"] for r in comparisons}),
                "raw_summary_input_files": {k: str(v) for k,v in resolved_inputs["sample_files"].items()},
                "input_file_sha256": {k: _sha256(v) for k,v in resolved_inputs["sample_files"].items()},
                "R_script": script_inspection["script"], "R_script_sha256": _sha256(Path(script_inspection["script"])),
                "R_version": r_environment.get("R_version"),
                "relevant_R_package_versions": r_environment.get("packages", {}),
                "cpm_normalization_detected": True, "cpm_detection_evidence": script_inspection["cpm_detection_evidence"],
                "normalization_method": script_inspection["normalization_method"], "analysis_method": script_inspection["analysis_method"],
                "pseudocount": script_inspection["pseudocount"], "significance_thresholds": script_inspection["significance_thresholds"],
                "fold_change_thresholds": script_inspection["fold_change_thresholds"], "multiple_testing_method": None,
                "classifier_annotation_used": bool(classifier_annotation),
                "classifier_prediction_table": str(classifier_annotation["prediction_table"]) if classifier_annotation else None,
                "warnings": [], "timestamp": _now()}
    paths["metadata"].write_text(json.dumps(metadata, indent=2, default=str)+"\n", encoding="utf-8")
    return paths["metadata"]


def _r_environment(outputs):
    path = next((p for p in outputs if p.name == "analysis_runtime_environment.csv"), None)
    if not path: return {}
    rows = _read_csv(path)
    if not rows: return {}
    row = rows[0]
    return {"R_version": row.pop("R_version", None), "packages": row}


def _write_status(path, result):
    existing = {}
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle: existing = {r["analysis_id"]: r for r in csv.DictReader(handle)}
    existing[result["analysis_id"]] = {key: result.get(key, "") for key in STATUS_COLUMNS}
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer=csv.DictWriter(handle, fieldnames=STATUS_COLUMNS); writer.writeheader(); writer.writerows(existing[k] for k in sorted(existing))


def run_treated_vs_parent_project(project_config: Path, analysis_id="treated_vs_parent_raw_cpm",
                                  comparison_design=None, comparisons=None, input_mode="auto",
                                  classifier_target=None, annotate_classifier=False, force=False,
                                  dry_run=False, keep_going=False, extra_args=None):
    started = time.monotonic(); start_time = _now(); config = load_project_for_treated_vs_parent(project_config)
    aid = _safe_id(analysis_id); paths = _paths(config, aid)
    design = Path(comparison_design).resolve() if comparison_design else infer_comparison_design(config)
    if design != _paths(config)["design"]: shutil.copy2(design, _paths(config)["design"]); design = _paths(config)["design"]
    rows = load_comparison_design(design); validation = validate_comparison_design(config, rows)
    _write_design_issues(_paths(config)["issues"], validation["issues"])
    if not validation["valid"]: raise ValueError("Comparison design has unresolved issues; review " + str(_paths(config)["issues"]))
    selected = get_included_comparisons(rows, comparisons)
    resolved = resolve_treated_vs_parent_inputs(config, selected, input_mode)
    script = Path(config["_repo_root"]) / "scripts" / "ytab_treated_vs_parent_screen.R"
    inspection = verify_cpm_normalization(inspect_treated_vs_parent_script(script))
    errors = validate_treated_vs_parent_inputs(config, selected, resolved, inspection)
    if errors: raise ValueError("\n".join(errors))
    annotation = resolve_optional_classifier_annotation(config, classifier_target) if annotate_classifier else None
    command = build_treated_vs_parent_command(config, aid, selected, resolved, paths["output"], annotation, extra_args)
    hashes = {str(p): _sha256(p) for p in resolved["sample_files"].values()}
    cache_key = {"comparison_design_hash": _sha256(design), "input_hashes": hashes, "script_hash": _sha256(script),
                 "analysis_parameters": {"input_mode":"raw-summary", "command":command[1:]},
                 "classifier_annotation": str(annotation["prediction_table"]) if annotation else None}
    if not force and paths["manifest"].is_file() and paths["stable"].is_file():
        old=json.loads(paths["manifest"].read_text(encoding="utf-8"))
        if old.get("status") == "success" and all(old.get(k)==v for k,v in cache_key.items()):
            return {**old, "status":"skipped", "command_run":command, "resolved_inputs":resolved,
                    "script_inspection":inspection, "manifest_file":str(paths["manifest"])}
    paths["output"].mkdir(parents=True, exist_ok=True); paths["log"].parent.mkdir(parents=True, exist_ok=True); paths["manifest"].parent.mkdir(parents=True, exist_ok=True)
    manifest={"project_id":config.get("project_id"), "species":config.get("species"), "analysis_id":aid,
              "input_mode":"raw-summary", "comparison_design_file":str(design), "selected_comparisons":selected,
              "parent_samples":sorted({r["parent_sample"] for r in selected}), "treated_samples":sorted({r["treated_sample"] for r in selected}),
              "input_files":{k:str(v) for k,v in resolved["sample_files"].items()}, "output_dir":str(paths["output"]),
              "stable_result_table":str(paths["stable"]), "comparison_summary_file":str(paths["summary"]),
              "run_metadata_file":str(paths["metadata"]), "log_file":str(paths["log"]), "status":"failed",
              "start_time":start_time, "command_run":command, "warnings":[], **cache_key}
    try:
        if dry_run:
            manifest["status"]="dry-run"; manifest["detected_outputs"]=[]
        else:
            with paths["log"].open("w", encoding="utf-8") as log:
                completed=subprocess.run(command, cwd=config["_repo_root"], stdout=log, stderr=subprocess.STDOUT, text=True)
            if completed.returncode: raise RuntimeError(f"R analysis exited with status {completed.returncode}; see {paths['log']}")
            outputs=detect_treated_vs_parent_outputs(paths["output"])
            stable=create_stable_treated_vs_parent_outputs(config, aid, selected, outputs)
            annotation_counts=annotate_treated_vs_parent_results(paths["stable"], annotation) if annotation else {}
            write_treated_vs_parent_metadata(config, aid, selected, resolved, inspection, outputs, annotation)
            paths["export"].mkdir(parents=True, exist_ok=True)
            for item in (paths["stable"], paths["summary"], paths["metadata"]): shutil.copy2(item, paths["export"] / item.name)
            manifest.update(status="success", detected_outputs=[str(p) for p in outputs], **stable, **annotation_counts)
    except Exception as exc:
        manifest["error_message"]=str(exc)
    manifest["end_time"]=_now(); manifest["elapsed_seconds"]=round(time.monotonic()-started,3)
    paths["manifest"].write_text(json.dumps(manifest,indent=2,default=str)+"\n",encoding="utf-8")
    status={"analysis_id":aid,"status":manifest["status"],"input_mode":"raw-summary","comparison_count":len(selected),
            "parent_samples":";".join(manifest["parent_samples"]),"treated_samples":";".join(manifest["treated_samples"]),
            "stable_result_table":str(paths["stable"]),"output_dir":str(paths["output"]),"log_file":str(paths["log"]),
            "elapsed_seconds":manifest["elapsed_seconds"],"message":manifest.get("error_message","")}
    _write_status(paths["status"],status)
    manifest.update(resolved_inputs=resolved,script_inspection=inspection,comparison_design_file=str(design),manifest_file=str(paths["manifest"]))
    return manifest
