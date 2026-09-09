"""Python orchestration for YTAB GO enrichment.

The statistical work is delegated to an R worker because topGO and rrvgo are
Bioconductor packages.  This module keeps YTAB-side concerns in Python:
project paths, gene-list normalization, C. glabrata-to-S. cerevisiae lookup,
progress records, and stable output manifests.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml

from ytab.pipeline.progress_tracker import ProgressTracker, make_job_id


LOOKUP_RELATIVE = Path(
    "resources/comparative/orthology/20260610_Cgla_CAGL_to_Scer_annotation_lookup.csv"
)
REQUIRED_LOOKUP_COLUMNS = {
    "cagl_id",
    "gwk60_id_clean",
    "qng_id",
    "cgla_common_name",
    "cg_to_sc_relationship",
    "pre_WGD_Ancestor",
    "scer_gene_id",
    "scer_gene_name",
    "SGD_essentiality",
    "SGD_description",
}
ONTOLOGY_LABELS = {
    "BP": "Biological Process",
    "MF": "Molecular Function",
    "CC": "Cellular Component",
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_id(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value or "").strip()).strip("._-")
    return value or datetime.now(timezone.utc).strftime("go_%Y%m%dT%H%M%SZ")


def _norm(value: object) -> str:
    value = "" if value is None else str(value)
    return value.strip().lower()


def _split_orfs(value: object) -> list[str]:
    out: list[str] = []
    for part in re.split(r"[;,]\s*|\s+", str(value or "").strip()):
        part = part.strip()
        if re.match(r"^Y[A-Z]{2}\d{3}[A-Z](?:-[A-Z])?$", part, re.I):
            out.append(part.upper())
    return list(dict.fromkeys(out))


def _repo_root_from_config(config: dict, project_config: Path) -> Path:
    """Resolve the repository root from project config without hard-coding a project.

    Existing project YAML files may store either ``repo_root`` or ``_repo_root``.
    Relative values are interpreted from the config file directory first, because
    project configs usually live under output/projects/<project>/config/.
    """
    project_config = project_config.expanduser().resolve()
    raw = config.get("_repo_root") or config.get("repo_root") or ""
    candidates: list[Path] = []
    if raw:
        raw_path = Path(str(raw)).expanduser()
        candidates.extend(
            [raw_path.resolve()] if raw_path.is_absolute() else [
                (project_config.parent / raw_path).resolve(),
                (Path.cwd() / raw_path).resolve(),
            ]
        )
    candidates.append(Path.cwd().resolve())
    candidates.extend(parent for parent in project_config.parents)
    unique_candidates = list(dict.fromkeys(candidates))
    for candidate in unique_candidates:
        if (
            (candidate / "app" / "shiny" / "app.R").is_file()
            and (candidate / LOOKUP_RELATIVE).is_file()
        ):
            return candidate
    for candidate in unique_candidates:
        if (candidate / "app" / "shiny" / "app.R").is_file():
            return candidate
    return project_config.parents[3]


def _resolve_path(value: object, repo_root: Path) -> Path:
    path = Path(str(value)).expanduser()
    return path.resolve() if path.is_absolute() else (repo_root / path).resolve()


def load_project(project_config: Path) -> tuple[dict, Path, Path]:
    project_config = project_config.expanduser().resolve()
    with project_config.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    repo_root = _repo_root_from_config(config, project_config)
    project_root = _resolve_path(config.get("output_project_dir", project_config.parents[1]), repo_root)
    return config, repo_root, project_root


def load_lookup(repo_root: Path) -> pd.DataFrame:
    path = repo_root / LOOKUP_RELATIVE
    if not path.is_file():
        raise FileNotFoundError(f"GO analysis requires the Cgla/Scer lookup file: {path}")
    lookup = pd.read_csv(path, dtype=str).fillna("")
    missing = sorted(REQUIRED_LOOKUP_COLUMNS.difference(lookup.columns))
    if missing:
        raise ValueError(f"Lookup file is missing required columns: {', '.join(missing)}")
    for col in ("gwk60_id_clean", "cagl_id", "qng_id", "scer_gene_id", "scer_gene_name", "cgla_common_name"):
        lookup[f"{col}_key"] = lookup[col].map(_norm)
    return lookup


def _lookup_index(lookup: pd.DataFrame) -> dict[str, int]:
    keys: dict[str, int] = {}
    for i, row in lookup.iterrows():
        for col in (
            "gwk60_id_clean_key",
            "cagl_id_key",
            "qng_id_key",
            "scer_gene_id_key",
            "scer_gene_name_key",
            "cgla_common_name_key",
        ):
            key = str(row.get(col, ""))
            if key and key not in keys:
                keys[key] = int(i)
    return keys


def read_gene_file(path: Path) -> list[str]:
    if not path or not path.is_file():
        return []
    text = path.read_text(encoding="utf-8", errors="replace")
    if "," in text or "\t" in text:
        try:
            sep = "\t" if "\t" in text.splitlines()[0] else ","
            df = pd.read_csv(path, sep=sep, dtype=str).fillna("")
            candidate_cols = [
                "input_gene",
                "feature_id",
                "Original feature ID",
                "CAGL ID",
                "Gene name",
                "gene",
                "Gene",
            ]
            for col in candidate_cols:
                if col in df.columns:
                    return [x for x in df[col].astype(str).tolist() if x.strip()]
            if len(df.columns):
                return [x for x in df.iloc[:, 0].astype(str).tolist() if x.strip()]
        except Exception:
            pass
    return [x for x in re.split(r"[\s,;]+", text) if x.strip()]


def parse_gene_inputs(genes: str = "", gene_file: Path | None = None) -> list[str]:
    out: list[str] = []
    if genes:
        out.extend(x for x in re.split(r"[\s,;]+", genes) if x.strip())
    if gene_file:
        out.extend(read_gene_file(gene_file))
    return list(dict.fromkeys(x.strip() for x in out if x.strip()))


def map_genes_to_scer(gene_ids: Iterable[str], lookup: pd.DataFrame) -> pd.DataFrame:
    index = _lookup_index(lookup)
    rows: list[dict[str, str]] = []
    for gene in gene_ids:
        key = _norm(gene)
        row = None
        if re.match(r"^Y[A-Z]{2}\d{3}[A-Z](?:-[A-Z])?$", str(gene), re.I):
            row = {
                "input_gene": gene,
                "matched": "true",
                "gwk60_id_clean": "",
                "cagl_id": "",
                "gene_name": gene.upper(),
                "cg_to_sc_relationship": "",
                "scer_gene_id": gene.upper(),
            }
        elif key in index:
            source = lookup.iloc[index[key]]
            orfs = _split_orfs(source.get("scer_gene_id", ""))
            row = {
                "input_gene": gene,
                "matched": "true" if orfs else "false",
                "gwk60_id_clean": source.get("gwk60_id_clean", ""),
                "cagl_id": source.get("cagl_id", ""),
                "gene_name": source.get("scer_gene_name", "") or source.get("cagl_id", "") or gene,
                "cg_to_sc_relationship": source.get("cg_to_sc_relationship", ""),
                "scer_gene_id": ";".join(orfs),
            }
        else:
            row = {
                "input_gene": gene,
                "matched": "false",
                "gwk60_id_clean": "",
                "cagl_id": "",
                "gene_name": gene,
                "cg_to_sc_relationship": "",
                "scer_gene_id": "",
            }
        rows.append(row)
    return pd.DataFrame(rows)


def _read_feature_column(path: Path) -> list[str]:
    try:
        if path.suffix.lower() in {".csv"}:
            df = pd.read_csv(path, dtype=str, nrows=200000).fillna("")
        else:
            df = pd.read_csv(path, dtype=str, sep="\t", nrows=200000).fillna("")
    except Exception:
        return []
    for col in ("feature_id", "gene_id", "Original feature ID", "Gene", "gene", "standard_name"):
        if col in df.columns:
            return [x for x in df[col].astype(str).tolist() if x.strip()]
    return []


def discover_project_feature_ids(project_root: Path) -> list[str]:
    candidates: list[Path] = []
    for pattern in (
        "treated_vs_parent/*/treated_vs_parent_results.csv",
        "classifier/*/essentiality_predictions*.csv",
        "summary*/**/feature_table*.csv",
        "summary*/**/feature_table*.tsv",
        "summary*/**/combined_feature_table*.txt",
    ):
        candidates.extend(project_root.glob(pattern))
    out: list[str] = []
    for path in candidates:
        if path.is_file():
            out.extend(_read_feature_column(path))
        if len(out) > 1000:
            break
    return list(dict.fromkeys(out))


def build_universe(project_root: Path, lookup: pd.DataFrame) -> pd.DataFrame:
    feature_ids = discover_project_feature_ids(project_root)
    if feature_ids:
        universe = map_genes_to_scer(feature_ids, lookup)
    else:
        universe = pd.DataFrame(
            {
                "input_gene": lookup["gwk60_id_clean"],
                "matched": "true",
                "gwk60_id_clean": lookup["gwk60_id_clean"],
                "cagl_id": lookup["cagl_id"],
                "gene_name": lookup["scer_gene_name"].where(lookup["scer_gene_name"].astype(bool), lookup["cagl_id"]),
                "cg_to_sc_relationship": lookup["cg_to_sc_relationship"],
                "scer_gene_id": lookup["scer_gene_id"].map(lambda x: ";".join(_split_orfs(x))),
            }
        )
    universe = universe[universe["scer_gene_id"].astype(str).str.len() > 0].copy()
    universe = universe.drop_duplicates(subset=["scer_gene_id"])
    return universe


def write_csv(path: Path, rows: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(path, index=False)


def run_go_analysis(
    project_config: Path,
    analysis_id: str,
    genes: list[str],
    ontologies: list[str],
    node_size: int,
    p_cutoff: float,
    top_terms: int,
    semantic_reduction: bool,
    rrvgo_threshold: float,
    source_label: str,
    force: bool = False,
    job_id: str | None = None,
    progress_file: Path | None = None,
) -> dict:
    config, repo_root, project_root = load_project(project_config)
    lookup = load_lookup(repo_root)
    analysis_id = _safe_id(analysis_id)
    out_dir = project_root / "go_analysis" / analysis_id
    manifest = out_dir / "manifest.json"
    if manifest.is_file() and not force:
        with manifest.open("r", encoding="utf-8") as handle:
            prior = json.load(handle)
        if prior.get("status") == "success":
            prior["status"] = "cached"
            return prior

    if not genes:
        raise ValueError("No target genes were provided for GO analysis.")
    ontologies = [x.upper() for x in ontologies if x.upper() in ONTOLOGY_LABELS]
    if not ontologies:
        raise ValueError("Select at least one GO ontology: BP, MF, or CC.")

    tracker = None
    if progress_file:
        tracker = ProgressTracker.create(
            progress_file,
            job_id or make_job_id("go_analysis"),
            str(config.get("project_id", project_root.name)),
            "go_analysis",
            ["prepare_inputs", "topgo", "semantic_reduction" if semantic_reduction else "summarize"],
        )
        tracker.start("Preparing GO analysis inputs")
        tracker.start_item("prepare_inputs", 1, "prepare_inputs")

    target = map_genes_to_scer(genes, lookup)
    target = target[target["scer_gene_id"].astype(str).str.len() > 0].copy()
    if target.empty:
        raise ValueError("None of the provided target genes could be mapped to S. cerevisiae ORF IDs.")
    universe = build_universe(project_root, lookup)
    universe_orfs = set()
    for value in universe["scer_gene_id"].astype(str):
        universe_orfs.update(_split_orfs(value))
    target["in_project_universe"] = target["scer_gene_id"].map(lambda v: any(orf in universe_orfs for orf in _split_orfs(v)))
    inputs_dir = out_dir / "inputs"
    tables_dir = out_dir / "tables"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    target_file = inputs_dir / "target_genes.csv"
    universe_file = inputs_dir / "universe_genes.csv"
    write_csv(target_file, target)
    write_csv(universe_file, universe)
    if tracker:
        tracker.finish_item("prepare_inputs", "success", [str(target_file), str(universe_file)], "GO inputs prepared")
        tracker.start_item("topgo", 2, "topgo")

    r_worker = repo_root / "scripts/local/ytab_topgo_rrvgo_worker.R"
    cmd = [
        "Rscript",
        str(r_worker),
        "--targets",
        str(target_file),
        "--universe",
        str(universe_file),
        "--out-dir",
        str(out_dir),
        "--ontologies",
        ",".join(ontologies),
        "--node-size",
        str(node_size),
        "--p-cutoff",
        str(p_cutoff),
        "--top-terms",
        str(top_terms),
    ]
    if semantic_reduction:
        cmd.extend(["--semantic-reduction", "--rrvgo-threshold", str(rrvgo_threshold)])

    start = _now()
    result = subprocess.run(cmd, cwd=repo_root, text=True, capture_output=True, check=False)
    if tracker:
        tracker.finish_item(
            "topgo",
            "success" if result.returncode == 0 else "failed",
            [str(tables_dir / "topgo_results.csv")],
            "topGO completed" if result.returncode == 0 else "topGO failed",
            result.stderr[-2000:] if result.returncode else "",
        )
    if result.returncode != 0:
        if tracker:
            tracker.finalize("failed", "GO analysis failed")
        raise RuntimeError(result.stderr or result.stdout or "R topGO worker failed.")
    if tracker:
        final_item = "semantic_reduction" if semantic_reduction else "summarize"
        tracker.start_item(final_item, 3, final_item)
        tracker.finish_item(final_item, "success", [str(tables_dir)], "GO result tables written")

    job_manifest_dir = project_root / "manifests" / "go_analysis"
    job_manifest_dir.mkdir(parents=True, exist_ok=True)
    job_manifest = job_manifest_dir / f"{analysis_id}.go_analysis_manifest.json"

    payload = {
        "schema_version": 1,
        "analysis_id": analysis_id,
        "project_id": str(config.get("project_id", project_root.name)),
        "source_label": source_label,
        "status": "success",
        "started_at": start,
        "completed_at": _now(),
        "lookup_file": str(repo_root / LOOKUP_RELATIVE),
        "target_gene_count": int(len(genes)),
        "mapped_target_rows": int(len(target)),
        "universe_gene_rows": int(len(universe)),
        "ontologies": ontologies,
        "node_size": int(node_size),
        "p_cutoff": float(p_cutoff),
        "top_terms": int(top_terms),
        "semantic_reduction": bool(semantic_reduction),
        "rrvgo_threshold": float(rrvgo_threshold),
        "output_dir": str(out_dir),
        "manifest": str(manifest),
        "job_manifest": str(job_manifest),
        "target_file": str(target_file),
        "universe_file": str(universe_file),
        "topgo_results": str(tables_dir / "topgo_results.csv"),
        "significant_results": str(tables_dir / "topgo_significant.csv"),
        "semantic_heatmap": str(tables_dir / "semantic_reduced_heatmap.csv"),
        "command_run": cmd,
        "stdout": result.stdout[-5000:],
        "stderr": result.stderr[-5000:],
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    with job_manifest.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    if tracker:
        tracker.finalize("success", "GO analysis complete")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run YTAB GO enrichment through topGO and optional rrvgo.")
    parser.add_argument("--project-config", type=Path, required=True)
    parser.add_argument("--analysis-id", default="")
    parser.add_argument("--genes", default="")
    parser.add_argument("--gene-file", type=Path)
    parser.add_argument("--ontologies", default="BP,MF,CC")
    parser.add_argument("--node-size", type=int, default=10)
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--top-terms", type=int, default=15)
    parser.add_argument("--semantic-reduction", action="store_true")
    parser.add_argument("--rrvgo-threshold", type=float, default=0.70)
    parser.add_argument("--source-label", default="Selected genes")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--job-id")
    parser.add_argument("--progress-file", type=Path)
    args = parser.parse_args(argv)
    genes = parse_gene_inputs(args.genes, args.gene_file)
    ontologies = [x.strip().upper() for x in re.split(r"[,; ]+", args.ontologies) if x.strip()]
    try:
        result = run_go_analysis(
            args.project_config,
            args.analysis_id,
            genes,
            ontologies,
            args.node_size,
            args.p_cutoff,
            args.top_terms,
            args.semantic_reduction,
            args.rrvgo_threshold,
            args.source_label,
            args.force,
            args.job_id,
            args.progress_file,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    writer = csv.writer(sys.stdout)
    writer.writerow(["Status", result.get("status", "")])
    writer.writerow(["Analysis ID", result.get("analysis_id", "")])
    writer.writerow(["Output directory", result.get("output_dir", "")])
    writer.writerow(["Target genes", result.get("target_gene_count", "")])
    writer.writerow(["Mapped target rows", result.get("mapped_target_rows", "")])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
