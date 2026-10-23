"""Output contract of the treated-vs-control z-score: the tables and manifest fields every result has, whatever the
normalization, library or species. The app only reads these files; if a plot needs something, it is added here and
written by the pipeline (ytab.fitness.zscore), never computed in the app.

`validate(result_dir)` lists what is missing; the z-score CLI runs it after every analysis and fails on any problem.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

VERSION = 1

MANIFEST = "zscore_manifest.json"
MANIFEST_FIELDS = ("status", "contract", "tables", "method", "normalization", "normalization_run", "value_unit", "value_scale",
                   "n_features", "library_sizes", "settings", "call_rule", "comparisons", "backgrounds", "calls")

HIT_COLUMNS = ("feature_id", "standard_name", "common_name", "call", "log2FC", "z", "mean_abundance", "n_pairs_valid",
               "n_pairs_depleted_support", "n_pairs_enriched_support", "n_pairs_depleted_called", "n_pairs_enriched_called",
               "pools_depleted_support", "pools_enriched_support", "passes_abundance", "z_strength_depleted", "rank_depleted",
               "z_strength_enriched", "rank_enriched", "z_strength_both", "rank_both")
# table -> required columns. Values are in the manifest's value_unit (CPM, or the normalization's reads per million).
TABLES = {
    "tables/treated_vs_parent.summary_by_feature.csv":  # one row per gene: calls, effects and plot abundances
        ("feature_id", "standard_name", "common_name", "final_call", "mean_log2FC", "mean_z", "combined_log2FC", "combined_z",
         "n_pairs", "n_enriched_pairs", "n_depleted_pairs", "call_log2FC", "call_z", "mean_abundance", "control_mean", "treated_mean"),
    "tables/treated_vs_parent.by_pool.log2fc_z.csv":    # one row per gene and control-treated pair
        ("feature_id", "contrast", "background", "pool", "parent_sample", "treated_sample", "parent_cpm", "treated_cpm",
         "A_exp", "log2FC", "sd_local", "z", "Zthr", "call", "standard_name", "common_name"),
    "tables/treated_vs_parent.combined.log2fc_z.csv":   # one row per gene and background x treated condition
        ("feature_id", "background", "treated_condition", "parent_mean", "treated_mean", "A_exp", "log2FC", "z", "Zthr", "call"),
    "tables/gene_values.long.csv":                      # the gene values the z-score used, per library
        ("feature_id", "sample", "condition", "background", "pool", "value"),
    "tables/control_control_z_scores.csv":              # the noise model's control-control z (call threshold)
        ("background", "parent_comparison", "feature_id", "A_ctrl", "M_ctrl", "z"),
    "tables/background_specific_z_thresholds.csv": ("background", "z_quantile", "z_threshold", "parent_comparison"),
    "tables/background_noise_model_coefficients.csv": ("background", "parameter", "value"),
    # top hits shown in the app (labels, top-hit table, heatmap): evidence and the order per direction
    "tables/top_hit_evidence.csv": HIT_COLUMNS,
    "tables/top_hit_evidence.by_pair.csv": ("contrast",)+HIT_COLUMNS,
}
# Written whenever raw reads are known (SummaryTables or the normalization run's raw matrix): manifest library_sizes.
LIBRARY_SIZES = ("tables/library_sizes.feature_reads.csv",
                 ("sample", "condition", "background", "pool", "total_feature_reads", "features_with_reads", "value_total"))


def validate(result_dir: Path) -> list[str]:
    """Problems with a z-score result directory (empty when it meets the contract)."""
    out, problems = Path(result_dir), []
    path = out/MANIFEST
    if not path.is_file():return [f"missing {MANIFEST}"]
    manifest = json.loads(path.read_text())
    problems += [f"{MANIFEST}: missing field {k}" for k in MANIFEST_FIELDS if k not in manifest]
    if manifest.get("contract") != VERSION:problems.append(f"{MANIFEST}: contract {manifest.get('contract')} is not {VERSION}")
    tables = dict(TABLES) | (dict([LIBRARY_SIZES]) if manifest.get("library_sizes") else {})
    for name, columns in tables.items():
        if not (out/name).is_file():problems.append(f"missing {name}"); continue
        header = pd.read_csv(out/name, nrows=0).columns
        problems += [f"{name}: missing column {c}" for c in columns if c not in header]
    listed = set(manifest.get("tables") or [])   # the tables this run wrote; anything else in tables/ is not from it
    problems += [f"{name} is not listed in {MANIFEST}" for name in tables if listed and name not in listed]
    problems += [f"{name} is listed but missing" for name in listed if not (out/name).is_file()]
    summary = out/"tables/treated_vs_parent.summary_by_feature.csv"
    if summary.is_file() and "n_features" in manifest:
        n = sum(1 for _ in open(summary))-1
        if n != manifest["n_features"]:problems.append(f"summary has {n} genes, manifest says {manifest['n_features']}")
    return problems


# Resampling and control-conditioned Gumbel (scripts/local/ytab_treated_vs_control.py): one folder per run, one subfolder
# per background x treated condition group.
TVC_MANIFEST_FIELDS = ("status", "contract", "kind", "method", "method_label", "normalization", "normalization_run", "fdr", "groups")
TVC_GENES = {"resampling": ("group", "background", "condition", "ORF", "Name", "status", "in_family", "mean_control", "mean_treatment",
                            "mean_abundance", "log2FC", "p_value", "q_bh", "call"),
             "gumbel": ("group", "background", "condition", "ORF", "Name", "status", "in_family", "n_eligible", "n_depleted", "r",
                        "null_mean_run", "log2_scaled_ratio", "p_value", "q_bh", "call")}


def validate_tvc(result_dir: Path) -> list[str]:
    """Problems with a resampling or Gumbel result folder (empty when it meets the contract)."""
    out = Path(result_dir); path = out/"manifest.json"
    if not path.is_file():return ["missing manifest.json"]
    manifest = json.loads(path.read_text()); problems = [f"manifest.json: missing field {k}" for k in TVC_MANIFEST_FIELDS if k not in manifest]
    if manifest.get("contract") != VERSION:problems.append(f"manifest.json: contract {manifest.get('contract')} is not {VERSION}")
    ok = [k for k, g in (manifest.get("groups") or {}).items() if g.get("status") == "success"]
    if ok:
        if not (out/"genes.tsv").is_file():return problems+["missing genes.tsv"]
        header = pd.read_csv(out/"genes.tsv", sep="\t", nrows=0).columns
        problems += [f"genes.tsv: missing column {c}" for c in TVC_GENES.get(manifest.get("method"), ()) if c not in header]
    problems += [f"missing {k}/sites.tsv.gz" for k in ok if not (out/k/"sites.tsv.gz").is_file()]
    return problems
