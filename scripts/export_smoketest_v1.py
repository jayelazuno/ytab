
#!/usr/bin/env python3
"""
export_smoketest_v1.py

Build an app-facing export bundle from YTAB smoke-test outputs.

Design:
- Reads only from output/smoketests/
- Writes only to output/exports/smoke_test_v1/
- Preserves raw per-sample files under qc/raw/
- Writes app-friendly aggregated tables under qc/tables/
- Copies browser tracks and classifier artifacts into stable locations
- Writes manifest.json for the Streamlit app

Current default scope:
- parent smoke test sample(s), starting with yH298-parent-pool1

Update SAMPLES later as more parent pools are added.
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Iterable, List

import pandas as pd

# Config

YTAB_ROOT = Path("/Users/jayelazuno/workspace/ytab")

SMOKE_ROOT = YTAB_ROOT / "output" / "smoketests"
EXPORT_ROOT = YTAB_ROOT / "output" / "exports" / "smoke_test_v1"

SAMPLES = [
    "yH298-parent-pool1",
]

WARNINGS: List[str] = []


# Helpers

def warn(msg: str) -> None:
    WARNINGS.append(msg)
    print(f"[WARN] {msg}")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def safe_copy(src: Path, dst: Path) -> bool:
    if not src.exists():
        warn(f"Missing file: {src}")
        return False
    ensure_dir(dst.parent)
    shutil.copy2(src, dst)
    return True


def copy_tree_contents(src_dir: Path, dst_dir: Path) -> None:
    if not src_dir.exists():
        warn(f"Missing directory: {src_dir}")
        return
    ensure_dir(dst_dir)
    for item in src_dir.iterdir():
        target = dst_dir / item.name
        if item.is_dir():
            shutil.copytree(item, target, dirs_exist_ok=True)
        else:
            shutil.copy2(item, target)


def find_first(root: Path, pattern: str) -> Path | None:
    matches = sorted(root.rglob(pattern))
    if not matches:
        return None
    return matches[0]


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported table type: {path}")


def add_sample_column(df: pd.DataFrame, sample: str) -> pd.DataFrame:
    if "sample" not in df.columns:
        df.insert(0, "sample", sample)
    return df


def concat_tables(paths_with_samples: Iterable[tuple[Path, str]], out_path: Path) -> None:
    frames = []
    for path, sample in paths_with_samples:
        if not path.exists():
            warn(f"Missing table: {path}")
            continue
        try:
            df = read_table(path)
            df = add_sample_column(df, sample)
            frames.append(df)
        except Exception as exc:
            warn(f"Failed reading {path}: {exc}")

    if not frames:
        warn(f"No data written for {out_path.name}")
        return

    combined = pd.concat(frames, ignore_index=True)
    ensure_dir(out_path.parent)

    if out_path.suffix.lower() == ".tsv":
        combined.to_csv(out_path, sep="\t", index=False)
    else:
        combined.to_csv(out_path, index=False)


def copy_text_files(paths_with_samples: Iterable[tuple[Path, str]], out_dir: Path) -> None:
    ensure_dir(out_dir)
    for path, sample in paths_with_samples:
        if not path.exists():
            warn(f"Missing text file: {path}")
            continue
        dst = out_dir / f"{sample}.{path.name}"
        shutil.copy2(path, dst)


def relative_to_export(path: Path) -> str:
    return str(path.relative_to(EXPORT_ROOT))


# Export steps

def export_browser() -> dict:
    browser_dir = EXPORT_ROOT / "browser"
    ensure_dir(browser_dir)

    exported = []

    for sample in SAMPLES:
        src_dir = SMOKE_ROOT / "create_hit_file" / sample

        bedgraph = src_dir / f"{sample}.insertions.bedgraph"
        wig = src_dir / f"{sample}.insertions.wig"

        if safe_copy(bedgraph, browser_dir / bedgraph.name):
            exported.append(browser_dir / bedgraph.name)
        if safe_copy(wig, browser_dir / wig.name):
            exported.append(browser_dir / wig.name)

    return {
        "browser_dir": relative_to_export(browser_dir),
        "tracks": [relative_to_export(p) for p in exported],
    }


def export_classifier() -> dict:
    classifier_root = SMOKE_ROOT / "classifier"

    export_tables = EXPORT_ROOT / "classifier" / "tables"
    export_images = EXPORT_ROOT / "classifier" / "images"
    ensure_dir(export_tables)
    ensure_dir(export_images)

    out = {}

    combined_tsv = find_first(classifier_root, "combined_glabrata_RF_G4.tsv")
    auc_png = find_first(classifier_root, "trial_AUC*.png")
    tables_xlsx = find_first(classifier_root, "tables.xlsx")

    if combined_tsv is not None:
        dst = export_tables / "combined_glabrata_RF_G4.tsv"
        safe_copy(combined_tsv, dst)
        out["classifier_table"] = relative_to_export(dst)

    if tables_xlsx is not None:
        dst = export_tables / "tables.xlsx"
        safe_copy(tables_xlsx, dst)
        out["classifier_workbook"] = relative_to_export(dst)

    if auc_png is not None:
        dst = export_images / "trial_auc_rf_g4.png"
        safe_copy(auc_png, dst)
        out["classifier_plot"] = relative_to_export(dst)

    return out


def export_qc_raw() -> None:
    raw_root = EXPORT_ROOT / "qc" / "raw"

    for sample in SAMPLES:
        # mapfastq raw
        map_src = SMOKE_ROOT / "mapfastq" / sample / f"{sample}.mapping_stats.csv"
        map_dst = raw_root / "mapfastq" / f"{sample}.mapping_stats.csv"
        safe_copy(map_src, map_dst)

        # diagnostics raw
        diag_src_dir = SMOKE_ROOT / "library_diagnostics" / sample
        diag_dst_dir = raw_root / "diagnostics" / sample
        copy_tree_contents(diag_src_dir, diag_dst_dir)

        # summary raw
        sum_src_dir = SMOKE_ROOT / "summary" / sample
        sum_dst_dir = raw_root / "summary" / sample
        copy_tree_contents(sum_src_dir, sum_dst_dir)


def export_qc_tables_and_images() -> dict:
    qc_tables = EXPORT_ROOT / "qc" / "tables"
    qc_images = EXPORT_ROOT / "qc" / "images"
    ensure_dir(qc_tables)
    ensure_dir(qc_images)

    manifest_paths = {}

    # -----------------------------
    # Mapping stats
    # -----------------------------
    mapping_inputs = [
        (SMOKE_ROOT / "mapfastq" / s / f"{s}.mapping_stats.csv", s)
        for s in SAMPLES
    ]
    mapping_out = qc_tables / "parent_mapping_stats.csv"
    concat_tables(mapping_inputs, mapping_out)
    manifest_paths["mapping_stats"] = relative_to_export(mapping_out)

    # -----------------------------
    # Library diagnostics
    # -----------------------------
    diag_root_summary = SMOKE_ROOT / "library_diagnostics" / "library_diagnostics.summary.csv"
    if diag_root_summary.exists():
        dst = qc_tables / "library_diagnostics_summary.csv"
        shutil.copy2(diag_root_summary, dst)
        manifest_paths["library_diagnostics_summary"] = relative_to_export(dst)
    else:
        warn(f"Missing diagnostics summary: {diag_root_summary}")

    diag_sample_summary_inputs = [
        (SMOKE_ROOT / "library_diagnostics" / s / f"{s}.summary.csv", s)
        for s in SAMPLES
    ]
    diag_sample_summary_out = qc_tables / "library_diagnostics_per_sample.csv"
    concat_tables(diag_sample_summary_inputs, diag_sample_summary_out)
    manifest_paths["library_diagnostics_per_sample"] = relative_to_export(diag_sample_summary_out)

    midlc_inputs = [
        (SMOKE_ROOT / "library_diagnostics" / s / f"{s}.midlc.csv", s)
        for s in SAMPLES
    ]
    midlc_out = qc_tables / "midlc_long.csv"
    concat_tables(midlc_inputs, midlc_out)
    manifest_paths["midlc_long"] = relative_to_export(midlc_out)

    seqbias_inputs = [
        (SMOKE_ROOT / "library_diagnostics" / s / f"{s}.seqbias_2_7.tsv", s)
        for s in SAMPLES
    ]
    seqbias_out = qc_tables / "seqbias_long.tsv"
    concat_tables(seqbias_inputs, seqbias_out)
    manifest_paths["seqbias_long"] = relative_to_export(seqbias_out)

    centromere_inputs = [
        (SMOKE_ROOT / "library_diagnostics" / s / f"{s}.centromere_bins.tsv", s)
        for s in SAMPLES
    ]
    centromere_out = qc_tables / "centromere_bins_long.tsv"
    concat_tables(centromere_inputs, centromere_out)
    manifest_paths["centromere_bins_long"] = relative_to_export(centromere_out)

    tss_inputs = [
        (SMOKE_ROOT / "library_diagnostics" / s / f"{s}.tss_metaplot.tsv", s)
        for s in SAMPLES
    ]
    tss_out = qc_tables / "tss_metaplot_long.tsv"
    concat_tables(tss_inputs, tss_out)
    manifest_paths["tss_metaplot_long"] = relative_to_export(tss_out)

    tts_inputs = [
        (SMOKE_ROOT / "library_diagnostics" / s / f"{s}.tts_metaplot.tsv", s)
        for s in SAMPLES
    ]
    tts_out = qc_tables / "tts_metaplot_long.tsv"
    concat_tables(tts_inputs, tts_out)
    manifest_paths["tts_metaplot_long"] = relative_to_export(tts_out)

    trna_inputs = [
        (SMOKE_ROOT / "library_diagnostics" / s / f"{s}.trna_metaplot.tsv", s)
        for s in SAMPLES
    ]
    trna_out = qc_tables / "trna_metaplot_long.tsv"
    concat_tables(trna_inputs, trna_out)
    manifest_paths["trna_metaplot_long"] = relative_to_export(trna_out)

    # Copy QC images
    for sample in SAMPLES:
        cent_png = SMOKE_ROOT / "library_diagnostics" / sample / f"{sample}.centromere_bias.png"
        meta_png = SMOKE_ROOT / "library_diagnostics" / sample / f"{sample}.metaplots.png"

        cent_dst = qc_images / f"{sample}.centromere_bias.png"
        meta_dst = qc_images / f"{sample}.metaplots.png"

        if safe_copy(cent_png, cent_dst):
            manifest_paths.setdefault("qc_images", []).append(relative_to_export(cent_dst))
        if safe_copy(meta_png, meta_dst):
            manifest_paths.setdefault("qc_images", []).append(relative_to_export(meta_dst))

    # Summary outputs
    summary_stats_inputs = [
        (SMOKE_ROOT / "summary" / s / "stats.csv", s)
        for s in SAMPLES
    ]
    summary_stats_out = qc_tables / "summary_stats.csv"
    concat_tables(summary_stats_inputs, summary_stats_out)
    manifest_paths["summary_stats"] = relative_to_export(summary_stats_out)

    hit_summary_inputs = [
        (SMOKE_ROOT / "summary" / s / "hit_summary.RDF_1.csv", s)
        for s in SAMPLES
    ]
    hit_summary_out = qc_tables / "hit_summary_long.csv"
    concat_tables(hit_summary_inputs, hit_summary_out)
    manifest_paths["hit_summary_long"] = relative_to_export(hit_summary_out)

    binned_hits_inputs = [
        (SMOKE_ROOT / "summary" / s / "binned_hits.RDF_1.csv", s)
        for s in SAMPLES
    ]
    binned_hits_out = qc_tables / "binned_hits_long.csv"
    concat_tables(binned_hits_inputs, binned_hits_out)
    manifest_paths["binned_hits_long"] = relative_to_export(binned_hits_out)

    feature_table_inputs = [
        (SMOKE_ROOT / "summary" / s / f"{s}.feature_table.RDF_1.csv", s)
        for s in SAMPLES
    ]
    feature_table_out = qc_tables / "feature_table_long.csv"
    concat_tables(feature_table_inputs, feature_table_out)
    manifest_paths["feature_table_long"] = relative_to_export(feature_table_out)

    all_hits_inputs = [
        (SMOKE_ROOT / "summary" / s / f"{s}.all_hits.csv", s)
        for s in SAMPLES
    ]
    all_hits_out = qc_tables / "all_hits_long.csv"
    concat_tables(all_hits_inputs, all_hits_out)
    manifest_paths["all_hits_long"] = relative_to_export(all_hits_out)

    analysis_inputs = [
        (SMOKE_ROOT / "summary" / s / f"{s}_analysis.csv", s)
        for s in SAMPLES
    ]
    analysis_out = qc_tables / "analysis_long.csv"
    concat_tables(analysis_inputs, analysis_out)
    manifest_paths["analysis_long"] = relative_to_export(analysis_out)

    corr_inputs = [
        (SMOKE_ROOT / "summary" / s / "insertion_vs_neighborhood_correlations.txt", s)
        for s in SAMPLES
    ]
    copy_text_files(corr_inputs, qc_tables)
    manifest_paths["correlation_text_dir"] = relative_to_export(qc_tables)

    return manifest_paths


def write_manifest(browser_info: dict, classifier_info: dict, qc_info: dict) -> Path:
    manifest = {
        "project": "YTAB smoke test",
        "version": "smoke_test_v1",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "samples": SAMPLES,
        "groups": {
            "parent": SAMPLES,
        },
        "paths": {
            **browser_info,
            **classifier_info,
            **qc_info,
        },
        "warnings": WARNINGS,
    }

    manifest_path = EXPORT_ROOT / "manifest.json"
    ensure_dir(manifest_path.parent)
    manifest_path.write_text(json.dumps(manifest, indent=2))
    return manifest_path


# Main

def main() -> None:
    ensure_dir(EXPORT_ROOT)

    print("============================================================")
    print("Building YTAB smoke_test_v1 export bundle")
    print(f"YTAB root   : {YTAB_ROOT}")
    print(f"Smoke root  : {SMOKE_ROOT}")
    print(f"Export root : {EXPORT_ROOT}")
    print(f"Samples     : {', '.join(SAMPLES)}")
    print("============================================================")

    export_qc_raw()
    browser_info = export_browser()
    classifier_info = export_classifier()
    qc_info = export_qc_tables_and_images()
    manifest_path = write_manifest(browser_info, classifier_info, qc_info)

    print("\nDone.")
    print(f"Manifest: {manifest_path}")
    if WARNINGS:
        print("\nCompleted with warnings:")
        for msg in WARNINGS:
            print(f" - {msg}")


if __name__ == "__main__":
    main()
