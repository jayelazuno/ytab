
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import pandas as pd
import streamlit as st


# Paths and loading helpers

APP_FILE = Path(__file__).resolve()
REPO_ROOT = APP_FILE.parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "output" / "exports" / "smoke_test_v1" / "manifest.json"


@st.cache_data(show_spinner=False)
def load_manifest(manifest_path: str) -> dict:
    path = Path(manifest_path)
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


@st.cache_data(show_spinner=False)
def load_table(path_str: str) -> pd.DataFrame:
    path = Path(path_str)
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


@st.cache_data(show_spinner=False)
def load_excel_sheet_names(path_str: str) -> list[str]:
    xls = pd.ExcelFile(path_str)
    return xls.sheet_names


@st.cache_data(show_spinner=False)
def load_excel_sheet(path_str: str, sheet_name: str) -> pd.DataFrame:
    return pd.read_excel(path_str, sheet_name=sheet_name)


def resolve_export_path(manifest_path: str, rel_path: str) -> Path:
    return Path(manifest_path).resolve().parent / rel_path


def maybe_load_df(manifest_path: str, rel_path: Optional[str]) -> Optional[pd.DataFrame]:
    if not rel_path:
        return None
    path = resolve_export_path(manifest_path, rel_path)
    if not path.exists():
        return None
    return load_table(str(path))


def metric_from_df(df: Optional[pd.DataFrame], column: str) -> Optional[float]:
    if df is None or df.empty or column not in df.columns:
        return None
    try:
        return float(df.iloc[0][column])
    except Exception:
        return None


def format_number(value: Optional[float], digits: int = 2) -> str:
    if value is None or pd.isna(value):
        return "NA"
    if abs(value) >= 1000:
        return f"{value:,.0f}"
    return f"{value:.{digits}f}"


# Small UI helpers


def show_table_preview(df: Optional[pd.DataFrame], title: str, n: int = 10) -> None:
    st.subheader(title)
    if df is None or df.empty:
        st.info("No data available.")
        return
    st.dataframe(df.head(n), use_container_width=True)


def show_image_if_exists(path: Path, caption: str) -> None:
    if path.exists():
        st.image(str(path), caption=caption, use_container_width=True)
    else:
        st.info(f"Missing image: {path.name}")


# App

st.set_page_config(page_title="YTAB", page_icon="🧬", layout="wide")
st.title("YTAB")
st.caption("Yeast Transposon Analysis Browser")

with st.sidebar:
    st.header("Data source")
    manifest_input = st.text_input("Manifest path", value=str(DEFAULT_MANIFEST))
    load_clicked = st.button("Load manifest", type="primary")

if "manifest_path" not in st.session_state:
    st.session_state.manifest_path = str(DEFAULT_MANIFEST)

if load_clicked:
    st.session_state.manifest_path = manifest_input

manifest_path = st.session_state.manifest_path
manifest_file = Path(manifest_path)

if not manifest_file.exists():
    st.error(f"Manifest not found: {manifest_file}")
    st.stop()

try:
    manifest = load_manifest(str(manifest_file))
except Exception as exc:
    st.error(f"Failed to load manifest: {exc}")
    st.stop()

samples = manifest.get("samples", [])
groups = manifest.get("groups", {})
paths = manifest.get("paths", {})
warnings = manifest.get("warnings", [])

with st.sidebar:
    st.success("Manifest loaded")
    st.write(f"**Project:** {manifest.get('project', 'YTAB')} ")
    st.write(f"**Version:** {manifest.get('version', 'NA')} ")
    st.write(f"**Samples:** {len(samples)}")
    if samples:
        selected_sample = st.selectbox("Sample", samples, index=0)
    else:
        selected_sample = None
else_block = None

# fallback if sidebar branch not entered as expected
if "selected_sample" not in locals():
    selected_sample = samples[0] if samples else None

mapping_df = maybe_load_df(manifest_path, paths.get("mapping_stats"))
libdiag_summary_df = maybe_load_df(manifest_path, paths.get("library_diagnostics_summary"))
libdiag_sample_df = maybe_load_df(manifest_path, paths.get("library_diagnostics_per_sample"))
midlc_df = maybe_load_df(manifest_path, paths.get("midlc_long"))
seqbias_df = maybe_load_df(manifest_path, paths.get("seqbias_long"))
centromere_df = maybe_load_df(manifest_path, paths.get("centromere_bins_long"))
tss_df = maybe_load_df(manifest_path, paths.get("tss_metaplot_long"))
tts_df = maybe_load_df(manifest_path, paths.get("tts_metaplot_long"))
trna_df = maybe_load_df(manifest_path, paths.get("trna_metaplot_long"))
summary_stats_df = maybe_load_df(manifest_path, paths.get("summary_stats"))
hit_summary_df = maybe_load_df(manifest_path, paths.get("hit_summary_long"))
feature_table_df = maybe_load_df(manifest_path, paths.get("feature_table_long"))
analysis_df = maybe_load_df(manifest_path, paths.get("analysis_long"))
classifier_df = maybe_load_df(manifest_path, paths.get("classifier_table"))


if selected_sample and mapping_df is not None and "sample" in mapping_df.columns:
    mapping_sample_df = mapping_df[mapping_df["sample"] == selected_sample].copy()
else:
    mapping_sample_df = mapping_df

if selected_sample and libdiag_sample_df is not None and "sample" in libdiag_sample_df.columns:
    libdiag_selected_df = libdiag_sample_df[libdiag_sample_df["sample"] == selected_sample].copy()
else:
    libdiag_selected_df = libdiag_sample_df

if selected_sample and midlc_df is not None and "sample" in midlc_df.columns:
    midlc_selected_df = midlc_df[midlc_df["sample"] == selected_sample].copy()
else:
    midlc_selected_df = midlc_df

if selected_sample and seqbias_df is not None and "sample" in seqbias_df.columns:
    seqbias_selected_df = seqbias_df[seqbias_df["sample"] == selected_sample].copy()
else:
    seqbias_selected_df = seqbias_df

if selected_sample and centromere_df is not None and "sample" in centromere_df.columns:
    centromere_selected_df = centromere_df[centromere_df["sample"] == selected_sample].copy()
else:
    centromere_selected_df = centromere_df

if selected_sample and tss_df is not None and "sample" in tss_df.columns:
    tss_selected_df = tss_df[tss_df["sample"] == selected_sample].copy()
else:
    tss_selected_df = tss_df

if selected_sample and tts_df is not None and "sample" in tts_df.columns:
    tts_selected_df = tts_df[tts_df["sample"] == selected_sample].copy()
else:
    tts_selected_df = tts_df

if selected_sample and trna_df is not None and "sample" in trna_df.columns:
    trna_selected_df = trna_df[trna_df["sample"] == selected_sample].copy()
else:
    trna_selected_df = trna_df

if selected_sample and summary_stats_df is not None and "sample" in summary_stats_df.columns:
    summary_stats_selected_df = summary_stats_df[summary_stats_df["sample"] == selected_sample].copy()
else:
    summary_stats_selected_df = summary_stats_df

if selected_sample and hit_summary_df is not None and "sample" in hit_summary_df.columns:
    hit_summary_selected_df = hit_summary_df[hit_summary_df["sample"] == selected_sample].copy()
else:
    hit_summary_selected_df = hit_summary_df

if selected_sample and feature_table_df is not None and "sample" in feature_table_df.columns:
    feature_table_selected_df = feature_table_df[feature_table_df["sample"] == selected_sample].copy()
else:
    feature_table_selected_df = feature_table_df

if selected_sample and analysis_df is not None and "sample" in analysis_df.columns:
    analysis_selected_df = analysis_df[analysis_df["sample"] == selected_sample].copy()
else:
    analysis_selected_df = analysis_df


tab_overview, tab_mapping, tab_diagnostics, tab_browser, tab_classifier = st.tabs([
    "Overview",
    "Mapping QC",
    "Library Diagnostics",
    "Insertion Browser",
    "Classifier",
])


with tab_overview:
    st.subheader("Smoke test overview")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Selected sample", selected_sample or "NA")
    c2.metric("Samples loaded", len(samples))
    c3.metric("Group", next((k for k, v in groups.items() if selected_sample in v), "NA") if selected_sample else "NA")
    c4.metric("Warnings", len(warnings))

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Total records", format_number(metric_from_df(mapping_sample_df, "total_records"), 0))
    m2.metric("Mapped %", format_number(metric_from_df(mapping_sample_df, "percent_mapped"), 2))
    m3.metric("MAPQ ≥20 %", format_number(metric_from_df(mapping_sample_df, "percent_mapq_ge20"), 2))
    m4.metric("MidLC est", format_number(metric_from_df(libdiag_selected_df, "midlc_est"), 2))

    st.markdown("### Export bundle")
    st.code(str(Path(manifest_path).resolve().parent), language="text")

    if warnings:
        with st.expander("Warnings"):
            for msg in warnings:
                st.write(f"- {msg}")

    col_a, col_b = st.columns(2)
    with col_a:
        show_table_preview(mapping_sample_df, "Mapping stats", n=5)
    with col_b:
        show_table_preview(libdiag_selected_df, "Library diagnostics summary", n=5)

    st.markdown("### Summary outputs")
    show_table_preview(summary_stats_selected_df, "Summary stats", n=20)


with tab_mapping:
    st.subheader("Mapping QC")

    if mapping_sample_df is None or mapping_sample_df.empty:
        st.info("No mapping stats found.")
    else:
        st.dataframe(mapping_sample_df, use_container_width=True)

        chart_cols = [
            c for c in [
                "total_records",
                "primary_mapped",
                "primary_unmapped",
                "mapq_ge20",
            ] if c in mapping_sample_df.columns
        ]
        if chart_cols:
            plot_df = mapping_sample_df.set_index("sample")[chart_cols]
            st.bar_chart(plot_df)

        pct_cols = [
            c for c in ["percent_mapped", "percent_duplicates", "percent_mapq_ge20"]
            if c in mapping_sample_df.columns
        ]
        if pct_cols:
            st.markdown("### Percentage metrics")
            st.bar_chart(mapping_sample_df.set_index("sample")[pct_cols])


with tab_diagnostics:
    st.subheader("Library diagnostics")

    diag_subtab1, diag_subtab2, diag_subtab3, diag_subtab4 = st.tabs([
        "MidLC",
        "Bias",
        "Metaplots",
        "Summary tables",
    ])

    with diag_subtab1:
        if midlc_selected_df is None or midlc_selected_df.empty:
            st.info("No MidLC table found.")
        else:
            st.dataframe(midlc_selected_df.head(25), use_container_width=True)
            numeric_cols = [c for c in midlc_selected_df.columns if c != "sample"]
            candidate_x = next((c for c in numeric_cols if "depth" in c.lower()), None)
            candidate_y = next((c for c in numeric_cols if "unique" in c.lower() or "site" in c.lower()), None)
            if candidate_x and candidate_y:
                st.line_chart(midlc_selected_df.set_index(candidate_x)[candidate_y])

    with diag_subtab2:
        left, right = st.columns(2)
        with left:
            st.markdown("#### Sequence bias")
            if seqbias_selected_df is None or seqbias_selected_df.empty:
                st.info("No sequence bias table found.")
            else:
                st.dataframe(seqbias_selected_df, use_container_width=True)

        with right:
            st.markdown("#### Centromere bias")
            if centromere_selected_df is None or centromere_selected_df.empty:
                st.info("No centromere bias table found.")
            else:
                st.dataframe(centromere_selected_df.head(30), use_container_width=True)

            qc_images = paths.get("qc_images", [])
            centromere_image = None
            for rel in qc_images:
                if selected_sample and rel.endswith(f"{selected_sample}.centromere_bias.png"):
                    centromere_image = resolve_export_path(manifest_path, rel)
                    break
            if centromere_image is not None:
                show_image_if_exists(centromere_image, f"{selected_sample} centromere bias")

    with diag_subtab3:
        st.markdown("#### TSS metaplot")
        if tss_selected_df is not None and not tss_selected_df.empty:
            st.dataframe(tss_selected_df.head(20), use_container_width=True)
        else:
            st.info("No TSS metaplot table found.")

        st.markdown("#### TTS metaplot")
        if tts_selected_df is not None and not tts_selected_df.empty:
            st.dataframe(tts_selected_df.head(20), use_container_width=True)
        else:
            st.info("No TTS metaplot table found.")

        st.markdown("#### tRNA metaplot")
        if trna_selected_df is not None and not trna_selected_df.empty:
            st.dataframe(trna_selected_df.head(20), use_container_width=True)
        else:
            st.info("No tRNA metaplot table found.")

        qc_images = paths.get("qc_images", [])
        metaplot_image = None
        for rel in qc_images:
            if selected_sample and rel.endswith(f"{selected_sample}.metaplots.png"):
                metaplot_image = resolve_export_path(manifest_path, rel)
                break
        if metaplot_image is not None:
            show_image_if_exists(metaplot_image, f"{selected_sample} metaplots")

    with diag_subtab4:
        show_table_preview(libdiag_selected_df, "Per-sample diagnostics summary", n=10)
        show_table_preview(summary_stats_selected_df, "Summary stats", n=20)
        show_table_preview(hit_summary_selected_df, "Hit summary", n=20)


with tab_browser:
    st.subheader("Insertion browser")
    st.info("This first app scaffold exposes the browser tracks and sample-level hit tables. IGV.js embedding can be added next.")

    tracks = paths.get("tracks", [])
    if tracks:
        track_table = pd.DataFrame({
            "track_file": [Path(t).name for t in tracks],
            "relative_path": tracks,
            "absolute_path": [str(resolve_export_path(manifest_path, t)) for t in tracks],
        })
        st.dataframe(track_table, use_container_width=True)
    else:
        st.info("No browser tracks listed in manifest.")

    st.markdown("### All hits")
    show_table_preview(analysis_selected_df, "Analysis table", n=20)
    show_table_preview(feature_table_selected_df, "Feature table", n=20)


with tab_classifier:
    st.subheader("Essentiality classifier")

    plot_rel = paths.get("classifier_plot")
    workbook_rel = paths.get("classifier_workbook")

    left, right = st.columns([1.2, 1])

    with left:
        if classifier_df is not None and not classifier_df.empty:
            st.markdown("#### Classifier table")
            st.dataframe(classifier_df, use_container_width=True)
        else:
            st.info("No classifier TSV found.")

    with right:
        if plot_rel:
            show_image_if_exists(resolve_export_path(manifest_path, plot_rel), "Classifier AUC plot")
        else:
            st.info("No classifier plot found.")

    if workbook_rel:
        workbook_path = resolve_export_path(manifest_path, workbook_rel)
        if workbook_path.exists():
            st.markdown("#### Workbook sheets")
            sheets = load_excel_sheet_names(str(workbook_path))
            selected_sheet = st.selectbox("Excel sheet", sheets)
            sheet_df = load_excel_sheet(str(workbook_path), selected_sheet)
            st.dataframe(sheet_df, use_container_width=True)

st.caption("YTAB app scaffold: smoke-test parent pool build")

