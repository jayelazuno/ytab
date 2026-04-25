from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import pandas as pd
import streamlit as st

APP_FILE = Path(__file__).resolve()
REPO_ROOT = APP_FILE.parents[1] if len(APP_FILE.parents) > 1 else Path.cwd()
DEFAULT_MANIFEST = REPO_ROOT / "output" / "exports" / "smoke_test_v1" / "manifest.json"


@st.cache_data(show_spinner=False)
def load_manifest(manifest_path: str) -> dict:
    with Path(manifest_path).open("r", encoding="utf-8") as fh:
        return json.load(fh)


@st.cache_data(show_spinner=False)
def load_table(path_str: str) -> pd.DataFrame:
    path = Path(path_str)
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    return pd.read_csv(path)


def resolve_export_path(manifest_path: str, rel_path: str | None) -> Optional[Path]:
    if not rel_path:
        return None
    return Path(manifest_path).resolve().parent / rel_path


def maybe_load_df(manifest_path: str, rel_path: Optional[str]) -> Optional[pd.DataFrame]:
    path = resolve_export_path(manifest_path, rel_path)
    if path is None or not path.exists():
        return None
    try:
        return load_table(str(path))
    except Exception:
        return None


def show_image_if_exists(path: Optional[Path], caption: str) -> None:
    if path is not None and path.exists():
        st.image(str(path), caption=caption, use_container_width=True)
    else:
        st.info(f"Missing image: {caption}")


def show_download_table(manifest_path: str, rel_path: Optional[str], label: str) -> None:
    path = resolve_export_path(manifest_path, rel_path)
    if path is None or not path.exists():
        return
    with path.open("rb") as fh:
        st.download_button(
            label=label,
            data=fh.read(),
            file_name=path.name,
            mime="text/csv" if path.suffix.lower() != ".xlsx" else "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            use_container_width=True,
        )


def panel_dict(manifest: dict) -> dict:
    return manifest.get("qc_panels", {})


def panels_for_tab(manifest: dict, tab_name: str) -> dict:
    return {
        k: v for k, v in panel_dict(manifest).items()
        if v.get("tab") == tab_name
    }


def render_panel_gallery(manifest_path: str, manifest: dict, tab_name: str) -> None:
    panels = panels_for_tab(manifest, tab_name)
    if not panels:
        st.info("No panels registered for this tab in manifest.")
        return

    panel_labels = {k: v.get("title", k) for k, v in panels.items()}
    selected_key = st.radio(
        "Choose plot",
        options=list(panel_labels.keys()),
        format_func=lambda k: panel_labels[k],
        horizontal=True,
    )
    selected = panels[selected_key]

    image_path = resolve_export_path(manifest_path, selected.get("image"))
    show_image_if_exists(image_path, selected.get("title", selected_key))

    with st.expander("Download / inspect source table"):
        table_rel = selected.get("table")
        show_download_table(manifest_path, table_rel, "Download source table")
        df = maybe_load_df(manifest_path, table_rel)
        if df is not None and not df.empty:
            st.dataframe(df.head(30), use_container_width=True)
        else:
            st.info("No previewable table found for this panel.")


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

manifest = load_manifest(str(manifest_file))
samples = manifest.get("samples", [])
groups = manifest.get("groups", {})
paths = manifest.get("paths", {})
warnings = manifest.get("warnings", [])

with st.sidebar:
    st.success("Manifest loaded")
    st.write(f"**Project:** {manifest.get('project', 'YTAB')}")
    st.write(f"**Version:** {manifest.get('version', 'NA')}")
    if samples:
        selected_sample = st.selectbox("Sample", samples, index=0)
    else:
        selected_sample = None

if "selected_sample" not in locals():
    selected_sample = samples[0] if samples else None

mapping_df = maybe_load_df(manifest_path, paths.get("mapping_stats"))
libdiag_summary_df = maybe_load_df(manifest_path, paths.get("library_diagnostics_summary"))
summary_stats_df = maybe_load_df(manifest_path, paths.get("summary_stats"))

selected_group = next((k for k, v in groups.items() if selected_sample in v), "NA") if selected_sample else "NA"

tab_overview, tab_mapping, tab_diagnostics, tab_feature, tab_genome, tab_browser, tab_classifier = st.tabs([
    "Overview",
    "Mapping QC",
    "Library Diagnostics",
    "Feature Summary",
    "Genome-wide Distribution",
    "Insertion Browser",
    "Classifier",
])

with tab_overview:
    st.subheader("Export overview")
    c1, c2, c3 = st.columns(3)
    c1.metric("Selected sample", selected_sample or "NA")
    c2.metric("Group", selected_group)
    c3.metric("Warnings", len(warnings))

    if warnings:
        with st.expander("Warnings"):
            for msg in warnings:
                st.write(f"- {msg}")

    st.markdown("### Available QC image panels")
    panels_df = pd.DataFrame([
        {
            "panel_key": key,
            "title": value.get("title", key),
            "tab": value.get("tab", "NA"),
            "image": value.get("image", "")
        }
        for key, value in panel_dict(manifest).items()
    ])
    if not panels_df.empty:
        st.dataframe(panels_df, use_container_width=True)

    if mapping_df is not None and not mapping_df.empty:
        st.markdown("### Mapping stats preview")
        st.dataframe(mapping_df.head(20), use_container_width=True)

    if libdiag_summary_df is not None and not libdiag_summary_df.empty:
        st.markdown("### Library diagnostics summary preview")
        st.dataframe(libdiag_summary_df.head(20), use_container_width=True)

with tab_mapping:
    st.subheader("Mapping QC")
    render_panel_gallery(manifest_path, manifest, "Mapping QC")

with tab_diagnostics:
    st.subheader("Library Diagnostics")
    render_panel_gallery(manifest_path, manifest, "Library Diagnostics")

with tab_feature:
    st.subheader("Feature Summary")
    render_panel_gallery(manifest_path, manifest, "Feature Summary")

with tab_genome:
    st.subheader("Genome-wide Distribution")
    render_panel_gallery(manifest_path, manifest, "Genome-wide Distribution")

with tab_browser:
    st.subheader("Insertion Browser")
    tracks = paths.get("tracks", [])
    if tracks:
        browser_df = pd.DataFrame({
            "track_file": [Path(t).name for t in tracks],
            "relative_path": tracks,
            "absolute_path": [str(resolve_export_path(manifest_path, t)) for t in tracks],
        })
        st.dataframe(browser_df, use_container_width=True)
    else:
        st.info("No browser tracks listed in manifest yet.")

    show_download_table(manifest_path, paths.get("all_hits_long"), "Download all hits table")
    if summary_stats_df is not None and not summary_stats_df.empty:
        st.markdown("### Summary stats preview")
        st.dataframe(summary_stats_df.head(20), use_container_width=True)

with tab_classifier:
    st.subheader("Classifier")
    workbook_rel = paths.get("classifier_workbook")
    show_download_table(manifest_path, workbook_rel, "Download classifier workbook")
    st.info("Classifier table and image can be wired here after the next export update.")

st.caption("YTAB image-first scaffold")
