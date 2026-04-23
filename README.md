# YTAB

**YTAB** = **Yeast Transposon Analysis Browser**

A reproducible platform for exploring yeast transposon-seq data, including:

- library and sample QC
- insertion mapping summaries
- gene-level abundance and fitness summaries
- essentiality and screen analysis
- reproducible workflow outputs for downstream visualization

## Planned structure

- `app/` — browser / UI
- `pipeline/` — Nextflow workflow
- `src/` — core analysis code
- `resources/` — reference files and annotations
- `configs/` — parameter files
- `containers/` — Dockerfiles
- `docs/` — architecture and migration notes
- `examples/` — small demo inputs
- `tests/` — unit and integration tests

## Current status

Repository scaffold created. Migration from legacy Hermes/Tn-seq repo in progress.

## NOTE (move this later)
- I had problems installing the conda env that I exported from the hpc because of the mismatch between linux hpc env and Mac Os linux. 
- I ended up installing a lighter version 
    python3 -m venv ~/venvs/ytab-app
    source ~/venvs/ytab-app/bin/activate
    python -m pip install --upgrade pip setuptools wheel
    pip install streamlit pandas numpy matplotlib plotly pyyaml