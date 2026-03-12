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
