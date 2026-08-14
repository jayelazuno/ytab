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

## Local quick start

```bash
cd "<repository-root>"
mamba activate ytab-local
./scripts/local/ytab_launch_app.sh \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml
```

The app binds to `127.0.0.1:3838` by default. Another host or port must be requested explicitly.
Use `python scripts/local/ytab_project_status.py --project-config <project.yaml> --show-next`
to inspect resumable stage state without running analysis.

YTAB has two separate scientific branches: parent-only MidLC normalization feeds the essentiality
classifier, while treated-versus-parent fitness uses raw per-sample SummaryTable output and performs
CPM normalization inside R. The orchestrator preserves this separation and reuses stage caches.

The Fitness Screen exposes generated comparison designs, temporary interactive subsets, explicit Preview and Run modes, and design-aware cache reuse. Classifier labels are optional and do not change fitness calculations. Matching, historical, and legacy results remain distinct; existing calls can be searched, filtered, and viewed as a call distribution or effect-versus-z plot. No P values or volcano plot are introduced.

## NOTE (move this later)
- I had problems installing the conda env that I exported from the hpc because of the mismatch between linux hpc env and Mac Os linux. 
- I ended up installing a lighter version 
    python3 -m venv ~/venvs/ytab-app
    source ~/venvs/ytab-app/bin/activate
    python -m pip install --upgrade pip setuptools wheel
    pip install streamlit pandas numpy matplotlib plotly pyyaml
# YTAB application workflow

YTAB now opens on a local landing page. For low-memory computers or large FASTQ inputs, preprocess first with two threads and then open the app:

```bash
mamba activate ytab-local
python scripts/local/ytab_run_pipeline.py \
  --project-config output/projects/<PROJECT_ID>/config/project.yaml \
  --profile core --threads 2 --keep-going
./scripts/local/ytab_launch_app.sh \
  --project-config output/projects/<PROJECT_ID>/config/project.yaml
```

Alternatively, launch `./scripts/local/ytab_launch_app.sh`, choose **Create new project**, select a local FASTQ directory, initialize it, and continue preprocessing. FASTQs stay in place and are not uploaded. Mapping runs one sample at a time, can take substantial time, and can be resumed later; completed cached stages are skipped. Raw SummaryTable completion makes a project analysis-ready. MidLC belongs only to the essentiality-classifier branch, while treated-versus-parent fitness uses raw SummaryTables and performs CPM normalization inside its R analysis.

The workspace keeps an active-job banner visible across tabs. Persistent progress files record the current sample, confirmed completions, elapsed time, recent state, and an approximate ETA after the first non-cached item finishes. Estimates vary with storage and system load. Cancelling preserves completed outputs so later cache-aware runs can resume without repeating successful samples.

Quality Control presents compact Mapping QC and raw Summary QC tables; secondary summary metrics and project-relative source paths remain available under **Detailed metrics**. Library Diagnostics has an independent searchable sample selection with all-eligible, parent, treated, and custom presets. **Preview command** validates inputs without scientific execution; **Run diagnostics** performs or reuses an exact sample-set-aware cached run, while **Force rerun** only bypasses a matching cache. Concise results remain separate from collapsed technical logs and paths.
