# Local YTAB setup

YTAB uses a Shiny user interface with a Python computational backend. The local environment includes the required Python and R packages plus Bowtie2 and SAMtools.

## Install the environment

Install [Miniforge](https://github.com/conda-forge/miniforge) if `mamba` or `conda` is not already available. From the YTAB repository root, run:

```bash
bash scripts/setup_ytab_local.sh
mamba activate ytab-local
```

The setup script prefers `mamba` and falls back to `conda`. It creates the `ytab-local` environment on the first run and updates it on later runs.

## Validate the installation

With the environment active, run:

```bash
python scripts/local/ytab_check_local_env.py
python scripts/local/ytab_check_resources.py
python scripts/local/ytab_smoke_reference_discovery.py
```

The environment checker verifies the Python, R, Bowtie2, and SAMtools dependencies. The resource checks inspect filenames and metadata only; they do not run alignment.

## Reference resources and indexing

Each reference belongs under:

```text
resources/species/<species>/reference_genome/
```

A runnable reference needs a FASTA, at least one GFF, GTF, or feature table, and a complete Bowtie2 index. Existing complete six-file Bowtie2 indexes are detected and reused, so indexing is skipped. Bowtie2 indexes are built from the FASTA only; GFF, GTF, and feature tables are used for gene and feature annotation.

To explicitly prepare a reference without running alignment, use:

```bash
python scripts/local/ytab_prepare_reference.py --species glabrata --threads 4
```

This creates a FASTA `.fai` when needed, builds a missing Bowtie2 index, and writes `reference_prepare_manifest.json` in the species reference directory. If the FASTA is missing, indexing cannot proceed: restore, download, or provide the reference FASTA first.

Future species can be added dynamically under:

```text
resources/species/<species>/reference_genome/
```

A species becomes selectable when this directory contains at least one recognized reference file. Kornmann is legacy/resource data and is not a selectable species.

The placeholder resource manager reports how to supply resources until a release manifest is available:

```bash
python scripts/local/ytab_download_resources.py
```

## Start the Shiny app

From the repository root, run:

```bash
Rscript app/shiny/run_app.R
```

The current app initializes and validates a local project. It does not run alignment yet.

## Run local mapping

Activate the local environment and first inspect the commands without running alignment:

```bash
mamba activate ytab-local

python scripts/local/ytab_run_mapfastq.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --threads 2 \
  --dry-run
```

Run a one-sample smoke mapping:

```bash
python scripts/local/ytab_smoke_mapfastq_one_sample.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --sample yH298-parent-pool1 \
  --threads 2
```

Run all included samples while continuing past individual failures:

```bash
python scripts/local/ytab_run_mapfastq.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --threads 2 \
  --keep-going
```

Mapping processes one sample at a time to keep memory use low. Samples with successful manifests and expected outputs are skipped unless `--force` is supplied. A complete Bowtie2 index must already exist; otherwise run reference preparation first.

Per-sample outputs are written under `output/projects/<PROJECT_ID>/mapfastq/`, with logs and restartable manifests in the project directories. Step 2 does not create hit files or run downstream analysis.

## Create hit files

Activate the environment and inspect all CreateHitFile commands without executing them:

```bash
mamba activate ytab-local

python scripts/local/ytab_run_create_hit_file.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --threads 2 \
  --dry-run
```

Run a one-sample smoke test:

```bash
python scripts/local/ytab_smoke_create_hit_file_one_sample.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --sample yH298-parent-pool1 \
  --threads 2
```

Run all included samples while continuing past individual failures:

```bash
python scripts/local/ytab_run_create_hit_file.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --threads 2 \
  --keep-going
```

CreateHitFile processes one sample at a time. Its input BAMs come from `output/projects/<PROJECT_ID>/mapfastq/`, and hit files are written under `output/projects/<PROJECT_ID>/create_hit_file/`. Successful cached samples are skipped unless `--force` is supplied. Step 3 does not run SummaryTable or any later analysis.

## Build summary tables

Activate the environment and inspect all SummaryTable commands without executing them:

```bash
mamba activate ytab-local

python scripts/local/ytab_run_summary_table.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --threads 2 \
  --dry-run
```

Run a one-sample smoke test:

```bash
python scripts/local/ytab_smoke_summary_table_one_sample.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --sample yH298-parent-pool1 \
  --threads 2
```

Run all included samples while continuing past individual failures:

```bash
python scripts/local/ytab_run_summary_table.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --threads 2 \
  --keep-going
```

SummaryTable processes one sample at a time. Input hit files come from `output/projects/<PROJECT_ID>/create_hit_file/`, and feature-level tables are written under `output/projects/<PROJECT_ID>/summary/`. Successful cached samples are skipped unless `--force` is supplied. Step 4 does not run LibraryDiagnostics or later analysis.

## Run library diagnostics

Activate the environment and inspect the cross-library diagnostics command without running it:

```bash
mamba activate ytab-local

python scripts/local/ytab_run_library_diagnostics.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --dry-run
```

Run diagnostics on a two-sample smoke subset:

```bash
python scripts/local/ytab_smoke_library_diagnostics.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --samples yH298-parent-pool1,yH298-parent-pool2
```

Run diagnostics across all included samples:

```bash
python scripts/local/ytab_run_library_diagnostics.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml
```

LibraryDiagnostics runs once across the selected hit files. Inputs come from `output/projects/<PROJECT_ID>/create_hit_file/`, diagnostics are written under `output/projects/<PROJECT_ID>/library_diagnostics/`, and useful QC files are copied to `output/exports/<PROJECT_ID>/qc/diagnostics/`. A prior successful run for the same selected samples is skipped unless `--force` is supplied. Missing optional reference resources are reported and their corresponding diagnostics are omitted.

This step produces raw library QC only. It does not normalize hit files or run downstream analysis.

## Explore MidLC normalization

Activate the environment and inspect the adaptive parent-only normalization command:

```bash
mamba activate ytab-local

python scripts/local/ytab_run_sample_normalization.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --targets auto \
  --sample-mode parents \
  --min-site-retention 0.95 \
  --threads 2 \
  --dry-run
```

Run the automatic parent-only smoke workflow:

```bash
python scripts/local/ytab_smoke_sample_normalization.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --targets auto \
  --sample-mode parents
```

Verify a decimal target manually:

```bash
python scripts/local/ytab_smoke_sample_normalization.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --targets 59.7 \
  --sample-mode parents
```

Run an optional exploratory sweep:

```bash
python scripts/local/ytab_run_sample_normalization.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --targets 20,30,40,50,60,80,100 \
  --sample-mode parents \
  --threads 2 \
  --keep-going
```

Normalization uses parent samples by default because their normalized hit files feed the later Gale-style classifier workflow. The default `auto` mode recommends a conservative target based on insertion-site retention; it is not a fixed target list. Any finite positive integer or decimal target can also be requested manually. Existing successful target outputs are skipped unless `--force` is used.

Inputs come from `output/projects/<PROJECT_ID>/create_hit_file/`. Normalized hit files and comparison artifacts are written under `output/projects/<PROJECT_ID>/sample_normalization/`, with app-facing copies under `output/exports/<PROJECT_ID>/qc/sample_normalization/`.

Feature-level retention requires SummaryTable on normalized hits and will be confirmed in a later step. Step 6 does not run SummaryTable on normalized hits, combine parent pools, or run the classifier or treated-vs-parent analysis.

## Evaluate normalized hit files with SummaryTable

Inspect the recommended-target commands without running SummaryTable:

```bash
mamba activate ytab-local

python scripts/local/ytab_run_summary_normalized.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --targets recommended \
  --sample-mode parents \
  --threads 2 \
  --dry-run
```

Run the recommended-target parent smoke workflow:

```bash
python scripts/local/ytab_smoke_summary_normalized.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --targets recommended \
  --sample-mode parents
```

Evaluate a specific decimal-tag target:

```bash
python scripts/local/ytab_smoke_summary_normalized.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --targets T059p7 \
  --sample-mode parents
```

Optionally process every available normalization target:

```bash
python scripts/local/ytab_run_summary_normalized.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --targets all \
  --sample-mode parents \
  --threads 2 \
  --keep-going
```

This step runs SummaryTable independently on each existing normalized hit file, fills feature-level retention metrics, and refines the normalization recommendation without overwriting the Step 6 site-retention recommendation. Outputs are written under `output/projects/<PROJECT_ID>/summary_normalized/`; evaluation and recommendation artifacts remain under `output/projects/<PROJECT_ID>/sample_normalization/` and the corresponding QC export directories.

Step 7 does not rerun normalization, combine parent pools, run SummaryTable on combined hits, run the classifier, or perform treated-vs-parent analysis.

## Combine normalized parent hit files

Inspect the recommended-target combine command without creating a combined file:

```bash
mamba activate ytab-local

python scripts/local/ytab_run_combine_hits.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --target recommended \
  --dry-run
```

Combine all included parent samples using the recommended target:

```bash
python scripts/local/ytab_smoke_combine_hits.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --target recommended
```

Verify a manually selected decimal target:

```bash
python scripts/local/ytab_smoke_combine_hits.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --target T059p7
```

Step 8 combines normalized parent hit files for one selected target, using parent samples by default. The combined file remains compatible with SummaryTable and is the input for the next step.

This step does not rerun normalization or normalized SummaryTable, run SummaryTable on the combined parent file, run the classifier, or perform treated-vs-parent analysis.

## Summarize the combined normalized parent hit file

Inspect the recommended-target SummaryTable command without running it:

```bash
mamba activate ytab-local

python scripts/local/ytab_run_summary_combined.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --target recommended \
  --threads 2 \
  --dry-run
```

Run the recommended-target smoke workflow:

```bash
python scripts/local/ytab_smoke_summary_combined.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --target recommended
```

Verify a manually selected decimal target:

```bash
python scripts/local/ytab_smoke_summary_combined.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --target T059p7
```

Step 9 runs SummaryTable once on the selected target's combined normalized parent hit file. It writes the stable feature table `combined_feature_table.<TAG>.txt`, copies SummaryTable statistics when available, and records target-keyed manifest and status files under `output/projects/<PROJECT_ID>/summary_combined/` and `manifests/summary_combined/`.

This step does not run the classifier or perform treated-vs-parent analysis.

## Run the essentiality classifier

The classifier consumes one stable combined parent feature table created by Steps 8 and 9. Inspect the recommended-target command without running the classifier:

```bash
mamba activate ytab-local

python scripts/local/ytab_run_classifier.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --target recommended \
  --dry-run
```

Run the recommended-target and manual-target smoke workflows:

```bash
python scripts/local/ytab_smoke_classifier.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --target recommended

python scripts/local/ytab_smoke_classifier.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --target T059p7
```

The automatically recommended target is not saved as the final biological selection. After reviewing normalization retention and library QC, explicitly persist a reviewed target with:

```bash
python scripts/local/ytab_run_classifier.py \
  --project-config output/projects/<PROJECT_ID>/config/project.yaml \
  --target <TARGET_TAG> \
  --save-final-target
```

Use `--save-final-target` only after biological review. Smoke-subset calls validate software behavior and must not be interpreted as final essentiality calls. Step 10 does not compare treated and parent libraries or perform treated-versus-parent fitness analysis.

## Run treated-versus-parent fitness analysis

Step 11 uses the raw per-sample feature tables created by Step 4. CPM normalization occurs inside
`scripts/ytab_treated_vs_parent_screen.R`; MidLC normalization is classifier-specific and is not used
here. Pairings live in the user-editable `output/projects/<PROJECT_ID>/config/comparison_design.csv`
and should be reviewed before analysis. Classifier predictions are optional post-analysis annotations.

```bash
mamba activate ytab-local
python scripts/local/ytab_init_comparison_design.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml --print-design
python scripts/local/ytab_run_treated_vs_parent.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --analysis-id H2O2_vs_parent --input-mode auto --dry-run
python scripts/local/ytab_smoke_treated_vs_parent.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --comparison-id <VALID_COMPARISON_ID> --input-mode auto
python scripts/local/ytab_run_treated_vs_parent.py \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml \
  --analysis-id H2O2_vs_parent --input-mode auto --keep-going
```

Optional classifier annotation:

```bash
python scripts/local/ytab_run_treated_vs_parent.py \
  --project-config output/projects/<PROJECT_ID>/config/project.yaml \
  --analysis-id <ANALYSIS_ID> --input-mode raw-summary \
  --annotate-classifier --classifier-target <TARGET_TAG>
```

Results from the 10,000-read smoke subset validate software behavior only and must not be interpreted as a final H2O2 fitness screen.
