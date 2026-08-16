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

## Gene & Domain Insertion Explorer

The Shiny app has a **Gene & Domain Insertion Explorer** tab for plotting
existing insertion hits across a queried gene. Search by gene name, systematic
name, locus tag, gene ID, or any alias available in the current project
annotation. Select one, multiple, or all available insertion tracks from raw
CreateHitFile outputs or from the combined parent hit file when present.

Domains are shown only when real domain annotations exist in the YTAB reference
resources. If no domain annotation is available for the selected gene, the gene
model and insertion tracks still plot and the app reports that no domain
annotation is available.

The CLI wrapper can also generate the same outputs:

```bash
python scripts/local/ytab_draw_gene_domain_insertions.py \
  --project-config output/projects/<PROJECT_ID>/config/project.yaml \
  --gene <GENE_OR_SYSTEMATIC_NAME> \
  --track-source raw \
  --samples all \
  --flank-bp 1000
```

Generated PNG figures, insertion tables, and manifests are saved under:

```text
output/projects/<PROJECT_ID>/gene_domain_explorer/
```

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

## Launch and operate the complete local application

```bash
cd "<repository-root>"
mamba activate ytab-local
./scripts/local/ytab_launch_app.sh \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml
```

The launcher defaults to the locally secure address `127.0.0.1:3838`. Use `--host` or `--port`
explicitly to change them. Manual launch is also supported:

```bash
Rscript app/shiny/run_app.R --host 127.0.0.1 --port 3838 \
  --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml
```

Inspect status and preview or resume the full cached workflow:

```bash
python scripts/local/ytab_project_status.py --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml --show-next
python scripts/local/ytab_run_pipeline.py --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml --profile all --threads 2 --dry-run --print-plan
python scripts/local/ytab_run_pipeline.py --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml --profile all --threads 2 --keep-going
```

Profiles `core`, `classifier`, `fitness`, and `all` run sequentially. The classifier profile uses the
parent MidLC branch; fitness uses raw SummaryTable inputs and CPM in R. Successful stage caches are
reused unless `--force` is supplied. Cancellation leaves partial outputs intact for later resumption.

Build presentation and reproducibility artifacts:

```bash
python scripts/local/ytab_build_project_report.py --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml
python scripts/local/ytab_export_project.py --project-config output/projects/local_glabrata_smoke_v1/config/project.yaml
```

Compact exports exclude FASTQ, BAM, SAM, raw hit, and normalized hit files by default. Outputs are
under `output/exports/<PROJECT_ID>/report/` and `output/exports/<PROJECT_ID>/bundles/`. Smoke-project
classifier and fitness calls validate software behavior and are not biological conclusions.
# Recommended application entry patterns

For large FASTQs or a computer with approximately 4–8 GB RAM, preprocess first:

```bash
mamba activate ytab-local
python scripts/local/ytab_run_pipeline.py \
  --project-config output/projects/<PROJECT_ID>/config/project.yaml \
  --profile core --threads 2 --keep-going
./scripts/local/ytab_launch_app.sh \
  --project-config output/projects/<PROJECT_ID>/config/project.yaml
```

To preprocess in the app, run `./scripts/local/ytab_launch_app.sh`, select **Create new project**, choose the FASTQ directory, initialize, and use the Preprocessing tab. The local directory chooser stores paths; it does not upload or copy FASTQs. Mapping is the longest and most memory-intensive stage, processes one sample at a time, and may be resumed after closing the app. Cached successful stages are skipped.

For mounted or network storage, pasting the FASTQ directory into **FASTQ directory path** and clicking **Use this path** is usually faster than browsing. **Browse local folders** remains available with focused local roots. Directory selection and scanning are separate: click **Scan FASTQs** to perform a nonrecursive, top-level filename and metadata scan. Paths containing spaces and pasted paths wrapped in matching quotes are supported. Project initialization does not start mapping automatically.

Raw SummaryTable completion defines analysis readiness. MidLC normalization is exclusive to the classifier branch. The treated-versus-parent branch uses raw per-sample SummaryTables and performs CPM normalization inside the R analysis.

The Preprocessing workspace uses a shared interactive sample selector. Project inclusion comes from `sample_sheet.csv`; temporary stage selection controls only the next run and does not edit inclusion. Presets select all included, parents, treated, incomplete, failed, or none. Mapping, hit-file creation, raw SummaryTable, and LibraryDiagnostics offer Run selected actions, while stage-specific Run all actions choose only dependency-valid samples. Dry run remains enabled by default. The status matrix reports FASTQ, mapping, BAM index, hit-file, and raw-summary state per sample. Built-in references are resolved relative to the repository root and complete Bowtie2 indexes are reused without rebuilding. Mapping continues one sample at a time by default.

Submitted preprocessing jobs write atomic progress records under `output/projects/<PROJECT_ID>/manifests/jobs/`; `manifests/orchestrator/current_job.json` points to the latest active record. The global workspace banner follows an active job across tabs and shows the current sample, confirmed processed count, elapsed time, and approximate ETA. The progress bar counts only finished, skipped, or failed items; an active sample remains separate until it finishes. ETA appears after one non-cached item completes and may vary with compression, storage speed, and system load. Cancelled jobs preserve completed outputs, and a later cache-aware resume skips successful samples.

**Change sample selection** returns to Preprocessing → Samples without clearing the temporary selection. **Switch project** returns to the YTAB landing page; internal navigation does not use browser history. A completed dry run is labeled **Dry run complete** and means commands were validated, not that scientific mapping finished. Completed job summaries can be dismissed without deleting progress files or logs. A one-sample mapping job generally finishes before it can provide a useful duration-based ETA; use two or more samples when observing ETA behavior is important.

## Quality Control and Library Diagnostics

Mapping QC and Summary QC show compact primary metrics without absolute paths. Summary QC secondary metrics are available under **Detailed metrics**. In Library Diagnostics, choose eligible samples with the searchable selector or the all, parents, treated, and clear presets. **Preview command** adds `--dry-run` and never records scientific completion. **Run diagnostics** performs real analysis unless an exact cache matches; **Force rerun** is optional and only bypasses that matching cache.

Cache identity includes the sorted sample set, hit-file hashes, relevant reference hashes, LibraryDiagnostics script hash, species, and runner/scientific parameters. Results are stored under `library_diagnostics/runs/<RUN_ID>/` and exported under `qc/diagnostics/<RUN_ID>/`. The concise result panel hides commands and absolute paths; expand **Technical details** when those are needed. **Diagnostic Files** inventories tables, plots, logs, and manifests with project-relative paths and lists existing plots without regenerating them.

The Library Diagnostics page keeps **Current selection** separate from stored results. A signed result is labeled matching only after its sample set and recorded input, reference, script, species, and parameter hashes validate; other completed results are historical. Unsigned older outputs are labeled **Legacy diagnostics result** and never claim to match the live selector. The eligibility table is collapsed under **Review eligible samples**, and **Technical details** remains collapsed unless troubleshooting is required.

## Fitness Screen

Generate or refresh parent–treated matches on the Design tab, review invalid or missing inputs, and choose a temporary comparison subset without rewriting the design CSV. Preview validates the design, raw SummaryTable inputs, and command without running R. Run executes or reuses an exact design-aware cache; Force rerun only bypasses that cache, and Keep going applies to multi-comparison failures. The fixed workflow is raw per-sample SummaryTable → CPM inside R. Optional classifier annotation adds labels without changing effects, z-scores, or calls. Historical and unsigned legacy results remain viewable but are not claimed as current matches. Smoke-project calls are software-validation output, not final biological conclusions.

Diagnostic Files shows either **File table** or **Plot gallery**, never both simultaneously. The table hides paths by default. The image-only gallery supports sample, plot-type, sample-set, and run filters, renders eight lazy thumbnails per page, and opens existing plots in a local modal. Viewing and filtering do not regenerate or modify plots.
