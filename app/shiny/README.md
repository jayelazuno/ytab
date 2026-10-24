# YTAB Shiny shell

`run_app.R` resolves the repository relative to its own location. The landing page supports new-project initialization through `shinyFiles`, dynamic existing-project discovery, and an optional demonstration project. Every workspace command uses one validator-backed active `project.yaml` path.

The workspace uses horizontal application tabs and nested stage tabs. Opening a project or changing a tab never runs scientific work. Raw SummaryTable success is the analysis-ready boundary; MidLC is in Normalization → MidLC (Classifier) and raw-summary CPM remains in Condition Comparisons → Fitness.

## Modules

| Module (top tab) | Subtabs | Code |
|---|---|---|
| Project | Samples · Logs & Provenance · Help | `R/ui_workspace.R` |
| Run | Samples · Reference · Mapping (mapping, hit files and SummaryTable in one run: `ytab_run_preprocessing.py`) · Run Profiles | `R/ui_preprocessing.R` |
| Quality Control | Overview (flagged QC table, `ytab qc`) · Reads & Mapping · Insertions & Genes · Saturation & Bias · Library Agreement | `R/qc_module.R`, `R/ui_qc.R`, `R/qc_library_correlation.R` |
| Normalization | Site Normalization (`ytab normalize-project`) · MidLC | `R/normalization_module.R` |
| Essentiality | Run (HMM; hybrid when labels exist, its gene features built from the chosen normalization: `ytab essentiality`, `ytab hybrid-inputs`) · Results | `R/essentiality_module.R` |
| Condition Comparisons | Design · Fitness: Run · Fitness: Results · Multi-condition · Genetic Interactions | `R/comparisons_module.R`, `R/ui_fitness.R` |
| Gene Explorer | Genome Browser · Gene & Domain Insertions · GO Enrichment · Across Species (incl. Essentiality across species: `ytab cross-species`) | `R/ui_workspace.R`, `R/cross_species.R` |
| Screen Planning | Coverage Planner · Mutant Identification | `R/screen_planning.R`, `R/mutant_screen.R` |

**Reports & exports** is a header button. The new analysis panels share `R/analysis_runs.R`:
- each run is a folder under the project (`library_correlation/`, `normalization/`, `essentiality/`, `multicondition/`, `genetic_interactions/`) with a manifest;
- results are shown as interactive plotly plots, with PNG/PDF and TSV downloads; downloads are publication quality
  (vector PDF through cairo, 600-dpi PNG);
- colours are defined once (`YTAB_*_COLOURS` in `R/analysis_runs.R`): species in the Okabe–Ito palette, the method under test in teal and comparators in
  grey, depleted navy and enriched brick red, QC flags green/yellow/red, one colour per normalization method;
- plots drawn from tables the commands write: Lorenz curve (QC › Saturation & Bias,
  `ytab qc` `plots/lorenz.tsv.gz`), gene-level replicate agreement (Normalization › Replicate agreement,
  `plots/pair_gene_agreement.tsv`, `pair_genes.tsv.gz`), insertion-target preference and ROC against curated labels
  (Essentiality › Results, `plots/target_model.tsv`, `plots/roc.tsv`; `ytab essentiality --refresh-plots RUN` adds
  them to an earlier run), and essentiality across species (Gene Explorer › Across Species, `ytab cross-species`:
  shared and species-specific essential sets, Cohen's κ between species, the ortholog table);
- the methods and validation for each panel are in a collapsed **Methods** section.

Multi-condition and genetic interactions read a Normalization run (BGC suggested). Navigation is defined once in `R/navigation.R` (`navigation_modules`, `navigation_inner_navsets` for GO and Across Species). Retired tab values (e.g. `preprocessing`, `go_analysis`, `genome_browser`) are translated, so old links still open the right place. The top-level menus are in `www/ytab_nav_menus.js`; bump its `?v=` cache-buster in `app.R` when they change.

The app must run from the `ytab-local` conda environment, which provides Python (yaml, scipy) and R with all packages.

FASTQs are never uploaded. Pasting a directory path is the quickest option for large or network-mounted storage; local folder browsing is also available. Choosing a directory does not scan it. **Scan FASTQs** explicitly performs one nonrecursive top-level scan, limited to 5,000 preview entries. Matching outer quotes are removed from pasted paths, while internal spaces are preserved. Execution uses raw `processx` argument vectors, and quoting is applied only to human-readable command previews.

Preprocessing sample selection is interactive and persists across nested tabs. Project-level `include` values remain read-only; stage selection is temporary. Presets cover included, parent, treated, incomplete, failed, and none. Each core stage presents selected/completed/remaining/failed counts and a per-sample status table. Raw SummaryTable success for every included sample defines analysis readiness. A stage whose own outputs are absent is still counted as success when a downstream output exists: a hit file or mapping stats imply the mapping ran (for example on the HPC, with BAMs left there), and a SummaryTable implies its hit file was made. Repository-relative built-in references are resolved from the repository root, and a complete existing Bowtie2 index is never rebuilt merely by opening or refreshing the workspace.

MapFastq, CreateHitFile, and SummaryTable use item-level persistent progress; LibraryDiagnostics reports reliable phases. Progress JSON is written atomically beneath `manifests/jobs`, and the current-job pointer supports progress recovery after the app is reopened. Confirmed item completions determine the progress bar. Approximate ETA uses the median current-run duration per input byte when sizes are available, otherwise median item duration; it is unavailable until the first non-cached completion. The global active-job banner polls every 1–1.5 seconds, while log panels retain only recent stdout/stderr lines. Cancellation preserves completed work and permits cache-aware resume.

Workspace tabs and nested preprocessing tabs use explicit stable values. **Change sample selection** navigates to the Samples tab while retaining selection; **Continue preprocessing** chooses the next dependency-valid tab; **Enter analysis browser** opens Quality Control → Mapping QC; and **Switch project** returns to landing without cancelling a running external job. Dry-run progress is explicitly validation-only: validated items may reach 100%, while executed and scientifically completed remain zero. Finished summaries can be dismissed from the session without deleting durable progress or logs.

Quality Control is implemented in `R/ui_qc.R`. Mapping and raw Summary QC use compact principal tables, with file paths and secondary summary metrics in collapsed details. Library Diagnostics owns a separate persistent QC selection and offers all eligible, parents, treated, clear, and custom selection. Preview and real execution are explicit; Force is a cache-bypass control only. Cache manifests are sample-set and input-hash aware, run-specific outputs are retained, and the shared job manager presents phase-level progress without fabricated item ETA. Concise result summaries keep full commands, stdout/stderr, hashes, and paths inside collapsed **Technical details**.

`R/qc_result_state.R` separates the live selector, a hash-validated matching result, and the latest historical result. Legacy unsigned manifests remain discoverable but cannot be labeled matching. Library Diagnostics uses three compact workflow sections, with sample review and Technical details collapsed by default.

Condition Comparisons keeps the fitness Design, Run and Results separate. Generated valid comparisons are searchable and temporary selection persists between those tabs without changing `comparison_design.csv`. Analysis IDs are generated and editable only under Advanced settings. Preview, real execution, cached reuse, Force cache bypass, Keep-going behavior, and optional classifier annotation are explicit. Fitness always uses raw SummaryTable inputs and CPM inside the established R analysis. Results distinguish matching, historical, and legacy output and offer existing-call summaries, filtering, stable-file downloads, a call chart, and an effect-versus-z plot. Technical provenance remains collapsed; no P values or volcano plot are added.
