# Command-line reference

Every analysis in the app is a YTAB command, and the app only displays the tables the commands write. This page lists
every command with the options it accepts, taken from the command's own `--help`.

## Running YTAB

Activate the environment (see Setup), then from the repository folder:

```sh
python scripts/local/ytab.py <command> [options]     # or: PYTHONPATH=src python -m ytab <command> [options]
python scripts/local/ytab.py modules                 # list the commands
python scripts/local/ytab.py <command> --help        # options of one command
```

**What a command needs.** Most commands act on one project and take `--project-config output/projects/<project>/config/project.yaml`.
The project file names the species reference and the sample sheet (`config/sample_sheet.csv`: one row per library with its
role (parent, control, treated), condition, background, pool and transposon). Commands that combine projects take
`--project NAME` (repeatable). Results go inside the project; commands that list `--output` as required need a folder named.

**Options shared by many commands.**

| option | meaning |
|---|---|
| `--project-config FILE` | the project's `config/project.yaml` |
| `--threads N` | worker threads (mapping, hit files, summaries) |
| `--force` | redo a step whose results already exist |
| `--output DIR` | where to write the results (required by some commands, see their options) |
| `--help` | every option of the command, with its default |

**From FASTQ to results, one project:**

```sh
python scripts/local/ytab.py init --project-id P --fastq-dir FASTQ/ --species glabrata    # make the project from local FASTQ files
python scripts/local/ytab.py preprocess --project-config P/config/project.yaml --threads 4
python scripts/local/ytab.py qc --project-config P/config/project.yaml
python scripts/local/ytab.py normalize-project --project-config P/config/project.yaml --method betageom --method ttr
python scripts/local/ytab.py essentiality --project-config P/config/project.yaml --libraries LIB1 LIB2 LIB3
python scripts/local/ytab.py comparison-design --project-config P/config/project.yaml
python scripts/local/ytab.py treated-vs-control --project-config P/config/project.yaml --method resampling     --normalization-run P/normalization/RUN --normalization betageom --output P/treated_vs_control/RUN
python scripts/local/ytab.py browser-tracks --project-config P/config/project.yaml
python scripts/local/ytab.py app --project-config P/config/project.yaml   # open the results
```

`ytab run --workflow WORKFLOW.yaml` runs an ordered list of these commands (see Command line); `ytab run --project-config
P/config/project.yaml --profile core|classifier|fitness|all` runs the standard pipeline.

## Commands by stage

### Start a project

| command | what it does |
|---|---|
| `init` | Initialize a project from local FASTQs |
| `import-hits` | Initialize a project from existing hit files |
| `reference` | Prepare species reference resources |
| `status` | Inspect resumable project stage status |

### Run everything

| command | what it does |
|---|---|
| `pipeline` | Existing core/classifier/fitness/all project pipeline |

### Sequencing and mapping

| command | what it does |
|---|---|
| `preprocess` | Mapping, hit files and SummaryTable for selected samples in one run (skips stages already done) |
| `mapfastq` | FASTQ clipping and reference alignment |
| `hit-files` | Create unnormalized insertion hit files |
| `merge-runs` | Merge sequencing runs of one library into one sample (summed hit file; runs excluded) so analyses count it once |
| `summary` | Raw per-sample gene summaries |

### Quality control

| command | what it does |
|---|---|
| `qc` | Library QC table for a project: all QC metrics, flagged (rule of thumb / project-relative) and commented, plus plot data |
| `library-qc` | Library diagnostics and insertion architecture |
| `tnseq-stats` | Per-library Tn-seq count statistics plus saved Hermes QC |
| `library-correlation` | Correlate Hermes libraries across pools, strains, conditions and species (gene insertion density) |

### Normalization

| command | what it does |
|---|---|
| `normalize-project` | Normalize a project's sites with chosen methods in one step (catalogue, experiment matrix, plot summaries) |
| `normalize` | Explicit site-normalization choices and method guidance |
| `normalize-experiment` | Evaluate normalization on design-defined controls and treatments |
| `normalize-genes` | Gene-level raw/TTR/CPM validation without changing fitness defaults |
| `ttr` | Empirical site catalogue and TTR normalization evaluation |
| `bgc-audit` | BGC numerical, seed and catalogue-sensitivity diagnostics |
| `abgc-audit` | Adaptive BGC search, tail and simulation sensitivity diagnostics |
| `normalization` | Legacy parent-only MidLC thinning to a fixed target (superseded by the midlc method of normalize-project) |
| `summary-normalized` | Summaries of normalized parent libraries |
| `combine-hits` | Combine independent parent insertion tracks |
| `summary-combined` | Combined-parent gene summaries |

### Essentiality

| command | what it does |
|---|---|
| `essentiality` | Essentiality for a project in one step: Hermes HMM (label-free), current RF scores and the hybrid when labels are given |
| `hermes-hmm` | Label-free Hermes window HMM essentiality for any yeast (haploid reference FASTA + GFF + hit files) |
| `hybrid-inputs` | Hybrid essentiality gene features from parents normalized with any method (default BGC; reuses earlier results) |
| `classifier` | Existing essentiality classifier |
| `cross-species` | Essentiality across species: HMM calls of one-to-one orthologs (S. cerevisiae hub) in two or more projects |

### Condition comparisons

| command | what it does |
|---|---|
| `comparison-design` | Create explicit fitness comparisons |
| `fitness` | Treated vs control, local-abundance z-score for a project (CPM by default, or any site normalization) |
| `fitness-zscore` | Local-abundance z-score on explicit inputs (feature tables or a normalization run) |
| `treated-vs-control` | Treated vs control for a project by resampling or the control-conditioned Gumbel (normalization run + comparison design) |
| `resampling` | Hermes site-restricted resampling on a ttr/totreads/nzmean normalization output |
| `control-gumbel` | Experimental control-conditioned count-aware Hermes Gumbel on a ttr/totreads/nzmean normalization output |
| `multicondition-project` | Multi-condition tests and their plots in one step for a project (reference from the project config) |
| `multicondition` | Hermes multi-condition tests (permutation ANOVA and all-pairs resampling) on a normalization output |
| `multicondition-plots` | Response correlation (corrplot) and response heatmap for a multicondition run |
| `genetic-interactions` | Hermes genetic interactions: strain x treatment response differences |

### Gene Explorer

| command | what it does |
|---|---|
| `browser-tracks` | Genome browser files for a project: reference with display chromosome names, gene BED, insertion tracks and an IGV session |
| `gene-insertions` | Gene/domain insertion plots and queries |
| `go` | GO enrichment and optional semantic reduction |

### Screen planning and mutants

| command | what it does |
|---|---|
| `coverage` | Analytical Poisson / equal-gene screen benchmark |
| `mutants` | Universal paired Sanger/ABI mutant identification |

### Datasets, results and sharing

| command | what it does |
|---|---|
| `datasets` | Catalogue of the reprocessed public datasets (study, strain, transposon, accession, QC, results in the app) for the Datasets page |
| `report` | Local HTML report of existing project outputs |
| `export` | Portable project-result export |
| `resource-bundle` | Bundle reprocessed projects (results, no FASTQ/BAM) with a checksum manifest, to open them in the app elsewhere |
| `app` | Open the Shiny app on existing results |

### Pre-registered benchmarks

| command | what it does |
|---|---|
| `tvc-benchmark` | Pre-registered treated-vs-control benchmark: Gale-calibrated simulation, method x normalization arms, paired tests |
| `tvc-phase2` | Treated-vs-control benchmark phase 2: every arm on real Gale FLU, a swap null, H2O2/Zn, and Gale's published scores |

## Options of each command

### init

Initialize a project from local FASTQs. Stage: Start a project.

```text
usage: ytab_init_local_project.py [-h] --project-id PROJECT_ID
                                  --fastq-dir FASTQ_DIR --species SPECIES
                                  [--threads THREADS] [--repo-root REPO_ROOT]
                                  [--force]

Initialize (but do not run) a local YTAB project.

options:
  -h, --help            show this help message and exit
  --project-id PROJECT_ID
  --fastq-dir FASTQ_DIR
  --species SPECIES
  --threads THREADS
  --repo-root REPO_ROOT
  --force               Overwrite existing config files; never delete outputs.
```

### import-hits

Initialize a project from existing hit files. Stage: Start a project.

```text
usage: ytab_import_existing_hit_project.py [-h] --project-id PROJECT_ID
                                           --species SPECIES --hit-dir HIT_DIR
                                           [--threads THREADS]
                                           [--sample-metadata SAMPLE_METADATA]
                                           [--repo-root REPO_ROOT] [--force]

Initialize a YTAB project from existing CreateHitFile outputs.

options:
  -h, --help            show this help message and exit
  --project-id PROJECT_ID
  --species SPECIES
  --hit-dir HIT_DIR
  --threads THREADS
  --sample-metadata SAMPLE_METADATA
  --repo-root REPO_ROOT
  --force               Overwrite config/manifests only; never delete hit
                        outputs.
```

### reference

Prepare species reference resources. Stage: Start a project.

```text
usage: ytab_prepare_reference.py [-h] --species SPECIES [--threads THREADS]
                                 [--repo-root REPO_ROOT] [--force]

Prepare FASTA and Bowtie2 indexes for one local YTAB reference.

options:
  -h, --help            show this help message and exit
  --species SPECIES
  --threads THREADS
  --repo-root REPO_ROOT
  --force
```

### status

Inspect resumable project stage status. Stage: Start a project. Takes `--project-config`.

```text
usage: ytab_project_status.py [-h] --project-config PROJECT_CONFIG [--json]
                              [--refresh]
                              [--profile {core,classifier,fitness,all}]
                              [--show-next]

Inspect YTAB project stage status without running analysis.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --json
  --refresh
  --profile {core,classifier,fitness,all}
  --show-next
```

### pipeline

Existing core/classifier/fitness/all project pipeline. Stage: Run everything. Takes `--project-config`.

```text
usage: ytab_run_pipeline.py [-h] --project-config PROJECT_CONFIG
                            [--profile {core,classifier,fitness,all}]
                            [--target TARGET] [--analysis-id ANALYSIS_ID]
                            [--stages STAGES] [--from-stage FROM_STAGE]
                            [--to-stage TO_STAGE] [--threads THREADS]
                            [--force] [--dry-run] [--keep-going]
                            [--print-plan]

Run existing YTAB stages sequentially with cache-aware runners.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --profile {core,classifier,fitness,all}
  --target TARGET
  --analysis-id ANALYSIS_ID
  --stages STAGES
  --from-stage FROM_STAGE
  --to-stage TO_STAGE
  --threads THREADS
  --force
  --dry-run
  --keep-going
  --print-plan
```

### preprocess

Mapping, hit files and SummaryTable for selected samples in one run (skips stages already done). Stage: Sequencing and mapping. Takes `--project-config`.

```text
usage: ytab_run_preprocessing.py [-h] --project-config PROJECT_CONFIG
                                 [--samples SAMPLES] [--threads THREADS]
                                 [--force] [--dry-run] [--keep-going]
                                 [--job-id JOB_ID]
                                 [--progress-file PROGRESS_FILE]
                                 [--junction-mate {R1,R2}]
                                 [--require-mate-alignment | --no-require-mate-alignment]

Map reads, create hit files and build SummaryTables for selected samples in one run (the app's Run > Mapping).

Each stage uses its own cache-aware runner (ytab_run_mapfastq.py, ytab_run_create_hit_file.py,
ytab_run_summary_table.py), so finished samples are skipped. A sample whose later outputs already exist is not sent
to the earlier stages: a hit file means it was mapped (e.g. on the HPC, with the BAM left there), and a SummaryTable
means its hit file was made; a sample with a SummaryTable is not run at all. A sample that fails a stage does not
continue to the next one.

  python scripts/local/ytab.py preprocess --project-config output/projects/P/config/project.yaml --samples s1,s2

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --samples SAMPLES     Comma-separated included sample names (default: all
                        included)
  --threads THREADS
  --force               Rerun every stage for the selected samples
  --dry-run
  --keep-going
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
  --junction-mate {R1,R2}
  --require-mate-alignment, --no-require-mate-alignment
```

### mapfastq

FASTQ clipping and reference alignment. Stage: Sequencing and mapping. Takes `--project-config`.

```text
usage: ytab_run_mapfastq.py [-h] --project-config PROJECT_CONFIG
                            [--samples SAMPLES] [--threads THREADS] [--force]
                            [--dry-run] [--keep-going] [--job-id JOB_ID]
                            [--progress-file PROGRESS_FILE]
                            [--junction-mate {R1,R2}]
                            [--require-mate-alignment | --no-require-mate-alignment]
                            [--hermes-cleanup | --no-hermes-cleanup]
                            [--transposon TRANSPOSON]
                            [--library-design LIBRARY_DESIGN]

Run MapFastq locally, one included project sample at a time.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --samples SAMPLES     Comma-separated included sample names
  --threads THREADS
  --force
  --dry-run
  --keep-going
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
  --junction-mate {R1,R2}
                        Override the designated Hermes/junction read for
                        selected samples (default R1)
  --require-mate-alignment, --no-require-mate-alignment
                        Require concordant pairs; default uses only the
                        designated junction mate
  --hermes-cleanup, --no-hermes-cleanup
                        Automatically trim observed Hermes/adaptor sequence
                        (default enabled)
  --transposon TRANSPOSON
                        Transposon profile in resources/transposons (default:
                        the project's, else hermes)
  --library-design LIBRARY_DESIGN
                        Library design within the transposon profile (default:
                        the profile's default)
```

### hit-files

Create unnormalized insertion hit files. Stage: Sequencing and mapping. Takes `--project-config`.

```text
usage: ytab_run_create_hit_file.py [-h] --project-config PROJECT_CONFIG
                                   [--samples SAMPLES] [--threads THREADS]
                                   [--force] [--dry-run] [--keep-going]
                                   [--job-id JOB_ID]
                                   [--progress-file PROGRESS_FILE]

Run CreateHitFile locally, one included project sample at a time.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --samples SAMPLES     Comma-separated included sample names
  --threads THREADS
  --force
  --dry-run
  --keep-going
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
```

### merge-runs

Merge sequencing runs of one library into one sample (summed hit file; runs excluded) so analyses count it once. Stage: Sequencing and mapping. Takes `--project-config`.

```text
usage: ytab_merge_runs.py [-h] --project-config PROJECT_CONFIG
                          --library LIBRARY --runs RUNS [RUNS ...] [--force]

Merge sequencing runs of one library into one project sample, so analyses count the library once.

A library sequenced on several lanes or runs has one row per run in the sample sheet. Comparisons and essentiality
treat each included sample as an independent library, so the runs must be summed first: this writes
<project>/create_hit_file/<library>/<library>_hits.txt (read counts summed per chromosome, strand and position), adds the
library to the sample sheet (role, condition, strain and pool taken from its runs, which must agree) and excludes the
runs, keeping their rows and QC. Then build the library's gene summary with `ytab summary --samples LIBRARY`.

  ytab merge-runs --project-config output/projects/P/config/project.yaml --library P_lib1 --runs RUN1 RUN2

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --library LIBRARY     Name of the merged sample (letters, numbers, - and _)
  --runs RUNS [RUNS ...]
                        Sample IDs of the runs of this library
  --force               Replace an existing merged sample of this name
```

### summary

Raw per-sample gene summaries. Stage: Sequencing and mapping. Takes `--project-config`.

```text
usage: ytab_run_summary_table.py [-h] --project-config PROJECT_CONFIG
                                 [--samples SAMPLES] [--threads THREADS]
                                 [--force] [--dry-run] [--keep-going]
                                 [--job-id JOB_ID]
                                 [--progress-file PROGRESS_FILE]

Run SummaryTable locally, one included project sample at a time.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --samples SAMPLES     Comma-separated included sample names
  --threads THREADS
  --force
  --dry-run
  --keep-going
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
```

### qc

Library QC table for a project: all QC metrics, flagged (rule of thumb / project-relative) and commented, plus plot data. Stage: Quality control. Takes `--project-config`.

```text
usage: ytab_qc.py [-h] --project-config PROJECT_CONFIG --output OUTPUT

Library QC table for a project: every per-library QC metric in one flagged, commented table plus plot data.

  python scripts/local/ytab.py qc --project-config output/projects/P/config/project.yaml --output output/projects/P/qc/RUN

Writes qc_table.tsv (one row per library), qc_flags.tsv (one row per library x metric: value, flag, reason, source,
comment, normalization hint), qc_table.xlsx (coloured, with the reasons as cell notes), metrics.tsv and plots/
(count histogram, QQ against a geometric fit and ranked counts). See ytab.qc.library_qc for the rules.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --output OUTPUT
```

### library-qc

Library diagnostics and insertion architecture. Stage: Quality control. Takes `--project-config`.

```text
usage: ytab_run_library_diagnostics.py [-h] --project-config PROJECT_CONFIG
                                       [--samples SAMPLES] [--threads THREADS]
                                       [--force] [--dry-run] [--job-id JOB_ID]
                                       [--progress-file PROGRESS_FILE]

Run raw LibraryDiagnostics across selected local project hit files.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --samples SAMPLES     Comma-separated included sample names
  --threads THREADS     Reserved runner resource setting
  --force
  --dry-run
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
```

### tnseq-stats

Per-library Tn-seq count statistics plus saved Hermes QC. Stage: Quality control.

```text
usage: ytab_tnseq_stats.py [-h] --catalogue CATALOGUE --output OUTPUT
                           [--diagnostics PROJECT=DIR]

Export per-library Tn-seq count statistics (adapted from TRANSIT's
tnseq_stats) plus saved Hermes QC (Excel-readable TSV).

options:
  -h, --help            show this help message and exit
  --catalogue CATALOGUE
                        Existing raw coordinate matrix and metadata directory
  --output OUTPUT       New output directory
  --diagnostics PROJECT=DIR
                        Repeat for existing library_diagnostics directories;
                        exact metadata project key
```

### library-correlation

Correlate Hermes libraries across pools, strains, conditions and species (gene insertion density). Stage: Quality control.

```text
usage: ytab_library_correlation.py [-h] [--project PROJECT]
                                   [--samples SAMPLES] [--orthologs ORTHOLOGS]
                                   [--metric {density,reads_per_kb}]
                                   [--min-reads MIN_READS]
                                   [--bootstrap BOOTSTRAP] [--seed SEED]
                                   --output OUTPUT

Correlate Hermes libraries across pools, strains, conditions and species, generalized from TRANSIT's corrplot (DeJesus et al. 2015).

Each library is summarized per gene (insertion density by default) and libraries are compared by Spearman
correlation across genes. Different species are matched through 1:1 ortholog groups (--orthologs).
Each pair is labelled same pool (shared clones), same strain/different pool, different strain or
different species; same-pool pairs measure library structure, not biology.

Libraries come from YTAB projects (--project NAME: every included sample, with condition, background
(strain) and pool from the sample sheet) and/or an explicit table (--samples TSV with columns
name, hits, species, gff and optionally strain, pool, condition). Hit files: YTAB *_hits.txt or any
CSV/TSV with chromosome, position and read-count columns.

Example (three C. glabrata projects and S. pombe Hermes data):
  ytab library-correlation --project Gale_2020 --project H2O2_screen_v1 --project Zn_toxicity_screen \
      --samples pombe.tsv --orthologs resources/species/orthologs/glabrata_cerevisiae_pombe_1to1.tsv \
      --output output/library_correlation/RUN

options:
  -h, --help            show this help message and exit
  --project PROJECT     YTAB project name under output/projects; repeat
  --samples SAMPLES     TSV: name, hits, species, gff[, strain, pool,
                        condition]; repeat
  --orthologs ORTHOLOGS
                        TSV of 1:1 ortholog groups (group, species, gene);
                        needed for several species
  --metric {density,reads_per_kb}
                        density (insertion sites per kb; default, comparable
                        across pools) or reads_per_kb (pool-specific clone
                        abundance)
  --min-reads MIN_READS
                        Reads for a site to count as an insertion (default 1)
  --bootstrap BOOTSTRAP
  --seed SEED
  --output OUTPUT
```

### normalize-project

Normalize a project's sites with chosen methods in one step (catalogue, experiment matrix, plot summaries). Stage: Normalization.

```text
usage: ytab_normalize_project.py [-h] [--project-config PROJECT_CONFIG]
                                 [--method {ttr,totreads,nzmean,midlc,quantile,betageom,aBGC,zinfnb}]
                                 [--output OUTPUT] [--trim TRIM]
                                 [--target TARGET] [--refresh-summaries RUN]
                                 [--normalize-together] [--sample SAMPLE]

Normalize one project's insertion-site counts with chosen methods, in one step (used by the app's Normalization tab).

Runs the existing, separately validated stages in order and adds plot-ready summaries:
  1. parent site catalogue (`ytab ttr`): every included parent/control sample; library = project::background::pool;
  2. experiment matrix and methods (`ytab normalize-experiment`): samples from config/comparison_design.csv, each
     method applied to the same raw matrix, scaled within each library; unnormalized counts (nonorm) are always
     included as the baseline;
  3. summaries for plotting: per sample and method, the change from its own raw counts (median and IQR of
     log2(normalized / raw) at its insertion sites, sites kept), positive-count distributions and matched-pair diagnostics;
     and each method's native before/after view (native()): TRANSIT-style per-sample statistics, the statistic a scaling
     method equalizes, quantiles against a fitted geometric and ranked counts (BGC, aBGC, unnormalized), and saturation
     curves with the thinning target (MidLC); and, for each matched pair, gene-level agreement under each method
     (gene_agreement(): per-gene log2 ratios, their median and the rank correlation of gene sums).
Without --method every available method is run, so any of them can be shown. --refresh-summaries RUN rewrites the plot
summaries of an earlier run (no renormalization).

The normalized output for downstream tests (multicondition, genetic-interactions, resampling) is OUTPUT/experiment/methods.

Example:
  ytab normalize-project --project-config output/projects/Zn_toxicity_screen/config/project.yaml \
      --method betageom --method ttr --output output/projects/Zn_toxicity_screen/normalization/run1

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --method {ttr,totreads,nzmean,midlc,quantile,betageom,aBGC,zinfnb}
                        Repeat; default every available method. Unnormalized
                        counts are always added
  --output OUTPUT       New output directory
  --trim TRIM
  --target TARGET
  --refresh-summaries RUN
                        Rewrite the plot summaries of an earlier run and exit
  --normalize-together  With --sample: normalize all samples as one group
                        (independently made libraries, e.g. one per strain)
  --sample SAMPLE       Normalize exactly these samples (repeat; roles from
                        the sample sheet) instead of the comparison design:
                        for multi-condition tests between strains or
                        conditions that are not pool-matched pairs
```

### normalize

Explicit site-normalization choices and method guidance. Stage: Normalization.

```text
usage: ytab_normalize.py [-h] [--list-methods]
                         [--describe {nonorm,ttr,totreads,nzmean,cpm,midlc,quantile,betageom,aBGC,zinfnb,emphist,QMgeom}]
                         [--catalogue CATALOGUE] [--output OUTPUT]
                         [--method {nonorm,ttr,totreads,nzmean,cpm,midlc,quantile,betageom,aBGC,zinfnb,emphist,QMgeom}]
                         [--group-by {library_id,project,dataset,cohort}]
                         [--observation OBSERVATION] [--trim TRIM]
                         [--target TARGET] [--bgc-cdf {analytic,simulation}]
                         [--bgc-seed BGC_SEED] [--bgc-draws BGC_DRAWS]
                         [--bgc-no-ttr] [--abgc-cdf {analytic,simulation}]
                         [--abgc-seed ABGC_SEED] [--abgc-draws ABGC_DRAWS]
                         [--abgc-tail {endpoint,midpoint,transit}]
                         [--abgc-fractions ABGC_FRACTIONS] [--abgc-no-ttr]

Explicit normalization choices and method explanations; no automatic
replacement of existing workflows.

options:
  -h, --help            show this help message and exit
  --list-methods
  --describe {nonorm,ttr,totreads,nzmean,cpm,midlc,quantile,betageom,aBGC,zinfnb,emphist,QMgeom}
  --catalogue CATALOGUE
                        Existing empirical catalogue directory containing
                        raw_counts.npz
  --output OUTPUT       New, separate output directory
  --method {nonorm,ttr,totreads,nzmean,cpm,midlc,quantile,betageom,aBGC,zinfnb,emphist,QMgeom}
                        Repeat to compare methods; no default method
  --group-by {library_id,project,dataset,cohort}
                        Required explicit scaling groups; library_id preserves
                        pool stratification
  --observation OBSERVATION
                        Optional explicit project::sample selection; repeat as
                        needed
  --trim TRIM
  --target TARGET
  --bgc-cdf {analytic,simulation}
  --bgc-seed BGC_SEED
  --bgc-draws BGC_DRAWS
  --bgc-no-ttr          Disable BGC's explicitly documented final TTR scaling
  --abgc-cdf {analytic,simulation}
  --abgc-seed ABGC_SEED
  --abgc-draws ABGC_DRAWS
  --abgc-tail {endpoint,midpoint,transit}
                        Explicit geometric scoring quantiles: endpoint (r/n;
                        'transit' is an accepted earlier name) or finite
                        midpoint sensitivity
  --abgc-fractions ABGC_FRACTIONS
                        Optional comma-separated trim fractions; default
                        0,0.02,...,0.30
  --abgc-no-ttr
```

### normalize-experiment

Evaluate normalization on design-defined controls and treatments. Stage: Normalization.

```text
usage: ytab_normalize_experiment.py [-h] --project-config PROJECT_CONFIG
                                    --catalogue CATALOGUE --output OUTPUT
                                    --method METHOD [--trim TRIM]
                                    [--target TARGET] [--no-plots]
                                    [--normalize-together] [--sample SAMPLE]

Evaluate explicit site-normalization choices on design-defined treatments and
controls.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --catalogue CATALOGUE
  --output OUTPUT
  --method METHOD
  --trim TRIM
  --target TARGET
  --no-plots
  --normalize-together  With --sample: one normalization group for all samples
                        (independent libraries)
  --sample SAMPLE       Normalize exactly these samples (repeat; roles from
                        the sample sheet) instead of the comparison design
```

### normalize-genes

Gene-level raw/TTR/CPM validation without changing fitness defaults. Stage: Normalization.

```text
usage: ytab_normalize_genes.py [-h] --validation VALIDATION --output OUTPUT
                               [--min-count MIN_COUNT]
                               [--fold-threshold FOLD_THRESHOLD]
                               [--existing-fitness EXISTING_FITNESS]
                               [--no-plots]
                               [--site-normalization SITE_NORMALIZATION]
                               [--compare-independent-pools]

Compare gene-level raw, TTR and CPM results without changing the fitness
workflow.

options:
  -h, --help            show this help message and exit
  --validation VALIDATION
                        Completed normalize-experiment output directory
  --output OUTPUT
  --min-count MIN_COUNT
  --fold-threshold FOLD_THRESHOLD
  --existing-fitness EXISTING_FITNESS
                        Optional existing fitness directory for exact default
                        CPM parity checks
  --no-plots
  --site-normalization SITE_NORMALIZATION
                        Optional normalize output with zero-preserving
                        quantile results from the identical raw matrix,
                        ordered samples and library_id groups
  --compare-independent-pools
                        Also compare independent pools sharing an explicit
                        background and comparison_group; between-pool controls
                        are not technical noise
```

### ttr

Empirical site catalogue and TTR normalization evaluation. Stage: Normalization.

```text
usage: ytab_ttr.py [-h] --manifest MANIFEST --output OUTPUT [--trim TRIM]
                   [--target TARGET] [--no-plots]

Build a parent site universe and evaluate site-level TTR, without changing
other analyses.

options:
  -h, --help           show this help message and exit
  --manifest MANIFEST  Explicit source projects, parent observations and
                       independent library IDs
  --output OUTPUT      New, versioned output directory
  --trim TRIM          Fraction trimmed from each tail of positive counts
                       (default .05)
  --target TARGET
  --no-plots
```

### bgc-audit

BGC numerical, seed and catalogue-sensitivity diagnostics. Stage: Normalization.

```text
usage: ytab_bgc_audit.py [-h] --site-run SITE_RUN --output OUTPUT
                         [--seed SEED] [--draws DRAWS]

Audit completed analytic BGC site results: simulation, fit and catalogue
sensitivity.

options:
  -h, --help           show this help message and exit
  --site-run SITE_RUN
  --output OUTPUT
  --seed SEED
  --draws DRAWS
```

### abgc-audit

Adaptive BGC search, tail and simulation sensitivity diagnostics. Stage: Normalization.

```text
usage: ytab_abgc_audit.py [-h] --site-run SITE_RUN --output OUTPUT
                          [--seed SEED] [--draws DRAWS]

Audit adaptive BGC endpoint and simulation sensitivity on an existing site
run.

options:
  -h, --help           show this help message and exit
  --site-run SITE_RUN
  --output OUTPUT
  --seed SEED
  --draws DRAWS
```

### normalization

Legacy parent-only MidLC thinning to a fixed target (superseded by the midlc method of normalize-project). Stage: Normalization. Takes `--project-config`.

```text
usage: ytab_run_sample_normalization.py [-h] --project-config PROJECT_CONFIG
                                        [--targets TARGETS]
                                        [--sample-mode {parents,all,treated}]
                                        [--samples SAMPLES]
                                        [--threads THREADS] [--force]
                                        [--dry-run] [--keep-going]
                                        [--min-site-retention MIN_SITE_RETENTION]
                                        [--auto-min-target AUTO_MIN_TARGET]
                                        [--auto-max-target AUTO_MAX_TARGET]
                                        [--auto-step AUTO_STEP]
                                        [--job-id JOB_ID]
                                        [--progress-file PROGRESS_FILE]

Explore local MidLC normalization targets and insertion-site retention.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --targets TARGETS     auto or comma-separated positive numeric values
  --sample-mode {parents,all,treated}
  --samples SAMPLES     Comma-separated included sample names; overrides
                        sample mode
  --threads THREADS
  --force
  --dry-run
  --keep-going
  --min-site-retention MIN_SITE_RETENTION
  --auto-min-target AUTO_MIN_TARGET
  --auto-max-target AUTO_MAX_TARGET
  --auto-step AUTO_STEP
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
```

### summary-normalized

Summaries of normalized parent libraries. Stage: Normalization. Takes `--project-config`.

```text
usage: ytab_run_summary_normalized.py [-h] --project-config PROJECT_CONFIG
                                      [--targets TARGETS]
                                      [--sample-mode {parents,all,treated}]
                                      [--samples SAMPLES] [--threads THREADS]
                                      [--force] [--dry-run] [--keep-going]
                                      [--min-feature-retention MIN_FEATURE_RETENTION]
                                      [--job-id JOB_ID]
                                      [--progress-file PROGRESS_FILE]

Run SummaryTable on normalized hit files and evaluate normalization targets.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --targets TARGETS
  --sample-mode {parents,all,treated}
  --samples SAMPLES     Comma-separated included sample names; overrides
                        sample mode
  --threads THREADS
  --force
  --dry-run
  --keep-going
  --min-feature-retention MIN_FEATURE_RETENTION
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
```

### combine-hits

Combine independent parent insertion tracks. Stage: Normalization. Takes `--project-config`.

```text
usage: ytab_run_combine_hits.py [-h] --project-config PROJECT_CONFIG
                                [--target TARGET] [--samples SAMPLES]
                                [--force] [--dry-run] [--job-id JOB_ID]
                                [--progress-file PROGRESS_FILE]

Combine normalized parent hit files for one selected normalization target.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --target TARGET
  --samples SAMPLES     Comma-separated included parent sample names
  --force
  --dry-run
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
```

### summary-combined

Combined-parent gene summaries. Stage: Normalization. Takes `--project-config`.

```text
usage: ytab_run_summary_combined.py [-h] --project-config PROJECT_CONFIG
                                    [--target TARGET] [--threads THREADS]
                                    [--force] [--dry-run] [--job-id JOB_ID]
                                    [--progress-file PROGRESS_FILE]

Run SummaryTable on one combined normalized parent hit file.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --target TARGET
  --threads THREADS
  --force
  --dry-run
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
```

### essentiality

Essentiality for a project in one step: Hermes HMM (label-free), current RF scores and the hybrid when labels are given. Stage: Essentiality.

```text
usage: ytab_essentiality.py [-h] --project-config PROJECT_CONFIG
                            [--libraries LIBRARIES [LIBRARIES ...]]
                            [--hmm-run HMM_RUN] [--features FEATURES]
                            [--rf RF] [--with-classifier] [--labels LABELS]
                            [--normalization NORMALIZATION]
                            [--normalization-run NORMALIZATION_RUN]
                            [--train-tiers TRAIN_TIERS] [--window WINDOW]
                            [--min-reads MIN_READS]
                            [--sequence-covariate {auto,motif,learned}]
                            [--output OUTPUT] [--refresh-plots RUN]

Essentiality for one project in one step (used by the app's Essentiality tab): Hermes window HMM, then the model table.

  1. HMM (label-free, any yeast; `ytab hermes-hmm`) on independent parent libraries of one strain, or reuse an
     earlier HMM fit (--hmm-run);
  2. model table (`ytab_essentiality_model_menu.py`): HMM calls, plus the current random forest's scores when its
     prediction table is given (--rf), plus the hybrid (random forest on gene features + HMM posteriors) when
     curated labels are also given (--labels);
  3. plot summaries: window posteriors and insertion counts along each chromosome (plots/windows.tsv.gz); the
     insertion-target preference (plots/target_model.tsv); with --labels, ROC curves of each model against the curated
     labels (plots/roc.tsv). --refresh-plots RUN writes the last two for an earlier run.

  ytab essentiality --project-config output/projects/Gale_2020/config/project.yaml \
      --libraries SRR11561508 SRR11561507 SRR11561506 --output output/projects/Gale_2020/essentiality/run1
  ... --labels resources/species/glabrata/essentiality_labels/labels.tsv            (adds the hybrid)

With --labels and no --rf, the hybrid's gene features are built first (`ytab hybrid-inputs`) from the same parent
libraries as the HMM, normalized with --normalization (default betageom; any method, or midlc). --rf instead uses an
existing classifier prediction table.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --libraries LIBRARIES [LIBRARIES ...]
                        Independent parent libraries of one strain (sample
                        IDs)
  --hmm-run HMM_RUN     Reuse an earlier HMM fit (a folder with genes.tsv,
                        fit.json, windows.npz)
  --features FEATURES   Gene features for the hybrid: a combined SummaryTable
                        feature table (default: built from the HMM's parents
                        with --normalization)
  --rf RF               Legacy classifier predictions
                        (essentiality_predictions.<target>.csv): adds its
                        score for comparison
  --with-classifier     Also score the built gene features with the legacy
                        classifier (comparison only)
  --labels LABELS       Curated labels (gene_id or cglab_id, label, tier),
                        e.g. resources/species/<species>/essentiality_labels/l
                        abels.tsv: adds the hybrid
  --normalization NORMALIZATION
                        Normalization for the hybrid's gene features when
                        --features is not given (method id, e.g. betageom, ttr
                        or midlc)
  --normalization-run NORMALIZATION_RUN
                        Site Normalization run for --normalization (default:
                        newest with the parents)
  --train-tiers TRAIN_TIERS
  --window WINDOW
  --min-reads MIN_READS
  --sequence-covariate {auto,motif,learned}
                        HMM sequence covariate: motif (Hermes nTnnnnAn, v2),
                        learned (target preference from the libraries' own
                        sites, v3), auto (learned when a library's sample-
                        sheet transposon is not Hermes, else motif)
  --output OUTPUT
  --refresh-plots RUN   Write the target-preference and ROC summaries of an
                        earlier run and exit
```

### hermes-hmm

Label-free Hermes window HMM essentiality for any yeast (haploid reference FASTA + GFF + hit files). Stage: Essentiality.

```text
usage: ytab_hermes_hmm.py [-h] [--fasta FASTA] [--gff GFF] [--hits HITS]
                          [--project PROJECT]
                          [--libraries LIBRARIES [LIBRARIES ...]]
                          [--window WINDOW] [--min-reads MIN_READS] [--each]
                          [--allow-high-mask]
                          [--sequence-covariate {motif,learned}] [--rescore]
                          output

Run the Hermes window HMM (label-free essentiality) on Hermes libraries of any yeast species.

Inputs: the genome FASTA the reads were aligned to, its GFF annotation, and one or more Hermes hit files
(YTAB *_hits.txt, or any CSV/TSV with chromosome, position and read-count columns). Libraries are fitted
jointly; use independent parent libraries of one strain.

  ytab hermes-hmm OUTPUT --fasta genome.fna --gff genes.gff --hits lib1_hits.txt --hits lib2_hits.txt
  ytab hermes-hmm OUTPUT --project PROJECT --libraries SAMPLE1 SAMPLE2      (genome and GFF from the project config)
With no inputs: Gale 2020 parents Cg1-Cg3 on the C. glabrata reference (the validated configuration).

Checks for other species:
- the FASTA must be the haploid reference used for alignment. A reference with two near-identical
  haplotypes (e.g. full diploid C. albicans) makes almost every 50-mer occur twice, so almost every window
  is masked. The run stops when > 50% of windows are masked (--allow-high-mask to override);
- hit-file chromosome names must match the FASTA (stops if > 10% of insertion sites are on unknown
  chromosomes);
- genes on chromosomes absent from the FASTA (e.g. the C. albicans B haplotype) are excluded and counted.
Ploidy: in a diploid, a heterozygous insertion in an essential gene is viable, so essential genes keep
insertions. Use haploid strains (or haploid derivatives) for essentiality.
Sequence covariate: `--sequence-covariate motif` (default, hermes_hmm_v2) counts the Hermes nTnnnnAn target;
`learned` (hermes_hmm_v3) learns the target-site preference from the libraries' own insertion sites (YTAB hit files
with their W/C strand), so any transposon (Hermes, Ac/Ds, ...) can be modelled; it also writes target_model.tsv.
Writes genes.tsv, windows.npz, fit.json. Reads no essentiality labels.
  ytab hermes-hmm FIT_DIR --rescore     rewrite FIT_DIR/genes.tsv with the current gene-scoring defaults (no refit)

positional arguments:
  output

options:
  -h, --help            show this help message and exit
  --fasta FASTA         Haploid genome FASTA the reads were aligned to
  --gff GFF             GFF annotation of that genome
  --hits HITS           Hermes hit file; repeat for each library
  --project PROJECT     YTAB project whose create_hit_file outputs and
                        reference to use
  --libraries LIBRARIES [LIBRARIES ...]
                        Independent parent libraries of one strain in
                        --project (sample IDs)
  --window WINDOW
  --min-reads MIN_READS
  --each                also fit each library alone
  --allow-high-mask     run even if > 50% of windows are unmappable
  --sequence-covariate {motif,learned}
                        motif: Hermes nTnnnnAn count (v2); learned: target
                        preference learned from the libraries' sites (v3, any
                        transposon)
  --rescore             OUTPUT is a saved fit (fit.json, windows.npz): rewrite
                        its genes.tsv with the current gene-scoring defaults,
                        without refitting
```

### hybrid-inputs

Hybrid essentiality gene features from parents normalized with any method (default BGC; reuses earlier results). Stage: Essentiality. Takes `--project-config`.

```text
usage: ytab_hybrid_inputs.py [-h] --project-config PROJECT_CONFIG
                             [--normalization NORMALIZATION]
                             [--normalization-run NORMALIZATION_RUN]
                             --samples SAMPLES [--threads THREADS] [--force]
                             [--with-classifier]

Build the hybrid essentiality model's gene features from the parents, normalized with any method.

The parents are normalized with the chosen method, written as hit files, combined and summarized (SummaryTable); the
result (summary_combined/<tag>/combined_feature_table.<tag>.txt) is what `ytab essentiality --features` takes.
--with-classifier also scores it with the legacy classifier (classifier/<tag>/essentiality_predictions.<tag>.csv),
for comparison only. Every method, MidLC common saturation included, is computed from the newest Site
Normalization run whose parent catalogue holds the parents (see ytab.pipeline.hybrid_inputs for the whole-number
rules). Results already built for the same method, run and parents are reused without being touched.

  python scripts/local/ytab.py hybrid-inputs --project-config output/projects/P/config/project.yaml \
      --normalization betageom --samples s1,s2,s3

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --normalization NORMALIZATION
                        Normalization method id (see ytab normalize --list-
                        methods), e.g. betageom, ttr or midlc
  --normalization-run NORMALIZATION_RUN
                        Site Normalization run to take the parent catalogue
                        from (default: newest that has the parents)
  --samples SAMPLES     Comma-separated parent libraries (one per independent
                        pool)
  --threads THREADS
  --force
  --with-classifier     Also score the features with the legacy classifier
                        (comparison only)
```

### classifier

Existing essentiality classifier. Stage: Essentiality. Takes `--project-config`.

```text
usage: ytab_run_classifier.py [-h] --project-config PROJECT_CONFIG
                              [--target TARGET] [--seed SEED] [--force]
                              [--dry-run] [--save-final-target]
                              [--save-final-target-only] [--job-id JOB_ID]
                              [--progress-file PROGRESS_FILE]

Run essentiality classification on a combined parent feature table.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --target TARGET
  --seed SEED
  --force
  --dry-run
  --save-final-target
  --save-final-target-only
                        Record an existing reviewed successful result without
                        rerunning classification
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
```

### cross-species

Essentiality across species: HMM calls of one-to-one orthologs (S. cerevisiae hub) in two or more projects. Stage: Essentiality.

```text
usage: ytab_cross_species.py [-h] --project-config PROJECT_CONFIG
                             [--essentiality-run ESSENTIALITY_RUN]
                             --output OUTPUT

Essentiality across species: label-free HMM calls of one-to-one orthologs in two or more projects (used by the app's
Gene Explorer > Across Species > Essentiality across species).

Genes are matched through one hub, S. cerevisiae systematic names (YAL001C ...): each project's reference config names
an ortholog table (`reference.orthology_file`) that maps its genes to S. cerevisiae; an S. cerevisiae project is the
hub itself. Only one-to-one orthologs are kept (a gene mapped to exactly one S. cerevisiae gene and back, per species),
and only hub genes with an ortholog scored by the HMM in every species, so each species is compared on the same genes.

Ortholog tables: tab-separated, '#' comments allowed. With a header that has a `scer_id` column, the first column is the
species' gene and `scer_id` its ortholog; without one, each row's first field is the species' gene and the first later
field that is a single S. cerevisiae systematic name is its ortholog (rows naming several genes, e.g. "YBR265W+YDR302W",
are ambiguous and skipped). This covers CGD, PomBase and YTAB-curated tables; a new species needs one such table.

Calls come from each project's essentiality run (newest unless --essentiality-run is given): P_ES_run and the HMM call
(ES or not); genes the HMM did not score (unmappable) are left out.

Outputs in OUTPUT/:
  orthologs.tsv      hub gene and, per species, the gene, its name, P_ES_run and ES (1/0)
  intersections.tsv  genes by the set of species in which they are essential (shared and species-specific sets)
  pairwise.tsv       per species pair: genes, both essential, either essential, Cohen's kappa of the calls (Cohen
                     1960) and the Spearman correlation of P_ES_run
  manifest.json

  ytab cross-species --project-config output/projects/Gale_2020/config/project.yaml \
      --project-config output/projects/Pombe_Hermes/config/project.yaml --output output/cross_species/run1

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
                        Repeat: one project per species (at least two)
  --essentiality-run ESSENTIALITY_RUN
                        Repeat in the same order as --project-config (default:
                        each project's newest run)
  --output OUTPUT
```

### comparison-design

Create explicit fitness comparisons. Stage: Condition comparisons. Takes `--project-config`.

```text
usage: ytab_init_comparison_design.py [-h] --project-config PROJECT_CONFIG
                                      [--overwrite] [--print-design]
                                      [--controls {one,all}]

Create and validate a parent-treated comparison design.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --overwrite
  --print-design
  --controls {one,all}  one: each treated library with the single parent of
                        its pool; all: with every control of its background
                        and pool (or its parents), one comparison per treated
                        condition (dose series, several conditions from one
                        pool)
```

### fitness

Treated vs control, local-abundance z-score for a project (CPM by default, or any site normalization). Stage: Condition comparisons. Takes `--project-config`.

```text
usage: ytab_run_treated_vs_parent.py [-h] --project-config PROJECT_CONFIG
                                     [--analysis-id ANALYSIS_ID]
                                     [--comparison-design COMPARISON_DESIGN]
                                     [--comparisons COMPARISONS]
                                     [--input-mode {auto,raw-summary}]
                                     [--annotate-classifier]
                                     [--classifier-target CLASSIFIER_TARGET]
                                     [--force] [--dry-run] [--keep-going]
                                     [--extra-r-arg EXTRA_R_ARG]
                                     [--normalization NORMALIZATION]
                                     [--normalization-run NORMALIZATION_RUN]
                                     [--job-id JOB_ID]
                                     [--progress-file PROGRESS_FILE]

Run raw-SummaryTable CPM treated-versus-parent fitness analysis.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --analysis-id ANALYSIS_ID
  --comparison-design COMPARISON_DESIGN
  --comparisons COMPARISONS
                        Comma-separated comparison or pair IDs
  --input-mode {auto,raw-summary}
  --annotate-classifier
  --classifier-target CLASSIFIER_TARGET
  --force
  --dry-run
  --keep-going
  --extra-r-arg EXTRA_R_ARG
  --normalization NORMALIZATION
                        cpm (default, from raw SummaryTables) or a Site
                        Normalization method, e.g. betageom
  --normalization-run NORMALIZATION_RUN
                        Site Normalization run directory (with --normalization
                        other than cpm)
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
```

### fitness-zscore

Local-abundance z-score on explicit inputs (feature tables or a normalization run). Stage: Condition comparisons.

```text
usage: ytab_fitness_zscore.py [-h] [--project-id PROJECT_ID]
                              [--analysis-id ANALYSIS_ID]
                              --comparison-design COMPARISON_DESIGN
                              [--comparisons COMPARISONS]
                              --output-dir OUTPUT_DIR
                              [--sample-table SAMPLE_TABLE]
                              [--normalized NORMALIZED]
                              [--normalization NORMALIZATION] [--gff GFF]
                              [--pseudocount PSEUDOCOUNT]
                              [--z-quantile Z_QUANTILE]
                              [--combine {replicates,pairs}]
                              [--recenter {none,trimmed-sum}]

Treated-vs-control local-abundance z-score (ytab.fitness.zscore): depleted and enriched genes, replicates averaged.

Default input: raw per-sample SummaryTable feature tables, as CPM (--sample-table SAMPLE=FILE for every design sample).
Any site normalization instead: --normalized RUN --normalization METHOD --gff GFF (gene totals of the method's
normalized site counts from a Site Normalization run).

  python scripts/local/ytab.py fitness-zscore --comparison-design output/projects/P/config/comparison_design.csv \
      --output-dir OUT --sample-table s1=output/projects/P/summary/s1/s1.feature_table.RDF_1.csv ...

options:
  -h, --help            show this help message and exit
  --project-id PROJECT_ID
  --analysis-id ANALYSIS_ID
  --comparison-design COMPARISON_DESIGN
  --comparisons COMPARISONS
                        Comma-separated comparison IDs
  --output-dir OUTPUT_DIR
  --sample-table SAMPLE_TABLE
                        SAMPLE=raw SummaryTable feature table (CPM input)
  --normalized NORMALIZED
                        Site Normalization run directory (instead of CPM)
  --normalization NORMALIZATION
                        Method in that run, e.g. betageom or ttr
  --gff GFF             GFF of the reference (with --normalized)
  --pseudocount PSEUDOCOUNT
  --z-quantile Z_QUANTILE
  --combine {replicates,pairs}
                        replicates: average each condition's replicates, one z
                        per gene (default); pairs: per-pair calls
  --recenter {none,trimmed-sum}
                        Scalar recentering of the treated values
```

### treated-vs-control

Treated vs control for a project by resampling or the control-conditioned Gumbel (normalization run + comparison design). Stage: Condition comparisons.

```text
usage: ytab_treated_vs_control.py [-h] --project-config PROJECT_CONFIG
                                  --method {resampling,gumbel}
                                  --normalization-run NORMALIZATION_RUN
                                  [--normalization NORMALIZATION]
                                  [--comparisons COMPARISONS] [--fdr FDR]
                                  --output OUTPUT

Treated-vs-control gene calls by resampling or the control-conditioned Gumbel (used by the app's Treated vs Control tab).

Reads a Site Normalization run (<project>/normalization/<run>) and the comparison design. The selected comparisons are
grouped by background and treated condition; each group is one test with the group's controls and treated libraries
(resampling permutes within each library when the group spans several pools). Genome FASTA and GFF come from the project.

  resampling   site-restricted permutation test, depleted and enriched; any normalization (BGC with trimmed-sum
               recentering by default). Validation status: docs/manuscript (pre-registered benchmark).
  gumbel       control-conditioned Gumbel, depletion runs only; scalar normalizations (TTR default, total reads,
               nonzero mean); needs >= 2 controls per group.

Output (ytab.fitness.contract.validate_tvc): manifest.json, genes.tsv (all groups) and <group>/sites.tsv.gz.

  ytab treated-vs-control --project-config output/projects/P/config/project.yaml --method resampling \
      --normalization-run output/projects/P/normalization/RUN --normalization betageom --output OUT

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --method {resampling,gumbel}
  --normalization-run NORMALIZATION_RUN
  --normalization NORMALIZATION
  --comparisons COMPARISONS
                        Comma-separated comparison IDs (default: all included)
  --fdr FDR
  --output OUTPUT
```

### resampling

Hermes site-restricted resampling on a ttr/totreads/nzmean normalization output. Stage: Condition comparisons.

```text
usage: ytab_resampling.py [-h] --normalized NORMALIZED
                          [--normalization {ttr,totreads,nzmean,betageom,aBGC}]
                          --gff GFF --reference REFERENCE --output OUTPUT
                          --control CONTROL --treatment TREATMENT
                          [--inference {exact,monte-carlo}]
                          [--opportunities {any-support,control-support}]
                          [--samples SAMPLES] [--grid-bins GRID_BINS]
                          [--min-sites MIN_SITES] [--seed SEED]
                          [--alternative {two-sided,depleted,enriched}]
                          [--recenter {none,median-gene,trimmed-sum}]
                          [--winsorize] [--pseudocount PSEUDOCOUNT]
                          [--fdr FDR] [--nterm NTERM] [--cterm CTERM]
                          [--recenter-trim RECENTER_TRIM]
                          [--asymmetry-warning ASYMMETRY_WARNING]
                          [--stratify-libraries] [--include-stop]

Hermes site-restricted resampling, adapted from TRANSIT's resampling test
(DeJesus et al. 2015) to empirical Hermes sites. Input is a YTAB `normalize`
output: ttr (default), totreads, nzmean, or the beta-geometric corrections
betageom / aBGC (BGC; recommended with --recenter trimmed-sum when replicate
noise differs between conditions). No normalization is run here. A warning is
printed when replicate noise is unequal.

options:
  -h, --help            show this help message and exit
  --normalized NORMALIZED
                        YTAB normalize output directory, or one method folder
                        inside it (e.g. .../ttr)
  --normalization {ttr,totreads,nzmean,betageom,aBGC}
                        Method folder to use when --normalized is a run
                        directory (default ttr)
  --gff GFF
  --reference REFERENCE
  --output OUTPUT
  --control CONTROL     Control observation ID; repeat for each replicate
  --treatment TREATMENT
                        Treatment observation ID; repeat for each replicate
  --inference {exact,monte-carlo}
  --opportunities {any-support,control-support}
  --samples SAMPLES
  --grid-bins GRID_BINS
  --min-sites MIN_SITES
  --seed SEED
  --alternative {two-sided,depleted,enriched}
  --recenter {none,median-gene,trimmed-sum}
                        Rescale the treatment once so typical genes have no
                        change (guards against normalization shifts under
                        widespread effects)
  --winsorize           Label-invariant winsorization: the largest site x
                        sample value of each gene is replaced by the second
                        largest
  --pseudocount PSEUDOCOUNT
  --fdr FDR
  --nterm NTERM
  --cterm CTERM
  --recenter-trim RECENTER_TRIM
  --asymmetry-warning ASYMMETRY_WARNING
  --stratify-libraries  Several independent libraries (e.g. paired pools):
                        permute labels only within each library_id
  --include-stop
```

### control-gumbel

Experimental control-conditioned count-aware Hermes Gumbel on a ttr/totreads/nzmean normalization output. Stage: Condition comparisons.

```text
usage: ytab_control_conditioned_gumbel.py [-h] --normalized NORMALIZED
                                          [--normalization {ttr,totreads,nzmean}]
                                          --gff GFF --reference REFERENCE
                                          --output OUTPUT --control CONTROL
                                          --treatment TREATMENT
                                          [--inference {exact,gumbel-mc}]
                                          [--calibration-sims CALIBRATION_SIMS]
                                          [--validation-sims VALIDATION_SIMS]
                                          [--min-control-reads MIN_CONTROL_READS]
                                          [--min-sites MIN_SITES]
                                          [--seed SEED] [--nodes NODES]
                                          [--alpha-site ALPHA_SITE]
                                          [--fdr FDR]
                                          [--tail-tolerance TAIL_TOLERANCE]
                                          [--ks-confidence KS_CONFIDENCE]
                                          [--nterm NTERM] [--cterm CTERM]
                                          [--eligibility {control-support,detectable}]
                                          [--treatment-dispersion {replicates,control}]
                                          [--dispersion-model {linear-quadratic,quadratic}]
                                          [--include-stop]

Experimental control-conditioned Hermes Gumbel on combined replicate tracks.
Input is a YTAB `normalize` output (ttr default, totreads or nzmean); no
normalization is run here.

options:
  -h, --help            show this help message and exit
  --normalized NORMALIZED
                        YTAB normalize output directory, or one method folder
                        inside it (e.g. .../ttr)
  --normalization {ttr,totreads,nzmean}
                        Method folder to use when --normalized is a run
                        directory (default ttr); must match when it is a
                        method folder
  --gff GFF
  --reference REFERENCE
  --output OUTPUT
  --control CONTROL     Matched control observation ID; repeat for each
                        control (>=2); combined into the null universe
  --treatment TREATMENT
                        Treatment observation ID; repeat for each replicate;
                        combined for the calls
  --inference {exact,gumbel-mc}
                        exact (default): deterministic exact conditional tail
                        p, Gumbel reported as a diagnostic; gumbel-mc: fitted
                        Gumbel p with Monte Carlo-aware adequacy gates (20000
                        simulations by default)
  --calibration-sims CALIBRATION_SIMS
                        default 1000 for exact (diagnostic only), 20000 for
                        gumbel-mc
  --validation-sims VALIDATION_SIMS
                        default 1000 for exact (diagnostic only), 20000 for
                        gumbel-mc
  --min-control-reads MIN_CONTROL_READS
  --min-sites MIN_SITES
  --seed SEED
  --nodes NODES
  --alpha-site ALPHA_SITE
  --fdr FDR
  --tail-tolerance TAIL_TOLERANCE
  --ks-confidence KS_CONFIDENCE
  --nterm NTERM
  --cterm CTERM
  --eligibility {control-support,detectable}
  --treatment-dispersion {replicates,control}
                        replicates: fit from treatment replicate disagreement
                        (needs >=2); control: assume control dispersion
                        explicitly
  --dispersion-model {linear-quadratic,quadratic}
                        linear-quadratic: Var=mu*(1+linear)+phi*mu^2 per
                        condition; quadratic: linear term fixed at 0
  --include-stop
```

### multicondition-project

Multi-condition tests and their plots in one step for a project (reference from the project config). Stage: Condition comparisons.

```text
usage: ytab_multicondition_project.py [-h] --project-config PROJECT_CONFIG
                                      [--normalized NORMALIZED]
                                      [--normalization NORMALIZATION]
                                      [--sample SAMPLE]
                                      [--reference-condition REFERENCE_CONDITION]
                                      [--adopt ADOPT]
                                      [--adopt-plots ADOPT_PLOTS]
                                      [--test {both,anova,all-pairs}]
                                      [--stratify-libraries]
                                      [--independent-libraries] [--fdr FDR]
                                      --output OUTPUT

Multi-condition tests and their plots in one step (used by the app's Condition Comparisons tab).

Runs `ytab multicondition` (permutation ANOVA and all-pairs resampling) on a normalization run, then
`ytab multicondition-plots` (response correlation and heatmap). The genome FASTA and GFF come from the project.
A plotting failure (e.g. fewer than two non-reference conditions for a correlation) is reported, not fatal.

  ytab multicondition-project --project-config output/projects/P/config/project.yaml \
      --normalized output/projects/P/normalization/RUN/experiment/methods --normalization betageom \
      --sample P::ctrl1=control --sample P::flu1=FLU ... --reference-condition control --output OUT

--adopt TESTS [--adopt-plots PLOTS] turns the outputs of an earlier `ytab multicondition` (and
`ytab multicondition-plots`) run into a run of this command (folders linked, not copied; manifest as for a new run),
so results made before this command existed show in the app without being recomputed.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --normalized NORMALIZED
  --normalization NORMALIZATION
  --sample SAMPLE       OBSERVATION_ID=CONDITION; repeat
  --reference-condition REFERENCE_CONDITION
  --adopt ADOPT         Folder of an earlier `ytab multicondition` run
                        (run_manifest.json, anova/pairs tables)
  --adopt-plots ADOPT_PLOTS
                        With --adopt: its `ytab multicondition-plots` folder
  --test {both,anova,all-pairs}
  --stratify-libraries
  --independent-libraries
                        Conditions are separate, independently made libraries
                        (e.g. one per strain)
  --fdr FDR
  --output OUTPUT
```

### multicondition

Hermes multi-condition tests (permutation ANOVA and all-pairs resampling) on a normalization output. Stage: Condition comparisons.

```text
usage: ytab_multicondition.py [-h] --normalized NORMALIZED
                              [--normalization {ttr,totreads,nzmean,betageom,aBGC}]
                              --gff GFF --reference REFERENCE --output OUTPUT
                              --sample SAMPLE
                              [--reference-condition REFERENCE_CONDITION]
                              [--test {both,anova,all-pairs}]
                              [--stratify-libraries] [--independent-libraries]
                              [--recenter {iterative-median,trimmed-sum,none}]
                              [--min-sites MIN_SITES] [--seed SEED]
                              [--fdr FDR] [--pseudocount PSEUDOCOUNT]
                              [--nterm NTERM] [--cterm CTERM]
                              [--asymmetry-warning ASYMMETRY_WARNING]
                              [--include-stop]

Hermes multi-condition tests: site-restricted permutation ANOVA and all-pairs resampling.

For one library (pool) split into several conditions (e.g. control, FLU, Zn, other drugs), or several
libraries each split into the same conditions (--stratify-libraries). Input is a YTAB `normalize` output;
BGC (betageom) is recommended. Assign every sample to a condition with --sample OBSERVATION_ID=CONDITION.

  anova      does a gene respond to any condition?  (q <= FDR -> "varies"; top_condition = largest |log2FC|)
  all-pairs  which conditions differ, and in which direction? (resampling v3 for every pair, one FDR family)

If the ANOVA replicate-noise warning is printed, use the all-pairs results: ANOVA p-values are
anti-conservative when one condition's replicates are much noisier (pre-registered evaluation:
output/essentiality/multicondition/benchmark_v2/RESULT.md).

Example:
  ytab multicondition --normalized output/normalization/RUN --normalization betageom \
      --gff genes.gff --reference genome.fna --output output/multicondition/RUN \
      --sample P::ctrl1=control --sample P::ctrl2=control --sample P::flu1=FLU --sample P::flu2=FLU \
      --sample P::zn1=Zn --sample P::zn2=Zn --reference-condition control

options:
  -h, --help            show this help message and exit
  --normalized NORMALIZED
                        YTAB normalize output directory, or one method folder
                        inside it
  --normalization {ttr,totreads,nzmean,betageom,aBGC}
                        Method folder when --normalized is a run directory
                        (betageom recommended)
  --gff GFF
  --reference REFERENCE
  --output OUTPUT
  --sample SAMPLE       OBSERVATION_ID=CONDITION; repeat for every sample
  --reference-condition REFERENCE_CONDITION
                        Condition used for recentering and log2FC (default:
                        the first condition listed)
  --test {both,anova,all-pairs}
  --stratify-libraries  Several libraries (pools), each split into the
                        conditions: permute only within library_id
  --independent-libraries
                        Conditions are separate, independently made libraries
                        (e.g. one per strain): permute among all samples at
                        each site, as TRANSIT
  --recenter {iterative-median,trimmed-sum,none}
                        Per-condition rescaling to the reference (iterative-
                        median is the validated default)
  --min-sites MIN_SITES
  --seed SEED
  --fdr FDR
  --pseudocount PSEUDOCOUNT
  --nterm NTERM
  --cterm CTERM
  --asymmetry-warning ASYMMETRY_WARNING
  --include-stop
```

### multicondition-plots

Response correlation (corrplot) and response heatmap for a multicondition run. Stage: Condition comparisons.

```text
usage: ytab_multicondition_plots.py [-h] --run RUN --output OUTPUT
                                    [--min-sites MIN_SITES] [--top TOP]
                                    [--pseudocount PSEUDOCOUNT]
                                    [--bootstrap BOOTSTRAP] [--seed SEED]

Response correlation and response heatmap for a `ytab multicondition` run, adapted from TRANSIT's corrplot and heatmap
(DeJesus et al. 2015).

  response_correlation.png/.tsv  Spearman correlation between conditions' responses (log2FC vs the reference),
                                 split-reference (unbiased) when every library has >= 2 reference replicates,
                                 with gene-bootstrap 95% CIs; on all tested genes and (table only) responsive genes
  heatmap.png / heatmap_matrix.tsv  log2FC vs reference for responsive genes (ANOVA, or all-pairs when ANOVA's
                                 noise warning fired), zero-centred colours, genes clustered on their profiles

Example:
  ytab multicondition-plots --run output/multicondition/RUN --output output/multicondition/RUN_plots

options:
  -h, --help            show this help message and exit
  --run RUN             Output directory of `ytab multicondition`
  --output OUTPUT       New output directory
  --min-sites MIN_SITES
                        Heatmap: minimum eligible sites per gene (default 3)
  --top TOP             Heatmap: only the top N responsive genes by ANOVA q
  --pseudocount PSEUDOCOUNT
  --bootstrap BOOTSTRAP
                        Gene-bootstrap resamples for correlation CIs (0 =
                        none)
  --seed SEED
```

### genetic-interactions

Hermes genetic interactions: strain x treatment response differences. Stage: Condition comparisons.

```text
usage: ytab_genetic_interactions.py [-h] --normalized NORMALIZED
                                    [--normalization {ttr,totreads,nzmean,betageom,aBGC}]
                                    --gff GFF --reference REFERENCE
                                    --output OUTPUT --sample SAMPLE
                                    --reference-strain REFERENCE_STRAIN
                                    --query-strain QUERY_STRAIN
                                    [--pseudocount PSEUDOCOUNT] [--fdr FDR]
                                    [--min-delta MIN_DELTA]
                                    [--clone-warning CLONE_WARNING]
                                    [--nterm NTERM] [--cterm CTERM]
                                    [--min-sites MIN_SITES] [--empirical-null]
                                    [--rope ROPE] [--correction {bh,holm}]
                                    [--confidence CONFIDENCE] [--include-stop]

Hermes genetic interactions, adapted from TRANSIT's GI method (DeJesus et al. 2017): does a gene's treatment
response differ between two strains?

delta = log2FC(query) - log2FC(reference); interaction classes Aggravating / Alleviating / Suppressive.
Each strain needs >= 2 control/treated pairs: replicate cultures split from one pool (paired in order),
or several pools each split into control and treated. Input is a YTAB `normalize` output (BGC recommended)
whose catalogue gives each sample's library (pool). Validated in a pre-registered simulation
(output/essentiality/gi/benchmark/RESULT.md).

Example:
  ytab genetic-interactions --normalized output/normalization/RUN --normalization betageom \
      --gff genes.gff --reference genome.fna --output output/gi/RUN \
      --reference-strain WT --query-strain mutA \
      --sample P::wt_c1=WT:control --sample P::wt_t1=WT:treated --sample P::wt_c2=WT:control --sample P::wt_t2=WT:treated \
      --sample P::mu_c1=mutA:control --sample P::mu_t1=mutA:treated --sample P::mu_c2=mutA:control --sample P::mu_t2=mutA:treated

options:
  -h, --help            show this help message and exit
  --normalized NORMALIZED
  --normalization {ttr,totreads,nzmean,betageom,aBGC}
  --gff GFF
  --reference REFERENCE
  --output OUTPUT
  --sample SAMPLE       OBSERVATION_ID=STRAIN:control|treated; repeat
  --reference-strain REFERENCE_STRAIN
                        Strain A (e.g. wild type)
  --query-strain QUERY_STRAIN
                        Strain B (e.g. the mutant)
  --pseudocount PSEUDOCOUNT
  --fdr FDR
  --min-delta MIN_DELTA
  --clone-warning CLONE_WARNING
  --nterm NTERM
  --cterm CTERM
  --min-sites MIN_SITES
  --empirical-null      Efron empirical null (inflate-only); off by default
                        (not used in validation)
  --rope ROPE           Call only interactions with |delta| demonstrably >
                        ROPE (TREAT test, McCarthy & Smyth 2009; 0.5 is a
                        common choice). 0 = test delta != 0
  --correction {bh,holm}
                        bh: false discovery rate (default, validated); holm:
                        family-wise error rate (NOT validated: 8.3% family-
                        wise error in simulation)
  --confidence CONFIDENCE
                        Level of the delta interval (delta_lower/upper)
  --include-stop
```

### browser-tracks

Genome browser files for a project: reference with display chromosome names, gene BED, insertion tracks and an IGV session. Stage: Gene Explorer.

```text
usage: ytab_browser_tracks.py [-h] --project-config PROJECT_CONFIG [--force]

Genome browser files for a project (Gene Explorer > Genome browser reads only these): reference with display
chromosome names, aliases, gene annotation BED and one insertion track per sample with a CreateHitFile hit table.

  ytab browser-tracks --project-config output/projects/P/config/project.yaml [--force]

Writes <project>/genome_browser/ (see ytab.domain_explorer.browser_tracks); rebuilds only when an input is newer.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --force
```

### gene-insertions

Gene/domain insertion plots and queries. Stage: Gene Explorer. Takes `--project-config`.

```text
usage: ytab_draw_gene_domain_insertions.py [-h]
                                           --project-config PROJECT_CONFIG
                                           [--gene GENE] [--query QUERY]
                                           [--list-tracks]
                                           [--track-source {raw,combined,combined_parent,all}]
                                           [--track-preset {all,parents,treated,matched_pairs,pool1_pair,pool2_pair,pool3_pair,pool4_pair,custom}]
                                           [--samples SAMPLES]
                                           [--flank-bp FLANK_BP]
                                           [--output OUTPUT]
                                           [--width-px WIDTH_PX] [--dpi DPI]
                                           [--label-mode {full,compact}]
                                           [--hide-site-counts]
                                           [--hide-domains] [--hide-direction]
                                           [--force] [--json]

Draw/query gene-domain insertion figures from existing YTAB project outputs.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --gene GENE           Gene identifier/name to resolve and plot
  --query QUERY         Query only; print candidate genes
  --list-tracks         List available insertion tracks
  --track-source {raw,combined,combined_parent,all}
  --track-preset {all,parents,treated,matched_pairs,pool1_pair,pool2_pair,pool3_pair,pool4_pair,custom}
                        Preset track selection; use custom with --samples for
                        explicit tracks
  --samples SAMPLES     Comma-separated sample/track names, or all
  --flank-bp FLANK_BP
  --output OUTPUT
  --width-px WIDTH_PX
  --dpi DPI
  --label-mode {full,compact}
  --hide-site-counts
  --hide-domains
  --hide-direction
  --force
  --json                Print machine-readable JSON
```

### go

GO enrichment and optional semantic reduction. Stage: Gene Explorer. Takes `--project-config`.

```text
usage: ytab_run_go_analysis.py [-h] --project-config PROJECT_CONFIG
                               [--analysis-id ANALYSIS_ID] [--genes GENES]
                               [--gene-file GENE_FILE]
                               [--ontologies ONTOLOGIES]
                               [--node-size NODE_SIZE] [--p-cutoff P_CUTOFF]
                               [--top-terms TOP_TERMS] [--semantic-reduction]
                               [--rrvgo-threshold RRVGO_THRESHOLD]
                               [--source-label SOURCE_LABEL] [--force]
                               [--job-id JOB_ID]
                               [--progress-file PROGRESS_FILE]

Run YTAB GO enrichment through topGO and optional rrvgo.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --analysis-id ANALYSIS_ID
  --genes GENES
  --gene-file GENE_FILE
  --ontologies ONTOLOGIES
  --node-size NODE_SIZE
  --p-cutoff P_CUTOFF
  --top-terms TOP_TERMS
  --semantic-reduction
  --rrvgo-threshold RRVGO_THRESHOLD
  --source-label SOURCE_LABEL
  --force
  --job-id JOB_ID
  --progress-file PROGRESS_FILE
```

### coverage

Analytical Poisson / equal-gene screen benchmark. Stage: Screen planning and mutants.

```text
usage: ytab_poisson_coverage.py [-h] --genes GENES [--target TARGET]
                                [--output OUTPUT] [--no-plot] [--force]
                                [--json]

Analytical equal-gene Poisson benchmark: expected recovery, not confidence
that every gene is recovered.

options:
  -h, --help       show this help message and exit
  --genes GENES    Number of genes in the selected benchmark universe N
  --target TARGET  Expected recovery percentage, strictly between 0 and 100
                   (default: 95)
  --output OUTPUT
  --no-plot
  --force
  --json           Print machine-readable result summary
```

### mutants

Universal paired Sanger/ABI mutant identification. Stage: Screen planning and mutants. Takes `--reference-config`.

```text
usage: ytab_identify_mutants.py [-h] --reference-config PROJECT_CONFIG
                                --input INPUT [INPUT ...] [--output OUTPUT]
                                [--hermes HERMES] [--threads THREADS]
                                [--min-length MIN_LENGTH]
                                [--min-mapq MIN_MAPQ] [--quality QUALITY]
                                [--overlap OVERLAP] [--window WINDOW]
                                [--junction-min-overlap JUNCTION_MIN_OVERLAP]
                                [--junction-gene-window JUNCTION_GENE_WINDOW]
                                [--junction-gene-min-bases JUNCTION_GENE_MIN_BASES]
                                [--min-identity MIN_IDENTITY]
                                [--error-rate ERROR_RATE]
                                [--fitness-result FITNESS_RESULT]
                                [--link-project-config LINK_PROJECT_CONFIG]
                                [--replace-current]
                                [--junction-read-end {start,end,require}]
                                [--condition CONDITION] [--pool POOL]

Identify individual Hermes mutants from local Sanger files or archives.

options:
  -h, --help            show this help message and exit
  --reference-config, --project-config PROJECT_CONFIG
                        Configuration selecting the alignment reference; does
                        not link screen results
  --input INPUT [INPUT ...]
  --output OUTPUT       Default: shared current identification loaded by
                        Shiny; specify another directory for standalone
                        results
  --hermes HERMES
  --threads THREADS
  --min-length MIN_LENGTH
  --min-mapq MIN_MAPQ
  --quality QUALITY
  --overlap OVERLAP
  --window WINDOW
  --junction-min-overlap JUNCTION_MIN_OVERLAP
  --junction-gene-window JUNCTION_GENE_WINDOW
  --junction-gene-min-bases JUNCTION_GENE_MIN_BASES
  --min-identity MIN_IDENTITY
  --error-rate ERROR_RATE
  --fitness-result FITNESS_RESULT
  --link-project-config LINK_PROJECT_CONFIG
  --replace-current
  --junction-read-end {start,end,require}
  --condition CONDITION
  --pool POOL
```

### datasets

Catalogue of the reprocessed public datasets (study, strain, transposon, accession, QC, results in the app) for the Datasets page. Stage: Datasets, results and sharing.

```text
usage: ytab_datasets.py [-h] [--output OUTPUT] [--designs DESIGNS]
                        [--find TEXT [TEXT ...]]

Catalogue of the public Tn-seq datasets YTAB has reprocessed, for the app's Datasets page and its search.

Each study is a row of resources/datasets/studies.tsv (citation, DOI, strain, transposon, accession, design); its runs
are the rows of pipeline/hpc/datasets/<dataset>/samples.csv whose `study` column matches. For every study and project
the catalogue adds what the runs are (roles, strains, conditions), their QC (newest `ytab qc` table) and which results
the app can open without computing anything (QC, essentiality, normalization, condition comparisons, genome browser,
library correlation, cross-species).

  ytab datasets                      # writes resources/datasets/catalog.tsv, libraries.tsv and docs/experimental_designs.md
  ytab datasets --find SRR29743984   # which study, library and analyses a run belongs to (also: a strain, drug or condition)

libraries.tsv gives each run its library (runs merged by `ytab merge-runs` share one) and every analysis that uses it, with
its role there (control, treated, reference condition, ...); experimental_designs.md lists the same per study, with the
commands that reproduce each analysis.

A new dataset appears once it has a dataset folder and a row in studies.tsv; nothing here names a study.

options:
  -h, --help            show this help message and exit
  --output OUTPUT       Folder for catalog.tsv and libraries.tsv
  --designs DESIGNS     Experimental-design page (searchable in the app)
  --find TEXT [TEXT ...]
                        Print the studies and runs matching TEXT (run ID,
                        strain, drug, condition) and the analyses using them;
                        writes nothing
```

### report

Local HTML report of existing project outputs. Stage: Datasets, results and sharing. Takes `--project-config`.

```text
usage: ytab_build_project_report.py [-h] --project-config PROJECT_CONFIG
                                    [--output OUTPUT] [--force]
                                    [--open-browser]

Build a lightweight HTML summary from existing YTAB outputs.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --output OUTPUT
  --force
  --open-browser
```

### export

Portable project-result export. Stage: Datasets, results and sharing. Takes `--project-config`.

```text
usage: ytab_export_project.py [-h] --project-config PROJECT_CONFIG
                              [--output OUTPUT] [--include-logs]
                              [--include-browser-tracks] [--include-hit-files]
                              [--include-bams] [--force] [--dry-run]

Create a compact reproducible YTAB project export.

options:
  -h, --help            show this help message and exit
  --project-config PROJECT_CONFIG
  --output OUTPUT
  --include-logs
  --include-browser-tracks
  --include-hit-files
  --include-bams
  --force
  --dry-run
```

### resource-bundle

Bundle reprocessed projects (results, no FASTQ/BAM) with a checksum manifest, to open them in the app elsewhere. Stage: Datasets, results and sharing.

```text
usage: ytab_resource_bundle.py [-h] --project PROJECT [--output OUTPUT]
                               [--list]

Bundle reprocessed projects so others can open them in the YTAB app without rerunning anything.

A bundle holds each project folder (configuration, hit files, QC, normalization, condition comparisons, essentiality,
gene explorer and screen-planning results; no FASTQ or BAM files), the cross-species runs that use only these
projects, and MANIFEST.tsv with the size and SHA-256 of every file. Unpacked at the repository root, the projects
appear in the app's project list.

  python scripts/local/ytab_resource_bundle.py --project Gale_2020 --project Segal_2018_albicans ... --output bundle.tar.gz
  python scripts/local/ytab_resource_bundle.py ... --list     # sizes only, writes nothing

options:
  -h, --help         show this help message and exit
  --project PROJECT  Project ID under output/projects (repeat)
  --output OUTPUT
  --list
```

### app

Open the Shiny app on existing results. Stage: Datasets, results and sharing. Takes `--project-config`.

```text
Usage: ytab app [--project-config FILE] [--host HOST] [--port PORT] [--no-browser] [--skip-env-check]
```

### tvc-benchmark

Pre-registered treated-vs-control benchmark: Gale-calibrated simulation, method x normalization arms, paired tests. Stage: Pre-registered benchmarks.

```text
usage: ytab_tvc_benchmark.py [-h] --output OUTPUT [--experiments EXPERIMENTS]
                             [--null-experiments NULL_EXPERIMENTS]
                             [--workers WORKERS] [--dev-seed DEV_SEED]
                             {calibrate,dev,blind,amend,report}

Treated-vs-control benchmark (docs/manuscript/preregistration_treated_vs_control.md).

  ytab tvc-benchmark calibrate --output DIR                calibration from the base project (default: Gale 2020 FLU)
  ytab tvc-benchmark dev --output DIR [--experiments 2]     development run on seeds 7001-7100 (never the benchmark seed)
  ytab tvc-benchmark blind --output DIR                     the pre-registered run; refuses unless the pre-registration is frozen
  ytab tvc-benchmark amend --output DIR                     amendment 1: rerun the amended z-score arms on the same experiments
  ytab tvc-benchmark report --output DIR                    metrics, paired tests, calibration gates -> report.md

The base project needs one Site Normalization run with two control and two treated libraries in its experiment matrix.

positional arguments:
  {calibrate,dev,blind,amend,report}

options:
  -h, --help            show this help message and exit
  --output OUTPUT
  --experiments EXPERIMENTS
                        per signal scenario
  --null-experiments NULL_EXPERIMENTS
  --workers WORKERS
  --dev-seed DEV_SEED
```

### tvc-phase2

Treated-vs-control benchmark phase 2: every arm on real Gale FLU, a swap null, H2O2/Zn, and Gale's published scores. Stage: Pre-registered benchmarks.

```text
usage: ytab_tvc_phase2.py [-h] --output OUTPUT [--workers WORKERS]

Treated-vs-control benchmark, phase 2 (pre-registration, "Real data, after the simulation"; descriptive).

  ytab tvc-phase2 --output DIR [--workers 6]

Runs every pre-registered arm (ytab.benchmarks.treated_vs_control.ARMS, the code of the amended blind run) on real
libraries, each design with 2 control and 2 treated libraries:
  gale_flu        Gale 2020 Cg1: no-drug controls 1-2 vs fluconazole 1-2
  gale_swap_null  control 1 + FLU 1 vs control 2 + FLU 2 (each side has one control and one FLU culture: no calls expected)
  h2o2_<bg>, zn_<bg>   H2O2 and Zn screens per background (illustration only)
and writes:
  calls.tsv            per design, arm and gene: score, call (+1/-1/0), p
  call_summary.tsv     calls per design and arm
  overlap.tsv          Jaccard overlap of depleted and of enriched calls between arms (gale_flu)
  publication.tsv      gale_flu arms vs Gale's published fluconazole scores
  control_control.tsv  the z-score's control-control comparison (fraction beyond its own call threshold)

options:
  -h, --help         show this help message and exit
  --output OUTPUT
  --workers WORKERS
```

