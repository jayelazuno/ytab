# Command-line analyses and master workflows

An opt-in [empirical site-catalogue and TTR evaluation](ttr_normalization.md) is
available as `python scripts/local/ytab.py ttr --help`. It uses an explicit
multi-project parent manifest and does not replace existing normalization stages.
The shared [normalization choices](normalization_methods.md) interface is
`python scripts/local/ytab.py normalize --help`; method descriptions and availability
are exposed with `--list-methods` and `--describe METHOD`.
Use `normalize-experiment` for design-selected controls/treatments, then
`normalize-genes` for separate raw/TTR/CPM gene-level validation. Both are opt-in
diagnostics; see the shared normalization guide for commands and interpretation.
`normalize-project` includes `midlc` (MidLC common saturation): every library in a
normalization group is thinned by binomial sampling to the group's smallest
depth-to-MidLC ratio, never upsampled. The `normalization` module (parent-only
MidLC thinning used by the `classifier` profile) uses the same rule: its automatic
target is the common saturation of the selected parents, and MidLC is computed
exactly (half the occupied sites expected under thinning) by both it and `library-qc`.
Earlier automatic targets (e.g. `T409600`) were read counts passed where a multiple
of MidLC is expected, so those runs did not thin; such targets are now rejected.

`tnseq-stats` exports [per-library Tn-seq count statistics plus saved Hermes QC](tnseq_stats.md).
`control-gumbel` runs the [experimental control-conditioned Gumbel](control_conditioned_gumbel.md)
on a `normalize` output (ttr default, totreads or nzmean), comparing combined treatment replicates
against combined matched controls. It runs no normalization of its own and does not replace the ML
classifier or the z-score fitness analysis.

### Treated vs control

Three methods call genes depleted or enriched in treated libraries relative to their controls
(pairs from `config/comparison_design.csv`):

| Method | Command | Normalization | Status |
|---|---|---|---|
| Local-abundance z-score | `fitness` (project) or `fitness-zscore` (explicit inputs) | total reads (CPM) by default; any Site Normalization method with `--normalization METHOD --normalization-run RUN` | established |
| Resampling | `treated-vs-control --method resampling` | BGC (default, with trimmed-sum recentering), aBGC, TTR, total reads, nonzero mean | passed its pre-registered calibration gates (FDR 6%) |
| Control-conditioned Gumbel | `treated-vs-control --method gumbel` | TTR (default), total reads, nonzero mean | passed its pre-registered calibration gates (FDR 1%); depletion only |

The z-score runs in Python (`src/ytab/fitness/zscore.py`); the former R script is retired. Its
noise model (SD = m1 + m2·A^m3) is fitted to full least-squares convergence. The R version (nlsLM)
stopped slightly earlier, so a few borderline calls differ from runs made with it: none of 5,234
genes on the H2O2 screen, 111 on Gale 2020, 18 on the Zn screen; log2 fold changes agree to 1e-15.

Calls now follow Gale et al. 2020: the replicate libraries of each background × treated condition are
averaged and each gene gets one z (`--combine replicates`, the default; final calls depleted,
enriched or none). The earlier rule, per-pair calls combined as consistent / single pool / mixed,
is `--combine pairs`. Per-pair z-scores are written in both cases
(`treated_vs_parent.by_pool.log2fc_z.csv`); the combined ones are in
`treated_vs_parent.combined.log2fc_z.csv`.

With a Site Normalization method, the gene values (the method's normalized site counts summed per gene) are expressed
per million of the mean library total before the z-score. Every library is multiplied by the same constant, so the
method's own between-library correction is kept, and the pseudocount (0.5) means the same as with CPM (amendment 1 of
`docs/manuscript/preregistration_treated_vs_control.md`).

**Outputs the app reads (`src/ytab/fitness/contract.py`).** Every run writes the same tables, whatever the normalization,
library or species, and the command fails if one is missing. The app only displays them; it computes nothing itself.
- z-score (`<project>/treated_vs_parent/<analysis ID>/`):
  - `zscore_manifest.json`: normalization and run, `value_unit`, call rule, call counts, settings.
  - `tables/treated_vs_parent.summary_by_feature.csv`: one row per gene, with the call, the effect and z the call is
    based on (`call_log2FC`, `call_z`), and plot abundances (`mean_abundance`, `control_mean`, `treated_mean`).
  - Per-pair and combined statistics.
  - `gene_values.long.csv`: the values the z-score used.
  - `library_sizes.feature_reads.csv`: reads before normalization.
  - The control–control z-scores and the noise-model fit.
  - `top_hit_evidence.csv` (and `.by_pair.csv`): the top hits highlighted on the MA plot, scatter and heatmap and
    used for GO enrichment.
    - Per gene: the effect and z behind its call, the abundance, and the pairs that support or call each direction.
    - Ranks per direction, following published hit lists (Gale et al. 2020: z-score calls with at least a 2-fold
      change, reported by z). Candidates are genes called in that direction with |log2FC| ≥ 1
      (`hit_min_log2FC`), ordered by |z|, then |log2FC|.
    - There's no abundance floor by default (`hit_min_abundance` = 0): the noise model already gives low-count genes
      a lower z.
    - The app only applies the viewer's direction, minimum supporting pools and number of hits.
- Resampling and Gumbel (`<output>/`): `manifest.json`, `genes.tsv` (all groups, with `mean_abundance` for resampling)
  and `<group>/sites.tsv.gz`.

Gene names and library sizes come from the raw SummaryTables when they exist. Otherwise they come from the GFF and the
Site Normalization run's raw site counts, so a run needs only the normalization run and the reference.

```sh
python scripts/local/ytab.py treated-vs-control \
  --project-config output/projects/P/config/project.yaml --method resampling \
  --normalization-run output/projects/P/normalization/RUN --normalization betageom \
  --output output/projects/P/treated_vs_control/resampling_RUN
```

## Genome browser files

```sh
python scripts/local/ytab.py browser-tracks --project-config output/projects/P/config/project.yaml
```

Writes `<project>/genome_browser/`:
- the reference with display chromosome names (and its index) and chromosome aliases;
- `genes.bed.gz`, labelled by common name, then gene ID (S. cerevisiae names for C. glabrata, the species' own otherwise);
- one insertion-site BED per sample with a CreateHitFile hit table;
- `igv_session.xml`, which opens all of them in IGV desktop.

The app's Gene Explorer > Genome browser runs this command when it opens, and only reads its files. "Open in IGV" sends
the same files to a running IGV desktop (port 60151); "Download for IGV" zips them with the session.

## Sample identity and display names

Every project uses the same sample metadata contract. `sample` in its configured
sample sheet is the stable internal identifier used for files, joins and analysis
selections. Optional `display_name` supplies human-readable names; without one,
the app shows the identifier. Study-specific names and experimental roles belong
in metadata, never in project-specific branches of application code.

After a run, use **Project & Data → Samples → Edit sample names** to change display
names without renaming files or rerunning analysis. The editor changes only
`display_name`, rejects duplicate names and stale edits, and is unavailable while
an app job is running. Blank names fall back to sample IDs. Saved analysis tables
and previously generated report/plot files retain their original identifiers;
the app's metadata-aware displays use the updated names.

YTAB supports individual analysis modules, a sequential project pipeline, and an
extensible ordered workflow. These commands delegate to the existing backends used
by Shiny; they do not introduce different scientific methods. Modules are named in an
independent registry, `src/ytab/cli/registry.py`.

Activate the existing YTAB environment and run from the repository root:

```sh
conda activate ytab-local
python scripts/local/ytab.py --help
python scripts/local/ytab.py modules
python scripts/local/ytab.py fitness --help
```

The equivalent module entry point is `PYTHONPATH=src python -m ytab`. This is a
source-checkout interface, not a claim that an installed `ytab` executable exists.
Existing `scripts/local/ytab_*.py` commands remain supported. CLI analyses run
headlessly; opening the browser requires the explicit `app` command.

HPC is optional, not required for the desktop app. The repository provides
[personal Argon/SGE stage jobs and master submitter](../pipeline/hpc/README.md), which use this
same CLI and existing project/output routes. Other clusters can supply their own
scheduler wrappers. The Argon master supports the existing core/classifier/fitness/all
profiles, adds per-run failure gates to scheduler holds, and preserves the same
scientific methods. See its guide for new projects, nested FASTQs and metadata.

New project configurations use a checkout-relative repository root and standard
downstream routes. Sample FASTQ filenames are relative to the explicitly configured
`fastq_dir` (`fastq_paths_relative_to: fastq_dir`), so only that base changes when
the same reads move between desktop and HPC. The CLI and app resolve these paths
consistently. Legacy configurations without this marker retain their original
path semantics. Moving inputs does not infer or change experimental sample roles.

## Individual modules

Use `modules --json` for the machine-readable registry and `<module> --help` for
the backend's exact options. Available modules include reference preparation,
mapping, hit files, summaries, library QC, site normalization (MidLC included), combined
parents, essentiality classification, fitness, GO enrichment, gene/domain insertion
plots, Poisson coverage, mutant identification, project status, reports and exports.
Initialization and hit-file import are also available.

Hermes mapping defaults to automatic sequence cleanup and the designated R1 read;
R2 and explicit concordant-pair requirements are supported. See
[Read alignment](hermes_alignment.md) for shared CLI/Shiny settings and reports.

**Transposons and library designs.** Trimming comes from a transposon profile
(`resources/transposons/<transposon>/profile.json`): `hermes` (`splinkerette_pcr`, the default;
`end_ligated_adapter`; `terminal_10bp`) and `ac_ds` (`ds_primer_tail`, the default; `satay_restriction`).
Select one in the project config (`mapping: {transposon: ..., library_design: ...}`), per sample with `transposon`
and `library_design` sample-sheet columns, or for one run:

```sh
python scripts/local/ytab.py mapfastq --project-config P/config/project.yaml \
  --transposon ac_ds --library-design ds_primer_tail
```

| option | meaning |
|---|---|
| `--transposon NAME` | transposon profile (default: the project's, else `hermes`) |
| `--library-design NAME` | design within the profile (default: the profile's `default_design`) |
| `--junction-mate R1\|R2`, `--require-mate-alignment`, `--hermes-cleanup` | as before; explicit settings override the design |

The backend `src/ytab/mapping/MapFastq.py` gains `--trim-3p-seq SEQ` (repeatable; all 3′ sequences trimmed in one
cutadapt pass) next to `--trim-5p-seq`, `--trim-5p-overlap` and `--discard-untrimmed`; a 5′ sequence starting with
`^` is anchored to the read start. The sequences of every design are listed in
[Public datasets](public_datasets.md#read-structures-and-trimming). Mapping manifests record the transposon,
design and profile used.

```sh
python scripts/local/ytab.py fitness \
  --project-config output/projects/H2O2_screen_v1/config/project.yaml \
  --analysis-id treated_vs_parent_raw_cpm --input-mode raw-summary

python scripts/local/ytab.py gene-insertions \
  --project-config output/projects/H2O2_screen_v1/config/project.yaml \
  --gene FMP30 --output /tmp/fmp30.png

python scripts/local/ytab.py coverage --genes 4379 --target 95 \
  --output /tmp/coverage --json

python scripts/local/ytab.py mutants \
  --reference-config output/projects/H2O2_screen_v1/config/project.yaml \
  --input data/mutant_screen/pfasta-Oct13-9-18-48.tar \
          data/mutant_screen/trace-Oct13-13-38-20.tar
```

The numbers and input paths above are examples, not defaults or recommendations.
Coverage uses the unchanged analytical Poisson helpers, with an explicitly supplied
gene universe and target **expected recovery percentage**, not probability of
recovering every gene. It exports `coverage.json`, `summary.csv`, `curve.csv` and a
headless `coverage.png` (`--no-plot` omits the image). Identical exports are reused;
different parameters at the same output require `--force` or another directory.
Retired empirical coverage models are not reintroduced into the UI or this command.

Mutant identification accepts `.seq`, `.ab1`, folders and supported archives;
matching exports are paired by the existing quality-aware backend. Reference
selection does **not** imply a project screen association. To explicitly associate
screen results, add `--link-project-config PATH`. The default output is the shared
`output/mutant_screen/current`, also loaded by Shiny. The existing cache and atomic
replacement behavior apply. `--output PATH` selects standalone output instead.

## Standard master run

```sh
python scripts/local/ytab.py run \
  --project-config output/projects/H2O2_screen_v1/config/project.yaml \
  --profile fitness --threads 4 --dry-run
```

Remove `--dry-run` to execute. Profiles retain the existing stage definitions:

- `core`: reference, mapping, hit files, raw summaries and library diagnostics.
- `classifier`: core plus parent MidLC thinning (common saturation), normalized/combined
  summaries and the essentiality classifier, for any species (the app builds hybrid inputs
  with `hybrid-inputs`, which takes any normalization).
- `fitness`: core plus comparison design and treated-versus-control fitness.
- `all`: all existing project-pipeline stages above, **not** every optional analysis.

`--stages`, `--from-stage` and `--to-stage` select existing project-pipeline stages.
The original pipeline retains its per-stage cache, project status and logs. Opt-in
site-level normalization methods are separate modules, not automatically inserted
into these existing profiles.

## Essentiality HMM for any transposon

`ytab hermes-hmm OUTPUT --project P --libraries S1 S2 --sequence-covariate learned` learns the target-site preference
from the libraries' own insertion sites (needs YTAB hit files, which carry the strand) instead of counting the Hermes
nTnnnnAn motif (`motif`, the default). Use it for Ac/Ds and other transposons; on Hermes it matches the motif model
(see the HMM methods page, section 8.4). It also writes `target_model.tsv`.

`ytab essentiality` also writes plot summaries for the app: `plots/target_model.tsv` (the target preference the fit
used, or, for a motif fit, one learned from the same libraries as a diagnostic) and, with `--labels`, `plots/roc.tsv`
(each model against the curated labels per tier; the hybrid's training genes out-of-fold).
`ytab essentiality --project-config P --refresh-plots RUN` adds both to an earlier run.

## Essentiality across species

`ytab cross-species --project-config A --project-config B [...] --output DIR` compares the newest essentiality runs of
two or more projects (one per species) on one-to-one orthologs, matched through S. cerevisiae systematic names with
each project's `reference.orthology_file` (an S. cerevisiae project is the hub). Writes `orthologs.tsv` (each
species' gene, score and call), `intersections.tsv` (shared and species-specific essential sets) and `pairwise.tsv`
(Cohen's kappa of the calls, Spearman correlation of the scores). A new species needs one ortholog table: tab-separated,
the species' gene first and a `scer_id` column, or a S. cerevisiae systematic name in a later column. The app shows it
in Gene Explorer > Across Species > Essentiality across species of every project it includes (runs are kept in
`output/cross_species/`).

## Species resources

Every species folder has the same layout: `reference_genome/` (FASTA, GFF, Bowtie2 index, chromosome names,
centromere and tRNA tracks), `annotation/`, `orthologs/`, `essentiality_labels/` (with `raw/` sources) and
`tnseq_data/` (published insertion data); see `resources/species/README.md`. Published libraries are reprocessed
from raw reads with the dataset manifests in `pipeline/hpc/datasets/` ([Public datasets](public_datasets.md)).

`scripts/local/ytab_build_species_resources.py` builds, for C. albicans and S. cerevisiae,
the files glabrata already has: centromere and tRNA tracks, chromosome display names
(`reference_genome/chromosome_names.tsv`), a gene annotation table
(`annotation/gene_annotation.tsv`: name, description, S. cerevisiae ortholog and its
viability) and curated essentiality labels (`essentiality_labels/`). The S. cerevisiae
S288C reference (FASTA and GFF, release R64-5-1) comes from the SGD release archive;
C. albicans uses the CGD A22 haplotype-A reference already in the repository. Label rules
are in the script's docstring; transposon-screen calls are kept as a comparator, never as
labels. The app uses these files for any project of that species.

Any further species needs only `resources/species/<species>/reference_genome/` with a
FASTA and a GFF3: SummaryTable takes genes from the GFF (`--feature-db gff`, NCBI, SGD,
CGD or Ensembl style; one gene per protein-coding locus, sequences absent from the
mapped FASTA left out), and the classifier (trained on S. cerevisiae) scores any
species' feature tables, attaching S. cerevisiae orthologs from
`annotation/gene_annotation.tsv` when present (`Classifier.py --species <species>`).
For gene names in the app, `ytab_build_species_resources.py --species <species> --annotation-only` writes
`annotation/gene_annotation.tsv` from the GFF (names, descriptions) and the species' ortholog table to S. cerevisiae
(`orthologs/`); `ytab_comparative_add_species.py --species <species> --label "<Genus species>"` then adds the
species to Gene Explorer > Across Species through those orthologs.

## Custom extensible master workflow

Copy and edit [analysis_workflow.yaml](examples/analysis_workflow.yaml), then run:

```sh
python scripts/local/ytab.py run --workflow analysis_workflow.yaml --dry-run
python scripts/local/ytab.py run --workflow analysis_workflow.yaml --threads 4
```

Workflow format:

```yaml
version: 1
project_config: output/projects/my_screen/config/project.yaml
steps:
  - id: core
    module: pipeline
    args: ["--profile", "core"]
  - id: fitness
    module: fitness
    args: ["--analysis-id", "treated_vs_parent_raw_cpm", "--input-mode", "raw-summary"]
  - id: benchmark
    module: coverage
    depends_on: []
    args: ["--genes", "4379", "--target", "95", "--output", "results/coverage"]
```

Relative argument paths and YAML `project_config` are relative to the **workflow
file's directory**. A CLI `--project-config` override is resolved relative to the
calling directory. Explicit per-step project/reference options take precedence.
`args` must be a list of strings: commands are passed without a shell. Step IDs must
be unique safe identifiers. By default each step depends on the previous step;
`depends_on: []` marks an independent analysis, and explicit dependencies must refer
to earlier steps. Execution is sequential, never a hidden parallel scientific run.

The master stops after a failure unless `--keep-going` is given. With that flag,
dependent steps remain blocked while independent steps can continue. Child failure
produces a nonzero master exit code; Ctrl-C records cancellation. `--force` and
`--threads` are forwarded only to compatible registered modules. Backend output
options remain authoritative: master `--output` selects only workflow logs/manifests.
Dry-run prints the plan without executing modules or writing outputs.

Every executed master run records step commands, working directories, exit codes,
logs and status in `output/cli_workflows/<run>/workflow.json` by default. Scientific
results remain in their modules' normal output locations; repeated analyses use
their existing caches rather than a second implementation.

## Viewing results and adding analyses

Project analyses produce the canonical files that the app already discovers.
Mutant identification's default current result is shared with the app. Standalone
Poisson exports and gene insertion images can be inspected without Shiny; the
Poisson app uses the same benchmark inputs/math, but does not import the CLI exports.

```sh
python scripts/local/ytab.py app \
  --project-config output/projects/H2O2_screen_v1/config/project.yaml
```

For a new analysis, provide a thin CLI runner around the shared backend, register its
project/threads/force capabilities, test it, and add its module to a workflow as
needed. Do not add it automatically to an established scientific profile without
an explicit dependency/design decision. This keeps both modular and master runs
extensible while preserving backward compatibility.
