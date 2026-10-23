# Reusable YTAB jobs for Argon

YTAB is primarily a desktop Shiny application. Its analysis modules and master
workflow also have a scheduler-independent [CLI](../../docs/command_line.md).
**No HPC account or scheduler is required to use YTAB.** The standard pipeline
can optionally run on Argon, either stage by stage or through the master submitter.

The scripts in this directory are personal University of Iowa
**Argon / SGE (`qsub`)** jobs,
not a universal HPC interface. Queue `UI`, parallel environment `smp`, email,
Conda initialization and scratch paths are personal site settings. Users of
Slurm, PBS or other clusters should wrap the ordinary CLI in their own scripts.
Public datasets are downloaded by `ytab_download_public.txt` (below); downloading is not an analysis stage and
is not called by the analysis jobs.

## Public datasets: download, then submit

Each folder in `datasets/` describes one public dataset: `samples.csv` (one row per sequencing run, keyed by its
accession, with role, pool, transposon and library design) and `dataset.env` (project ID and species).

| dataset | species | runs | transposon / library design |
|---|---|---|---|
| `segal_2018_albicans` | albicans | 11 | `ac_ds` / `ds_primer_tail` |
| `pombe_hermes` | pombe | 26 (25 included) | `hermes` / `terminal_10bp` (Guo 2013), `end_ligated_adapter` (Grech 2019, Levin 2020) |
| `michel_2017_satay` | cerevisiae | 7 (WT ×2, rapamycin, 4 mutants) | `ac_ds` / `satay_restriction` |
| `gale_2020_glabrata` | glabrata | 7 | `hermes` / `splinkerette_pcr` (re-download only; the project exists) |

```bash
cd /Users/jayelazuno/workspace/ytab
# 1. One download job per dataset (NCBI SRA by default, 3 runs at once; resumable; SOURCE=ena for ENA):
bash pipeline/hpc/ytab_public_dataset.sh segal_2018_albicans download
bash pipeline/hpc/ytab_public_dataset.sh pombe_hermes download
bash pipeline/hpc/ytab_public_dataset.sh michel_2017_satay download
bash pipeline/hpc/ytab_public_dataset.sh pombe_hermes status       # files verified so far
# 2. When a download has finished, create its project and submit the stage chain
#    (reference, mapping, hit files, SummaryTable, library diagnostics, report, export):
bash pipeline/hpc/ytab_public_dataset.sh segal_2018_albicans submit --profile core --threads 8
```

Runs come from NCBI SRA (prefetch, vdb-validate, fasterq-dump, pigz); each is accepted only if its archive validates and its read count equals ENA's, and a run that fails there is fetched from ENA with MD5 checks. ENA alone (`SOURCE=ena`) is slow from Argon (about 1–2 MB/s per file). `PARALLEL` sets runs (SRA) or files (ENA) processed at once. FASTQs go to `/nfsscratch/jayelazuno/workspace/public_data/<dataset>/fastq` (`DATA_ROOT` overrides the base); SRA archives are kept in `<dataset>/sra`. The
export stage bundles each project's summaries, hit files and logs under `output/exports/<project>` for pulling back.
Transposon profiles are in `resources/transposons/<transposon>/profile.json`; see `docs/hermes_alignment.md`.

## One-command master submission

Run this launcher with **bash on the login node**, not with `qsub`:

```bash
cd /Users/jayelazuno/workspace/ytab

# New project: the same standard routes are derived from the project ID.
bash pipeline/hpc/submit_pipeline.sh \
  --project-id My_screen --species glabrata \
  --fastq-dir /nfsscratch/jayelazuno/workspace/my_reads \
  --profile core --threads 8

# Existing Gale project: preserve its configuration and sample metadata.
bash pipeline/hpc/submit_pipeline.sh \
  --project-config output/projects/Gale_2020/config/project.yaml \
  --profile core --threads 8

# Full existing pipeline, including classifier and fitness branches.
# Review parent/control roles, independent pools and comparisons first.
bash pipeline/hpc/submit_pipeline.sh \
  --project-config output/projects/My_screen/config/project.yaml \
  --profile all --threads 8 --stage-threads mapfastq=16
```

The launcher activates `tnseq`, initializes a new project if requested, and submits
one job per stage with `qsub -terse` and a serial `-hold_jid` chain. It returns after
submission; **the jobs keep running without the terminal remaining open**. Each
stage uses the existing scientific CLI and canonical project outputs; no new
normalization or classifier logic is introduced.

Profiles follow the local master:

- `core` (default): reference preparation, alignment, hit files, raw summaries and QC.
- `classifier`: core plus existing parent MidLC normalization, normalized/combined
  summaries and the existing random-forest essentiality classifier, unchanged.
- `fitness`: core plus comparison design and raw-summary/CPM fitness.
- `all`: both existing analysis branches.

Every profile finishes with report and compact export jobs. The export explicitly
includes hit files, browser tracks and logs but not FASTQs or BAMs. It remains a
compact bundle, not a copy of every raw/normalized feature table: transfer complete
`summary/` and analysis directories as well when needed for app browsing, or
regenerate raw summaries locally from the included hit files.

`--stages mapfastq,create_hit_file,summary,library_diagnostics,report,export`
selects a subset in standard order; excluded prerequisites must already exist.
`--stage-threads STAGE=N` may be repeated. The default is eight slots for each job;
reference/alignment need at least two, other stages permit one. Override queue,
parallel environment or email with `--queue`, `--parallel-env` and `--email`.
`--junction-mate R2` and `--require-mate-alignment` /
`--no-require-mate-alignment` affect mapping only; omitted options preserve the
project's existing defaults. `--target` selects the existing classifier target;
`--analysis-id` identifies fitness results. Fitness always uses raw-summary inputs.

### FASTQ discovery and experimental metadata

New projects accept flat or nested FASTQ folders, including SRA download trees.
Nested inputs get a small symlink-only flat view in `config/argon_fastq_inputs/`;
read files are never copied or modified. Sample IDs follow FASTQ basenames through
the ordinary YTAB discovery rules. Duplicate basenames are rejected rather than
silently overwriting or merging libraries. Accession-only filenames remain
accession IDs; the submitter does **not** infer biological identities from SRR IDs
or folder names. Use one consistently named R1/R2 pair per sample under the
existing discovery convention; prepare multi-lane inputs using the normal pipeline
conventions rather than silently treating lanes as independent biological pools.

For a new project, `--sample-sheet /path/metadata.csv` can supply reviewed metadata.
It must have a `sample` column matching every discovered sample exactly once;
sequencing paths/layout come from discovery, while conditions, roles, pool IDs,
background and inclusion flags come from the supplied sheet. If `include` is
absent all samples are included. For an existing project, edit its own sample
sheet first and use `--project-config`; initialization never overwrites it.

New configurations store the repository root relative to `project.yaml`, standard
downstream routes relative to the checkout, and FASTQ filenames relative to
`fastq_dir`. To relocate the same reads for an existing FASTQ project:

```bash
bash pipeline/hpc/submit_pipeline.sh \
  --project-config output/projects/My_screen/config/project.yaml \
  --fastq-dir /new/location/of/the/same/reads --profile core
```

This changes the FASTQ base without rediscovering samples or changing roles,
comparisons, references or output routes. Legacy absolute FASTQ paths are converted
only if they are inside the declared old base; unrelated paths are rejected.
Relative subdirectory structure must be preserved; keep an existing symlink view
or point to an equivalent view. A different experiment should get a new project.
Active master jobs must finish or be cancelled before changing the base.

Imported-hit projects such as the H2O2 release retain their entry mode and need no
FASTQs. Their default master profiles skip reference indexing, mapping and hit
creation and begin with existing hit files. A FASTQ-base override cannot convert
an imported-hit project into a sequencing project.

Only the FASTQ location normally changes for core processing on this installation.
However, **parent/treated roles, independent pool identities and comparisons are
experiment-specific**. Review these before selecting `classifier`, `fitness` or
`all`. The existing runners retain their validation and will fail on invalid
designs rather than inventing biological replication. The configured reference
FASTA/annotation must be available; missing reference indexes are prepared by the
reference job. There are no new scientific defaults for accession-only experiments.

### Preview, monitoring, failure handling and resubmission

Add `--dry-run` to either master command for a read-only preview: no project is
created, no files are written and no jobs are submitted. This is a submission
preview, not a guarantee that every downstream analysis has valid biological data.

Each actual submission has a unique directory under
`output/projects/<project>/manifests/argon/<run>/`, containing `plan.json`, job IDs
and per-stage JSON receipts. Scheduler logs are under the matching
`output/projects/<project>/logs/argon/<run>/`. The master prints both locations.

```bash
qstat -u jayelazuno
tail -f output/projects/My_screen/logs/argon/RUN_ID/mapfastq.log
qdel JOB_ID
```

`hold_jid` waits for completion, not success. Therefore each worker additionally
requires successful receipts from **this submission's** earlier stages. Failures,
missing receipts after `qdel`, and configuration/sample-sheet edits block downstream
work, even if old scientific outputs happen to exist. Where SGE delivers SIGTERM,
the worker stops its child process group and records cancellation; abrupt kills
may leave a running/missing receipt, which also blocks subsequent stages.
Queued downstream jobs are not automatically deleted: they start when their hold
releases, record a blocked result and exit without scientific execution. Cancel
the remaining job IDs in `plan.json` if you do not want those brief jobs to start.

Do not edit project configuration/sample sheets while jobs are queued/running.
Only one master chain per project is allowed; a submission lock prevents
simultaneous submissions. A still-queued/running or partially submitted chain
requires checking/cancelling its jobs before using `--replace-active`. This option
does not cancel anything itself. Partial `qsub` failures record already-submitted
IDs in the plan and do not silently delete jobs. The guard applies to master
submissions, so do not run separate standalone/app jobs on the same project at the
same time. After a completed or failed chain, submit a new run using
`--project-config`; existing backend caches apply. `--force` forwards the normal
force/overwrite options (including comparison-design regeneration); use it only
when intentionally refreshing the requested stages and their downstream results.

## Individual jobs for every standard stage

`ytab_stage.txt` is a reusable job template for every stage; no copied scientific
scripts are needed. Submit, for example:

```bash
qsub -N YTAB_Summary -v PROJECT_ID=Gale_2020,STAGE=summary pipeline/hpc/ytab_stage.txt
qsub -N YTAB_QC -v PROJECT_ID=Gale_2020,STAGE=library_diagnostics pipeline/hpc/ytab_stage.txt
qsub -N YTAB_Classifier -v PROJECT_ID=My_screen,STAGE=classifier pipeline/hpc/ytab_stage.txt
qsub -N YTAB_Export -v PROJECT_ID=My_screen,STAGE=export pipeline/hpc/ytab_stage.txt
```

Supported `STAGE` values are:

```text
reference               mapfastq                create_hit_file
summary                 library_diagnostics     sample_normalization
summary_normalized      combined_hits           summary_combined
classifier              comparison_design       treated_vs_parent
report                  export
```

These standalone jobs require existing upstream outputs; they do not create a
dependency chain. Use the master for success-gated chaining. `PROJECT_CONFIG` can
replace `PROJECT_ID`. `FORCE=1`, `TARGET`, `ANALYSIS_ID`, `JUNCTION_MATE` and
`REQUIRE_MATE_ALIGNMENT=0|1` provide the corresponding settings. The master uses
explicit plan options instead of inheriting these standalone scientific overrides.
Use master `--dry-run` for previews. The older standalone mapping/hit-file wrappers
below retain their existing `DRY_RUN`, `SAMPLES` and `KEEP_GOING` options.

Optional analyses needing separate inputs (GO enrichment, individual Sanger mutant
identification, gene-specific plots and Poisson calculations) are not automatically
part of a project profile. Their ordinary modular CLI/custom workflow interfaces
remain available; this is not a new scientific workflow definition.

## Same pipeline and same outputs

The job calls `python scripts/local/ytab.py mapfastq --project-config ...`.
It neither discovers FASTQs independently nor treats R1/R2 as separate samples.
The project's existing YAML/sample sheet/reference resolve inputs and final
outputs exactly as for local CLI/Shiny runs. There is no HPC-specific sample
manifest, duplicate analysis implementation, new project route or automatic
scratch-to-project copying. Prepare the project and reference index first.

The wrapper handles only environment activation, SGE CPU allocation, options and
temporary storage. Hermes cleanup and cache/provenance checks belong to the shared
runner. The default junction mate is R1; the other mate is not required unless
project settings or an explicit override request concordant paired alignment.

## Submit from a writable directory on Argon

From the YTAB repository root:

```bash
qsub -v PROJECT_ID=H2O2_screen_v1 pipeline/hpc/ytab_mapFastq.txt

# Replace this example ID with an existing configured Gale project.
qsub -v PROJECT_ID=Gale_2020,JUNCTION_MATE=R1 pipeline/hpc/ytab_mapFastq.txt

# Preview without aligning reads; the normal runner may write preview manifests.
qsub -v PROJECT_ID=Gale_2020,DRY_RUN=1 pipeline/hpc/ytab_mapFastq.txt

# Explicitly request R2 or paired agreement.
qsub -v PROJECT_ID=Gale_2020,JUNCTION_MATE=R2,REQUIRE_MATE_ALIGNMENT=0 pipeline/hpc/ytab_mapFastq.txt
qsub -v PROJECT_ID=Gale_2020,JUNCTION_MATE=R1,REQUIRE_MATE_ALIGNMENT=1 pipeline/hpc/ytab_mapFastq.txt

# Request more cores through SGE; NSLOTS is passed to the CLI.
qsub -pe smp 16 -v PROJECT_ID=Gale_2020 pipeline/hpc/ytab_mapFastq.txt
```

`PROJECT_CONFIG=/absolute/path/project.yaml` can replace `PROJECT_ID`; no project
is silently chosen by default. Configuration/reference paths must be valid on the
HPC. `YTAB` defaults to `/Users/jayelazuno/workspace/ytab` and can be overridden.
Blank mate/agreement options preserve per-project/per-sample settings. Optional
`FORCE=1` bypasses mapping caches and `KEEP_GOING=0` stops after the first sample
failure. `SAMPLES` passes the CLI's comma-separated selection; SGE's `-v` also uses
commas, so export it first and use `qsub -v SAMPLES,PROJECT_ID=... ...` when selecting
multiple samples.

Default environment loading is Conda `tnseq` under the user's Miniforge installation.
It must already contain the [YTAB dependencies](../../environment.local.yml),
including Python Cutadapt, Biopython, pysam, PyYAML and pandas, plus Bowtie2 and
samtools. The job installs nothing. `YTAB_CONDA_ENV`/`YTAB_CONDA_SH` override these
personal environment settings. `YTAB_ENV_MODE=modules` instead loads the
space-separated names in exported `YTAB_MODULES`; choose the site's actual module
names/versions and ensure all Python packages remain available. The SRA Toolkit
module used for downloading is not an alignment dependency.

## Storage, monitoring and reruns

Final sorted BAMs, indices, mapping QC and manifests stay at the project's existing
configured output paths. The wrapper does **not** relocate them. Each job gets a
unique `TMPDIR` under `/nfsscratch/jayelazuno/workspace/ytab_alignment`; Hermes
cleanup honors it. Other final/intermediate alignment files still use the mapping
output directory. Override `SCRATCH_ROOT` only if another temporary location is
needed. Allocate/check project output quota as well as scratch quota.

SGE writes `YTAB_MapFastq.log` in the submission directory; the job also writes
`alignment.log` inside its unique scratch directory. Failed/successful job logs
are retained. Site scratch-purge policies apply. Submit only one mapping job per
project at a time. Existing input FASTQs are never deleted.

```bash
qstat -u jayelazuno
tail -f YTAB_MapFastq.log
```

## Create hit files on Argon after successful alignment

The second wrapper calls the existing `ytab.py hit-files` runner. It reads each
included sample's `mapfastq/<sample>/<sample>.sorted.bam` and requires its BAM
index, the configured reference FASTA and annotation. It does not realign reads
or change BAMs. No FASTQ directory or Bowtie2 executable is needed for this stage.
The default environment is the same `tnseq` environment; Conda/module overrides,
`SAMPLES`, `FORCE`, `DRY_RUN` and `KEEP_GOING` work as for the mapping wrapper.

```bash
# Run after checking the mapping log for successful completion.
qsub -v PROJECT_ID=Gale_2020 pipeline/hpc/ytab_CreateHitFile.txt

# Optional preview; validates inputs and may write preview manifests.
qsub -v PROJECT_ID=Gale_2020,DRY_RUN=1 pipeline/hpc/ytab_CreateHitFile.txt

# If alignments were deliberately replaced, rebuild affected hit files.
qsub -v PROJECT_ID=Gale_2020,FORCE=1 pipeline/hpc/ytab_CreateHitFile.txt

tail -f YTAB_CreateHitFile.log
```

SGE can delay submission until a mapping job finishes using
`qsub -hold_jid MAPPING_JOB_ID -v PROJECT_ID=Gale_2020 pipeline/hpc/ytab_CreateHitFile.txt`.
Replace `MAPPING_JOB_ID` with the actual numeric job ID. **SGE holds wait for
completion, not success**; missing BAMs/indices fail validation, but a hold is not
a substitute for checking mapping logs and sample statuses.

Final outputs remain in the configured project's `create_hit_file/`,
`manifests/create_hit_file/` and `logs/create_hit_file/`, with tracks in
`output/exports/<project>/browser_tracks/create_hit_file/`. SGE writes
`YTAB_CreateHitFile.log`; a second `hit_files.log` is retained in a unique
temporary directory under `/nfsscratch/jayelazuno/workspace/ytab_hit_files`
(override with `SCRATCH_ROOT`). Submit only one hit-file job per project at a time.
Neither wrapper submits the other or any downstream analysis automatically.

## Continue locally after hit-file creation

Transfer the matching project's `create_hit_file/`, `manifests/create_hit_file/`
and `logs/create_hit_file/` directories, plus its exported browser tracks, back
to the same project on the PC. Keep sample IDs and reference assembly/annotation
consistent, and retain local configuration paths rather than overwriting them
with HPC paths. These hit files are the inputs for local SummaryTable generation;
the large BAMs are not required for that stage.

```bash
python scripts/local/ytab.py summary --project-config PATH/TO/project.yaml
```

Alternatively, open the project in Shiny and run **SummaryTable**. If hit files
were not generated on Argon, transfer alignments and run Hit Files locally first:

Transfer the matching project's `mapfastq/` directory (including sorted BAMs,
indices, mapping statistics and Hermes cleanup reports) and `manifests/mapfastq/`
back to the same project on the PC. The project's `logs/mapfastq/` and
`output/exports/<project>/qc/mapping/` are also useful to retain. Keep sample IDs
consistent and use the same reference assembly/annotation. Keep the PC's local
configuration/reference paths rather than replacing them with HPC absolute paths.

Open that project in Shiny and continue with **Hit Files**, or use the ordinary CLI:

```bash
python scripts/local/ytab.py hit-files --project-config PATH/TO/project.yaml
python scripts/local/ytab.py summary --project-config PATH/TO/project.yaml
```

Do not rerun the mapping stage merely to load HPC results; mapping cache keys can
legitimately differ between machines. No scheduler is needed for the local steps.

The old hardcoded multi-stage/Zn scripts remain retired and recoverable from Git.
The master and shared stage template above replace that submission logic without
copying the scientific pipeline. General metadata and specialized diagnostics in
`scripts/hpc` are separate utilities. After intentionally replacing alignments,
refresh affected downstream outputs rather than reusing old hit tables/results.
