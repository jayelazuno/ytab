# Personal Argon alignment job

YTAB is primarily a desktop Shiny application. Its analysis modules and master
workflow also have a scheduler-independent [CLI](../../docs/command_line.md).
**No HPC account or scheduler is required to use YTAB.** Alignment is the optional
offload represented here; this does not imply that all other analyses are cheap
for every dataset.

`ytab_mapFastq.txt` is Joshua's University of Iowa **Argon / SGE (`qsub`)** job,
not a universal HPC interface. Queue `UI`, parallel environment `smp`, email,
Conda initialization and scratch paths are personal site settings. Users of
Slurm, PBS or other clusters should wrap the ordinary CLI in their own scripts.
The separate personal `sra_access.txt` download utility is retained; it is not an
analysis stage and is not called by the alignment job.

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

## Continue locally after alignment

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

The older multi-stage jobs, project-specific Zn submission chain and HPC-only job
smoke wrappers were retired. Their tracked versions remain recoverable from Git;
no scientific CLI modules or completed results were removed. General metadata and
specialized diagnostic utilities in `scripts/hpc` are not a supported full-pipeline
scheduler interface. Run subsequent stages through the desktop app or ordinary
CLI/master runner. After intentionally replacing alignments, refresh the affected
downstream outputs rather than reusing old hit tables/results.
