# Command-line analyses and master workflows

YTAB supports individual analysis modules, a sequential project pipeline, and an
extensible ordered workflow. These commands delegate to the existing backends used
by Shiny; they do not introduce different scientific methods. Inspired by TRANSIT's
named-method dispatch, YTAB has an independent module registry in
`src/ytab/cli/registry.py`. No TRANSIT implementation is copied here.

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

## Individual modules

Use `modules --json` for the machine-readable registry and `<module> --help` for
the backend's exact options. Available modules include reference preparation,
mapping, hit files, summaries, library QC, existing MidLC normalization, combined
parents, essentiality classification, fitness, GO enrichment, gene/domain insertion
plots, Poisson coverage, mutant identification, project status, reports and exports.
Initialization and hit-file import are also available.

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
- `classifier`: core plus parent MidLC normalization, normalized/combined summaries
  and the essentiality classifier.
- `fitness`: core plus comparison design and treated-versus-control fitness.
- `all`: all existing project-pipeline stages above, **not** every optional analysis.

`--stages`, `--from-stage` and `--to-stage` select existing project-pipeline stages.
The original pipeline retains its per-stage cache, project status and logs. MidLC
normalization and fitness's internal CPM normalization remain distinct; this work
does not implement new TRANSIT normalization algorithms.

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
