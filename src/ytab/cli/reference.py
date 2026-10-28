"""Command-line reference (docs/cli_reference.md): every registered module with its description and its own --help.

  ytab modules --markdown > docs/cli_reference.md

Built from the registry and each module's argparse help, so the page always matches the code.
"""
import os
import subprocess
import sys
from .registry import MODULES, ROOT, command

# Order of the analysis, as the app's tabs follow it; a module missing here goes under "Other".
STAGES = [
    ("Start a project", ["init", "import-hits", "reference", "status"]),
    ("Run everything", ["pipeline"]),
    ("Sequencing and mapping", ["preprocess", "mapfastq", "hit-files", "merge-runs", "summary"]),
    ("Quality control", ["qc", "library-qc", "tnseq-stats", "library-correlation"]),
    ("Normalization", ["normalize-project", "normalize", "normalize-experiment", "normalize-genes", "ttr", "bgc-audit", "abgc-audit",
                       "normalization", "summary-normalized", "combine-hits", "summary-combined"]),
    ("Essentiality", ["essentiality", "hermes-hmm", "hybrid-inputs", "classifier", "cross-species"]),
    ("Condition comparisons", ["comparison-design", "fitness", "fitness-zscore", "treated-vs-control", "resampling", "control-gumbel",
                               "multicondition-project", "multicondition", "multicondition-plots", "genetic-interactions"]),
    ("Gene Explorer", ["browser-tracks", "gene-insertions", "go"]),
    ("Screen planning and mutants", ["coverage", "mutants"]),
    ("Datasets, results and sharing", ["datasets", "report", "export", "resource-bundle", "app"]),
    ("Pre-registered benchmarks", ["tvc-benchmark", "tvc-phase2"]),
]

INTRO = """# Command-line reference

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
python scripts/local/ytab.py treated-vs-control --project-config P/config/project.yaml --method resampling \
    --normalization-run P/normalization/RUN --normalization betageom --output P/treated_vs_control/RUN
python scripts/local/ytab.py browser-tracks --project-config P/config/project.yaml
python scripts/local/ytab.py app --project-config P/config/project.yaml   # open the results
```

`ytab run --workflow WORKFLOW.yaml` runs an ordered list of these commands (see Command line); `ytab run --project-config
P/config/project.yaml --profile core|classifier|fitness|all` runs the standard pipeline.

## Commands by stage
"""


def module_help(name):
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONPATH=os.pathsep.join(filter(None, [str(ROOT/"src"), os.environ.get("PYTHONPATH")])))
    try:
        out = subprocess.run(command(name, ["--help"]), capture_output=True, text=True, timeout=120, env=env)
        text = (out.stdout or out.stderr).strip()
    except (OSError, subprocess.TimeoutExpired) as e:
        text = f"(help unavailable: {e})"
    return text.replace("```", "'''")


def markdown():
    placed = {n for _, names in STAGES for n in names}
    stages = [(t, [n for n in names if n in MODULES]) for t, names in STAGES]
    stages.append(("Other", [n for n in MODULES if n not in placed]))
    stages = [(t, names) for t, names in stages if names]
    lines = [INTRO]
    for title, names in stages:
        lines += [f"### {title}", "", "| command | what it does |", "|---|---|"]
        lines += [f"| `{n}` | {MODULES[n].description} |" for n in names]
        lines.append("")
    lines += ["## Options of each command", ""]
    for title, names in stages:
        for n in names:
            m = MODULES[n]
            lines += [f"### {n}", "", f"{m.description}. Stage: {title}." + (f" Takes `{m.project_option}`." if m.project_option else ""), "",
                      "```text", module_help(n), "```", ""]
    return "\n".join(lines)


if __name__ == "__main__":sys.stdout.write(markdown())
