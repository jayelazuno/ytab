#!/usr/bin/env python3
"""Plan the donor-read counting for every dataset.

Counting donor reads does not need the libraries mapped again: `ytab donor-reads` puts the reads against the donor's
own index (a few kilobases) and only sends the reads that hit it, or that carry a recurring junction flank, to the
genome. One pass per library, no hit files or downstream output touched. A library being mapped for the first time
gets the same table for free, because `ytab mapfastq` writes it as it runs.

Prints, per dataset, the libraries, the donor attached (or the flank estimate), and the command. Nothing is submitted
from here: the cluster session is the user's.

  python scripts/local/ytab_donor_remap_plan.py              # the plan
  python scripts/local/ytab_donor_remap_plan.py --commands   # just the commands
"""
from __future__ import annotations
import argparse, csv
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[2]; DATASETS = ROOT/"pipeline/hpc/datasets"
LOCAL = ("H2O2_screen_v1", "Zn_toxicity_screen", "local_glabrata_smoke_v1")   # the lab's own data, on this machine


def project_of(dataset: Path) -> str:
    for line in (dataset/"dataset.env").read_text().splitlines():
        if line.startswith("PROJECT_ID="):return line.split("=", 1)[1].strip().strip('"')
    return ""


def donor_of(project: str) -> str:
    config = ROOT/"output/projects"/project/"config/project.yaml"
    if not config.is_file():return ""
    return ((yaml.safe_load(config.read_text()) or {}).get("mapping") or {}).get("donor_fasta", "")


def rows():
    out = []
    for d in sorted(DATASETS.iterdir()):
        if (d/"samples.csv").is_file():
            project = project_of(d)
            out.append((d.name, project, sum(1 for _ in csv.DictReader(open(d/"samples.csv"))), donor_of(project)))
    for project in LOCAL:
        sheet = ROOT/"output/projects"/project/"config/sample_sheet.csv"
        if sheet.is_file():
            n = sum(1 for r in csv.DictReader(open(sheet)) if str(r.get("include", "true")).strip().lower() in ("1", "true", "yes", "y"))
            out.append(("(local)", project, n, donor_of(project)))
    return out


def main(commands_only: bool):
    table = rows()
    if not commands_only:
        print(f"{'dataset':34} {'project':32} {'libraries':>9}  donor")
        for dataset, project, n, donor in table:
            print(f"{dataset:34} {project:32} {n:>9}  {Path(donor).parent.name if donor else 'flank estimate'}")
        print(f"\n{sum(n for _, _, n, _ in table)} libraries. One donor pass each; nothing is remapped to the genome,"
              f"\nso hit files, insertions and every downstream result stay as they are. Rerun `ytab qc` afterwards so"
              f"\nthe QC table picks up the donor metrics.\n\nCommands:")
    for dataset, project, _, _ in table:
        config = f"output/projects/{project}/config/project.yaml"
        if dataset == "(local)":print(f"ytab donor-reads --project-config {config} --threads 8")
        else:
            # the projects already exist on Argon, so the chain is submitted against the existing config rather than
            # through ytab_public_dataset.sh, which refuses to recreate a project
            print(f"# {dataset}\nbash pipeline/hpc/submit_pipeline.sh --project-config {config} --stages donor_reads --threads 8")


if __name__ == "__main__":
    a = argparse.ArgumentParser(); a.add_argument("--commands", action="store_true"); main(a.parse_args().commands)
