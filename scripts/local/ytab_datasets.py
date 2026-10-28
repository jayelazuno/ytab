#!/usr/bin/env python3
"""Catalogue of the public Tn-seq datasets YTAB has reprocessed, for the app's Datasets page and its search.

Each study is a row of resources/datasets/studies.tsv (citation, DOI, strain, transposon, accession, design); its runs
are the rows of pipeline/hpc/datasets/<dataset>/samples.csv whose `study` column matches. For every study and project
the catalogue adds what the runs are (roles, strains, conditions), their QC (newest `ytab qc` table) and which results
the app can open without computing anything (QC, essentiality, normalization, condition comparisons, genome browser,
library correlation, cross-species).

  ytab datasets                      # writes resources/datasets/catalog.tsv and libraries.tsv
  ytab datasets --output DIR

A new dataset appears once it has a dataset folder and a row in studies.tsv; nothing here names a study.
"""
import argparse
import csv
from collections import Counter
import json
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[2]
STUDIES = ROOT/"resources/datasets/studies.tsv"
DATASETS = ROOT/"pipeline/hpc/datasets"
PROJECTS = ROOT/"output/projects"

# result folder -> (app tab, file that marks a finished run)
RESULTS = [("essentiality", "Essentiality", "models/genes.tsv"), ("normalization", "Site normalization", "manifest.json"),
           ("treated_vs_parent", "Treated vs control (z-score)", "zscore_manifest.json"), ("treated_vs_control", "Treated vs control (resampling/Gumbel)", "manifest.json"),
           ("multicondition", "Multi-condition", "manifest.json"), ("genetic_interactions", "Genetic interactions", "run_manifest.json"),
           ("library_correlation", "Library correlation", "manifest.json")]


def read_tsv(path, delimiter="\t"):
    with open(path, newline="", encoding="utf-8") as h:return list(csv.DictReader(h, delimiter=delimiter))


def env(path):
    return dict(line.strip().split("=", 1) for line in open(path) if "=" in line and not line.startswith("#"))


def finished_runs(folder, marker):
    """Runs under a result folder that have their marker file (one level, or two for per-method folders)."""
    if not folder.is_dir():return []
    runs = [d for d in folder.iterdir() if d.is_dir() and ((d/marker).exists() or any((c/marker).exists() for c in d.iterdir() if c.is_dir()))]
    return sorted(r.name for r in runs)


def newest_qc(project):
    runs = sorted((PROJECTS/project/"qc").glob("qc_*/qc_table.tsv"))
    return {r["sample"]:r for r in read_tsv(runs[-1])} if runs else {}


def cross_species_runs():
    out = {}
    for m in sorted((ROOT/"output/cross_species").glob("*/manifest.json")):
        for p in json.loads(m.read_text()).get("projects", []):out.setdefault(p.get("project"), []).append(m.parent.name)
    return out


def number(value):
    try:return float(value)
    except (TypeError, ValueError):return None


def counts(values):
    return "; ".join(f"{k} {n}" for k, n in Counter(v for v in values if v).most_common())


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output", type=Path, default=ROOT/"resources/datasets", help="Folder for catalog.tsv and libraries.tsv")
    a = p.parse_args()
    cross = cross_species_runs(); catalog, libraries = [], []
    for study in read_tsv(STUDIES):
        for dataset in study["datasets"].split():
            folder = DATASETS/dataset; project = env(folder/"dataset.env")["PROJECT_ID"]
            runs = [r for r in read_tsv(folder/"samples.csv", ",") if r.get("study") == study["study"]]
            if not runs:continue
            qc = newest_qc(project); root = PROJECTS/project
            included = [r for r in runs if r.get("include", "True").strip().lower() not in ("false", "0", "no")]
            for r in runs:
                q = qc.get(r["sample"], {})
                libraries.append(dict(study=study["study"], citation=study["citation"], project=project, sample=r["sample"], display_name=r.get("display_name", ""),
                                      role=r["library_role"], condition=r["condition"], background=r["background"], pool=r["pool"], transposon=r["transposon"],
                                      library_design=r["library_design"], included=r in included, input_reads=q.get("input_reads", ""),
                                      percent_mapped=q.get("percent_mapped", ""), insertion_sites=q.get("unique_sites", ""),
                                      sites_per_kb=q.get("insertions_per_kb", ""), top_site_share=q.get("max_share", ""), top_site=q.get("max_site", ""),
                                      notes=r.get("notes", "")))
            density = [x for x in (number(qc.get(r["sample"], {}).get("insertions_per_kb")) for r in included) if x is not None]
            share = [(number(qc[r["sample"]].get("max_share")), qc[r["sample"]].get("max_site", "")) for r in included if number(qc.get(r["sample"], {}).get("max_share")) is not None]
            top = max(share, default=(None, ""))
            results = {label:finished_runs(root/folder_name, marker) for folder_name, label, marker in RESULTS}
            results["Genome browser"] = ["genome_browser"] if (root/"genome_browser/manifest.json").exists() else []
            results["Cross-species"] = cross.get(project, [])
            tabs = ["Library QC"] if qc else []
            tabs += [f"{label} ({len(v)})" for label, v in results.items() if v]
            catalog.append(dict(study=study["study"], citation=study["citation"], title=study["title"], journal=study["journal"], year=study["year"],
                                doi=study["doi"], pmid=study["pmid"], species=study["species"], strain=study["strain"], transposon=study["transposon"],
                                read_design=study["read_design"], accession=study["accession"], design=study["design"], project=project,
                                runs=len(runs), runs_included=len(included), roles=counts(r["library_role"] for r in included),
                                backgrounds=counts(r["background"] for r in included), conditions=counts(r["condition"] for r in included),
                                median_sites_per_kb=round(median(density), 1) if density else "", max_top_site_share=round(top[0], 3) if top[0] is not None else "",
                                max_top_site=top[1], in_app="; ".join(tabs),
                                essentiality_runs=" ".join(results["Essentiality"]), project_config=f"output/projects/{project}/config/project.yaml"))
    a.output.mkdir(parents=True, exist_ok=True)
    for name, rows in (("catalog.tsv", catalog), ("libraries.tsv", libraries)):
        with open(a.output/name, "w", newline="", encoding="utf-8") as h:
            w = csv.DictWriter(h, list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(json.dumps(dict(studies=len({r["study"] for r in catalog}), projects=len({r["project"] for r in catalog}), runs=len(libraries),
                          output=str(a.output)), indent=1))


if __name__ == "__main__":main()
