#!/usr/bin/env python3
"""Correlate Hermes libraries across pools, strains, conditions and species, generalized from TRANSIT's corrplot (DeJesus et al. 2015).

Each library is summarized per gene (insertion density by default) and libraries are compared by Spearman
correlation across genes. Different species are matched through 1:1 ortholog groups (--orthologs).
Each pair is labelled same pool (shared clones), same strain/different pool, different strain or
different species; same-pool pairs measure library structure, not biology.

Libraries come from YTAB projects (--project NAME: every included sample, with condition, background
(strain) and pool from the sample sheet) and/or an explicit table (--samples TSV with columns
name, hits, species, gff and optionally strain, pool, condition). Hit files: YTAB *_hits.txt or any
CSV/TSV with chromosome, position and read-count columns.

Example (three C. glabrata projects and S. pombe Hermes data):
  ytab library-correlation --project Gale_2020 --project H2O2_screen_v1 --project Zn_toxicity_screen \\
      --samples pombe.tsv --orthologs resources/species/orthologs/glabrata_cerevisiae_pombe_1to1.tsv \\
      --output output/library_correlation/RUN
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality.library_correlation import Sample, run


def project_samples(name):
    base = ROOT/"output/projects"/name; config = yaml.safe_load((base/"config/project.yaml").read_text())
    species = config.get("species") or config["reference"]["species"]; gff = ROOT/config["reference"]["gff"]
    out = []
    with (base/"config/sample_sheet.csv").open() as handle:
        for r in csv.DictReader(handle):
            if str(r.get("include", "true")).lower() not in ("true", "1", "yes"):continue
            hits = ROOT/r["hit_file"] if r.get("hit_file") else base/"create_hit_file"/r["sample"]/f"{r['sample']}_hits.txt"
            if not hits.is_file():raise ValueError(f"Missing hit file for {name}::{r['sample']}: {hits}")
            out.append(Sample(name=f"{name}::{r['sample']}", hits=str(hits), species=species, gff=str(gff),
                              strain=r.get("background") or r.get("guessed_background") or "", pool=r.get("pool") or r.get("guessed_pool") or "",
                              condition=r.get("condition") or r.get("guessed_condition") or "", source=name))
    return out


def table_samples(path):
    out = []
    with Path(path).open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            resolve = lambda p:str(Path(p) if Path(p).is_absolute() else (Path(path).parent/p if (Path(path).parent/p).exists() else ROOT/p))
            out.append(Sample(name=r["name"], hits=resolve(r["hits"]), species=r["species"], gff=resolve(r["gff"]), strain=r.get("strain", ""),
                              pool=r.get("pool", ""), condition=r.get("condition", ""), source=str(path)))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--project", action="append", default=[], help="YTAB project name under output/projects; repeat")
    parser.add_argument("--samples", action="append", default=[], help="TSV: name, hits, species, gff[, strain, pool, condition]; repeat")
    parser.add_argument("--orthologs", help="TSV of 1:1 ortholog groups (group, species, gene); needed for several species")
    parser.add_argument("--metric", choices=["density", "reads_per_kb"], default="density",
                        help="density (insertion sites per kb; default, comparable across pools) or reads_per_kb (pool-specific clone abundance)")
    parser.add_argument("--min-reads", type=int, default=1, help="Reads for a site to count as an insertion (default 1)")
    parser.add_argument("--bootstrap", type=int, default=500); parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        samples = [s for p in args.project for s in project_samples(p)]+[s for t in args.samples for s in table_samples(t)]
        if len(samples) < 2:raise ValueError("Give at least two libraries (--project and/or --samples)")
        manifest = run(samples, args.output, args.metric, args.orthologs, args.bootstrap, args.seed, args.min_reads)
    except (ValueError, KeyError, OSError) as error:parser.exit(1, f"ERROR: {error}\n")
    print(json.dumps({k:manifest[k] for k in ("status", "metric", "samples", "by_pair_type", "files")}, indent=2))


if __name__ == "__main__":main()
