#!/usr/bin/env python3
"""Insertion density of parent (untreated) libraries: our H2O2 parents vs Gale 2020 parents (Cg1-3).

Per library, at full depth and thinned to a common depth (the smallest library; seed 2026): distinct
insertion sites (>= 1 read), genes hit (>= 1 site in the CDS, N/C-terminal trimming), median sites per gene,
and genes with >= 3 sites; the same with sites requiring >= 5 reads. Also the union over each project's
parents. Hit files only. Writes output/essentiality/pooled_gumbel/library_complexity/parent_coverage.tsv.
"""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.normalization.ttr import write_table

OUT = ROOT/"output/essentiality/pooled_gumbel/library_complexity"
GFF = ROOT/"resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff"
PARENTS = {"H2O2_screen_v1":[f"yH29{8 if p < 3 else 9}-parent-pool{p}" for p in range(1, 5)],
           "Gale_2020":["SRR11561508", "SRR11561507", "SRR11561506"]}   # Gale Cg1, Cg2, Cg3 parents
NAMES = {"SRR11561508":"Gale Cg1 parent", "SRR11561507":"Gale Cg2 parent", "SRR11561506":"Gale Cg3 parent"}
SEED = 2026


def load(project, sample):
    path = ROOT/"output/projects"/project/"create_hit_file"/sample/f"{sample}_hits.txt"
    chrom = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=(0,), dtype=str, ndmin=1)
    pos, count = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=(11, 12), dtype=np.int64, ndmin=2).T
    keep = count > 0
    out = {}
    for c, p, n in zip(chrom[keep], pos[keep], count[keep]):out[(c, int(p))] = out.get((c, int(p)), 0)+int(n)   # both strands of one position
    return out


def coverage(sites, genes):
    idx = ccg.gene_sites(genes, sorted(sites), ccg.Settings()); n = np.array([len(i) for i in idx])
    return dict(sites=len(sites), genes_hit=float(np.mean(n > 0)), median_sites_per_gene=float(np.median(n)), genes_with_3plus=float(np.mean(n >= 3)))


def main():
    rng = np.random.default_rng(SEED); genes = ccg.annotations(GFF)
    libs = {(p, s):load(p, s) for p, samples in PARENTS.items() for s in samples}
    common = min(sum(v.values()) for v in libs.values())
    thinned = {}
    for key, v in libs.items():
        k = list(v); c = np.array([v[x] for x in k]); t = rng.binomial(c, common/c.sum())
        thinned[key] = {x:int(n) for x, n in zip(k, t) if n > 0}
    rows = []
    for depth_label, data in (("full depth", libs), (f"thinned to {common/1e6:.1f}M reads", thinned)):
        for min_reads in (1, 5):
            sets = {key:{s for s, n in v.items() if n >= min_reads} for key, v in data.items()}
            for (project, sample), s in sets.items():
                rows.append(dict(depth=depth_label, min_reads=min_reads, library=NAMES.get(sample, f"H2O2 {sample}"),
                                 reads=sum(data[(project, sample)].values())) | coverage(s, genes))
            for project in PARENTS:
                union = set().union(*(s for (p, _), s in sets.items() if p == project))
                rows.append(dict(depth=depth_label, min_reads=min_reads, library=f"UNION {project} parents ({len(PARENTS[project])})",
                                 reads=sum(sum(data[k].values()) for k in data if k[0] == project)) | coverage(union, genes))
    write_table(OUT/"parent_coverage.tsv", list(rows[0]), rows)
    for r in rows:print("\t".join(str(round(v, 3)) if isinstance(v, float) else str(v) for v in r.values()))


if __name__ == "__main__":main()
