#!/usr/bin/env python3
"""Library complexity comparison: our libraries (zinc, H2O2) vs Gale 2020.

Per library: distinct sites at matched depths (binomial thinning, seed 2026); saturation as the share of
the library's final sites already seen at 1M reads, and the new sites gained by the second half of the
sequencing (per million reads; ~0 = sequencing exhausted the library); clones (>= 1 read per million);
gene coverage (genes with >= 1 clone among genes with >= 1 catalogue site, N/C-terminal trimming). A
Michaelis-Menten asymptote was tried and dropped: it pins at the observed count because the curves keep a
slowly growing singleton tail. Hit files only. Writes
output/essentiality/pooled_gumbel/library_complexity/comparison.tsv and comparison.png.
"""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.normalization.ttr import write_table

OUT = ROOT/"output/essentiality/pooled_gumbel/library_complexity"
PROJECTS = ("Zn_toxicity_screen", "H2O2_screen_v1", "Gale_2020")
GFF = ROOT/"resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff"
DEPTHS = np.array([25_000, 50_000, 100_000, 250_000, 500_000, 1_000_000, 2_000_000, 5_000_000])
SEED = 2026


def load(path):
    chrom = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=(0,), dtype=str, ndmin=1)
    pos, count = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=(11, 12), dtype=np.int64, ndmin=2).T
    keep = count > 0
    return chrom[keep], pos[keep], count[keep]


def main():
    rng = np.random.default_rng(SEED); genes = ccg.annotations(GFF); rows, curves = [], {}
    for project in PROJECTS:
        for path in sorted((ROOT/"output/projects"/project/"create_hit_file").glob("*/*_hits.txt")):
            chrom, pos, c = load(path); total = int(c.sum())
            depths = DEPTHS[DEPTHS < total]; sites = [int((rng.binomial(c, d/total) > 0).sum()) for d in depths]
            x = np.r_[depths, total].astype(float); y = np.r_[sites, len(c)].astype(float)
            half = int((rng.binomial(c, .5) > 0).sum())
            clones = c >= total/1e6
            idx = ccg.gene_sites(genes, list(zip(chrom, pos)), ccg.Settings())
            covered = np.array([clones[i].any() for i in idx]); per_gene = np.array([clones[i].sum() for i in idx])
            at = lambda d:int((rng.binomial(c, d/total) > 0).sum()) if d <= total else None
            rows.append(dict(project=project, library=path.parent.name, reads=total, sites=len(c), clones=int(clones.sum()),
                             sites_at_1M=at(1_000_000), sites_at_5M=at(5_000_000), share_of_sites_seen_by_1M=at(1_000_000)/len(c),
                             new_sites_per_M_reads_last_half=(len(c)-half)/(total/2e6), reads_per_site=total/len(c),
                             reads_per_clone=float(c[clones].sum()/clones.sum()), genes_with_clone=float(covered.mean()),
                             median_clones_per_gene=float(np.median(per_gene))))
            curves[(project, path.parent.name)] = (x, y)
    write_table(OUT/"comparison.tsv", list(rows[0]), rows)
    summary = []
    for project in PROJECTS:
        r = [x for x in rows if x["project"] == project]
        summary.append(dict(project=project, libraries=len(r), **{k:f"{np.median([x[k] for x in r]):.4g} ({min(x[k] for x in r):.3g}-{max(x[k] for x in r):.3g})"
                       for k in ("reads", "sites", "clones", "sites_at_1M", "share_of_sites_seen_by_1M", "new_sites_per_M_reads_last_half", "reads_per_clone",
                                 "genes_with_clone", "median_clones_per_gene")}))
    write_table(OUT/"comparison_summary.tsv", list(summary[0]), summary)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    colours = {"Zn_toxicity_screen":"#c0392b", "H2O2_screen_v1":"#e69f00", "Gale_2020":"#2c7fb8"}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for (project, name), (x, y) in curves.items():
        axes[0].plot(x, y, marker="o", ms=3, color=colours[project], alpha=.7, lw=1)
    axes[0].set_xscale("log"); axes[0].set_xlabel("reads (thinned; last point = sequenced depth)"); axes[0].set_ylabel("distinct sites")
    axes[0].set_title("Rarefaction"); handles = [plt.Line2D([], [], color=c, marker="o", label=l) for l, c in zip(["zinc", "H2O2", "Gale"], colours.values())]
    axes[0].legend(handles=handles, frameon=False)
    for j, (key, label) in enumerate((("clones", "clones (>= 1 read per million)"), ("genes_with_clone", "fraction of genes with >= 1 clone")), 1):
        for i, project in enumerate(PROJECTS):
            v = [r[key] for r in rows if r["project"] == project]
            axes[j].scatter(np.full(len(v), i)+np.linspace(-.15, .15, len(v)), v, color=colours[project], s=30)
        axes[j].set_xticks(range(3)); axes[j].set_xticklabels(["zinc", "H2O2", "Gale"]); axes[j].set_ylabel(label)
        if key == "clones":axes[j].set_yscale("log")
    for ax in axes:ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(OUT/"comparison.png", dpi=150)


if __name__ == "__main__":main()
