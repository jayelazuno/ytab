#!/usr/bin/env python3
"""Library complexity of Hermes hit files at equal sequencing depth (zinc, H2O2 and Gale).

Reads are thinned binomially (seeded) to common depths; at each depth we count distinct insertion
sites and 'clones' (sites with >= 1 read per million of the thinned library), and the read share of
the 1,000 largest sites. Hit files carry no fragment ends, so PCR duplicates cannot be separated from
molecules here; see scripts/hpc/ytab_fragment_complexity.py for that. Writes
output/essentiality/pooled_gumbel/library_complexity/.
"""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.normalization.ttr import write_table

OUT = ROOT/"output/essentiality/pooled_gumbel/library_complexity"
DEPTHS = (250_000, 500_000, 1_000_000, 2_000_000, 5_000_000, 9_000_000)
SEED = 2026


def libraries():
    for project in ("Zn_toxicity_screen", "H2O2_screen_v1", "Gale_2020"):
        for path in sorted((ROOT/"output/projects"/project/"create_hit_file").glob("*/*_hits.txt")):
            yield project, path.parent.name, path


def counts(path):
    values = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=(12,), dtype=np.int64, ndmin=1)
    return values[values > 0]


def summary(c):
    s = np.sort(c)[::-1]; clones = c >= c.sum()/1e6
    return dict(reads=int(c.sum()), sites=int(len(c)), clones=int(clones.sum()), reads_per_site=float(c.sum()/len(c)),
                reads_per_clone=float(c[clones].sum()/max(clones.sum(), 1)), singleton_sites=float(np.mean(c == 1)),
                top1000_read_share=float(s[:1000].sum()/c.sum()))


def main():
    OUT.mkdir(parents=True, exist_ok=True); rng = np.random.default_rng(SEED); full, curves = [], []
    for project, sample, path in libraries():
        c = counts(path); full.append(dict(project=project, sample=sample) | summary(c))
        for depth in DEPTHS:
            if depth > c.sum():continue
            thinned = rng.binomial(c, depth/c.sum()); thinned = thinned[thinned > 0]
            curves.append(dict(project=project, sample=sample, depth=depth) | summary(thinned))
        print(project, sample, full[-1]["reads"], full[-1]["sites"], full[-1]["clones"], flush=True)
    write_table(OUT/"libraries.tsv", list(full[0]), full); write_table(OUT/"rarefaction.tsv", list(curves[0]), curves)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for key, ax in zip(("sites", "clones"), axes):
        for project, colour in (("Zn_toxicity_screen", "#c0392b"), ("H2O2_screen_v1", "#e69f00"), ("Gale_2020", "#2c7fb8")):
            for sample in sorted({r["sample"] for r in curves if r["project"] == project}):
                pts = [(r["depth"], r[key]) for r in curves if r["sample"] == sample and r["project"] == project]
                ax.plot(*zip(*pts), marker="o", ms=3, color=colour, alpha=.8, label=project if sample == sorted({r["sample"] for r in curves if r["project"] == project})[0] else None)
        ax.set_xscale("log"); ax.set_xlabel("reads (thinned)"); ax.set_ylabel(key if key == "sites" else "clones (>= 1 read per million)")
        ax.grid(alpha=.3)
    axes[0].legend(frameon=False); fig.suptitle("Library complexity at equal depth: zinc, H2O2, Gale"); fig.tight_layout()
    fig.savefig(OUT/"rarefaction.png", dpi=150)


if __name__ == "__main__":main()
