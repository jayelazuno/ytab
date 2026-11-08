#!/usr/bin/env python3
"""Do the zinc cultures carry the clones of their own parent pool? (H2O2 thawed-parent libraries as reference)

Clones: sites with >= 1 read per million of the library. For every ordered pair of libraries (A, B):
the fraction of A's clones seen in B (>= 1 read), the share of A's reads at those clones, and the clone
Jaccard. If pools are distinct and libraries capture their pool well, a zinc culture's clones should be
found far more often in its own parent pool than in the other pools. Hit-file coordinates of the two
projects were checked to agree (exact-position peak, Hermes +-7 bp duplication pattern as within zinc).
Writes output/essentiality/pooled_gumbel/parent_pool_comparison/.
"""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.normalization.ttr import write_table

OUT = ROOT/"output/essentiality/pooled_gumbel/parent_pool_comparison"
LIBS = ([("H2O2_screen_v1", f"yH29{8 if p < 3 else 9}-parent-pool{p}", f"parent{p}", p) for p in range(1, 5)]+
        [("Zn_toxicity_screen", f"yH29{8 if p < 3 else 9}-parent-pool{p}-mock", f"Zn-mock{p}", p) for p in range(1, 5)]+
        [("Zn_toxicity_screen", f"yH29{8 if p < 3 else 9}-parent-pool{p}-1_5mM-Zn-treated", f"Zn-treated{p}", p) for p in range(1, 5)]+
        [("H2O2_screen_v1", f"yH29{8 if p < 3 else 9}-H2O2-treated-facs-pool{p}", f"H2O2-treated{p}", p) for p in range(1, 5)])


def load(project, sample):
    path = ROOT/"output/projects"/project/"create_hit_file"/sample/f"{sample}_hits.txt"
    chrom = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=(0,), dtype=str, ndmin=1)
    pos, count = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=(11, 12), dtype=np.int64, ndmin=2).T
    keep = count > 0
    return {(c, int(p)):int(n) for c, p, n in zip(chrom[keep], pos[keep], count[keep])}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = {label:load(project, sample) for project, sample, label, _ in LIBS}; pools = {label:p for *_, label, p in LIBS}
    clones = {k:{s for s, n in v.items() if n >= sum(v.values())/1e6} for k, v in data.items()}
    rows = []
    for a in data:
        total = sum(data[a][s] for s in clones[a])
        for b in data:
            if a == b:continue
            found = {s for s in clones[a] if s in data[b]}
            rows.append(dict(a=a, b=b, same_pool=pools[a] == pools[b], a_clones=len(clones[a]), clones_found_in_b=len(found)/len(clones[a]),
                             read_share_found_in_b=sum(data[a][s] for s in found)/total,
                             clone_jaccard=len(clones[a] & clones[b])/len(clones[a] | clones[b])))
    write_table(OUT/"pairs.tsv", list(rows[0]), rows)
    group = lambda r:(r["a"].rstrip("1234"), r["b"].rstrip("1234"), "same pool" if r["same_pool"] else "other pools")
    summary = {}
    for r in rows:summary.setdefault(group(r), []).append(r)
    out = [dict(a=k[0], b=k[1], pools=k[2], pairs=len(v), clones_found_in_b=float(np.mean([r["clones_found_in_b"] for r in v])),
                read_share_found_in_b=float(np.mean([r["read_share_found_in_b"] for r in v])), clone_jaccard=float(np.mean([r["clone_jaccard"] for r in v])))
           for k, v in sorted(summary.items())]
    write_table(OUT/"summary.tsv", list(out[0]), out)
    for r in out:print("\t".join(str(round(v, 3)) if isinstance(v, float) else str(v) for v in r.values()))


if __name__ == "__main__":main()
