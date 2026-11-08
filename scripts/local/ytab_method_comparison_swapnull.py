#!/usr/bin/env python3
"""Method-comparison study: the swap null on Gale 2020 (control 1 + FLU 1 vs control 2 + FLU 2; drug effects cancel, so
calls are false calls) for every method x normalization arm of the pre-registered benchmark (ytab.benchmarks.
treated_vs_control.ALL_ARMS), run in-process on the count matrix exactly as the phase-2 benchmark did (the catalogue's
role check does not apply: each side deliberately mixes a control and a treated culture).
Output: output/comparison/Gale_2020/sweep/swap_null_arms.tsv (per arm: genes tested, calls, depleted, enriched)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT/"src")); sys.path.append(str(ROOT/"scripts/local"))
from ytab.benchmarks import treated_vs_control as tvc  # noqa: E402
from ytab.essentiality import control_conditioned_gumbel as ccg, resampling  # noqa: E402
import ytab_tvc_phase2 as phase2  # noqa: E402

run, samples = phase2.DESIGNS["gale_swap_null"]
raw, sites = phase2.load(run, samples)
genes = ccg.annotations(phase2.GFF)
gene_r = ccg.gene_sites(genes, sites, resampling.Settings()); gene_g = ccg.gene_sites(genes, sites, ccg.Settings())
positions = np.array([p for _, p in sites], dtype=np.int64)
arms = {k: v for k, v in tvc.ALL_ARMS.items() if not k.endswith("_pairs")}
rows = []
for arm, spec in arms.items():
    r = tvc.run_arms(raw, None, genes, gene_r, gene_g, positions, arms={arm: spec})[arm]
    c = r["call"].to_numpy(); rows.append(dict(arm=arm, method=spec[0], normalization=spec[1], genes=int(r["score"].notna().sum()),
                                               calls=int((c != 0).sum()), depleted=int((c < 0).sum()), enriched=int((c > 0).sum())))
    print(rows[-1], flush=True)
out = ROOT/"output/comparison/Gale_2020/sweep"; out.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(out/"swap_null_arms.tsv", sep="\t", index=False)
