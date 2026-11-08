#!/usr/bin/env python3
"""Simulation check of the response correlation (multicondition_plots): bias and CI coverage.

Uses the multi-condition benchmark generator (scripts/local/ytab_hermes_anova_benchmark_v2.py) with a new
seed (818), 8 experiments per scenario, BGC + iterative-median recentering as in `ytab multicondition`.
Per pair of non-reference conditions: split-reference (replicate-matched) rho and its gene-bootstrap 95% CI,
and the shared-reference rho, on all tested genes and on ANOVA-responsive genes.
Expected sign of the true response correlation: 0 (null; response in one condition only), < 0 (opposite
effects), > 0 (shared depletion). Writes output/essentiality/multicondition/plots_check/.
"""
from concurrent.futures import ProcessPoolExecutor
import importlib.util
from pathlib import Path
import sys

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import multicondition as mc
from ytab.essentiality import multicondition_plots as mp
from ytab.normalization.methods import normalize
from ytab.normalization.ttr import write_table

spec = importlib.util.spec_from_file_location("bench", ROOT/"scripts/local/ytab_hermes_anova_benchmark_v2.py")
bench = importlib.util.module_from_spec(spec); spec.loader.exec_module(bench)
OUT = ROOT/"output/essentiality/multicondition/plots_check"
SEED, EXPERIMENTS = 818, 8
SCENARIOS = {"null_3x2":"0", "null_4x2":"0", "depletion_4x_one_condition":"0", "opposite_effects":"<0", "shared_depletion_4x":">0"}


def task(item):
    name, k = item
    raw, conds, libs, genes, idx, truth, names = bench.simulate(name, k, SEED)
    values = normalize(sparse.csc_matrix(raw), "betageom", libs)[0].toarray()
    code = np.array([names.index(c) for c in conds])
    centred, _ = mc.recenter(values, code, 0, idx, "iterative-median")
    S, N = mc.gene_sample_sums(centred, idx)
    tested = N.max(1) > 0
    rows, _ = mc.anova(values, conds, genes, idx, settings=mc.Settings(seed=SEED))
    responsive = np.array([r["call"] == "varies" for r in rows])
    out = []
    for subset, mask in (("all tested", tested), ("responsive", responsive)):
        for r in mp.response_correlation(S, N, conds, libs, "C", mask, boot=300, seed=k):
            out.append(dict(scenario=name, experiment=k, gene_set=subset, pair=f"{r['condition_a']}-{r['condition_b']}", genes=r["genes"],
                            rho=r["rho"], ci_low=r["ci_low"], ci_high=r["ci_high"], rho_shared=r["rho_shared_reference"]))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(4) as pool:rows = [r for chunk in pool.map(task, [(s, k) for s in SCENARIOS for k in range(EXPERIMENTS)]) for r in chunk]
    write_table(OUT/"experiments.tsv", list(rows[0]), rows)
    summary = []
    for s, expected in SCENARIOS.items():
        for subset in ("all tested", "responsive"):
            for pair in sorted({r["pair"] for r in rows if r["scenario"] == s}):
                sel = [r for r in rows if r["scenario"] == s and r["gene_set"] == subset and r["pair"] == pair and r["genes"] >= 10]
                if not sel:continue
                rho = np.array([r["rho"] for r in sel]); shared = np.array([r["rho_shared"] for r in sel])
                cover = sum(r["ci_low"] <= 0 <= r["ci_high"] for r in sel)
                summary.append(dict(scenario=s, expected=expected, gene_set=subset, pair=pair, experiments=len(sel), median_genes=int(np.median([r["genes"] for r in sel])),
                                    rho_split_mean=round(float(rho.mean()), 3), rho_split_sd=round(float(rho.std()), 3), rho_shared_mean=round(float(shared.mean()), 3),
                                    ci_covers_zero=f"{cover}/{len(sel)}"))
    write_table(OUT/"summary.tsv", list(summary[0]), summary)
    for r in summary:print(r)


if __name__ == "__main__":main()
