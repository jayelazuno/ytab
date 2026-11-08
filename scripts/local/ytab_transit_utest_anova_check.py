#!/usr/bin/env python3
"""Fit check: do TRANSIT's Mann-Whitney U-test and ANOVA suit Hermes two-condition screens?

Exploratory (not an adoption benchmark): faithful ports of TRANSIT utest.py and anova.py
(docs/codex/transit, commit 50e23d6) and one Hermes-adapted rank test, run next to resampling v3 on
- the resampling simulator (scripts/local/ytab_resampling_benchmark.py scenarios) with a fresh seed, and
- Gale 2020 FLU vs control (TTR): a balanced label-swap null and the real comparison.
Tests (per gene; sites = gene sites with reads in any compared sample, as TRANSIT's default):
  utest        TRANSIT: Mann-Whitney U on control vs treatment site x replicate counts pooled (two-sided)
  anova        TRANSIT: one-way F on the same pooled counts, MSE moderated by alpha = 1000
  signed_rank  adapted: Wilcoxon signed-rank on per-site log2((mean T + 1)/(mean C + 1)) (site-paired)
  resampling   YTAB resampling v3 (BGC + trimmed-sum recentering, exact site-restricted)
BH over genes; a call is q <= 0.05 with the direction of the mean change.
Writes output/essentiality/transit_tests/fit_check/.
"""
from concurrent.futures import ProcessPoolExecutor
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.essentiality import resampling
from ytab.normalization.methods import normalize
from ytab.normalization.ttr import write_table

spec = importlib.util.spec_from_file_location("v1", ROOT/"scripts/local/ytab_resampling_benchmark.py")
v1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(v1)
OUT = ROOT/"output/essentiality/transit_tests/fit_check"
SEED, EXPERIMENTS = 8128, 8
v1.SEED = SEED
v1.SCENARIOS.update({"null_noisy_control":dict(phi=.5, treatment_phi=.02)})
SCENARIOS = ["null_equal_depth", "null_shallow_treatment", "null_high_dispersion", "null_jackpot_prior", "null_noisy_treatment",
             "null_noisy_control", "null_bottleneck_treatment", "complete_loss_block", "depletion_4x_block", "depletion_2x_block",
             "uniform_depletion_2x_whole_gene", "enrichment_4x_block"]
ALPHA = 1000.   # TRANSIT anova default moderation
FDR = .05


def bh(p):
    p = np.asarray(p, float); n = len(p); order = np.argsort(p); q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order]*n/np.arange(1, n+1))[::-1])[::-1]
    return np.minimum(q, 1.)


def gene_tests(control, treatment, idx):
    """p-values and directions for the three rank/ANOVA tests, one row per gene."""
    rows = []
    for ii in idx:
        c, t = control[ii], treatment[ii]; keep = (c.sum(1)+t.sum(1)) > 0; c, t = c[keep], t[keep]
        row = dict(sites=int(keep.sum()), direction=float(np.sign(t.mean()-c.mean())) if keep.any() else 0.)
        if keep.sum() == 0:rows.append(row | dict(utest=1., anova=1., signed_rank=1.)); continue
        d1, d2 = c.ravel(), t.ravel()
        try:row["utest"] = float(stats.mannwhitneyu(d1, d2, alternative="two-sided").pvalue)
        except ValueError:row["utest"] = 1.
        if keep.sum() <= 1:row["anova"] = 1.
        else:
            groups = [d1, d2]; allc = np.concatenate(groups); grand = allc.mean(); k, n = 2, len(allc)
            msr = sum(g.size*(g.mean()-grand)**2 for g in groups)/(k-1)
            mse = sum(((g-g.mean())**2).sum() for g in groups)/(n-k)+ALPHA
            row["anova"] = float(stats.f.sf(msr/mse, k-1, n-k))
        diff = np.log2((t.mean(1)+1)/(c.mean(1)+1))
        row["signed_rank"] = float(stats.wilcoxon(diff).pvalue) if np.any(diff != 0) else 1.
        rows.append(row)
    for test in ("utest", "anova", "signed_rank"):
        q = bh([r[test] for r in rows])
        for r, qq in zip(rows, q):r[f"{test}_call"] = ("depleted" if r["direction"] < 0 else "enriched") if qq <= FDR else "not significant"
    return rows


def score(calls, pvals, truth):
    calls, pvals = np.asarray(calls), np.asarray(pvals); null = truth == 0
    called = calls != "not significant"; correct = (calls == np.where(truth < 0, "depleted", "enriched")) & ~null
    return dict(unaffected=int(null.sum()), p05=int((pvals[null] <= .05).sum()), calls=int(called.sum()),
                false=int((called & ~correct).sum()), affected=int((~null).sum()), hits=int(correct.sum()))


def synthetic_task(item):
    name, k = item
    raw, scaled, _, genes, idx, pos, truth = v1.make_data(name, k)
    rows = gene_tests(scaled[:, :2], scaled[:, 2:], idx)
    out = {t:score([r[f"{t}_call"] for r in rows], [r[t] for r in rows], truth) for t in ("utest", "anova", "signed_rank")}
    bgc = normalize(sparse.csc_matrix(raw), "betageom", ["synthetic"]*4)[0].toarray()
    rr, _ = resampling.analyze(bgc[:, :2], bgc[:, 2:], genes, idx, pos, resampling.Settings(recenter="trimmed-sum"))
    out["resampling"] = score([r.get("call", "not significant") if r.get("call") in ("depleted", "enriched") else "not significant" for r in rr],
                              [r.get("p_value", 1.) if r.get("in_family") else 1. for r in rr], truth)
    return name, k, out


def gale():
    design, genes, idx, _ = v1.gale_arrays("ttr", OUT/"_scratch"); s = design["scaled"]
    real = gene_tests(s[:, :2], s[:, 2:], idx)
    swap = gene_tests(s[:, [0, 2]], s[:, [1, 3]], idx)
    table = pd.DataFrame([dict(ORF=g["id"], Name=g["name"]) | r for g, r in zip(genes, real)])
    return table, {t:int(sum(r[f"{t}_call"] != "not significant" for r in swap)) for t in ("utest", "anova", "signed_rank")}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(4) as pool:res = list(pool.map(synthetic_task, [(s, k) for s in SCENARIOS for k in range(EXPERIMENTS)]))
    agg = {}
    for name, _, out in res:
        for test, sc in out.items():
            a = agg.setdefault((name, test), dict.fromkeys(sc, 0))
            for key, v in sc.items():a[key] += v
    rows = []
    for (name, test), a in agg.items():
        rows.append(dict(scenario=name, test=test, type1_p05=round(a["p05"]/a["unaffected"], 4),
                         fdr=round(a["false"]/a["calls"], 3) if a["calls"] else 0., calls=a["calls"],
                         sensitivity=round(a["hits"]/a["affected"], 3) if a["affected"] else None))
    write_table(OUT/"synthetic.tsv", list(rows[0]), rows)
    table, swap = gale(); table.to_csv(OUT/"gale_ttr_genes.tsv", sep="\t", index=False)
    v3 = pd.read_csv(ROOT/"output/essentiality/resampling/Gale_2020_bgc_v3/genes.tsv", sep="\t")
    strict = pd.read_csv(ROOT/"output/essentiality/resampling/v3_gale_prioritization/strict_published_genes.tsv", sep="\t")
    m = table.merge(v3[["ORF", "call"]].rename(columns={"call":"resampling_call"}), on="ORF", how="left")
    gale_summary = dict(label_swap_null_calls=swap | dict(resampling=6),
                        real_calls={t:m[f"{t}_call"].value_counts().to_dict() for t in ("utest", "anova", "signed_rank")} |
                                   dict(resampling=m.resampling_call.value_counts().to_dict()),
                        strict_published_depleted_recovered={t:int(m[m.ORF.isin(strict.ORF)][f"{t}_call"].eq("depleted").sum())
                                                             for t in ("utest", "anova", "signed_rank")} | dict(resampling=18, of=len(strict)))
    (OUT/"gale.json").write_text(json.dumps(gale_summary, indent=2)+"\n")
    print(pd.DataFrame(rows).pivot(index="scenario", columns="test", values=["type1_p05", "fdr", "sensitivity"]).to_string())
    print(json.dumps(gale_summary, indent=1))


if __name__ == "__main__":main()
