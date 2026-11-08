#!/usr/bin/env python3
"""DEVELOPMENT ONLY (v4 candidates; not acceptance). Seeds never used for evaluation.

Candidates on the v1 generator + the v3 held-out scenarios:
  bgc_rc        BGC (per sample) + trimmed-sum recentering, exact mean-difference permutation (v3 reference)
  bgcj_rc       joint BGC per condition (one beta-geometric fit over a condition's replicates) + recentering
  bgc_rc_stud   BGC + recentering, studentized site-restricted permutation (Welch-type statistic with per-site
                within-condition variances), Monte Carlo p = (b+1)/(S+1)
Usage: ytab_resampling_v4_dev.py TAG SEED EXPERIMENTS [--samples S] [--workers N]
"""
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import importlib.util
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
from scipy.stats import geom

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import resampling
from ytab.essentiality.control_conditioned_gumbel import benjamini_hochberg
from ytab.normalization.betageom import beta_geometric_cdf, fit_parameters
from ytab.normalization.methods import normalize
from ytab.normalization.ttr import ttr_parameters, write_table

spec = importlib.util.spec_from_file_location("v1", ROOT/"scripts/local/ytab_resampling_benchmark.py")
v1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(v1)
v1.SCENARIOS.update({"null_noisy_control":dict(phi=.5, treatment_phi=.02), "enrichment_4x_noisy_treatment":dict(fold=4., phi=.02, treatment_phi=.5)})
OUT = ROOT/"output/essentiality/resampling/v4_development"
SAMPLES = int(sys.argv[sys.argv.index("--samples")+1]) if "--samples" in sys.argv else 2000


def joint_bgc(raw, conditions, trim=.05, target=100):
    """One beta-geometric fit per condition (pooled positive counts of its replicates); each sample mapped through it, then TTR."""
    out = np.zeros(raw.shape, dtype=float); means = raw.sum(0)/(raw > 0).sum(0); target_mean = float(means.mean())
    for c in set(conditions):
        cols = [j for j, x in enumerate(conditions) if x == c]
        fit = fit_parameters(np.concatenate([raw[raw[:, j] > 0, j] for j in cols]).astype(float), raw.shape[0]*len(cols))
        for j in cols:
            pos = raw[:, j] > 0; prob = np.clip(beta_geometric_cdf(raw[pos, j], fit["alpha"], fit["beta"]), np.nextafter(0., 1.), np.nextafter(1., 0.))
            values = np.maximum(geom.ppf(prob, 1/target_mean), 1.)
            out[pos, j] = values*ttr_parameters(values, raw.shape[0], trim, target)["factor"]
    return out


SPLITS = list(combinations(range(4), 2))   # treatment label pairs among 4 samples; observed (2,3) is last


def studentized_p(A, B, rng, samples):
    """Site-restricted permutation of z = (mean_B - mean_A)/sqrt(sum_i (vB_i/2 + vA_i/2))/n, two-sided, Monte Carlo."""
    V = np.concatenate([A, B], 1); n = len(V)
    sums = np.stack([V[:, list(s)].sum(1) for s in SPLITS], 1); tot = V.sum(1)[:, None]
    sq = np.stack([((V[:, list(s)]-V[:, list(s)].mean(1, keepdims=True))**2).sum(1) for s in SPLITS], 1)
    sqA = np.stack([((V[:, [k for k in range(4) if k not in s]]-V[:, [k for k in range(4) if k not in s]].mean(1, keepdims=True))**2).sum(1) for s in SPLITS], 1)
    def z(idx):
        rows = np.arange(n)   # idx: draws x sites, the chosen split per site
        b = sums[rows, idx].sum(-1); a = tot.sum()-b
        var = (sq[rows, idx].sum(-1)+sqA[rows, idx].sum(-1))/2/n**2   # per-site sample variances (1 df each) / 2 reps
        t = (b-a)/(2*n); return np.where(var > 0, t/np.sqrt(np.maximum(var, 1e-300)), 0.)
    observed = float(z(np.full((1, n), len(SPLITS)-1))[0])
    draws = z(rng.integers(0, len(SPLITS), size=(samples, n)))
    return (np.sum(np.abs(draws) >= abs(observed)-1e-12)+1)/(samples+1), observed


def task(item):
    name, k, seed = item; v1.SEED = seed
    raw, scaled, factors, genes, idx, pos, truth = v1.make_data(name, k)
    bgc = normalize(sparse.csc_matrix(raw), "betageom", ["synthetic"]*4)[0].toarray(); bgcj = joint_bgc(raw, ["c", "c", "t", "t"])
    out = {}
    for vid, m in (("bgc_rc", bgc), ("bgcj_rc", bgcj)):
        rows, _ = resampling.analyze(m[:, :2], m[:, 2:], genes, idx, pos, resampling.Settings(recenter="trimmed-sum"))
        out[vid] = v1.score(rows, truth)
    rows, diag = resampling.analyze(bgc[:, :2], bgc[:, 2:], genes, idx, pos, resampling.Settings(recenter="trimmed-sum"))
    A, B = bgc[:, :2], bgc[:, 2:]*diag["treatment_rescale"]; rng = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(k, 99)))
    family = [(r, i) for r, i in zip(rows, idx) if r["in_family"]]
    stud = []
    for r, i in family:
        i = np.asarray(i); i = i[(A[i].sum(1)+B[i].sum(1)) > 0]
        p, zobs = studentized_p(A[i], B[i], rng, SAMPLES); stud.append((p, zobs))
    q = benjamini_hochberg([p for p, _ in stud]); srows = []
    for (r, _), (p, zobs), qq in zip(family, stud, q):
        srows.append(dict(r, p_value=p, q_bh=qq, delta_mean=zobs, call="not significant" if qq > .05 else ("depleted" if zobs < 0 else "enriched")))
    fam_ids = {r["ORF"] for r, _ in family}; srows += [r for r in rows if r["ORF"] not in fam_ids]
    order = {g["id"]:n for n, g in enumerate(genes)}; srows.sort(key=lambda r:order[r["ORF"]])
    out["bgc_rc_stud"] = v1.score(srows, truth)
    return name, k, out


def main():
    tag, seed, experiments = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    workers = int(sys.argv[sys.argv.index("--workers")+1]) if "--workers" in sys.argv else 6
    with ProcessPoolExecutor(workers) as pool:results = list(pool.map(task, [(n, k, seed) for n in v1.SCENARIOS for k in range(experiments)]))
    rows = []
    variants = list(results[0][2])
    for vid in variants:
        for n in v1.SCENARIOS:
            m = v1.aggregate([r[2][vid] for r in results if r[0] == n])
            rows.append(dict(variant=vid, scenario=n, p05=round(m["p05_rate"], 4), ci_low=round(m["p05_ci"][0], 4), calls=m["calls"], false=m["false"],
                             exp_false=m["experiments_with_false_call"], sens=None if m["sensitivity"] is None else round(m["sensitivity"], 3),
                             fdr=None if m["fdr"] is None else round(m["fdr"], 3)))
    OUT.mkdir(parents=True, exist_ok=True); write_table(OUT/f"{tag}.tsv", list(rows[0]), rows)
    nulls = [n for n in v1.SCENARIOS if n.startswith("null_")]; signals = [n for n in v1.SCENARIOS if n not in nulls]
    for vid in variants:
        r = [x for x in rows if x["variant"] == vid]; worst = max((x for x in r if x["scenario"] in nulls), key=lambda x:x["p05"])
        fdr = sum(x["false"] for x in r if x["scenario"] in signals)/max(sum(x["calls"] for x in r if x["scenario"] in signals), 1)
        print(f"{vid:14s} worst-null p05 {worst['p05']:.3f} ({worst['scenario']}, CI low {worst['ci_low']:.3f})  pooled FDR {fdr:.3f}  "
              f"null exps w/ false {sum(x['exp_false'] for x in r if x['scenario'] in nulls)}/{len(nulls)*experiments}  "
              f"mean sens {np.mean([x['sens'] or 0 for x in r if x['scenario'] in signals]):.3f}")


if __name__ == "__main__":main()
