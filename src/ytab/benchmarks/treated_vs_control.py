"""Treated-vs-control benchmark: Gale-calibrated simulation, method x normalization arms, metrics and paired tests.

Design: docs/manuscript/preregistration_treated_vs_control.md. Steps:
- calibrate(): site abundances, library exposures and abundance-dependent negative-binomial dispersions from real
  control and treated libraries (one project, two control/treated pairs);
- simulate(): one experiment (2 controls, 2 treated) with planted gene effects and gene-level truth;
- run_arms(): every method x normalization arm on one experiment -> per-gene ranking score and call;
- score(): AUC, empirical FDR, power, precision at K, null calibration;
- paired_tests(): sign-flip permutation and bootstrap interval of paired arm differences, Holm within families.
Species-agnostic: any project with two control/treated pairs and a GFF can be the calibration base.
"""
from __future__ import annotations

from dataclasses import dataclass
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import roc_auc_score

from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.essentiality import resampling
from ytab.fitness import zscore
from ytab.normalization.methods import normalize

SCENARIOS = {
    "null": dict(fraction=0.),
    "depletion": dict(fraction=.05, kinds=("dep4", "dep2")),
    "enrichment": dict(fraction=.05, kinds=("enr4", "enr2")),
    "mixed": dict(fraction=.05, kinds=("dep4", "dep2", "enr4", "enr2")),
    "complete_loss": dict(fraction=.03, kinds=("loss",)),
    "widespread": dict(fraction=.20, kinds=("dep4", "dep2", "enr4", "enr2")),
    "widespread_asymmetric": dict(fraction=.20, kinds=("enr4", "enr2")*4+("dep4", "dep2")),   # 80% enriched: composition shift
    "treated_jackpots": dict(fraction=.05, kinds=("dep4", "dep2", "enr4", "enr2"), jackpots=(.001, 50.)),
}
FOLDS = {"dep4":.25, "dep2":.5, "enr4":4., "enr2":2., "loss":0.}
ARMS = {  # arm -> (method, normalization, settings)
    "zscore_cpm": ("zscore", "cpm", {}), "zscore_ttr": ("zscore", "ttr", {}), "zscore_bgc": ("zscore", "betageom", {}),
    "resampling_totreads": ("resampling", "totreads", {}), "resampling_ttr": ("resampling", "ttr", {}),
    "resampling_bgc": ("resampling", "betageom", dict(recenter="trimmed-sum")),
    "gumbel_totreads": ("gumbel", "totreads", {}), "gumbel_ttr": ("gumbel", "ttr", {}),
    # secondary (descriptive): the earlier per-pair z-score call rule
    "zscore_cpm_pairs": ("zscore", "cpm", dict(combine="pairs")), "zscore_ttr_pairs": ("zscore", "ttr", dict(combine="pairs")),
    "zscore_bgc_pairs": ("zscore", "betageom", dict(combine="pairs")),
}
BINS = 20


@dataclass
class Calibration:
    sites: list            # (chromosome, position) per site
    abundance: np.ndarray  # lambda_s (mean depth-scaled control count)
    exposures: np.ndarray  # controls then treated, relative to the control mean
    control_phi: tuple     # (log mean bin centres, log phi) for interpolation
    treated_phi: tuple


def _dispersion(x1, x2):
    """Abundance-dependent NB dispersion phi(mu) from two replicate tracks (depth-scaled), method of moments in
    BINS abundance quantile bins: Var = mu + phi mu^2."""
    mu = (x1+x2)/2; keep = mu > 0; mu, var = mu[keep], ((x1-x2)**2/2)[keep]
    edges = np.unique(np.quantile(mu, np.linspace(0, 1, BINS+1)))
    centres, phis = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        b = (mu >= lo) & (mu <= hi)
        if b.sum() < 50:continue
        m, v = mu[b].mean(), var[b].mean()
        centres.append(np.log(m)); phis.append(np.log(max((v-m)/m**2, 1e-4)))
    return np.array(centres), np.array(phis)


def calibrate(raw: sparse.spmatrix, sites: list, control_cols: list[int], treated_cols: list[int]) -> Calibration:
    raw = sparse.csc_matrix(raw); cols = control_cols+treated_cols
    counts = {j:np.asarray(raw[:, j].todense()).ravel().astype(float) for j in cols}
    totals = np.array([counts[j].sum() for j in cols]); scale = totals[:len(control_cols)].mean()
    scaled = {j:counts[j]/t*scale for j, t in zip(cols, totals)}
    abundance = np.mean([scaled[j] for j in control_cols], axis=0)
    keep = abundance > 0
    c1, c2 = control_cols[:2]; t1, t2 = treated_cols[:2]
    return Calibration(sites=[s for s, k in zip(sites, keep) if k], abundance=abundance[keep], exposures=totals/scale,
                       control_phi=_dispersion(scaled[c1][keep], scaled[c2][keep]),
                       treated_phi=_dispersion(scaled[t1][keep], scaled[t2][keep]))


def _phi(curve, mu):
    x, y = curve
    return np.exp(np.interp(np.log(np.maximum(mu, 1e-9)), x, y))


def _nb(rng, mu, phi):
    size = 1/phi
    return rng.negative_binomial(size, size/(size+mu))


def simulate(cal: Calibration, gene_indices: list, scenario: str, rng: np.random.Generator):
    """Raw counts (sites x 4: control 1, control 2, treated 1, treated 2) and gene truth (+1, -1, 0)."""
    spec = SCENARIOS[scenario]; n_genes = len(gene_indices)
    effect = np.ones(len(cal.abundance)); truth = np.zeros(n_genes, dtype=int)
    eligible = [g for g, idx in enumerate(gene_indices) if len(idx) >= 3]
    k = int(round(spec["fraction"]*len(eligible)))
    if k:
        chosen = rng.choice(eligible, k, replace=False); kinds = spec["kinds"]
        for i, g in enumerate(chosen):
            kind = kinds[i % len(kinds)]; idx = np.asarray(sorted(gene_indices[g]))
            if i % 2:                      # contiguous block of half the gene's sites
                length = max(2, len(idx)//2); start = rng.integers(0, len(idx)-length+1); idx = idx[start:start+length]
            effect[idx] = FOLDS[kind]; truth[g] = 1 if FOLDS[kind] > 1 else -1
    columns = []
    for j, e in enumerate(cal.exposures[:4]):
        treated = j >= 2
        mu = cal.abundance*e*(effect if treated else 1.)
        draw = _nb(rng, mu, _phi(cal.treated_phi if treated else cal.control_phi, np.maximum(mu, 1e-9))) if mu.any() else np.zeros_like(mu)
        if treated and spec.get("jackpots"):
            share, factor = spec["jackpots"]; hit = rng.random(len(mu)) < share
            draw = draw + np.where(hit, _nb(rng, mu*(factor-1), _phi(cal.treated_phi, np.maximum(mu*(factor-1), 1e-9))), 0)
        columns.append(draw)
    return np.column_stack(columns).astype(np.int64), truth


def _gene_totals(matrix, gene_indices, names):
    m = np.asarray(matrix)
    return pd.DataFrame({c:[m[idx, j].sum() if len(idx) else 0. for idx in gene_indices] for j, c in enumerate(["c1", "c2", "t1", "t2"])}, index=names)


DESIGN = ("comparison_id,pair_id,parent_sample,treated_sample,background,pool,parent_condition,treated_condition,include,notes\n"
          "p1,p1,c1,t1,sim,1,parent,treated,true,\np2,p2,c2,t2,sim,2,parent,treated,true,\n")


def _zscore(values: pd.DataFrame, settings: dict | None = None):
    """Ranking score: the combined z (replicates averaged) or the mean per-pair z, matching the call rule."""
    st = zscore.Settings(**(settings or {}))
    with tempfile.TemporaryDirectory() as tmp:
        design = Path(tmp)/"design.csv"; design.write_text(DESIGN)
        values = values.rename_axis("feature_id")
        zscore.run(design, Path(tmp)/"out", normalized=dict(genes=values, label="benchmark"), settings=st)
        s = pd.read_csv(Path(tmp)/"out/tables/treated_vs_parent.summary_by_feature.csv").set_index("feature_id")
    call = s.final_call.astype(str)
    direction = np.where(call.str.endswith("enriched"), 1, np.where(call.str.endswith("depleted"), -1, 0))
    score = s.combined_z if st.combine == "replicates" else s.mean_z
    return pd.DataFrame(dict(score=score.fillna(0.), call=direction, p=np.nan), index=s.index)


def run_arms(raw, cal, genes, gene_indices_r, gene_indices_g, positions, arms=ARMS):
    """Per arm: DataFrame indexed by gene ID with signed ranking score, signed call (+1/-1/0) and p-value."""
    names = [g["id"] for g in genes]; out = {}; cache = {}
    groups = ["sim"]*4
    for arm, (method, norm, settings) in arms.items():
        if method == "zscore":
            if norm == "cpm":
                totals = _gene_totals(raw, gene_indices_r, names); values = totals/totals.sum()*1e6
            else:
                if norm not in cache:cache[norm] = normalize(sparse.csc_matrix(raw), norm, groups)
                values = _gene_totals(cache[norm][0].toarray(), gene_indices_r, names)
            out[arm] = _zscore(values, settings).reindex(names).fillna(dict(score=0., call=0))
        elif method == "resampling":
            if norm not in cache:cache[norm] = normalize(sparse.csc_matrix(raw), norm, groups)
            m = cache[norm][0].toarray()
            rows, _ = resampling.analyze(m[:, :2], m[:, 2:], genes, gene_indices_r, positions, resampling.Settings(**settings))
            out[arm] = _rows(rows, signed=True)
        else:
            if norm not in cache:cache[norm] = normalize(sparse.csc_matrix(raw), norm, groups)
            factors = np.array([r["factor"] for r in cache[norm][1]], dtype=float)
            inverse = 1/factors; exposures = inverse/inverse[:2].mean()
            rows, _, _ = ccg.analyze(raw[:, :2], raw[:, 2:], exposures, genes, gene_indices_g, positions, ccg.Settings(**settings))
            out[arm] = _rows(rows, signed=False)
    return out


def _rows(rows, signed):
    d = pd.DataFrame(rows).set_index("ORF")
    p = pd.to_numeric(d.p_value, errors="coerce")
    tested = d.in_family.astype(bool) & p.notna()
    lp = -np.log10(np.clip(p.where(tested, 1.), 1e-300, 1))
    if signed:
        sign = np.sign(pd.to_numeric(d.log2FC, errors="coerce").fillna(0.)); score = lp*sign
        call = np.where(d.call == "enriched", 1, np.where(d.call == "depleted", -1, 0))
    else:
        score = -lp; call = np.where(d.call == "depleted run", -1, 0)   # depletion only: lower = more depleted
    return pd.DataFrame(dict(score=score.astype(float), call=call, p=p.where(tested)), index=d.index)


def score(result: pd.DataFrame, truth: pd.Series, depletion_only=False):
    """Metrics of one arm on one experiment. truth: +1/-1/0 per gene (same index)."""
    r = result.reindex(truth.index).fillna(dict(score=0., call=0)); t = truth.to_numpy(); s = r.score.to_numpy(); c = r.call.to_numpy()
    out = {}
    aucs = []
    for direction in ((-1,) if depletion_only else (-1, 1)):
        keep = (t == direction) | (t == 0)
        if (t[keep] == direction).any() and (t[keep] == 0).any():
            auc = roc_auc_score(t[keep] == direction, s[keep]*direction); out[f"auc_{'depletion' if direction < 0 else 'enrichment'}"] = auc; aucs.append(auc)
    out["auc"] = float(np.mean(aucs)) if aucs else np.nan
    affected = (t != 0) if not depletion_only else (t == -1)
    calls = c != 0; correct = calls & (c == t)
    out |= dict(calls=int(calls.sum()), false_calls=int((calls & ~correct).sum()), affected=int(affected.sum()),
                correct_calls=int(correct.sum()), power=correct.sum()/affected.sum() if affected.any() else np.nan,
                fdp=(calls & ~correct).sum()/max(calls.sum(), 1))   # false discovery proportion; 0 without calls
    k = int(affected.sum())
    if k:
        rank_score = -s if depletion_only else np.abs(s)
        top = np.argsort(-rank_score, kind="stable")[:k]
        hit = (t[top] != 0) & (np.sign(s[top]) == np.sign(t[top])) if not depletion_only else (t[top] == -1)
        out["precision_at_k"] = float(hit.mean())
    p = r.p.to_numpy(dtype=float); unaffected = (t == 0) & np.isfinite(p)
    out |= dict(unaffected_tested=int(unaffected.sum()), unaffected_p05=int((p[unaffected] <= .05).sum()))
    return out


def sign_flip(diff: np.ndarray, flips=100_000, seed=0):
    """Two-sided paired sign-flip permutation p-value of mean(diff)."""
    rng = np.random.default_rng(seed); diff = diff[np.isfinite(diff)]
    obs = abs(diff.mean()); signs = rng.choice((-1., 1.), size=(flips, len(diff)))
    return float(((np.abs((signs*diff).mean(axis=1)) >= obs-1e-15).sum()+1)/(flips+1))


def bootstrap_ci(diff: np.ndarray, draws=10_000, seed=0):
    rng = np.random.default_rng(seed); diff = diff[np.isfinite(diff)]
    means = diff[rng.integers(0, len(diff), (draws, len(diff)))].mean(axis=1)
    return tuple(np.percentile(means, [2.5, 97.5]))


def holm(pvalues: dict) -> dict:
    order = sorted(pvalues, key=pvalues.get); m = len(order); out, running = {}, 0.
    for i, k in enumerate(order):
        running = max(running, min(1., (m-i)*pvalues[k])); out[k] = running
    return out
