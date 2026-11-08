#!/usr/bin/env python3
"""DEVELOPMENT ONLY (not acceptance): does TRANSIT's genetic-interaction (GI) method suit Hermes, and is a
site-paired alternative feasible? Seeds used here (5150, 5151) are development seeds.

Design (TRANSIT GI): strain A and strain B, each with condition 1 (control) and condition 2 (treated).
Hermes: each strain is its own library (pool, its own insertion sites); control and treated are split from
that pool (shared clones), 2 replicate cultures each. GI for a gene: delta = log2FC(B) - log2FC(A).

Arms:
  transit_gi   faithful port of TRANSIT gi.py (commit 50e23d6): TTR; per group, pooled site x replicate
               counts; conjugate normal posterior of the mean, truncated at 0, priors from 1%-trimmed means
               of gene means and variances; 10,000 draws; GI when the 95% HDI of delta excludes ROPE +-0.5
               (TRANSIT default "HDI"); also "prob" (P(delta in ROPE) < 0.05).
  site_paired  candidate: BGC per pool; within each strain, treated recentered to control (iterative median);
               per replicate pair r, gene response R_r = mean over the gene's sites of
               log2((treated_ir + 5)/(control_ir + 5)); R = mean_r R_r, s^2 = var_r R_r; s^2 moderated per strain
               towards a trend in log(sites) (limma-style empirical Bayes, Smyth 2004);
               t = (R_B - R_A)/sqrt(s~^2_A/m + s~^2_B/m), Welch-Satterthwaite df; BH; GI at q <= 0.05.
Writes output/essentiality/gi/fit_check/.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import sys

import numpy as np
from scipy import sparse, stats
from scipy.special import digamma, polygamma

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import multicondition as mc
from ytab.essentiality.control_conditioned_gumbel import benjamini_hochberg
from ytab.normalization.methods import normalize
from ytab.normalization.ttr import write_table

OUT = ROOT/"output/essentiality/gi/fit_check"
SEEDS, EXPERIMENTS, GENES = (5150, 5151), 4, 300
BASE = dict(m=float(np.log(30)), tau=1.5, phi=.2, phi_B=None, depth_B=1.)
SCENARIOS = {   # effects: per gene class, (A treated fold, B treated fold, B baseline fold in both conditions)
    "null_no_effects":dict(),
    "null_shared_response":dict(shared=.25),                         # 10% genes depleted 4x by treatment in both strains: no GI
    "null_baseline_difference":dict(baseline_B=.25),                 # 10% genes 4x lower in B in both conditions: no GI
    "null_shared_response_unequal_depth":dict(shared=.25, depth_B=.3),
    "null_shared_response_noisy_B":dict(shared=.25, phi=.02, phi_B=.5),
    "gi_aggravating":dict(gi_B=.25),                                 # 10% genes depleted 4x by treatment in B only
    "gi_alleviating":dict(gi_A=.25),                                 # 10% genes depleted 4x by treatment in A only
    "gi_mixed":dict(shared=.25, gi_B=.25),                           # 10% shared (no GI) + 10% B-only (GI)
}


def simulate(name, k, seed):
    cfg = BASE | SCENARIOS[name]
    rng = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(sorted(SCENARIOS).index(name), k)))
    classes = rng.permutation(np.repeat([0, 1, 2], [GENES-2*(GENES//10), GENES//10, GENES//10]))   # 0 none, 1 first effect, 2 second
    strains = {}
    for strain in "AB":
        sites = np.clip(np.round(np.exp(rng.uniform(np.log(3), np.log(120), GENES))), 3, 120).astype(int)
        lam = np.exp(cfg["m"]+cfg["tau"]*rng.standard_normal(sites.sum())); gene_of = np.repeat(np.arange(GENES), sites)
        base = np.ones(GENES); treat = np.ones(GENES)
        first = classes == 1; second = classes == 2
        if "shared" in cfg:treat[first] = cfg["shared"]
        if strain == "B" and "baseline_B" in cfg:base[first] = cfg["baseline_B"]
        if strain == "A" and "gi_A" in cfg:treat[first] = cfg["gi_A"]
        if strain == "B" and "gi_B" in cfg:(treat.__setitem__(second, cfg["gi_B"]) if "shared" in cfg else treat.__setitem__(first, cfg["gi_B"]))
        # block effect: half of each affected gene's sites, from a quarter in
        offs = np.concatenate([np.arange(n) for n in sites]); block = (offs >= sites[gene_of]//4) & (offs < sites[gene_of]//4+np.maximum(2, sites[gene_of]//2))
        tf = np.where(block, treat[gene_of], 1.); bf = np.where(block, base[gene_of], 1.)
        phi = cfg["phi_B"] if strain == "B" and cfg["phi_B"] else cfg["phi"]; size = 1/phi
        depth = cfg["depth_B"] if strain == "B" else 1.
        draw = lambda mean:rng.negative_binomial(size, size/(size+mean))
        cols = [draw(depth*e*lam*bf) for e in (1., 1.3)]+[draw(depth*e*lam*bf*tf) for e in (.8, 1.1)]
        idx, start = [], 0
        for n in sites:idx.append(np.arange(start, start+n)); start += n
        strains[strain] = dict(raw=np.column_stack(cols), idx=idx)
    gi = (classes == 2) if ("shared" in cfg and "gi_B" in cfg) else (classes == 1) if any(x in cfg for x in ("gi_A", "gi_B")) else np.zeros(GENES, bool)
    return strains, gi


def sample_post(data, S, mu0, s20, rng, k0=1., nu0=1.):
    n = len(data); s2 = np.var(data, ddof=1); ybar = np.mean(data); kn, nun = k0+n, nu0+n
    mun = (k0*mu0+n*ybar)/kn; s2n = (nu0*s20+(n-1)*s2+(k0*n/kn)*(ybar-mu0)**2)/nun
    s2_post = 1/stats.gamma.rvs(nun/2, scale=2/(s2n*nun), size=S, random_state=rng)
    sd = np.sqrt(s2_post/kn); a, b = (0-mun)/sd, (1e6-mun)/sd
    return stats.truncnorm.rvs(a, b, loc=mun, scale=sd, size=S, random_state=rng)


def hdi(x, mass=.95):
    x = np.sort(x); m = int(np.ceil(mass*len(x))); w = x[m:]-x[:len(x)-m]; i = int(np.argmin(w)); return x[i], x[i+m]


def transit_gi(strains, rng, S=10000, rope=.5):
    groups = {}
    for s in "AB":
        v = normalize(sparse.csc_matrix(strains[s]["raw"]), "ttr", [s]*4)[0].toarray()
        groups[s+"1"] = [v[i][:, :2].ravel() for i in strains[s]["idx"]]; groups[s+"2"] = [v[i][:, 2:].ravel() for i in strains[s]["idx"]]
    prior = {g:(stats.trim_mean([x.mean() for x in d], .01), stats.trim_mean([x.var() for x in d], .01)) for g, d in groups.items()}
    hdi_call, prob_call = [], []
    for j in range(GENES):
        post = {g:sample_post(groups[g][j], S, *prior[g], rng) for g in groups}
        delta = np.log2(post["B2"]/post["B1"])-np.log2(post["A2"]/post["A1"])
        lo, hi = hdi(delta); hdi_call.append(lo > rope or hi < -rope)
        prob_call.append(np.mean((delta >= -rope) & (delta <= rope)) < .05)
    return np.array(hdi_call), np.array(prob_call)


def moderate(s2, d, covariate):
    """limma-style moderated variances with a prior trend in `covariate` (Smyth 2004, method of moments)."""
    ok = s2 > 0; e = np.log(s2[ok])-digamma(d/2)+np.log(d/2)
    X = np.column_stack([np.ones(ok.sum()), covariate[ok]]); beta = np.linalg.lstsq(X, e, rcond=None)[0]
    resid = e-X@beta; excess = np.mean(resid**2)*ok.sum()/(ok.sum()-2)-polygamma(1, d/2)
    if excess <= 0:d0 = np.inf
    else:
        lo, hi = 1e-3, 1e6                                   # solve trigamma(d0/2) = excess
        for _ in range(200):
            mid = np.sqrt(lo*hi); (lo, hi) = (mid, hi) if polygamma(1, mid/2) > excess else (lo, mid)
        d0 = mid
    s0 = np.exp(np.column_stack([np.ones(len(s2)), covariate])@beta+(digamma(d0/2)-np.log(d0/2) if np.isfinite(d0) else 0))
    post = s0 if not np.isfinite(d0) else (d0*s0+d*s2)/(d0+d)
    return post, d0


def site_paired(strains, pc=5.):
    R, S2, nsites = {}, {}, {}
    for s in "AB":
        v = normalize(sparse.csc_matrix(strains[s]["raw"]), "betageom", [s]*4)[0].toarray(); idx = strains[s]["idx"]
        v, _ = mc.recenter(v, np.array([0, 0, 1, 1]), 0, idx, "iterative-median")
        rep = []
        for c, t in ((0, 2), (1, 3)):
            rep.append([np.mean(np.log2((v[i, t]+pc)/(v[i, c]+pc))[v[i].sum(1) > 0]) if (v[i].sum(1) > 0).any() else np.nan for i in idx])
        rep = np.array(rep); R[s] = rep.mean(0); S2[s] = rep.var(0, ddof=1); nsites[s] = np.array([(v[i].sum(1) > 0).sum() for i in idx])
    ok = np.isfinite(R["A"]) & np.isfinite(R["B"]); p = np.ones(GENES); t = np.zeros(GENES)
    var = {}
    for s in "AB":
        post, d0 = moderate(S2[s][ok], 1., np.log(np.maximum(nsites[s][ok], 1))); var[s] = (post/2, d0+1)
    se2 = var["A"][0]+var["B"][0]
    df = se2**2/(var["A"][0]**2/min(var["A"][1], 1e6)+var["B"][0]**2/min(var["B"][1], 1e6))
    t[ok] = (R["B"][ok]-R["A"][ok])/np.sqrt(se2); p[ok] = 2*stats.t.sf(np.abs(t[ok]), df)
    q = benjamini_hochberg(p); return p, q <= .05


def task(item):
    name, k, seed = item
    strains, gi = simulate(name, k, seed); rng = np.random.default_rng(seed+k)
    h, pr = transit_gi(strains, rng); p, call = site_paired(strains)
    score = lambda called:dict(null_genes=int((~gi).sum()), false=int((called & ~gi).sum()), calls=int(called.sum()), hits=int((called & gi).sum()), gi_genes=int(gi.sum()))
    return name, dict(transit_hdi=score(h), transit_prob=score(pr), site_paired=score(call) | dict(p05_null=int((p[~gi] <= .05).sum())))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(4) as pool:res = list(pool.map(task, [(s, k, seed) for seed in SEEDS for s in SCENARIOS for k in range(EXPERIMENTS)]))
    agg = {}
    for name, out in res:
        for arm, sc in out.items():
            a = agg.setdefault((name, arm), {})
            for key, v in sc.items():a[key] = a.get(key, 0)+v
    rows = [dict(scenario=n, arm=arm, null_call_rate=round(a["false"]/a["null_genes"], 4), type1_p05=round(a["p05_null"]/a["null_genes"], 4) if "p05_null" in a else None,
                 calls=a["calls"], fdr=round(a["false"]/a["calls"], 3) if a["calls"] else 0., sensitivity=round(a["hits"]/a["gi_genes"], 3) if a["gi_genes"] else None)
            for (n, arm), a in agg.items()]
    write_table(OUT/"dev_results.tsv", list(rows[0]), rows)
    for r in rows:print(r)


if __name__ == "__main__":main()
