#!/usr/bin/env python3
"""Pre-registered evaluation of the window HMM v3 (hermes_hmm_v3: learned target-site covariate) on Gale 2020 Cg1-Cg3.

blind:   label-free gates of benchmark_v3/PREREGISTRATION.md: (1) the learned target model rediscovers Hermes'
         nTnnnnAn (two most informative offsets 5 bp apart, T then A); (2) determinism; (3) EM monotone and
         converged; (4) synthetic recovery (v2's gates, seed 60331); (5) agreement with v2 on Gale (Spearman of
         P_ES_run >= 0.90, kappa of ES calls >= 0.80). Writes decision.json. Reads no labels.
unblind: only after a pass; non-inferiority to v2 against S. cerevisiae ortholog (silver) labels: the 95% paired
         bootstrap CI (2000 resamples, seed 60331) of AUC_v3 - AUC_v2 must lie above -0.01.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src")); sys.path.append(str(ROOT/"scripts/local"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.essentiality import hermes_hmm as hh
from ytab.normalization.ttr import sha256, write_table
import ytab_hermes_hmm as runner

BENCH = ROOT/"output/essentiality/hermes_hmm/benchmark_v3"
V2 = ROOT/"output/essentiality/hermes_hmm/benchmark_v2/gale_blinded/joint_run1/genes.tsv"
SEED = 60331
GENE_STATE_P = dict(ES=.20, GD=.05, NE=.73, GA=.02)      # as v2: generic, pre-declared; not fitted to any labels
SCENARIOS = {"well_specified":{}, "accessibility_domains":dict(access_sd=.5, block=10), "partial_essential":dict(partial=.5),
             "low_complexity":dict(scale=1/3), "lognormal_poisson":dict(lognormal=True)}
GATED = ("well_specified", "accessibility_domains", "partial_essential")


def check_preregistration():
    frozen = (BENCH/"PREREGISTRATION.sha256").read_text().split()
    if sha256(BENCH/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


_TARGET = {}


def target_model(settings):
    """The learned target-site model of the Gale parents (deterministic; learned once per process)."""
    if "m" not in _TARGET:_TARGET["m"] = runner.learn_target(runner.GALE, [None]*len(runner.GALE), settings, "Gale_2020", runner.FASTA)
    return _TARGET["m"]


def context(settings):
    genes = ccg.annotations(runner.GFF); feats = runner.features(settings, genes, target=target_model(settings))
    libs = [runner.hit_sites(s, settings.min_reads) for s in runner.GALE]
    chroms, scaling = runner.design(feats, libs, settings)
    return genes, chroms


def truth_windows(genes, chroms, rng, partial=None):
    """Gene states drawn independently; trimmed-CDS windows take the gene state (central fraction only if partial)."""
    names = [c["name"] for c in chroms]; states = {c["name"]:np.full(len(c["Y"]), 2, dtype=int) for c in chroms}
    p = np.array(list(GENE_STATE_P.values())); truth = {}
    for g in genes:
        if g["chrom"] not in states:continue
        s = int(rng.choice(4, p=p)); truth[g["id"]] = s
        ws = [w for w, _ in hh.gene_windows(g, settings_global.window, settings_global)]
        if partial and s in (0, 1) and len(ws) > 2:
            n = max(1, int(round(partial*len(ws)))); lo = (len(ws)-n)//2; ws = ws[lo:lo+n]
        for w in ws:states[g["chrom"]][w] = s
    return states, truth


def simulate(genes, chroms, fitted, scenario, rng):
    cfg = SCENARIOS[scenario]; states, truth = truth_windows(genes, chroms, rng, cfg.get("partial"))
    out = []
    for c in chroms:
        st = states[c["name"]]; logmu = fitted["a"][None, :]+(c["X"]@fitted["b"])[:, None]+np.log(fitted["r"][st])[:, None]
        if "scale" in cfg:logmu += np.log(cfg["scale"])
        if "access_sd" in cfg:
            n = len(st); blocks = rng.normal(0, cfg["access_sd"], n//cfg["block"]+1); logmu += np.repeat(blocks, cfg["block"])[:n, None]
        mu = np.exp(logmu); psi = fitted["psi"]
        if cfg.get("lognormal"):   # shared lognormal accessibility with the same variance instead of gamma
            u = np.exp(rng.normal(-.5*np.log1p(psi), np.sqrt(np.log1p(psi)), (len(mu), 1)))
        else:
            u = rng.gamma(1/psi, psi, (len(mu), 1))
        Y = rng.poisson(mu*u)
        out.append(dict(name=c["name"], X=c["X"], Y=Y.astype(float), keep=c["keep"]))
    return out, truth


def gene_scores(genes, chroms, fitres):
    G, offsets = fitres["gamma"], np.cumsum([0]+[len(c["Y"]) for c in chroms])
    gamma = {c["name"]:G[offsets[i]:offsets[i+1]] for i, c in enumerate(chroms)}; keep = {c["name"]:c["keep"] for c in chroms}
    chain = {c["name"]:fitres["chains"][i] for i, c in enumerate(chroms)}
    return hh.gene_table(genes, set(gamma), gamma, keep, settings_global, chain, fitres["A"])


def auc(score, label):
    score, label = np.asarray(score, float), np.asarray(label, bool)
    order = np.argsort(score); ranks = np.empty(len(score)); ranks[order] = np.arange(1, len(score)+1)
    # average ranks for ties
    _, inv, counts = np.unique(score, return_inverse=True, return_counts=True); sums = np.bincount(inv, ranks); ranks = (sums/counts)[inv]
    n1, n0 = label.sum(), (~label).sum(); return float((ranks[label].sum()-n1*(n1+1)/2)/(n1*n0))


def recovery(rows, truth):
    r = [x for x in rows if x["status"] == "called" and x["ORF"] in truth]
    es = np.array([truth[x["ORF"]] == 0 for x in r]); p = np.array([x["P_ES_run"] for x in r]); called = np.array([x["call"] == "ES" for x in r])
    hi = p >= .9
    return dict(genes=len(r), true_es=int(es.sum()), auc=auc(p, es), es_calls=int(called.sum()),
                precision=float(es[called].mean()) if called.any() else None, recall=float(called[es].mean()),
                high_confidence=int(hi.sum()), high_confidence_true=float(es[hi].mean()) if hi.any() else None,
                state_accuracy=float(np.mean([hh.STATES.index(x["call"]) == truth[x["ORF"]] for x in r])))


def synthetic_task(item):
    scenario, k, fitted = item
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(list(SCENARIOS).index(scenario), k)))
    genes, chroms = context(settings_global); sim, truth = simulate(genes, chroms, fitted, scenario, rng)
    fitres = hh.fit(sim, settings_global); rows = gene_scores(genes, sim, fitres)
    rec = recovery(rows, truth)
    return scenario, k, rec | dict(r=fitres["r"].tolist(), b=fitres["b"].tolist(), iterations=fitres["iterations"])


def kappa(a, b):
    a, b = np.asarray(a, bool), np.asarray(b, bool); po = np.mean(a == b); pe = a.mean()*b.mean()+(1-a.mean())*(1-b.mean())
    return float((po-pe)/(1-pe))


settings_global = hh.Settings()


def rediscovery(model):
    """Gate 1: the two most informative offsets are 5 bp apart, preferring T at the first and A at the second."""
    order = np.argsort(model["information"])[::-1][:2]; k1, k2 = sorted(order)
    j1, j2 = int(model["offsets"][k1]), int(model["offsets"][k2])
    b1, b2 = "ACGT"[int(np.argmax(model["logodds"][k1, :4]))], "ACGT"[int(np.argmax(model["logodds"][k2, :4]))]
    return dict(offsets=[j1, j2], bases=[b1, b2], information_bits=[float(model["information"][k1]), float(model["information"][k2])],
                passed=bool(j2-j1 == 5 and b1 == "T" and b2 == "A"))


def blind(workers):
    import pandas as pd
    from scipy.stats import spearmanr
    prereg = check_preregistration()
    if (BENCH/"decision.json").exists():raise SystemExit("decision.json exists; phase 1 already ran")
    model = target_model(settings_global); redis = rediscovery(model)
    runs = [runner.run(BENCH/"gale_blinded"/f"joint_run{i}", runner.GALE, settings_global, f"joint run {i}", sequence_covariate="learned", target=model)
            for i in (1, 2)]
    (s1, rows1), (s2, rows2) = runs
    identical = [(r["ORF"], r.get("P_ES_run"), r["call"]) for r in rows1] == [(r["ORF"], r.get("P_ES_run"), r["call"]) for r in rows2]
    fitted = dict(a=np.array(s1["a"]), b=np.array(s1["b"]), r=np.array(s1["r"]), psi=float(s1["psi"]))
    with ProcessPoolExecutor(workers) as pool:
        synth = list(pool.map(synthetic_task, [(s, 0, fitted) for s in SCENARIOS]))
    rec = {s:r for s, _, r in synth}
    write_table(BENCH/"synthetic.tsv", ["scenario"]+[k for k in rec["well_specified"] if k not in ("r", "b")],
                [dict(scenario=s) | {k:v for k, v in r.items() if k not in ("r", "b")} for s, r in rec.items()])
    g_syn = all(rec[s]["auc"] >= .95 and (rec[s]["precision"] or 0) >= .90 and rec[s]["recall"] >= .80 and
                (rec[s]["high_confidence_true"] or 0) >= .90 for s in GATED)
    v3 = pd.DataFrame(rows1); v2 = pd.read_csv(V2, sep="\t")
    both = v3[v3.status == "called"][["ORF", "P_ES_run", "call"]].merge(v2[v2.status == "called"][["ORF", "P_ES_run", "call"]], on="ORF", suffixes=("_v3", "_v2"))
    rho = float(spearmanr(both.P_ES_run_v3, both.P_ES_run_v2).correlation)
    k = kappa(both.call_v3 == "ES", both.call_v2 == "ES")
    agreement = dict(genes=len(both), spearman=rho, kappa=k, es_v3=int((both.call_v3 == "ES").sum()), es_v2=int((both.call_v2 == "ES").sum()))
    gates = dict(g1_rediscovery=redis["passed"], g2_determinism=identical, g3_em=bool(s1["monotone"] and s1["converged"]),
                 g4_synthetic=g_syn, g5_agreement_with_v2=bool(rho >= .90 and k >= .80))
    decision = dict(status="success", phase="blind", preregistration=prereg, passed=all(gates.values()), gates=gates,
                    rediscovery=redis, agreement_with_v2=agreement,
                    gale_joint=dict(gene_calls=s1["gene_calls"], r=s1["r"], b=s1["b"], psi=s1["psi"], A=s1["A"], masked_fraction=s1["masked_fraction"],
                                    iterations=s1["iterations"], covariates=s1["covariate_scaling"]["covariates"]),
                    synthetic=rec, blinding="no essentiality labels, RF output or Gale calls were read",
                    implementation_sha256=sha256(Path(hh.__file__)), runner_sha256=sha256(Path(runner.__file__)), script_sha256=sha256(Path(__file__)))
    (BENCH/"decision.json").write_text(json.dumps(hh._json_ready(decision), indent=2)+"\n")
    print(json.dumps(hh._json_ready({k:v for k, v in decision.items() if not k.endswith("sha256")}), indent=1))


def unblind():
    import pandas as pd
    prereg = check_preregistration(); decision = json.loads((BENCH/"decision.json").read_text())
    if not decision["passed"]:raise SystemExit("Phase 1 did not pass; phase 2 is not run")
    v3 = pd.read_csv(BENCH/"gale_blinded/joint_run1/genes.tsv", sep="\t"); v2 = pd.read_csv(V2, sep="\t")
    lab = pd.read_csv(ROOT/"resources/species/glabrata/essentiality_labels/labels.tsv", sep="\t", dtype=str).rename(columns={"cglab_id":"ORF"})
    d = (lab[["ORF", "label", "tier"]].merge(v3[["ORF", "status", "P_ES_run", "call"]], on="ORF")
         .merge(v2[["ORF", "status", "P_ES_run", "call"]], on="ORF", suffixes=("_v3", "_v2")))
    e = d[(d.tier == "silver") & d.label.isin(["essential", "nonessential"]) & (d.status_v3 == "called") & (d.status_v2 == "called")]
    e.to_csv(BENCH/"evaluation.tsv", sep="\t", index=False)
    y = (e.label == "essential").to_numpy(); p3, p2 = e.P_ES_run_v3.to_numpy(), e.P_ES_run_v2.to_numpy()
    rng = np.random.default_rng(SEED); boots = []
    for _ in range(2000):
        i = rng.integers(0, len(e), len(e)); yy = y[i]
        if yy.all() or not yy.any():continue
        boots.append(auc(p3[i], yy)-auc(p2[i], yy))
    ci = [float(np.quantile(boots, .025)), float(np.quantile(boots, .975))]
    out = dict(preregistration=prereg, labels="silver (S. cerevisiae ortholog)", genes=len(e), essential=int(y.sum()),
               auc_v3=auc(p3, y), auc_v2=auc(p2, y), auc_difference=auc(p3, y)-auc(p2, y), auc_difference_ci=ci,
               non_inferior=bool(ci[0] > -.01))
    for v in ("v3", "v2"):
        es = (e[f"call_{v}"] == "ES").to_numpy()
        out[v] = dict(es_calls=int(es.sum()), precision=float(y[es].mean()) if es.any() else None, recall=float(es[y].mean()))
    (BENCH/"unblinded.json").write_text(json.dumps(hh._json_ready(out), indent=2)+"\n"); print(json.dumps(hh._json_ready(out), indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("phase", choices=["blind", "unblind"]); p.add_argument("--workers", type=int, default=5)
    a = p.parse_args(); blind(a.workers) if a.phase == "blind" else unblind()


if __name__ == "__main__":main()
