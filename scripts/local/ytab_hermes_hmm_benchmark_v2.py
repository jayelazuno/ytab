#!/usr/bin/env python3
"""Pre-registered evaluation of the Hermes window HMM v2 (hermes_hmm_v2: P_ES_run gene score) on Gale 2020 parents Cg1-Cg3.

blind:   label-free. (1) Gale joint fit, run twice (determinism, EM monotone/converged); (2) each library fitted
         alone (replicate consistency); (3) synthetic recovery: genomes simulated from the fitted Gale model on the
         real windows/covariates/masks/genes, with known gene states, under pre-declared scenarios.
         Applies the gates in benchmark/PREREGISTRATION.md; writes decision.json. Reads no labels.
unblind: only after a pass; compares gene P(ES) and calls with the curated labels
         (resources/species/glabrata/essentiality_labels/labels.tsv), the current random forest (RF - G4 score),
         the TRANSIT HMM port and Gale 2020 (comparator).
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

BENCH = ROOT/"output/essentiality/hermes_hmm/benchmark_v2"
SEED = 52417
GENE_STATE_P = dict(ES=.20, GD=.05, NE=.73, GA=.02)      # generic, pre-declared; not fitted to any labels
SCENARIOS = {"well_specified":{}, "accessibility_domains":dict(access_sd=.5, block=10), "partial_essential":dict(partial=.5),
             "low_complexity":dict(scale=1/3), "lognormal_poisson":dict(lognormal=True)}
GATED = ("well_specified", "accessibility_domains", "partial_essential")


def check_preregistration():
    frozen = (BENCH/"PREREGISTRATION.sha256").read_text().split()
    if sha256(BENCH/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def context(settings):
    genes = ccg.annotations(runner.GFF); feats = runner.features(settings, genes)
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


def blind(workers):
    prereg = check_preregistration()
    if (BENCH/"decision.json").exists():raise SystemExit("decision.json exists; phase 1 already ran")
    runs = [runner.run(BENCH/"gale_blinded"/f"joint_run{i}", runner.GALE, settings_global, f"joint run {i}") for i in (1, 2)]
    each = {s:runner.run(BENCH/"gale_blinded"/s, [s], settings_global, s) for s in runner.GALE}
    (s1, rows1), (s2, rows2) = runs
    identical = [(r["ORF"], r.get("P_ES_run"), r["call"]) for r in rows1] == [(r["ORF"], r.get("P_ES_run"), r["call"]) for r in rows2]
    fitted = dict(a=np.array(s1["a"]), b=np.array(s1["b"]), r=np.array(s1["r"]), psi=float(s1["psi"]))
    calls = {s:{r["ORF"]:r["call"] for r in rows if r["status"] == "called"} for s, (_, rows) in each.items()}
    kappas = {}
    for i, a in enumerate(runner.GALE):
        for b in runner.GALE[i+1:]:
            common = sorted(set(calls[a]) & set(calls[b]))
            kappas[f"{a}~{b}"] = dict(genes=len(common), kappa=kappa([calls[a][g] == "ES" for g in common], [calls[b][g] == "ES" for g in common]))
    with ProcessPoolExecutor(workers) as pool:
        synth = list(pool.map(synthetic_task, [(s, 0, fitted) for s in SCENARIOS]))
    rec = {s:r for s, _, r in synth}
    write_table(BENCH/"synthetic.tsv", ["scenario"]+[k for k in rec["well_specified"] if k not in ("r", "b")],
                [dict(scenario=s) | {k:v for k, v in r.items() if k not in ("r", "b")} for s, r in rec.items()])
    g_syn = all(rec[s]["auc"] >= .95 and (rec[s]["precision"] or 0) >= .90 and rec[s]["recall"] >= .80 and
                (rec[s]["high_confidence_true"] or 0) >= .90 for s in GATED)
    gates = dict(g1_determinism=identical, g2_em=bool(s1["monotone"] and s1["converged"]), g3_synthetic=g_syn)
    decision = dict(status="success", phase="blind", preregistration=prereg, passed=all(gates.values()), gates=gates,
                    gale_joint=dict(gene_calls=s1["gene_calls"], r=s1["r"], b=s1["b"], psi=s1["psi"], A=s1["A"], masked_fraction=s1["masked_fraction"],
                                    iterations=s1["iterations"]),
                    each_library={s:dict(gene_calls=v[0]["gene_calls"], r=v[0]["r"], sites=v[0]["sites_per_library"]) for s, v in each.items()},
                    replicate_kappa=kappas, synthetic=rec, blinding="no essentiality labels, RF output or Gale calls were read",
                    implementation_sha256=sha256(Path(hh.__file__)), script_sha256=sha256(Path(__file__)))
    (BENCH/"decision.json").write_text(json.dumps(hh._json_ready(decision), indent=2)+"\n")
    print(json.dumps(hh._json_ready({k:v for k, v in decision.items() if k not in ("implementation_sha256", "script_sha256")}), indent=1))


def unblind():
    import pandas as pd
    prereg = check_preregistration(); decision = json.loads((BENCH/"decision.json").read_text())
    if not decision["passed"]:raise SystemExit("Phase 1 did not pass; phase 2 is not run")
    hmm = pd.read_csv(BENCH/"gale_blinded/joint_run1/genes.tsv", sep="\t")
    lab = pd.read_csv(ROOT/"resources/species/glabrata/essentiality_labels/labels.tsv", sep="\t", dtype=str).rename(columns={"cglab_id":"ORF"})
    rf = pd.read_csv(ROOT/"output/projects/Gale_2020/classifier/T409600/essentiality_predictions.T409600.csv").rename(columns={"Standard name":"ORF"})
    tp = pd.read_csv(ROOT/"output/essentiality/hermes_hmm/development/transit_port/genes.tsv", sep="\t")
    d = (lab[["ORF", "label", "tier", "gale2020_comparator", "cg_direct_call", "sc_call"]].merge(hmm, on="ORF", how="left")
         .merge(rf[["ORF", "RF - G4"]], on="ORF", how="left").merge(tp[["ORF", "call"]].rename(columns={"call":"transit_port_call"}), on="ORF", how="left"))
    d.to_csv(BENCH/"evaluation.tsv", sep="\t", index=False)
    rng = np.random.default_rng(SEED); out = dict(preregistration=prereg)
    for tier in ("silver", "gold", "all"):
        e = d[(d.label.isin(["essential", "nonessential"])) & (d.status == "called") & d["RF - G4"].notna() & ((d.tier == tier) if tier != "all" else True)]
        y = (e.label == "essential").to_numpy()
        if y.sum() == 0 or (~y).sum() == 0:out[tier] = dict(genes=len(e)); continue
        a_h, a_r = auc(e.P_ES_run, y), auc(e["RF - G4"], y)
        boots = []
        for _ in range(2000):
            i = rng.integers(0, len(e), len(e)); yy = y[i]
            if yy.all() or not yy.any():continue
            boots.append(auc(e.P_ES_run.to_numpy()[i], yy)-auc(e["RF - G4"].to_numpy()[i], yy))
        es = (e.call == "ES").to_numpy(); tes = (e.transit_port_call == "ES").to_numpy()
        out[tier] = dict(genes=len(e), essential=int(y.sum()), auc_hmm=a_h, auc_rf=a_r, auc_difference=a_h-a_r,
                         auc_difference_ci=[float(np.quantile(boots, .025)), float(np.quantile(boots, .975))],
                         hmm_es_calls=int(es.sum()), hmm_precision=float(y[es].mean()) if es.any() else None, hmm_recall=float(es[y].mean()),
                         transit_port_es_calls=int(tes.sum()), transit_port_recall=float(tes[y].mean()))
    g = d[d.gale2020_comparator.isin(["essential", "nonessential"]) & (d.status == "called")]
    out["gale_comparator"] = dict(genes=len(g), agree=float(((g.call == "ES") == (g.gale2020_comparator == "essential")).mean()),
                                  auc_vs_gale=auc(g.P_ES_run, g.gale2020_comparator == "essential"))
    disc = d[d.tier.eq("gold") & d.sc_call.isin(["essential", "nonessential"]) & (d.cg_direct_call != d.sc_call)]
    out["species_discordant_gold"] = disc[["ORF", "cg_direct_call", "sc_call", "P_ES_run", "call", "RF - G4"]].to_dict(orient="records")
    (BENCH/"unblinded.json").write_text(json.dumps(hh._json_ready(out), indent=2)+"\n"); print(json.dumps(hh._json_ready(out), indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("phase", choices=["blind", "unblind"]); p.add_argument("--workers", type=int, default=5)
    a = p.parse_args(); blind(a.workers) if a.phase == "blind" else unblind()


if __name__ == "__main__":main()
