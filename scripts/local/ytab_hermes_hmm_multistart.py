#!/usr/bin/env python3
"""Label-free robustness check of the Hermes window HMM v2: Baum-Welch from several starting points.

Start 0 is the default (Settings.init_r, init_stay, zero covariate effects, psi 0.5). Starts 1..N draw
r_ES ~ U(0.005, 0.1), r_GD ~ U(0.15, 0.6), r_GA ~ U(1.5, 5), stay ~ U(0.7, 0.98), covariate effects
~ N(0, 0.5), psi ~ U(0.1, 1), with seed 2718. Gale Cg1-Cg3, joint fit. Reports, per start: log-likelihood,
parameters, iterations, and agreement of gene P_ES_run / calls with start 0.
Writes output/essentiality/hermes_hmm/multistart/. Reads no labels.
"""
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src")); sys.path.append(str(ROOT/"scripts/local"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.essentiality import hermes_hmm as hh
from ytab.normalization.ttr import write_table
import ytab_hermes_hmm as runner

OUT = ROOT/"output/essentiality/hermes_hmm/multistart"
STARTS, SEED = 7, 2718


def start(i, L, P, Y):
    s = hh.Settings()
    if i == 0:return None, None
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(i,)))
    r = np.array([rng.uniform(.005, .1), rng.uniform(.15, .6), 1., rng.uniform(1.5, 5)])
    theta = hh.pack(np.log(r), np.log(np.maximum(Y.mean(0), 1e-3)), rng.normal(0, .5, P), np.log(rng.uniform(.1, 1)))
    stay = rng.uniform(.7, .98); A = np.full((4, 4), (1-stay)/3); np.fill_diagonal(A, stay)
    return theta, A


def task(i):
    settings = hh.Settings(); genes = ccg.annotations(runner.GFF); feats = runner.features(settings, genes)
    libs = [runner.hit_sites(s, settings.min_reads) for s in runner.GALE]; chroms, _ = runner.design(feats, libs, settings)
    Y = np.vstack([c["Y"] for c in chroms]); keep = np.concatenate([c["keep"] for c in chroms])
    theta, A = start(i, Y.shape[1], chroms[0]["X"].shape[1], Y[keep])
    f = hh.fit(chroms, settings, theta=theta, A=A)
    G, offsets = f["gamma"], np.cumsum([0]+[len(c["Y"]) for c in chroms])
    gamma = {c["name"]:G[offsets[k]:offsets[k+1]] for k, c in enumerate(chroms)}; kp = {c["name"]:c["keep"] for c in chroms}
    chain = {c["name"]:f["chains"][k] for k, c in enumerate(chroms)}
    rows = hh.gene_table(genes, set(gamma), gamma, kp, settings, chain, f["A"])
    return i, dict(loglik=f["loglik"], iterations=f["iterations"], converged=f["converged"], r=f["r"].tolist(), b=f["b"].tolist(),
                   psi=f["psi"], stay=np.diag(f["A"]).tolist(), init=None if theta is None else dict(r=np.exp(hh.unpack(theta, 3, 3)[0]).tolist(),
                   stay=float(A[0, 0]))), {r["ORF"]:(r.get("P_ES_run"), r["call"]) for r in rows if r["status"] == "called"}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(4) as pool:res = list(pool.map(task, range(STARTS+1)))
    base = res[0][2]; genes = sorted(base); rows = []
    for i, info, calls in res:
        common = [g for g in genes if g in calls]
        p0 = np.array([base[g][0] for g in common]); p1 = np.array([calls[g][0] for g in common])
        es0 = np.array([base[g][1] == "ES" for g in common]); es1 = np.array([calls[g][1] == "ES" for g in common])
        rows.append(dict(start=i, loglik=round(info["loglik"], 2), delta_loglik=round(info["loglik"]-res[0][1]["loglik"], 2), iterations=info["iterations"],
                         converged=info["converged"], r_ES=round(info["r"][0], 4), r_GD=round(info["r"][1], 4), r_GA=round(info["r"][3], 3), psi=round(info["psi"], 4),
                         b=";".join(f"{x:.3f}" for x in info["b"]), es_calls=int(es1.sum()), es_call_agreement=float(np.mean(es0 == es1)),
                         es_calls_changed=int((es0 != es1).sum()), p_es_run_spearman=float(spearmanr(p0, p1)[0]),
                         max_abs_p_diff=float(np.max(np.abs(p0-p1))), init=json.dumps(info["init"])))
    write_table(OUT/"starts.tsv", list(rows[0]), rows)
    for r in rows:print({k:v for k, v in r.items() if k != "init"})


if __name__ == "__main__":main()
