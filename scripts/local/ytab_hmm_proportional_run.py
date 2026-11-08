#!/usr/bin/env python3
"""Pre-registered test: an essential run proportional to gene length instead of a fixed 3 windows.

  python scripts/local/ytab_hmm_proportional_run.py freeze     # record the SHA-256 and UTC time of PREREGISTRATION.md
  python scripts/local/ytab_hmm_proportional_run.py score      # rescore the sensitivity fits; reads no labels
  python scripts/local/ytab_hmm_proportional_run.py evaluate   # refuses if PREREGISTRATION.md changed; writes result.json, arms.tsv

The fits are those of output/essentiality/hmm_sensitivity/fits/ (not refitted). Their posterior chains were not saved,
so `score` rebuilds them with one forward-backward pass from the fitted parameters in fit.json and checks that the
rebuilt window posteriors match the saved ones before scoring.
Design: output/essentiality/hmm_proportional_run/PREREGISTRATION.md.
"""
from dataclasses import replace
import csv
import datetime as dt
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src")); sys.path.append(str(ROOT/"scripts/local"))
from ytab.essentiality import hermes_hmm as hh
import ytab_hermes_hmm as runner
import ytab_hmm_sensitivity as sens

OUT = ROOT/"output/essentiality/hmm_proportional_run"
FITS = sens.OUT/"fits"
SEED, BOOTSTRAPS = 27182, 2000
FRACTIONS = {"prop_20":.2, "prop_33":1/3, "prop_50":.5}            # tested, on the default (100-bp) fits
REPORTED = {"window_50": ["prop_33"], "window_200": ["prop_33"]}    # reported only
MAX_GAMMA_DIFF = 1e-6


def freeze():
    stamp = OUT/"PREREGISTRATION.sha256"
    if stamp.exists():raise SystemExit("already frozen")
    stamp.write_text(f"{sens.fsb.sha256(OUT/'PREREGISTRATION.md')} {dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}\n")
    print(stamp.read_text().strip())


def check():
    frozen = (OUT/"PREREGISTRATION.sha256").read_text().split()
    if sens.fsb.sha256(OUT/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def rebuild(species, arm):
    """Genes, chromosomes, window posteriors and posterior chains of a saved fit, from its parameters (no refitting)."""
    fit_dir = FITS/species/arm; info = json.loads((fit_dir/"fit.json").read_text())
    settings = hh.Settings(**{k:tuple(v) if k == "init_r" else v for k, v in info["settings"].items()})
    project, libraries = sens.SPECIES[species]; fasta, gff = runner.project_reference(project)
    genes = runner.ccg.annotations(Path(gff).resolve())
    target = runner.learn_target(libraries, [None]*len(libraries), settings, project, fasta)
    feats = runner.features(settings, genes, fasta, target)
    libs = [runner.hit_sites(s, settings.min_reads, project) for s in libraries]
    genes, _ = runner.check_inputs(feats, genes, libs, libraries, settings)
    chroms, _ = runner.design(feats, libs, settings)
    theta = hh.pack(np.log(np.array(info["r"])), np.array(info["a"]), np.array(info["b"]), np.log(info["psi"]))
    A, pi = np.array(info["A"]), np.array(info["pi"])
    saved = np.load(fit_dir/"windows.npz"); gamma, chain, keep, diff = {}, {}, {}, 0.
    for c in chroms:
        g, _, _, B, beta = hh.forward_backward(hh.emission_loglik(c["Y"], c["X"], c["keep"], theta), A, pi, full=True)
        diff = max(diff, float(np.abs(g-saved[f"{c['name']}_gamma"]).max()))
        gamma[c["name"]], chain[c["name"]], keep[c["name"]] = g, (B, beta), c["keep"]
    return settings, genes, gamma, chain, keep, A, diff


def proportional_scores(genes, gamma, chain, keep, A, settings, fraction):
    """P_ES_run with the run length m = max(settings.min_es_windows, ceil(fraction x mappable windows)), capped at the
    gene's mappable windows; everything else as hh.gene_table."""
    rows = []
    for gene in genes:
        if gene["chrom"] not in gamma:continue
        ws = hh.gene_windows(gene, settings.window, settings)
        if not ws:continue
        idx = np.array([w for w, _ in ws]); ov = np.array([o for _, o in ws], dtype=float); kp = keep[gene["chrom"]][idx]
        if not kp.any() or ov[kp].sum()/ov.sum() < settings.min_gene_mappable:continue
        n = int(kp.sum()); m = max(1, min(max(settings.min_es_windows, math.ceil(fraction*n)), n))
        span = np.arange(idx.min(), idx.max()+1); B, beta = chain[gene["chrom"]]
        p = hh.run_probability(gamma[gene["chrom"]][span[0]], A, B[span], beta[span], np.isin(span, idx[kp]), m)
        rows.append(dict(ORF=gene["id"], status="called", windows=n, es_run_windows=m, P_ES_run=p))
    return rows


def score():
    check(); log = []
    for species in sens.SPECIES:
        for arm, names in [("default", list(FRACTIONS))]+list(REPORTED.items()):
            settings, genes, gamma, chain, keep, A, diff = rebuild(species, arm)
            line = f"{species} {arm}: max |rebuilt - saved posterior| = {diff:.2e}"; print(line, flush=True); log.append(line)
            if diff > MAX_GAMMA_DIFF:raise SystemExit(f"rebuilt posteriors differ from the saved fit ({diff:.2e}); not scoring")
            for name in names:
                rows = proportional_scores(genes, gamma, chain, keep, A, settings, FRACTIONS[name]); (OUT/"scores"/species).mkdir(parents=True, exist_ok=True)
                runner.write_table(OUT/"scores"/species/f"genes.{arm}.{name}.tsv", list(rows[0]), rows)
    (OUT/"score.log").write_text("\n".join(log)+"\n")


def evaluate():
    prereg = check(); rng = np.random.default_rng(SEED); rows = []
    for species in sens.SPECIES:
        lab = sens.labels(species); ref = sens.genes(FITS/species/"default"/"genes.tsv")
        for arm, names in [("default", list(FRACTIONS))]+list(REPORTED.items()):
            for name in names:
                alt = sens.genes(OUT/"scores"/species/f"genes.{arm}.{name}.tsv")
                common = sorted(set(ref) & set(alt) & set(lab)); y = np.array([lab[g] for g in common])
                a, b = np.array([alt[g] for g in common]), np.array([ref[g] for g in common])
                boots = []
                for _ in range(BOOTSTRAPS):
                    i = rng.integers(0, len(y), len(y))
                    if 0 < y[i].sum() < len(i):boots.append(sens.fsb.auc(a[i], y[i])-sens.fsb.auc(b[i], y[i]))
                lo, hi = np.percentile(boots, [2.5, 97.5]); d = sens.fsb.auc(a, y)-sens.fsb.auc(b, y)
                from scipy.stats import spearmanr
                rows.append(dict(species=species, fit=arm, arm=name, tested=arm == "default", genes=len(y), essential=int(y.sum()),
                                 auc_default=sens.fsb.auc(b, y), auc_arm=sens.fsb.auc(a, y), delta=d, delta_low=lo, delta_high=hi,
                                 direction="better" if lo > 0 else "worse" if hi < 0 else "neither",
                                 spearman=spearmanr(a, b).statistic, kappa_es=sens.kappa(a >= .5, b >= .5),
                                 es_calls_arm=int((a >= .5).sum()), es_calls_default=int((b >= .5).sum())))
    claims = {}
    for name in FRACTIONS:
        r = [x for x in rows if x["arm"] == name and x["tested"]]
        claims[name] = ("better" if all(x["direction"] == "better" for x in r) else
                        "worse" if any(x["direction"] == "worse" for x in r) else "no difference shown")
    better = [n for n, c in claims.items() if c == "better"]
    adopt = max(better, key=lambda n: np.mean([x["delta"] for x in rows if x["arm"] == n and x["tested"]])) if better else None
    result = dict(preregistration=prereg, claims=claims, adopt=adopt)
    with open(OUT/"arms.tsv", "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    (OUT/"result.json").write_text(json.dumps(hh._json_ready(result), indent=2))
    print(json.dumps(hh._json_ready(result), indent=2))


if __name__ == "__main__":
    {"freeze":freeze, "score":score, "evaluate":evaluate}[sys.argv[1]]()
