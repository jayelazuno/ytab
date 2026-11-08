#!/usr/bin/env python3
"""Pre-registered test: does library density decide whether a proportional essential run helps? (S. pombe, thinned)

  python scripts/local/ytab_hmm_run_density.py freeze     # record the SHA-256 and UTC time of PREREGISTRATION.md
  python scripts/local/ytab_hmm_run_density.py score      # rescore the saturation-benchmark fits; reads no labels
  python scripts/local/ytab_hmm_run_density.py evaluate   # refuses if PREREGISTRATION.md changed; writes result.json, arms.tsv

Fits: output/essentiality/saturation/runs/<species>/q<q>_s<seed>/hmm/joint (thinned libraries; not refitted). Their
posterior chains are rebuilt from fit.json as in ytab_hmm_proportional_run.py and checked against the saved posteriors.
Design: output/essentiality/hmm_run_density/PREREGISTRATION.md.
"""
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

OUT = ROOT/"output/essentiality/hmm_run_density"
RUNS = ROOT/"output/essentiality/saturation/runs"
SEED, BOOTSTRAPS, MAX_GAMMA_DIFF = 16180, 2000, 1e-6
FRACTIONS = {"prop_20":.2, "prop_33":1/3, "prop_50":.5}
FITS = [(q, s) for q, s in ((1, 0), (2, 1), (2, 2), (4, 1), (4, 2), (8, 1), (8, 2))]   # non-converged fits are skipped
TESTED = ("pombe", (4, 1), (4, 2))


def freeze():
    stamp = OUT/"PREREGISTRATION.sha256"
    if stamp.exists():raise SystemExit("already frozen")
    stamp.write_text(f"{sens.fsb.sha256(OUT/'PREREGISTRATION.md')} {dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}\n")
    print(stamp.read_text().strip())


def check():
    frozen = (OUT/"PREREGISTRATION.sha256").read_text().split()
    if sens.fsb.sha256(OUT/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def fit_dir(species, q, seed):
    return RUNS/species/f"q{q}_s{seed}"


def rebuild(species, q, seed):
    """Settings, genes, window posteriors, chains, mask and A of a saved thinned fit (no refitting)."""
    run = fit_dir(species, q, seed); hmm_dir = run/"hmm/joint"; info = json.loads((hmm_dir/"fit.json").read_text())
    settings = hh.Settings(**{k:tuple(v) if k == "init_r" else v for k, v in info["settings"].items()})
    project, libraries = sens.SPECIES[species]
    paths = [run/"hits"/f"{lib}_hits.txt" for lib in libraries]
    fasta, gff = ROOT/info["inputs"]["fasta"], ROOT/info["inputs"]["gff"]
    genes = runner.ccg.annotations(Path(gff).resolve())
    target = runner.learn_target(libraries, paths, settings, project, fasta)
    feats = runner.features(settings, genes, fasta, target)
    libs = [runner.hit_sites(s, settings.min_reads, project, p) for s, p in zip(libraries, paths)]
    genes, _ = runner.check_inputs(feats, genes, libs, libraries, settings)
    chroms, _ = runner.design(feats, libs, settings)
    theta = hh.pack(np.log(np.array(info["r"])), np.array(info["a"]), np.array(info["b"]), np.log(info["psi"]))
    A, pi = np.array(info["A"]), np.array(info["pi"])
    saved = np.load(hmm_dir/"windows.npz"); gamma, chain, keep, diff = {}, {}, {}, 0.
    for c in chroms:
        g, _, _, B, beta = hh.forward_backward(hh.emission_loglik(c["Y"], c["X"], c["keep"], theta), A, pi, full=True)
        diff = max(diff, float(np.abs(g-saved[f"{c['name']}_gamma"]).max()))
        gamma[c["name"]], chain[c["name"]], keep[c["name"]] = g, (B, beta), c["keep"]
    return settings, genes, gamma, chain, keep, A, diff, info


def scores(genes, gamma, chain, keep, A, settings, fraction):
    """P_ES_run per called gene with run m = max(min_es_windows, ceil(fraction x mappable windows)), capped (hh.gene_table)."""
    from dataclasses import replace
    rows = hh.gene_table(genes, set(gamma), gamma, keep, replace(settings, es_run_fraction=fraction), chain, A)
    return {r["ORF"].removeprefix("gene:"):r["P_ES_run"] for r in rows if r.get("status") == "called"}


def score():
    check(); log = []
    for species in sens.SPECIES:
        for q, seed in FITS:
            fit_json = fit_dir(species, q, seed)/"hmm/joint/fit.json"
            if not fit_json.exists():continue
            if not json.loads(fit_json.read_text())["converged"]:   # DEVIATIONS.md 1: a capped fit cannot be rebuilt exactly
                line = f"{species} q{q}_s{seed}: not converged, skipped"; print(line, flush=True); log.append(line); continue
            settings, genes, gamma, chain, keep, A, diff, info = rebuild(species, q, seed)
            line = f"{species} q{q}_s{seed}: converged={info['converged']} max |rebuilt - saved posterior| = {diff:.2e}"
            print(line, flush=True); log.append(line)
            if diff > MAX_GAMMA_DIFF:raise SystemExit(f"rebuilt posteriors differ from the saved fit ({diff:.2e}); not scoring")
            out = OUT/"scores"/species/f"q{q}_s{seed}"; out.mkdir(parents=True, exist_ok=True)
            for name, f in {"fixed":0., **FRACTIONS}.items():
                with open(out/f"{name}.tsv", "w", newline="") as handle:
                    w = csv.writer(handle, delimiter="\t", lineterminator="\n"); w.writerow(["ORF", "P_ES_run"])
                    w.writerows(sorted(scores(genes, gamma, chain, keep, A, settings, f).items()))
            (out/"fit.json").write_text(json.dumps(dict(converged=info["converged"], iterations=info["iterations"],
                                                        sites_per_coding_kb=json.loads((fit_dir(species, q, seed)/"run.json").read_text())["sites_per_coding_kb"])))
    (OUT/"score.log").write_text("\n".join(log)+"\n")


def read(path):
    with open(path) as handle:return {r["ORF"]:float(r["P_ES_run"]) for r in csv.DictReader(handle, delimiter="\t")}


def evaluate():
    prereg = check(); rng = np.random.default_rng(SEED); rows = []
    for species in sens.SPECIES:
        lab = sens.labels(species)
        for q, seed in FITS:
            d = OUT/"scores"/species/f"q{q}_s{seed}"
            if not d.exists():continue
            info = json.loads((d/"fit.json").read_text()); ref = read(d/"fixed.tsv")
            for name in FRACTIONS:
                alt = read(d/f"{name}.tsv"); common = sorted(set(ref) & set(alt) & set(lab)); y = np.array([lab[g] for g in common])
                a, b = np.array([alt[g] for g in common]), np.array([ref[g] for g in common]); boots = []
                for _ in range(BOOTSTRAPS):
                    i = rng.integers(0, len(y), len(y))
                    if 0 < y[i].sum() < len(i):boots.append(sens.fsb.auc(a[i], y[i])-sens.fsb.auc(b[i], y[i]))
                lo, hi = np.percentile(boots, [2.5, 97.5])
                rows.append(dict(species=species, q=q, seed=seed, sites_per_coding_kb=info["sites_per_coding_kb"], converged=info["converged"],
                                 tested=species == TESTED[0] and (q, seed) in TESTED[1:] and name == "prop_50", arm=name, genes=len(y),
                                 essential=int(y.sum()), auc_fixed=sens.fsb.auc(b, y), auc_arm=sens.fsb.auc(a, y), delta=sens.fsb.auc(a, y)-sens.fsb.auc(b, y),
                                 delta_low=lo, delta_high=hi, direction="better" if lo > 0 else "worse" if hi < 0 else "neither"))
    primary = [r for r in rows if r["tested"]]
    if len(primary) != 2 or not all(r["converged"] for r in primary):verdict = "not evaluable (a primary fit is missing or did not converge)"
    else:verdict = "supported" if all(r["direction"] == "better" for r in primary) else "not supported"
    pombe = {q:np.mean([r["delta"] for r in rows if r["species"] == "pombe" and r["arm"] == "prop_50" and r["q"] == q and r["converged"]])
             for q in (1, 2, 4)}
    result = dict(preregistration=prereg, primary=verdict, pombe_prop50_mean_delta_by_q=pombe,
                  secondary_monotone=bool(pombe[1] < pombe[2] < pombe[4]))
    with open(OUT/"arms.tsv", "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    (OUT/"result.json").write_text(json.dumps(hh._json_ready(result), indent=2))
    print(json.dumps(hh._json_ready(result), indent=2))


if __name__ == "__main__":
    {"freeze":freeze, "score":score, "evaluate":evaluate}[sys.argv[1]]()
