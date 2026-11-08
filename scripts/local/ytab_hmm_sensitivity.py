#!/usr/bin/env python3
"""Pre-registered parameter sensitivity of the window HMM v3 in four species.

  python scripts/local/ytab_hmm_sensitivity.py freeze     # record the SHA-256 and UTC time of PREREGISTRATION.md
  python scripts/local/ytab_hmm_sensitivity.py fit        # 16 label-free fits (2 at a time); reads no labels
  python scripts/local/ytab_hmm_sensitivity.py evaluate   # refuses if PREREGISTRATION.md changed; writes result.json, arms.tsv

Design: output/essentiality/hmm_sensitivity/PREREGISTRATION.md.
"""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import csv
import datetime as dt
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src")); sys.path.append(str(ROOT/"scripts/local"))
from ytab.essentiality import hermes_hmm as hh
import ytab_hermes_hmm as runner

OUT = ROOT/"output/essentiality/hmm_sensitivity"
SEED, BOOTSTRAPS, TOLERANCE = 31415, 2000, 0.02
SPECIES = {  # species: (project, libraries)
    "glabrata": ("Gale_2020", ["SRR11561508", "SRR11561507", "SRR11561506"]),
    "albicans": ("Segal_2018_albicans", ["SRR7824843", "SRR7824841", "SRR7824838"]),
    "cerevisiae": ("Michel_2017_SATAY", ["ERR1533147", "ERR1533148"]),
    "pombe": ("Pombe_Hermes", [f"ERR2641{n}" for n in range(199, 213)]),
}
FIT_ARMS = {"default":{}, "window_50":dict(window=50), "window_200":dict(window=200), "min_reads_2":dict(min_reads=2)}
RESCORE = {"default":{"run_2":dict(min_es_windows=2), "run_4":dict(min_es_windows=4), "cterm_0":dict(cterm=0.),
                      "cterm_20":dict(cterm=20.), "nterm_5":dict(nterm=5.)},
           "window_50":{"window_50_run_6":dict(min_es_windows=6)}, "window_200":{"window_200_run_2":dict(min_es_windows=2)}}
TESTED = ["window_50", "window_200", "min_reads_2", "run_2", "run_4", "cterm_0", "cterm_20", "nterm_5"]
REPORTED = ["window_50_run_6", "window_200_run_2"]

spec = importlib.util.spec_from_file_location("fsb", ROOT/"scripts/local/ytab_four_species_benchmark.py")
fsb = importlib.util.module_from_spec(spec); spec.loader.exec_module(fsb)


def freeze():
    stamp = OUT/"PREREGISTRATION.sha256"
    if stamp.exists():raise SystemExit("already frozen")
    stamp.write_text(f"{fsb.sha256(OUT/'PREREGISTRATION.md')} {dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}\n")
    print(stamp.read_text().strip())


def check():
    frozen = (OUT/"PREREGISTRATION.sha256").read_text().split()
    if fsb.sha256(OUT/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def fit_one(job):
    species, arm = job
    project, libraries = SPECIES[species]
    settings = replace(hh.Settings(), **FIT_ARMS[arm])
    rescore = {name:replace(settings, **change) for name, change in RESCORE.get(arm, {}).items()}
    fasta, gff = runner.project_reference(project)
    out = OUT/"fits"/species/arm
    if (out/"fit.json").exists():return species, arm, "exists"
    try:
        s, _ = runner.run(out, libraries, settings, f"{species} {arm}", fasta=fasta, gff=gff, project=project,
                          sequence_covariate="learned", rescore=rescore)
        return species, arm, f"converged={s['converged']} iterations={s['iterations']} seconds={s['seconds']:.0f}"
    except Exception as error:
        return species, arm, f"FAILED: {error}"


def fit():
    check()
    jobs = [(sp, arm) for sp in SPECIES for arm in FIT_ARMS]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for species, arm, status in pool.map(fit_one, jobs):
            print(species, arm, status, flush=True)


def labels(species):
    if species == "glabrata":
        rows = csv.DictReader(open(ROOT/"resources/species/glabrata/essentiality_labels/labels.tsv"), delimiter="\t")
        return {r["cglab_id"]:int(r["label"] == "essential") for r in rows if r["tier"] == "silver" and r["label"] in ("essential", "nonessential")}
    lab = fsb.labels(species); mask = fsb.masked(species)
    return {g:y for g, y in lab.items() if g not in mask}


def genes(path):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r.get("status") == "called" and r.get("P_ES_run") not in (None, "", "nan"):
            out[r["ORF"].removeprefix("gene:")] = float(r["P_ES_run"])
    return out


def kappa(a, b):
    a, b = np.asarray(a, bool), np.asarray(b, bool); po = np.mean(a == b)
    pe = a.mean()*b.mean()+(1-a.mean())*(1-b.mean())
    return float((po-pe)/(1-pe)) if pe < 1 else 1.0


def arm_path(species, arm):
    base = OUT/"fits"/species
    if arm in FIT_ARMS:return base/arm/"genes.tsv", base/arm/"fit.json"
    parent = next(p for p, names in RESCORE.items() if arm in names)
    return base/parent/f"genes.{arm}.tsv", base/parent/"fit.json"


def evaluate():
    prereg = check(); rng = np.random.default_rng(SEED); rows, result = [], dict(preregistration=prereg, species={}, claims={})
    for species in SPECIES:
        lab = labels(species); ref_path, ref_fit = arm_path(species, "default")
        ref = genes(ref_path); ref_info = json.loads(ref_fit.read_text())
        result["species"][species] = dict(default_converged=ref_info["converged"], default_iterations=ref_info["iterations"])
        for arm in TESTED+REPORTED:
            path, fit_path = arm_path(species, arm)
            info = json.loads(fit_path.read_text()); alt = genes(path)
            common = sorted(set(ref) & set(alt) & set(lab)); y = np.array([lab[g] for g in common])
            a, b = np.array([alt[g] for g in common]), np.array([ref[g] for g in common])
            boots = []
            for _ in range(BOOTSTRAPS):
                i = rng.integers(0, len(y), len(y))
                if 0 < y[i].sum() < len(i):boots.append(fsb.auc(a[i], y[i])-fsb.auc(b[i], y[i]))
            lo, hi = np.percentile(boots, [2.5, 97.5])
            evaluable = bool(info["converged"] and ref_info["converged"])
            from scipy.stats import spearmanr
            rows.append(dict(species=species, arm=arm, tested=arm in TESTED, genes=len(y), essential=int(y.sum()),
                             auc_default=fsb.auc(b, y), auc_arm=fsb.auc(a, y), delta=fsb.auc(a, y)-fsb.auc(b, y), delta_low=lo, delta_high=hi,
                             equivalent=bool(-TOLERANCE <= lo and hi <= TOLERANCE), evaluable=evaluable,
                             spearman=spearmanr(a, b).statistic, kappa_es=kappa(a >= .5, b >= .5),
                             es_calls_arm=int((a >= .5).sum()), es_calls_default=int((b >= .5).sum()),
                             converged=info["converged"], iterations=info["iterations"], seconds=round(info["seconds"])))
    for arm in TESTED:
        r = [x for x in rows if x["arm"] == arm]
        result["claims"][arm] = ("not made (a fit did not converge)" if not all(x["evaluable"] for x in r) else
                                 "robust" if all(x["equivalent"] for x in r) else "not robust")
    with open(OUT/"arms.tsv", "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    (OUT/"result.json").write_text(json.dumps(hh._json_ready(result), indent=2))
    print(json.dumps(hh._json_ready(result), indent=2))


if __name__ == "__main__":
    {"freeze":freeze, "fit":fit, "evaluate":evaluate}[sys.argv[1]]()
