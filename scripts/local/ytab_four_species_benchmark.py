#!/usr/bin/env python3
"""Four-species essentiality benchmark: label-free window HMM v3 vs Segal et al. 2018's random forest.

Implements output/essentiality/four_species/PREREGISTRATION.md (frozen; both phases refuse to run if it changed).
  fit:       HMM v3 fits on the pre-registered libraries (reads no labels); writes fits.json
  evaluate:  reads the independent labels and Segal's RF scores; writes evaluation.tsv and result.json
"""
import argparse
import csv
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.normalization.ttr import sha256, write_table

BENCH = ROOT/"output/essentiality/four_species"
SEED = 81723
MARGIN = -0.05
BOOTSTRAPS = 2000
SEGAL = ROOT/"resources/species/albicans/essentiality_labels/raw/segal_2018_figshare"

FITS = {  # name: (species, project, libraries, role)
    "albicans_primary": ("albicans", "Segal_2018_albicans", ["SRR7824843", "SRR7824841", "SRR7824838"], "primary"),
    "albicans_all": ("albicans", "Segal_2018_albicans", [f"SRR78248{n}" for n in range(38, 49)], "secondary"),
    "cerevisiae_primary": ("cerevisiae", "Michel_2017_SATAY", ["ERR1533147", "ERR1533148"], "primary"),
    "pombe_primary": ("pombe", "Pombe_Hermes", ["SRR327340"], "primary"),
    "pombe_grech": ("pombe", "Pombe_Hermes", [f"ERR2641{n}" for n in range(199, 213)], "secondary"),
}
SEGAL_SHEETS = {  # species: (Dataset 2 sheet, header row, score column, verdict column, Dataset 1 sheet, published CV AUC)
    "albicans": ("2A-Ca classification", 5, "CaTn", "Ca RF verdict CaTn", "Dataset 1A. Calb", 0.997),
    "cerevisiae": ("2B- Sc classification", 8, "ScTn", "Sc RF verdict ScTn", "Dataset 1B. Sc", 0.989),
    "pombe": ("2C- Sp classification", 6, "SpTn", "Sp RF verdict SpTn", "Dataset 1C. Sp", 0.966),
}
MASKED_NAMES = {"cerevisiae": ["ADE2", "URA3", "CYC1", "GAL1", "GAL10", "HIS3", "LEU2", "MET15", "LYS2"], "albicans": ["ADE2"], "pombe": []}


def check_preregistration():
    frozen = (BENCH/"PREREGISTRATION.sha256").read_text().split()
    if sha256(BENCH/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def fit_one(name):
    species, project, libraries, role = FITS[name]
    out = BENCH/"hmm"/name
    cmd = [sys.executable, str(ROOT/"scripts/local/ytab_hermes_hmm.py"), str(out), "--project", project, "--libraries", *libraries,
           "--sequence-covariate", "learned"]
    run = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    (out).mkdir(parents=True, exist_ok=True); (out/"stdout.txt").write_text(run.stdout); (out/"stderr.txt").write_text(run.stderr)
    record = dict(name=name, species=species, project=project, libraries=libraries, role=role, command=" ".join(cmd[1:]), returncode=run.returncode)
    fit = out/"joint"/"fit.json"
    if run.returncode == 0 and fit.exists():
        s = json.loads(fit.read_text())
        record.update(converged=s.get("converged"), iterations=s.get("iterations"), masked_fraction=s.get("masked_fraction"),
                      notes=s.get("notes"), gene_calls=s.get("gene_calls"), sites_per_library=s.get("sites_per_library"),
                      target_information=s.get("target_model", {}).get("information_bits"),
                      genes_sha256=sha256(out/"joint"/"genes.tsv"), passed=bool(s.get("converged")) and s.get("masked_fraction", 1) < 0.5)
    else:
        record.update(passed=False, error=run.stderr.strip().splitlines()[-1:] if run.stderr else ["no output"])
    return record


def fit(workers):
    prereg = check_preregistration()
    with ProcessPoolExecutor(max_workers=workers) as pool:records = list(pool.map(fit_one, FITS))
    (BENCH/"fits.json").write_text(json.dumps(dict(preregistration=prereg, script_sha256=sha256(Path(__file__)), fits=records), indent=1))
    for r in records:print(r["name"], "passed" if r["passed"] else "FAILED", r.get("iterations"), r.get("masked_fraction"), r.get("gene_calls"))


# ---------------------------------------------------------------- evaluation (reads labels)

def auc(score, label):
    """ROC AUC (Mann-Whitney; ties averaged). label: 1 essential, 0 nonessential."""
    from scipy.stats import rankdata
    score, label = np.asarray(score, float), np.asarray(label, int); n1 = label.sum(); n0 = len(label)-n1
    ranks = rankdata(score); return float((ranks[label == 1].sum()-n1*(n1+1)/2)/(n1*n0))


def calls(pred, label):
    pred, label = np.asarray(pred, bool), np.asarray(label, bool); tp = int((pred & label).sum())
    precision = tp/max(1, int(pred.sum())); recall = tp/max(1, int(label.sum()))
    return dict(called=int(pred.sum()), precision=precision, recall=recall, f1=2*precision*recall/(precision+recall) if precision+recall else 0.0)


def hmm_genes(name):
    rows = {}
    for r in csv.DictReader(open(BENCH/"hmm"/name/"joint"/"genes.tsv"), delimiter="\t"):
        gid = r["ORF"].removeprefix("gene:")
        if r["P_ES_run"] not in ("", "nan", "NA"):rows[gid] = dict(p=float(r["P_ES_run"]), call=r["call"])
    return rows


def sgd_maps():
    sgdid, by_name = {}, {}
    for f in csv.reader(open(ROOT/"resources/species/cerevisiae/reference_genome/SGD_features.tab"), delimiter="\t"):
        if len(f) > 4 and f[3] and f[1] == "ORF":
            sgdid[f[0]] = f[3]
            for alias in [f[4] or f[3]]+(f[5].split("|") if len(f) > 5 and f[5] else []):by_name.setdefault(alias, f[3])  # MET15 is an alias of MET17
    return sgdid, by_name


def labels(species):
    """Independent labels: 1 essential, 0 nonessential."""
    if species == "pombe":
        vocab = {"inviable":1, "viable":0}
        return {g:vocab[v] for g, v in (l.rstrip("\n").split("\t")[:2] for l in open(ROOT/"resources/species/pombe/essentiality_labels/raw/FYPOviability.tsv")) if v in vocab}
    out = {}
    for r in csv.DictReader(open(ROOT/f"resources/species/{species}/essentiality_labels/labels.tsv"), delimiter="\t"):
        if r["tier"] == "gold" and r["label"] in ("essential", "nonessential"):out[r["gene_id"]] = int(r["label"] == "essential")
    return out


def masked(species):
    names = MASKED_NAMES[species]
    if not names:return set()
    if species == "cerevisiae":
        _, by_name = sgd_maps(); missing = [n for n in names if n not in by_name]
        if missing:raise SystemExit(f"masked genes not found: {missing}")
        return {by_name[n] for n in names}
    found = {r["name"]:r["gene_id"] for r in csv.DictReader(open(ROOT/f"resources/species/{species}/essentiality_labels/labels.tsv"), delimiter="\t") if r["name"] in names}
    if set(found) != set(names):raise SystemExit(f"masked genes not found: {set(names)-set(found)}")
    return set(found.values())


def segal(species):
    import pandas as pd
    sheet, header, score, verdict, train_sheet, cv_auc = SEGAL_SHEETS[species]
    d = pd.read_excel(SEGAL/"Dataset 2 ORF classifications.xlsx", sheet_name=sheet, header=header)
    d = d[d["Standard name"].notna()]
    ids = d["Standard name"].astype(str).str.strip()
    if species == "cerevisiae":
        sgdid, _ = sgd_maps(); ids = ids.map(lambda x:sgdid.get(x, ""))
    # The sheets' two headers are swapped relative to their contents: the column headed "<sp> RF verdict <sp>Tn" holds
    # the RF score and the one headed "<sp>Tn" holds the verdict (Ess/NE). Read them by content, and check it.
    numeric = pd.to_numeric(d[verdict], errors="coerce"); words = d[score].astype(str).str.strip()
    if numeric.notna().mean() < 0.95 or not words[numeric.notna()].isin(["Ess", "NE"]).all():
        raise SystemExit(f"Unexpected layout of Segal Dataset 2 sheet {sheet!r}")
    rf = {}
    for gid, s, v in zip(ids, numeric, words):
        if gid and np.isfinite(s):rf[gid] = dict(score=float(s), call=v == "Ess")
    t = pd.read_excel(SEGAL/"Dataset 1 Training sets and exclusions.xlsx", sheet_name=train_sheet, header=7)
    t = t.iloc[:, :4]; t.columns = ["id", "sc", "sp", "usage"]; t = t[t["id"].notna()]
    tid = t["id"].astype(str).str.strip()
    if species == "cerevisiae":tid = tid.map(lambda x:sgdid.get(x, ""))
    usage = {"Essential":1, "Non essential":0}  # genes the RF was trained on; FN, FP and conflicting genes were held out
    training = {g:usage[str(u).strip()] for g, u in zip(tid, t["usage"]) if g and str(u).strip() in usage}
    return rf, training, cv_auc


def bootstrap_delta(a, b, y, rng):
    n = len(y); deltas = []
    for _ in range(BOOTSTRAPS):
        i = rng.integers(0, n, n)
        if 0 < y[i].sum() < n:deltas.append(auc(a[i], y[i])-auc(b[i], y[i]))
    return np.percentile(deltas, [2.5, 97.5]).tolist()


def evaluate():
    prereg = check_preregistration(); fits = json.loads((BENCH/"fits.json").read_text())
    passed = {f["name"]:f["passed"] for f in fits["fits"]}
    rng = np.random.default_rng(SEED); rows, result = [], dict(preregistration=prereg, script_sha256=sha256(Path(__file__)), species={})
    for species in ("albicans", "cerevisiae", "pombe"):
        primary = f"{species}_primary"; res = dict(fit_passed=passed[primary]); result["species"][species] = res
        if not passed[primary]:continue
        hmm, lab, mask = hmm_genes(primary), labels(species), masked(species)
        rf, training, cv_auc = segal(species)
        steps = dict(labelled=len(lab), with_rf=len(set(lab) & set(rf)), with_hmm=len(set(lab) & set(rf) & set(hmm)))
        genes = sorted(set(lab) & set(rf) & set(hmm) - mask)
        steps["after_mask"] = len(genes); genes = [g for g in genes if g not in training]; steps["after_training_removed"] = len(genes)
        y = np.array([lab[g] for g in genes]); h = np.array([hmm[g]["p"] for g in genes]); r = np.array([rf[g]["score"] for g in genes])
        a_h, a_r = auc(h, y), auc(r, y); ci = bootstrap_delta(h, r, y, rng)
        res.update(gene_steps=steps, genes=len(genes), essential=int(y.sum()), auc_hmm=a_h, auc_rf=a_r, delta=a_h-a_r, delta_ci95=ci,
                   non_inferior=ci[0] > MARGIN, superior=ci[0] > 0,
                   calls_hmm=calls([hmm[g]["call"] == "ES" for g in genes], y), calls_rf=calls([rf[g]["call"] for g in genes], y))
        tg = sorted(set(training) & set(hmm) - mask); ty = np.array([training[g] for g in tg])
        res["training_set"] = dict(genes=len(tg), essential=int(ty.sum()), auc_hmm=auc([hmm[g]["p"] for g in tg], ty), rf_published_cv_auc=cv_auc)
        secondary = {"albicans":"albicans_all", "pombe":"pombe_grech"}.get(species)
        if secondary and passed.get(secondary):
            hs = hmm_genes(secondary); common = [i for i, g in enumerate(genes) if g in hs]
            ys = y[common]; a_s = auc([hs[genes[i]]["p"] for i in common], ys); a_p = auc(h[common], ys)
            res["secondary_libraries"] = dict(fit=secondary, genes=len(common), auc_hmm=a_s, auc_hmm_primary_same_genes=a_p, delta=a_s-a_p)
        for g, yy, hh_, rr in zip(genes, y, h, r):
            rows.append(dict(species=species, gene=g, label=int(yy), hmm_P_ES_run=hh_, hmm_call=hmm[g]["call"], rf_score=rr, rf_call=int(rf[g]["call"])))
    result["claim_holds"] = all(v.get("non_inferior", False) for v in result["species"].values())
    write_table(BENCH/"evaluation.tsv", list(rows[0]), rows)
    (BENCH/"result.json").write_text(json.dumps(result, indent=1)); print(json.dumps(result, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("phase", choices=["fit", "evaluate"]); p.add_argument("--workers", type=int, default=5)
    a = p.parse_args(); fit(a.workers) if a.phase == "fit" else evaluate()


if __name__ == "__main__":main()
