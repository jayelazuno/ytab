#!/usr/bin/env python3
"""Pre-registered benchmark: hybrid essentiality model (existing RF features + Hermes HMM v2 posteriors)
against the individual models, on Gale 2020 Cg1-Cg3 with the curated C. glabrata labels.

Arms (see output/essentiality/hybrid/benchmark/PREREGISTRATION.md):
  hmm            HMM v2 P_ES_run (unsupervised; no training)
  rf_existing    current RF score "RF - G4" (trained on S. cerevisiae data; unchanged)
  rf_features    random forest retrained on C. glabrata silver labels, existing RF features only (control)
  hybrid         same learner, existing RF features + HMM features
  hmm_features   same learner, HMM features only (reported)
  *_logistic     logistic-regression versions of the three supervised arms (reported)
Silver: out-of-fold predictions, 5 folds grouped by chromosome. Gold: trained on all silver, scored on gold.
"""
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.normalization.ttr import sha256

BENCH = ROOT/"output/essentiality/hybrid/benchmark"
SEED = 90210
RF_FEATURES = ["Length", "Hits", "Reads", "Longest free interval", "Freedom index", "Neighborhood index", "Upstream hits 100"]
LOG_FEATURES = ["Length", "Hits", "Reads", "Longest free interval", "Upstream hits 100"]
HMM_FEATURES = ["P_ES_run", "P_ES", "P_GD", "P_NE", "P_GA", "mappable_fraction"]
HMM = ROOT/"output/essentiality/hermes_hmm/benchmark_v2/gale_blinded/joint_run1/genes.tsv"
RF = ROOT/"output/projects/Gale_2020/classifier/T409600/essentiality_predictions.T409600.csv"
LABELS = ROOT/"resources/species/glabrata/essentiality_labels/labels.tsv"
CHROM_FOLD = {f"CP0482{30+i}.1":i % 5 for i in range(13)}   # chromosomes A..M -> folds 0..4


def check_preregistration():
    frozen = (BENCH/"PREREGISTRATION.sha256").read_text().split()
    if sha256(BENCH/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def auc(score, label):
    score, label = np.asarray(score, float), np.asarray(label, bool)
    order = np.argsort(score, kind="mergesort"); ranks = np.empty(len(score)); ranks[order] = np.arange(1, len(score)+1)
    _, inv, counts = np.unique(score, return_inverse=True, return_counts=True); ranks = (np.bincount(inv, ranks)/counts)[inv]
    n1, n0 = label.sum(), (~label).sum(); return float((ranks[label].sum()-n1*(n1+1)/2)/(n1*n0))


def learner(kind):
    if kind == "rf":return RandomForestClassifier(n_estimators=500, min_samples_leaf=5, n_jobs=4, random_state=SEED)
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, C=1.))


def design(d, features, kind):
    X = d[features].astype(float).copy()
    if kind == "logistic":
        for f in features:
            if f in LOG_FEATURES:X[f] = np.log1p(X[f])
    return X.to_numpy()


def data():
    lab = pd.read_csv(LABELS, sep="\t", dtype=str).rename(columns={"cglab_id":"ORF"})
    hmm = pd.read_csv(HMM, sep="\t"); rf = pd.read_csv(RF).rename(columns={"Standard name":"ORF"})
    d = (lab[["ORF", "name", "chromosome", "label", "tier", "cg_direct_call", "sc_call", "gale2020_comparator"]]
         .merge(hmm[["ORF", "status", "call"]+HMM_FEATURES], on="ORF", how="left").merge(rf[["ORF", "RF - G4"]+RF_FEATURES], on="ORF", how="left"))
    usable = (d.status == "called") & d[RF_FEATURES+["RF - G4"]].notna().all(1) & d[HMM_FEATURES].notna().all(1)
    return d, usable


def paired_boot(y, a, b, rng, n=2000):
    out = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y)); yy = y[i]
        if yy.all() or not yy.any():continue
        out.append(auc(a[i], yy)-auc(b[i], yy))
    return [float(np.quantile(out, .025)), float(np.quantile(out, .975))]


def metrics(y, s, threshold=.5):
    called = s >= threshold
    return dict(auc=auc(s, y), precision=float(y[called].mean()) if called.any() else None, recall=float(called[y].mean()),
                calls=int(called.sum()), brier=float(np.mean((s-y)**2)) if np.all((s >= 0) & (s <= 1)) else None)


def main():
    prereg = check_preregistration()
    if (BENCH/"result.json").exists():raise SystemExit("result.json exists; the benchmark already ran")
    d, usable = data()
    sil = d[usable & (d.tier == "silver") & d.label.isin(["essential", "nonessential"])].reset_index(drop=True)
    gold = d[usable & (d.tier == "gold") & d.label.isin(["essential", "nonessential"])].reset_index(drop=True)
    ys, yg = (sil.label == "essential").to_numpy(), (gold.label == "essential").to_numpy()
    folds = sil.chromosome.map(CHROM_FOLD).to_numpy()
    arms = {"rf_features":RF_FEATURES, "hybrid":RF_FEATURES+HMM_FEATURES, "hmm_features":HMM_FEATURES}
    oof, gpred, importances = {}, {}, {}
    for kind in ("rf", "logistic"):
        for arm, feats in arms.items():
            name = arm if kind == "rf" else f"{arm}_logistic"; p = np.zeros(len(sil))
            for k in range(5):
                tr, te = folds != k, folds == k
                m = learner(kind).fit(design(sil[tr], feats, kind), ys[tr]); p[te] = m.predict_proba(design(sil[te], feats, kind))[:, 1]
            oof[name] = p
            m = learner(kind).fit(design(sil, feats, kind), ys); gpred[name] = m.predict_proba(design(gold, feats, kind))[:, 1]
            if kind == "rf":importances[name] = dict(zip(feats, m.feature_importances_.round(4).tolist()))
    oof["hmm"], gpred["hmm"] = sil.P_ES_run.to_numpy(), gold.P_ES_run.to_numpy()
    oof["rf_existing"], gpred["rf_existing"] = sil["RF - G4"].to_numpy(), gold["RF - G4"].to_numpy()
    rng = np.random.default_rng(SEED)
    result = dict(preregistration=prereg, silver_genes=len(sil), silver_essential=int(ys.sum()), gold_genes=len(gold), gold_essential=int(yg.sum()),
                  excluded_unusable_labelled=int((~usable & d.label.isin(["essential", "nonessential"])).sum()),
                  silver={a:metrics(ys, s) for a, s in oof.items()}, gold={a:metrics(yg, s) for a, s in gpred.items()},
                  primary=dict(comparison="hybrid - rf_features (silver, out-of-fold AUC)", difference=auc(oof["hybrid"], ys)-auc(oof["rf_features"], ys),
                               ci=paired_boot(ys, oof["hybrid"], oof["rf_features"], rng)),
                  secondary={f"hybrid - {b}":dict(difference=auc(oof["hybrid"], ys)-auc(oof[b], ys), ci=paired_boot(ys, oof["hybrid"], oof[b], rng))
                             for b in ("rf_existing", "hmm")},
                  gold_secondary={f"hybrid - {b}":dict(difference=auc(gpred["hybrid"], yg)-auc(gpred[b], yg), ci=paired_boot(yg, gpred["hybrid"], gpred[b], rng))
                                  for b in ("rf_features", "rf_existing", "hmm")},
                  feature_importance=importances)
    result["primary"]["improvement_claimed"] = result["primary"]["ci"][0] > 0
    disc = gold[gold.sc_call.isin(["essential", "nonessential"]) & gold.cg_direct_call.isin(["essential", "nonessential"]) & (gold.cg_direct_call != gold.sc_call)]
    result["species_discordant_gold"] = [dict(ORF=r.ORF, name=r.name, cg=r.cg_direct_call, sc=r.sc_call,
                                              **{a:round(float(gpred[a][i]), 3) for a in ("hmm", "rf_existing", "rf_features", "hybrid")})
                                         for i, r in gold.iterrows() if r.ORF in set(disc.ORF)]
    g = sil.gale2020_comparator.isin(["essential", "nonessential"])
    result["gale_comparator_agreement_silver"] = {a:float(np.mean((s[g] >= .5) == (sil.gale2020_comparator[g] == "essential"))) for a, s in oof.items()}
    pred = sil[["ORF", "name", "chromosome", "label"]].assign(**{f"oof_{a}":s for a, s in oof.items()})
    pred.to_csv(BENCH/"silver_out_of_fold.tsv", sep="\t", index=False)
    gold[["ORF", "name", "label", "cg_direct_call", "sc_call"]].assign(**{f"pred_{a}":s for a, s in gpred.items()}).to_csv(BENCH/"gold_predictions.tsv", sep="\t", index=False)
    result["script_sha256"] = sha256(Path(__file__))
    (BENCH/"result.json").write_text(json.dumps(result, indent=2)+"\n"); print(json.dumps(result, indent=1))


if __name__ == "__main__":main()
