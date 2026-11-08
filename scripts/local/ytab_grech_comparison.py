#!/usr/bin/env python3
"""Confirmatory follow-up of the YTAB HMM vs Grech et al. 2019 comparison (pre-registered).

  python scripts/local/ytab_grech_comparison.py freeze     # record the SHA-256 and UTC time of PREREGISTRATION.md
  python scripts/local/ytab_grech_comparison.py evaluate   # refuses if PREREGISTRATION.md changed; writes result.json

Design: output/essentiality/grech_comparison/PREREGISTRATION.md.
"""
import datetime as dt
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/"output/essentiality/grech_comparison"
GRECH = ROOT/"resources/species/pombe/tnseq_data/grech_2019/grech2019_supplementary_tables.xlsx"
YTAB = ROOT/"output/essentiality/four_species/hmm/pombe_grech/joint/genes.tsv"
SEED, BOOTSTRAPS, MARGIN, MIN_CLASS = 20261024, 2000, -0.05, 10

spec = importlib.util.spec_from_file_location("fsb", ROOT/"scripts/local/ytab_four_species_benchmark.py")
fsb = importlib.util.module_from_spec(spec); spec.loader.exec_module(fsb)
auc, sha256 = fsb.auc, fsb.sha256


def freeze():
    stamp = OUT/"PREREGISTRATION.sha256"
    if stamp.exists():raise SystemExit("already frozen")
    stamp.write_text(f"{sha256(OUT/'PREREGISTRATION.md')} {dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}\n")
    print(stamp.read_text().strip())


def check():
    frozen = (OUT/"PREREGISTRATION.sha256").read_text().split()
    if sha256(OUT/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def load():
    t = pd.read_excel(GRECH, sheet_name="Table 2", header=14)
    t = t[t["gene"].notna() & t["mean.HMM"].notna()].copy(); t["gene"] = t["gene"].astype(str).str.strip()
    y = pd.read_csv(YTAB, sep="\t"); y = y[y["status"] == "called"].copy()
    y["gene"] = y["ORF"].str.removeprefix("gene:")
    labels = pd.Series(fsb.labels("pombe"), name="pombase")
    d = t.merge(y[["gene", "P_ES_run", "P_ES", "P_GD"]], on="gene").merge(labels, left_on="gene", right_index=True)
    return d


def bootstrap(stat, n, rng):
    out = []
    for _ in range(BOOTSTRAPS):
        i = rng.integers(0, n, n); v = stat(i)
        if v is not None and np.isfinite(v):out.append(v)
    return np.array(out)


def h1(d, rng):
    s = d[d["essential"].astype(int) != d["pombase"]]
    y, a, b = s["pombase"].to_numpy(), s["P_ES_run"].to_numpy(), -s["mean.HMM"].to_numpy()
    res = dict(genes=len(s), essential=int(y.sum()), nonessential=int(len(y)-y.sum()),
               grech_training_essential_now_viable=int(((s["essential"] == 1) & (s["pombase"] == 0)).sum()),
               grech_training_nonessential_now_inviable=int(((s["essential"] == 0) & (s["pombase"] == 1)).sum()))
    if min(res["essential"], res["nonessential"]) < MIN_CLASS:return dict(res, testable=False)
    delta = auc(a, y)-auc(b, y)
    bs = bootstrap(lambda i:auc(a[i], y[i])-auc(b[i], y[i]) if 0 < y[i].sum() < len(i) else None, len(y), rng)
    # reported only: genes whose labels agree, and calls at the usual thresholds
    agree = d[d["essential"].astype(int) == d["pombase"]]; ya = agree["pombase"].to_numpy()
    return dict(res, testable=True, auc_ytab=auc(a, y), auc_grech=auc(b, y), delta=delta, ci95=np.percentile(bs, [2.5, 97.5]).tolist(),
                p_one_sided=float(np.mean(bs <= 0)),
                ytab_es_call_matches_pombase=float(np.mean((a >= .5) == (y == 1))),
                agreeing_genes=dict(genes=len(agree), auc_ytab=auc(agree["P_ES_run"].to_numpy(), ya), auc_grech=auc(-agree["mean.HMM"].to_numpy(), ya)))


def correlations(d, col, rng):
    s = d[(d["pombase"] == 0) & d[col].notna()]
    x = pd.to_numeric(s[col], errors="coerce"); s = s[x.notna()]; x = x[x.notna()].to_numpy()
    a, b = (s["P_ES"]+s["P_GD"]).to_numpy(), -s["mean.HMM"].to_numpy()
    r = lambda u, v:-spearmanr(u, v).statistic
    bs = bootstrap(lambda i:r(a[i], x[i])-r(b[i], x[i]), len(x), rng)
    return dict(genes=len(x), r_ytab=r(a, x), r_grech=r(b, x), delta=r(a, x)-r(b, x), ci95=np.percentile(bs, [2.5, 97.5]).tolist(),
                p_one_sided=float(np.mean(bs <= MARGIN)), superior_ci_above_0=bool(np.percentile(bs, 2.5) > 0),
                pearson_ytab=-pearsonr(a, x).statistic, pearson_grech=-pearsonr(b, x).statistic)


def evaluate():
    prereg = check(); rng = np.random.default_rng(SEED); d = load()
    r1 = h1(d, rng); r2 = correlations(d, "Malecki2016.KO.colony.size", rng); r2s = correlations(d, "Sideri2015.barseq", rng)
    tests = {"H1":r1["p_one_sided"]} if r1["testable"] else {}
    tests["H2"] = r2["p_one_sided"]
    supported, level = {}, [0.05/len(tests), 0.05]
    for k, (name, p) in enumerate(sorted(tests.items(), key=lambda kv:kv[1])):
        supported[name] = bool(p <= level[k] if len(tests) == 2 else p <= 0.05) and (k == 0 or all(supported.values()))
    result = dict(preregistration=prereg, script_sha256=sha256(Path(__file__)), genes_scored=len(d), H1=r1, H2=r2,
                  H2_secondary_sideri=r2s, holm=dict(p=tests, supported=supported))
    (OUT/"result.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    {"freeze":freeze, "evaluate":evaluate}[sys.argv[1]]()
