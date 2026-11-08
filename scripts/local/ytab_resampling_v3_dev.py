#!/usr/bin/env python3
"""DEVELOPMENT ONLY (not acceptance): resampling v3 variants on the v1 synthetic scenarios with development
seeds (never used for the final evaluation). Writes output/essentiality/resampling/v3_development/<tag>.tsv.
Usage: ytab_resampling_v3_dev.py TAG SEED EXPERIMENTS VARIANT_JSON [--workers N]
A variant may name "normalization" (any available YTAB method, applied to the raw counts; default ttr).
"""
from concurrent.futures import ProcessPoolExecutor
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from scipy import sparse

from ytab.essentiality import resampling
from ytab.normalization.methods import normalize
from ytab.normalization.ttr import write_table

spec = importlib.util.spec_from_file_location("v1", ROOT/"scripts/local/ytab_resampling_benchmark.py")
v1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(v1)
OUT = ROOT/"output/essentiality/resampling/v3_development"


def task(item):
    name, k, seed, variants = item
    v1.SEED = seed
    raw, scaled, factors, genes, idx, pos, truth = v1.make_data(name, k)
    cache = {"ttr":scaled}; out = {}
    for vid, spec in variants.items():
        spec = dict(spec); method = spec.pop("normalization", "ttr")   # optional per-variant normalization of the raw counts
        if method not in cache:cache[method] = normalize(sparse.csc_matrix(raw), method, ["synthetic"]*raw.shape[1])[0].toarray()
        m = cache[method]
        out[vid] = v1.score(resampling.analyze(m[:, :2], m[:, 2:], genes, idx, pos, resampling.Settings(**spec))[0], truth)
    return name, k, out


def main():
    tag, seed, experiments, variants = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), json.loads(sys.argv[4])
    workers = int(sys.argv[sys.argv.index("--workers")+1]) if "--workers" in sys.argv else 6
    with ProcessPoolExecutor(workers) as pool:
        results = list(pool.map(task, [(n, k, seed, variants) for n in v1.SCENARIOS for k in range(experiments)]))
    rows = []
    for vid in variants:
        for n in v1.SCENARIOS:
            m = v1.aggregate([r[2][vid] for r in results if r[0] == n])
            rows.append(dict(variant=vid, scenario=n, p05=round(m["p05_rate"], 4), ci_low=round(m["p05_ci"][0], 4), calls=m["calls"], false=m["false"],
                             exp_false=m["experiments_with_false_call"], sens=None if m["sensitivity"] is None else round(m["sensitivity"], 3),
                             fdr=None if m["fdr"] is None else round(m["fdr"], 3)))
    OUT.mkdir(parents=True, exist_ok=True); write_table(OUT/f"{tag}.tsv", list(rows[0]), rows)
    nulls = [n for n in v1.SCENARIOS if n.startswith("null_")]; signals = [n for n in v1.SCENARIOS if n not in nulls]
    for vid in variants:
        r = [x for x in rows if x["variant"] == vid]
        worst = max((x for x in r if x["scenario"] in nulls), key=lambda x:x["p05"])
        fdr = sum(x["false"] for x in r if x["scenario"] in signals)/max(sum(x["calls"] for x in r if x["scenario"] in signals), 1)
        sens = sum(x["sens"] or 0 for x in r if x["scenario"] in signals)/len(signals)
        print(f"{vid:28s} worst-null p05 {worst['p05']:.3f} ({worst['scenario']}, CI low {worst['ci_low']:.3f})  pooled FDR {fdr:.3f}  "
              f"null exps w/ false {sum(x['exp_false'] for x in r if x['scenario'] in nulls)}/{len(nulls)*experiments}  mean sens {sens:.3f}")


if __name__ == "__main__":main()
