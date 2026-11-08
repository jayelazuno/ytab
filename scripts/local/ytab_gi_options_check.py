#!/usr/bin/env python3
"""Pre-registered check of the GI options adapted from TRANSIT's GI method: delta interval, ROPE (TREAT) and Holm FWER.

Gates in output/essentiality/gi/options_check/PREREGISTRATION.md. Simulator and scenarios: ytab_gi_benchmark.py,
plus "gi_small" (query-only treatment fold 0.85: a real but small interaction, |delta| well below 0.5).
Fresh seed 94949, 8 experiments per scenario x design.
"""
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import genetic_interactions as gi
from ytab.normalization.methods import normalize
from ytab.normalization.ttr import sha256

spec = importlib.util.spec_from_file_location("bench", ROOT/"scripts/local/ytab_gi_benchmark.py")
bench = importlib.util.module_from_spec(spec); spec.loader.exec_module(bench)
bench.SCENARIOS["gi_small"] = dict(first=(1., .85, 1.))
OUT = ROOT/"output/essentiality/gi/options_check"; SEED, EXPERIMENTS, ROPE = 94949, 8, .5


def main():
    frozen = (OUT/"PREREGISTRATION.sha256").read_text().split()
    if sha256(OUT/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after freezing")
    if (OUT/"decision.json").exists():raise SystemExit("already ran")
    cover, fwer_exps, fwer_false, rope = {}, 0, 0, dict(small_genes=0, small_p05=0, calls=0, false=0, null_calls=0)
    for name in bench.SCENARIOS:
        if name in bench.REPORTED_ONLY:continue
        for design in ("split", "pools"):
            for k in range(EXPERIMENTS):
                raw, st, co, li, idx, truth = bench.simulate(name, design, k, SEED)
                v = normalize(sparse.csc_matrix(raw), "betageom", li)[0].toarray(); genes = [dict(id=f"g{j}", name="") for j in range(len(truth))]
                base, _ = gi.analyze(v, st, co, li, genes, idx)
                holm, _ = gi.analyze(v, st, co, li, genes, idx, settings=gi.Settings(correction="holm"))
                rp, _ = gi.analyze(v, st, co, li, genes, idx, settings=gi.Settings(rope=ROPE))
                c = cover.setdefault((name, design), [0, 0])
                for t, r in zip(truth, base):
                    if r["status"] == "tested" and not t:c[0] += 1; c[1] += r["delta_lower"] <= 0 <= r["delta_upper"]
                if name != "gi_small":
                    fwer_exps += 1; fwer_false += any(r["call"] in gi.CLASSES and not t for t, r in zip(truth, holm))
                for t, r in zip(truth, rp):
                    if r["status"] != "tested":continue
                    called = r["call"] in gi.CLASSES
                    if name == "gi_small" and t:rope["small_genes"] += 1; rope["small_p05"] += r["p_value"] <= .05
                    if name.startswith("gi_") and name != "gi_small":rope["calls"] += called; rope["false"] += called and not t
                    if name.startswith("null"):rope["null_calls"] += called
    coverage = {f"{n} {d}":v[1]/v[0] for (n, d), v in cover.items()}
    fwer = fwer_false/fwer_exps; ci = [float(beta.ppf(.025, fwer_false, fwer_exps-fwer_false+1)) if fwer_false else 0., float(beta.ppf(.975, fwer_false+1, fwer_exps-fwer_false))]
    small = rope["small_p05"]/rope["small_genes"]; rope_fdr = rope["false"]/rope["calls"] if rope["calls"] else 0.
    gates = {"interval coverage (non-GI genes) in [0.925, 0.975] everywhere":dict(value=[min(coverage.values()), max(coverage.values())],
                                                                                passed=all(.925 <= x <= .975 for x in coverage.values())),
             "Holm FWER (experiments with >= 1 false call) <= 0.075":dict(value=fwer, ci=ci, experiments=fwer_exps, passed=fwer <= .075),
             "ROPE: small interactions (|delta| < 0.5) with p <= 0.05 <= 0.075":dict(value=small, genes=rope["small_genes"], passed=small <= .075),
             "ROPE: FDR among calls in GI scenarios <= 0.10":dict(value=rope_fdr, calls=rope["calls"], passed=rope_fdr <= .10)}
    decision = dict(preregistration=dict(sha256=frozen[0], frozen_utc=frozen[1]), passed=all(g["passed"] for g in gates.values()), gates=gates,
                    coverage=coverage, rope_null_scenario_calls=rope["null_calls"], script_sha256=sha256(Path(__file__)), module_sha256=sha256(Path(gi.__file__)))
    (OUT/"decision.json").write_text(json.dumps(decision, indent=2)+"\n"); print(json.dumps(decision, indent=1))


if __name__ == "__main__":main()
