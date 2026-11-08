#!/usr/bin/env python3
"""R28: aBGC against BGC for the multi-condition tests (R11 simulator) and genetic interactions (R12a simulator).

Each simulated experiment is normalized twice, BGC and aBGC, and both arms are tested on the same counts (paired).
Simulators, scenarios, tests and gates are those of R11 (ytab_hermes_anova_benchmark_v2.py) and R12a
(ytab_gi_benchmark.py), imported unchanged; only the normalization differs. Gates and decision rule:
docs/manuscript/preregistration_abgc_conditions.md.

dev:    development seeds, 1 experiment per scenario; checks the code only (writes output/essentiality/abgc_r28/development/).
blind:  evaluation seeds, 8 experiments per scenario; writes output/essentiality/abgc_r28/{synthetic_*.tsv,decision.json}.
  python scripts/local/ytab_r28_abgc_benchmark.py blind [--workers 6]
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.normalization import methods as norm_methods  # noqa: E402
from ytab.normalization.ttr import sha256, write_table  # noqa: E402

OUT = ROOT/"output/essentiality/abgc_r28"
PREREG = ROOT/"docs/manuscript/preregistration_abgc_conditions.md"
SOURCES = dict(mc=ROOT/"scripts/local/ytab_hermes_anova_benchmark_v2.py", gi=ROOT/"scripts/local/ytab_gi_benchmark.py")
ARMS = ("betageom", "aBGC")
SEEDS = dict(dev=dict(mc=2801, gi=2802), blind=dict(mc=28028, gi=28029))   # R11 used 71717, R12a 93939
EXPERIMENTS = dict(dev=1, blind=8)
AUC_MARGIN, BOOTSTRAPS = -.01, 10000
_modules = {}


def module(family, arm):
    """The R11 or R12a benchmark module with its normalization fixed to `arm` (loaded once per process)."""
    key = (family, arm)
    if key not in _modules:
        spec = importlib.util.spec_from_file_location(f"r28_{family}_{arm}", SOURCES[family])
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
        mod.normalize = lambda matrix, _method, groups, *a, **k: norm_methods.normalize(matrix, arm, groups, *a, **k)
        _modules[key] = mod
    return _modules[key]


def auc(scores, truth):
    """Probability that an affected gene has a smaller p than an unaffected one (ties count half)."""
    s, t = np.asarray(scores, float), np.asarray(truth, bool); ok = np.isfinite(s); s, t = s[ok], t[ok]
    if t.all() or not t.any():return np.nan
    r = rankdata(s)
    return float((r[t].sum()-t.sum()*(t.sum()+1)/2)/(t.sum()*(~t).sum()))


def work(item):
    family, arm, args = item; mod = module(family, arm)
    if family == "mc":
        name, k, seed = args; result = mod.run_experiment(args)
        truth = mod.simulate(name, k, seed)[5]
        responsive = [len(set(truth[g].values())) > 1 for g in range(mod.GENES)]
        p = [result[3].get(f"g{g:04d}") for g in range(mod.GENES)]
    else:
        name, design, k, seed = args; result = mod.task(args)
        responsive = mod.simulate(name, design, k, seed)[5]; p = result[4]
    p = np.array([np.nan if x is None else float(x) for x in p])
    return family, arm, args, result, auc(-p, responsive) if any(responsive) else np.nan


def items(family, phase):
    seed, n = SEEDS[phase][family], EXPERIMENTS[phase]; mod = module(family, "betageom")
    base = [(s, k, seed) for s in mod.SCENARIOS for k in range(n)] if family == "mc" else mod.items(seed, n)
    return [(family, arm, b) for arm in ARMS for b in base]


def mc_gates(mod, results, again, first):
    rows = mod.aggregate(results); gates = {}
    for r in rows:
        if r["scenario"] in mod.REPORTED_ONLY:continue
        noisy = r["scenario"] == "null_3x2_one_noisy_condition"
        if r["test"] == "anova" and not noisy:gates[f"type1 anova {r['scenario']}"] = dict(value=r["type1_p05"], passed=.025 <= r["type1_p05"] <= .075)
        if r["test"] == "all_pairs":gates[f"type1 all_pairs {r['scenario']}"] = dict(value=r["type1_p05"], passed=.025 <= r["type1_p05"] <= .075)
        if not r["scenario"].startswith("null"):gates[f"fdr {r['test']} {r['scenario']}"] = dict(value=r["fdr"], passed=r["fdr"] <= .10)
    gates["determinism"] = dict(passed=first[3] == again[3])
    warn = {(n, k):out["anova"]["warnings"] for n, k, out, _ in results}
    single = [n for n in mod.SCENARIOS if n.startswith("null") and mod.DESIGNS[(mod.BASE | mod.SCENARIOS[n])["design"]]["pools"] == 1]
    noisy = {"null_3x2_one_noisy_condition", "null_3x2_very_noisy_condition"}
    gates["noise warning fires in every noisy-condition experiment"] = dict(passed=all(v for (n, _), v in warn.items() if n in noisy))
    gates["noise warning silent in every other single-pool null experiment"] = dict(passed=not any(v for (n, _), v in warn.items() if n in single and n not in noisy))
    return rows, gates


def gi_gates(mod, results, again, first):
    rows = mod.aggregate(results); gates = {}
    for r in rows:
        if r["scenario"] in mod.REPORTED_ONLY:continue
        gates[f"type1 {r['scenario']} {r['design']}"] = dict(value=r["type1_p05"], passed=.025 <= r["type1_p05"] <= .075)
        if r["scenario"].startswith("gi_"):gates[f"fdr {r['scenario']} {r['design']}"] = dict(value=r["fdr"], passed=r["fdr"] <= .10)
    gates["determinism"] = dict(passed=again[4] == first[4])
    return rows, gates


def auc_difference(by_arm, seed):
    """Paired aBGC - BGC AUC over experiments with an effect: mean and percentile bootstrap 95% CI (experiments resampled)."""
    keys = sorted(k for k in by_arm["aBGC"] if np.isfinite(by_arm["aBGC"][k]) and np.isfinite(by_arm["betageom"].get(k, np.nan)))
    d = np.array([by_arm["aBGC"][k]-by_arm["betageom"][k] for k in keys])
    rng = np.random.default_rng(seed); boot = d[rng.integers(0, len(d), (BOOTSTRAPS, len(d)))].mean(axis=1)
    return dict(experiments=len(d), auc_bgc=float(np.mean([by_arm["betageom"][k] for k in keys])), auc_abgc=float(np.mean([by_arm["aBGC"][k] for k in keys])),
                difference=float(d.mean()), ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))])


def evaluate(phase, workers):
    target = OUT if phase == "blind" else OUT/"development"; target.mkdir(parents=True, exist_ok=True)
    jobs = items("mc", phase)+items("gi", phase)
    with ProcessPoolExecutor(workers) as pool:done = list(pool.map(work, jobs, chunksize=1))
    decision = dict(phase=phase, seeds=SEEDS[phase], experiments_per_scenario=EXPERIMENTS[phase], families={})
    for family in ("mc", "gi"):
        fam = {}; by_arm = {a:{} for a in ARMS}
        for arm in ARMS:
            mod = module(family, arm); results = [r for f, a, _, r, _ in done if f == family and a == arm]
            for f, a, args, _, value in done:
                if f == family and a == arm:by_arm[arm][tuple(args[:-1])] = value
            first_args = ("null_3x2", 0, SEEDS[phase]["mc"]) if family == "mc" else ("gi_mixed", "split", 0, SEEDS[phase]["gi"])
            again = work((family, arm, first_args))[3]
            first = next(r for f, a, args, r, _ in done if f == family and a == arm and tuple(args) == first_args)
            rows, gates = (mc_gates if family == "mc" else gi_gates)(mod, results, again, first)
            write_table(target/f"synthetic_{family}_{arm}.tsv", list(rows[0]), [r | dict(type1_ci=";".join(f"{x:.4f}" for x in r["type1_ci"]) if r["type1_ci"][0] is not None else "") for r in rows])
            fam[arm] = dict(gates_passed=all(g["passed"] for g in gates.values()), failed=[k for k, g in gates.items() if not g["passed"]], gates=gates)
        fam["auc"] = auc_difference(by_arm, SEEDS[phase][family])
        fam["abgc_becomes_default"] = bool(fam["aBGC"]["gates_passed"] and fam["auc"]["ci95"][0] > AUC_MARGIN)
        decision["families"]["multicondition" if family == "mc" else "genetic_interactions"] = fam
    if phase == "blind":
        frozen = (PREREG.with_suffix(".sha256")).read_text().split()
        decision["preregistration"] = dict(sha256=frozen[0], frozen_utc=frozen[1])
    decision["script_sha256"] = sha256(Path(__file__))
    (target/"decision.json").write_text(json.dumps(decision, indent=2, default=float)+"\n")
    summary = {f:dict(bgc_gates=v["betageom"]["gates_passed"], abgc_gates=v["aBGC"]["gates_passed"], abgc_failed=v["aBGC"]["failed"], auc=v["auc"],
                      abgc_becomes_default=v["abgc_becomes_default"]) for f, v in decision["families"].items()}
    print(json.dumps(summary, indent=1, default=float))


def blind(workers):
    frozen = PREREG.with_suffix(".sha256")
    if not frozen.exists():raise SystemExit("Pre-registration not frozen; refusing to run")
    lines = [l.split() for l in frozen.read_text().splitlines() if l.strip()]
    if sha256(PREREG) != lines[0][0]:raise SystemExit("Pre-registration changed after freezing; refusing to run")
    changed = [path for digest, path in (l[:2] for l in lines[1:]) if sha256(ROOT/path) != digest]
    if changed:raise SystemExit(f"Frozen code changed: {', '.join(changed)}; refusing to run")
    if (OUT/"decision.json").exists():raise SystemExit("decision.json exists; the blind evaluation already ran")
    evaluate("blind", workers)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("phase", choices=["dev", "blind"]); ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    evaluate("dev", a.workers) if a.phase == "dev" else blind(a.workers)
