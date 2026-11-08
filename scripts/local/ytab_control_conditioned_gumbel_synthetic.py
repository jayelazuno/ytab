#!/usr/bin/env python3
"""Stored synthetic calibration and sensitivity evidence for the control-conditioned Gumbel extension.

Two control and two treatment replicate tracks per experiment, combined within condition as in the method.

Generic simulation settings; never tuned to Gale or any published labels.
Repeated-experiment scenarios redraw controls and treatment; the conditional scenario fixes
controls and true abundances and redraws only the treatment.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import beta

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"src"))
from ytab.essentiality.control_conditioned_gumbel import Settings, analyze, synthetic_experiment
from ytab.normalization.ttr import sha256, write_table

GENES = 300
BASE = dict(m=float(np.log(30)), tau=1.5, phi=.2, exposures=(1., 1.3, .8, 1.1))   # two controls, two treatments
SCENARIOS = {
    "null_equal_depth": dict(exposures=(1., 1., 1., 1.)),
    "null_unequal_depth_shallow_treatment": dict(exposures=(1., 1.3, .35, .3)),
    "null_unequal_depth_deep_treatment": dict(exposures=(1., .8, 2.5, 2.)),
    "null_low_read_controls": dict(m=float(np.log(3))),
    "null_high_dispersion": dict(phi=1.),
    "null_misspecified_heavy_tailed_abundance": dict(prior="jackpot"),
    "null_conditional_fixed_controls": dict(conditional=True),
    "null_noisy_treatment_replicates": dict(phi=.02, treatment_phi=.5),
    "null_noisy_treatment_assumed_control_dispersion": dict(phi=.02, treatment_phi=.5, treatment_dispersion="control"),
    "partial_depletion_0.25_noisy_treatment": dict(phi=.02, treatment_phi=.5, fraction=.1, fold=.25, block=.5),
    "null_bottleneck_treatment_replicates": dict(phi=.02, treatment_phi=.05, treatment_linear=4.),
    "null_bottleneck_treatment_quadratic_only": dict(phi=.02, treatment_phi=.05, treatment_linear=4., dispersion_model="quadratic"),
    "partial_depletion_0.25_bottleneck_treatment": dict(phi=.02, treatment_phi=.05, treatment_linear=4., fraction=.1, fold=.25, block=.5),
    "complete_loss_block": dict(fraction=.1, fold=0., block=.5),
    "partial_depletion_0.25_block": dict(fraction=.1, fold=.25, block=.5),
    "partial_depletion_0.5_block": dict(fraction=.1, fold=.5, block=.5),
    "isolated_complete_loss_scattered": dict(fraction=.1, fold=0., block=.5, scattered=True),
}


def interval(successes, trials):
    """Clopper-Pearson 95% interval."""
    if not trials:return None, None
    lo = beta.ppf(.025, successes, trials-successes+1) if successes else 0.
    hi = beta.ppf(.975, successes+1, trials-successes) if successes < trials else 1.
    return float(lo), float(hi)


def make_data(name, spec, seed):
    config = BASE | spec; rng = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(sorted(SCENARIOS).index(name),)))
    sites = np.clip(np.round(np.exp(rng.uniform(np.log(3), np.log(120), GENES))), 3, 120).astype(int)
    effects = np.ones(sites.sum()); affected = np.zeros(GENES, dtype=bool)
    if config.get("fraction"):
        affected[rng.choice(GENES, int(config["fraction"]*GENES), replace=False)] = True
        starts = np.concatenate([[0], np.cumsum(sites)[:-1]])
        for g in np.flatnonzero(affected):
            length = max(2, int(round(config["block"]*sites[g])))
            offsets = (np.arange(0, 2*length, 2) if config.get("scattered") else np.arange(length))+sites[g]//4
            effects[starts[g]+offsets[offsets < sites[g]]] = config["fold"]
    if config.get("prior") == "jackpot":
        # Abundance prior deliberately not lognormal: 5% jackpot-like sites at 20x the lognormal abundance.
        lam = np.exp(config["m"]+config["tau"]*rng.standard_normal(sites.sum()))*np.where(rng.random(sites.sum()) < .05, 20., 1.)
        size = 1/config["phi"]; e = np.asarray(config["exposures"])
        draw = lambda mean:rng.negative_binomial(size, size/(size+mean))
        _, _, genes, idx, pos = synthetic_experiment(rng, sites, e, treatments=2)
        controls = np.column_stack([draw(x*lam) for x in e[:2]]); treatment = np.column_stack([draw(x*lam) for x in e[2:]])
        return controls, treatment, genes, idx, pos, affected, e
    e = np.asarray(config["exposures"])
    return (*synthetic_experiment(rng, sites, e, config["m"], config["tau"], config["phi"], effects,
                                  treatments=2, treatment_phi=config.get("treatment_phi"),
                                  treatment_linear=config.get("treatment_linear")), affected, e)


def one(task):
    name, experiment, seed, settings = task
    spec = SCENARIOS[name]
    if spec.get("conditional"):
        # Same controls and true abundances every replicate; only treatment counts are redrawn.
        rng = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(99,)))
        sites = np.clip(np.round(np.exp(rng.uniform(np.log(3), np.log(120), GENES))), 3, 120).astype(int)
        e = np.asarray(BASE["exposures"]); lam = np.exp(BASE["m"]+BASE["tau"]*rng.standard_normal(sites.sum()))
        size = 1/BASE["phi"]
        controls = np.column_stack([rng.negative_binomial(size, size/(size+x*lam)) for x in e[:2]])
        _, _, genes, idx, pos = synthetic_experiment(np.random.default_rng(0), sites, e, treatments=2)
        treatment_rng = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(100, experiment)))
        treatment = np.column_stack([treatment_rng.negative_binomial(size, size/(size+x*lam)) for x in e[2:]]); affected = np.zeros(GENES, dtype=bool)
    else:
        controls, treatment, genes, idx, pos, affected, e = make_data(name, spec, seed*1000+experiment)
    for field in ("treatment_dispersion", "dispersion_model"):
        if spec.get(field):settings = replace(settings, **{field:spec[field]})
    rows, _, diagnostics = analyze(controls, treatment, e, genes, idx, pos, settings)
    return name, experiment, rows, affected, diagnostics


def summarize(name, results):
    tested = called = false = true_called = truth = 0; p05 = p01 = 0; experiments_with_false = 0
    inadequate = unvalidated = 0; strata = {}; observed = expected = 0.
    for _, _, rows, affected, diagnostics in results:
        observed += diagnostics["observed_depleted_sites"]; expected += diagnostics["expected_depleted_sites_null"]
        inadequate += diagnostics["family_null_fit_inadequate"]; unvalidated += diagnostics["family_tail_magnitude_unvalidated"]
        false_here = 0
        for row, hit in zip(rows, affected):
            if not row["in_family"]:continue
            positive = row["call"] == "depleted run"
            if hit:
                truth += 1; true_called += positive
                stratum = "3-19 sites" if row["n_eligible"] < 20 else "20-59 sites" if row["n_eligible"] < 60 else ">=60 sites"
                hits, total = strata.get(stratum, (0, 0)); strata[stratum] = (hits+positive, total+1)
            else:
                tested += 1; false_here += positive; p05 += row["p_family"] <= .05; p01 += row["p_family"] <= .01
        called += sum(r["call"] == "depleted run" for r in rows); false += false_here; experiments_with_false += false_here > 0
    out = dict(scenario=name, experiments=len(results), unaffected_tests=tested,
        unaffected_p_le_05=p05/max(tested, 1), unaffected_p_le_05_ci=interval(p05, tested),
        unaffected_p_le_01=p01/max(tested, 1), unaffected_p_le_01_ci=interval(p01, tested),
        false_calls=false, experiments_with_any_false_call=experiments_with_false,
        empirical_fdr=false/called if called else None, calls=called,
        site_flags_observed_over_expected=observed/expected if expected else None,
        null_fit_inadequate_genes=inadequate, tail_magnitude_unvalidated_genes=unvalidated)
    if truth:
        out |= dict(affected_tests=truth, sensitivity=true_called/truth, sensitivity_ci=interval(true_called, truth),
                    sensitivity_by_size={k:dict(sensitivity=h/t, ci=interval(h, t), genes=t) for k, (h, t) in sorted(strata.items())})
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--experiments", type=int, default=8)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--scenario", action="append", choices=sorted(SCENARIOS))
    parser.add_argument("--inference", choices=["exact", "gumbel-mc"], default=Settings().inference)
    parser.add_argument("--calibration-sims", type=int, default=None, help="default 1000 for exact, 20000 for gumbel-mc")
    parser.add_argument("--validation-sims", type=int, default=None)
    args = parser.parse_args()
    output = Path(args.output)
    if output.exists() and any(output.iterdir()):parser.exit(1, "ERROR: choose a new output directory\n")
    settings = Settings(inference=args.inference, calibration_sims=args.calibration_sims, validation_sims=args.validation_sims); names = args.scenario or list(SCENARIOS)
    tasks = [(name, k, args.seed, settings) for name in names for k in range(args.experiments)]
    with ProcessPoolExecutor(args.workers) as pool:results = list(pool.map(one, tasks))
    summaries = [summarize(name, [r for r in results if r[0] == name]) for name in names]
    output.mkdir(parents=True, exist_ok=True)
    flat = [{k:(json.dumps(v) if isinstance(v, (dict, tuple, list)) else v) for k, v in s.items()} for s in summaries]
    write_table(output/"scenarios.tsv", list(dict.fromkeys(k for s in flat for k in s)), flat)
    (output/"summary.json").write_text(json.dumps(dict(status="success", settings=asdict(settings), genes_per_experiment=GENES,
        base_model=BASE, scenarios=SCENARIOS, experiments=args.experiments, seed=args.seed, results=summaries,
        implementation_sha256=sha256(Path(sys.modules[analyze.__module__].__file__)), script_sha256=sha256(__file__),
        notes=["Generic generative settings; not fitted to or tuned on Gale or publication labels",
               "Repeated-experiment scenarios redraw controls; the conditional scenario redraws only treatment",
               "Rates carry Clopper-Pearson 95% intervals; small tails are not validated by these counts"]), indent=2)+"\n")
    for s in summaries:
        print(s["scenario"], f"p<=.05:{s['unaffected_p_le_05']:.4f}", f"falsecalls:{s['false_calls']}",
              f"sens:{s.get('sensitivity')}", f"obs/exp:{s['site_flags_observed_over_expected']:.3f}", f"inadequate:{s['null_fit_inadequate_genes']}")


if __name__ == "__main__":main()
