#!/usr/bin/env python3
"""Pre-registered blinded evaluation of stratified (paired independent-library) Hermes resampling.

blind:   zinc-like synthetic calibration + label-free zinc diagnostics; applies the gates in
         benchmark_v2_paired/PREREGISTRATION.md; writes decision.json. Reads no fitness, z-score,
         classifier, Gale or Gumbel output.
unblind: only after a pass; compares zinc calls with the existing zinc z-score fitness results.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.essentiality import resampling
from ytab.normalization.methods import normalize
from ytab.normalization.ttr import sha256, write_table

BENCH = ROOT/"output/essentiality/resampling/benchmark_v2_paired"
GENES = 300; EXPERIMENTS = 8; SEED = 3031
POOLS = [(1., 1.2), (.8, 1.), (1.1, .7), (.9, 1.3)]
BASE = dict(m=float(np.log(30)), tau=1.5, phi=.2)
SCENARIOS = {
    "null_equal_depth": dict(exposures=[(1., 1.)]*4),
    "null_shallow_treatment": dict(treated_scale=.3),
    "null_deep_treatment": dict(treated_scale=2.2),
    "null_low_read": dict(m=float(np.log(3))),
    "null_high_dispersion": dict(phi=1.),
    "null_jackpot_prior": dict(jackpot=True),
    "null_noisy_treatment": dict(phi=.02, treatment_phi=.5),
    "null_bottleneck_treatment": dict(phi=.02, treatment_phi=.05, treatment_linear=4.),
    "complete_loss_block": dict(fold=0.),
    "depletion_4x_block": dict(fold=.25),
    "depletion_2x_block": dict(fold=.5),
    "scattered_complete_loss": dict(fold=0., scattered=True),
    "depletion_4x_noisy_treatment": dict(fold=.25, phi=.02, treatment_phi=.5),
    "depletion_4x_bottleneck": dict(fold=.25, phi=.02, treatment_phi=.05, treatment_linear=4.),
    "uniform_depletion_4x_whole_gene": dict(fold=.25, whole=True),
    "uniform_depletion_2x_whole_gene": dict(fold=.5, whole=True),
    "enrichment_4x_block": dict(fold=4.),
}
SETTINGS = resampling.Settings(inference="exact", opportunities="any-support")
LIBS = ["L1", "L2", "L3", "L4"]
ZN = dict(controls=[f"Zn_toxicity_screen::yH29{8 if p < 3 else 9}-parent-pool{p}-mock" for p in range(1, 5)],
          treatments=[f"Zn_toxicity_screen::yH29{8 if p < 3 else 9}-parent-pool{p}-1_5mM-Zn-treated" for p in range(1, 5)],
          normalizations={"ttr":"output/normalization/zn_ttr_validation_v1/methods/ttr", "totreads":"output/normalization/zn_totreads_sites_v1/totreads",
                          "nzmean":"output/normalization/zn_nzmean_sites_v1/nzmean"},
          gff="resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff",
          reference="resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.fna")


def check_preregistration():
    frozen = (BENCH/"PREREGISTRATION.sha256").read_text().split()
    if sha256(BENCH/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def interval(k, n):
    return ((float(beta.ppf(.025, k, n-k+1)) if k else 0.), (float(beta.ppf(.975, k+1, n-k)) if k < n else 1.)) if n else (None, None)


def make_data(name, experiment):
    """Four independent pools sharing gene layout and biological effects; returns raw, TTR-normalized, oracle-normalized."""
    config = BASE | SCENARIOS[name]
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(sorted(SCENARIOS).index(name), experiment)))
    sites = np.clip(np.round(np.exp(rng.uniform(np.log(3), np.log(120), GENES))), 3, 120).astype(int)
    effects = np.ones(sites.sum()); truth = np.zeros(GENES, dtype=int)
    if "fold" in config:
        affected = rng.choice(GENES, GENES//10, replace=False); starts = np.concatenate([[0], np.cumsum(sites)[:-1]])
        for g in affected:
            if config.get("whole"):offsets = np.arange(sites[g])
            else:
                length = max(2, int(round(.5*sites[g])))
                offsets = (np.arange(0, 2*length, 2) if config.get("scattered") else np.arange(length))+sites[g]//4
            effects[starts[g]+offsets[offsets < sites[g]]] = config["fold"]
        truth[affected] = 1 if config["fold"] > 1 else -1
    pools = [(m, t*config.get("treated_scale", 1.)) for m, t in config.get("exposures", POOLS)]
    columns, true_depth = [], []
    for mock, treated in pools:
        e = np.array([mock, treated])
        if config.get("jackpot"):
            lam = np.exp(config["m"]+config["tau"]*rng.standard_normal(sites.sum()))*np.where(rng.random(sites.sum()) < .05, 20., 1.)
            size = 1/config["phi"]; draw = lambda mean:rng.negative_binomial(size, size/(size+mean))
            _, _, genes, idx, pos = ccg.synthetic_experiment(rng, sites, e, treatments=1)
            columns += [draw(mock*lam), draw(treated*lam)]
        else:
            C, T, genes, idx, pos = ccg.synthetic_experiment(rng, sites, e, config["m"], config["tau"], config["phi"], effects, treatments=1,
                                                             treatment_phi=config.get("treatment_phi"), treatment_linear=config.get("treatment_linear"))
            columns += [C[:, 0], T]
        true_depth += [mock, treated]
    raw = np.column_stack(columns)   # mock1, treated1, mock2, treated2, ...
    scaled, _ = normalize(sparse.csc_matrix(raw), "ttr", [p for p in LIBS for _ in range(2)])
    return raw, scaled.toarray(), raw/np.array(true_depth), genes, idx, pos, truth


def split(matrix):return matrix[:, 0::2], matrix[:, 1::2]


def score(rows, truth):
    out = dict(unaffected=0, p05=0, false=0, calls=0, affected=0, hits=0)
    for row, t in zip(rows, truth):
        if not row.get("in_family"):
            out["affected"] += t != 0; continue
        called = row["call"] in {"depleted", "enriched"}; direction = 1 if row["call"] == "enriched" else -1 if row["call"] == "depleted" else 0
        out["calls"] += called
        if t == 0:out["unaffected"] += 1; out["p05"] += row["p_value"] <= .05; out["false"] += called
        else:out["affected"] += 1; out["hits"] += called and direction == t; out["false"] += called and direction != t
    return out


def aggregate(runs):
    total = {k:sum(r[k] for r in runs) for k in ("unaffected", "p05", "false", "calls", "affected", "hits")}
    total |= dict(experiments_with_false_call=sum(r["false"] > 0 for r in runs), p05_rate=total["p05"]/max(total["unaffected"], 1),
                  p05_ci=interval(total["p05"], total["unaffected"]), sensitivity=total["hits"]/total["affected"] if total["affected"] else None,
                  fdr=total["false"]/total["calls"] if total["calls"] else None)
    return total


def synthetic_task(item):
    name, k = item
    raw, ttr, oracle, genes, idx, pos, truth = make_data(name, k)
    out = {}
    for arm, matrix in (("ttr", ttr), ("oracle", oracle)):
        control, treatment = split(matrix)
        out[arm] = score(resampling.analyze(control, treatment, genes, idx, pos, SETTINGS, LIBS, LIBS)[0], truth)
    return name, k, out


def zinc_arrays(normalization):
    design = ccg.load_design(ROOT/ZN["normalizations"][normalization], normalization, ROOT/ZN["reference"], BENCH/"_unused_output",
                             ZN["controls"], ZN["treatments"], min_controls=1, stratified=True)
    genes = ccg.annotations(ROOT/ZN["gff"]); idx = ccg.gene_sites(genes, design["sites"], SETTINGS)
    return design, genes, idx, np.array([p for _, p in design["sites"]], dtype=np.int64)


def zinc_task(item):
    normalization, seed, swap = item
    design, genes, idx, pos = zinc_arrays(normalization)
    s = design["scaled"]; control, treatment = s[:, :4].copy(), s[:, 4:].copy()
    if swap:   # flip labels in pools 2 and 4 (one per parental strain): true Zn effects cancel across pools
        for p in (1, 3):control[:, p], treatment[:, p] = s[:, 4+p], s[:, p]
    libs = design["libraries"]
    rows, diagnostics = resampling.analyze(control, treatment, genes, idx, pos, resampling.Settings(**{**SETTINGS.__dict__, "seed":seed}), libs[:4], libs[4:])
    return normalization, seed, swap, rows, diagnostics


def jaccard(a, b):return len(a & b)/len(a | b) if a | b else 1.


def blind(workers):
    prereg = check_preregistration()
    if (BENCH/"decision.json").exists():raise SystemExit("decision.json exists; phase 1 already ran")
    with ProcessPoolExecutor(workers) as pool:synthetic = list(pool.map(synthetic_task, [(n, k) for n in SCENARIOS for k in range(EXPERIMENTS)]))
    with ProcessPoolExecutor(workers) as pool:
        zinc = list(pool.map(zinc_task, [(n, s, False) for n in ZN["normalizations"] for s in (2026, 2027)]+[("ttr", 2026, True)]))
    nulls = [n for n in SCENARIOS if n.startswith("null_")]; signals = [n for n in SCENARIOS if n not in nulls]
    metrics = {arm:{n:aggregate([r[2][arm] for r in synthetic if r[0] == n]) for n in SCENARIOS} for arm in ("ttr", "oracle")}
    rows = [dict(arm=arm, scenario=n, unaffected_p_le_05=m["p05_rate"], ci_low=m["p05_ci"][0], ci_high=m["p05_ci"][1], false_calls=m["false"], calls=m["calls"],
                 experiments_with_false_call=m["experiments_with_false_call"], sensitivity=m["sensitivity"], fdr=m["fdr"])
            for arm in metrics for n, m in metrics[arm].items()]
    write_table(BENCH/"synthetic.tsv", list(rows[0]), rows)
    (BENCH/"zinc_blinded").mkdir(exist_ok=True); calls, family = {}, {}
    for norm, seed, swap, zrows, diagnostics in zinc:
        tag = f"{norm}_seed{seed}"+("_swapnull" if swap else "")
        calls[tag] = {r["ORF"] for r in zrows if r["call"] in {"depleted", "enriched"}}; family[tag] = diagnostics["family_size"]
        keys = ["ORF", "Name", "n_sites", "libraries_with_sites", "mean_control", "mean_treatment", "log2FC", "library_log2FC", "delta_mean", "p_value", "q_bh", "call", "status"]
        write_table(BENCH/"zinc_blinded"/f"{tag}.tsv", keys, [{k:r.get(k) for k in keys} for r in zrows])
    ttr = metrics["ttr"]; null_runs = [r[2]["ttr"] for r in synthetic if r[0] in nulls]
    fdr = sum(ttr[n]["false"] for n in signals)/max(sum(ttr[n]["calls"] for n in signals), 1)
    seed_j = {n:jaccard(calls[f"{n}_seed2026"], calls[f"{n}_seed2027"]) for n in ZN["normalizations"]}
    swap_fraction = len(calls["ttr_seed2026_swapnull"])/family["ttr_seed2026_swapnull"]
    gates = dict(g1a=all(ttr[n]["p05_ci"][0] <= .05 for n in nulls), g1b=sum(r["false"] > 0 for r in null_runs)/len(null_runs) <= .10,
                 g1c=fdr <= .10, g2=min(seed_j.values()) >= .99, g3=swap_fraction <= .005)
    decision = dict(status="success", phase="blind", preregistration=prereg, passed=all(gates.values()), gates=gates,
        failing_null_scenarios=[n for n in nulls if ttr[n]["p05_ci"][0] > .05], pooled_signal_fdr=fdr,
        null_experiments_with_false_call=sum(r["false"] > 0 for r in null_runs)/len(null_runs),
        mean_signal_sensitivity=float(np.mean([ttr[n]["sensitivity"] for n in signals])), seed_jaccard=seed_j,
        swap_null_calls=len(calls["ttr_seed2026_swapnull"]), swap_null_fraction=swap_fraction,
        reported_only=dict(zinc_calls={n:len(calls[f"{n}_seed2026"]) for n in ZN["normalizations"]},
            normalization_jaccard={f"{a}-{b}":jaccard(calls[f"{a}_seed2026"], calls[f"{b}_seed2026"]) for a, b in (("ttr", "totreads"), ("ttr", "nzmean"), ("totreads", "nzmean"))},
            oracle_pooled_signal_fdr=sum(metrics["oracle"][n]["false"] for n in signals)/max(sum(metrics["oracle"][n]["calls"] for n in signals), 1),
            oracle_failing_null_scenarios=[n for n in nulls if metrics["oracle"][n]["p05_ci"][0] > .05]),
        settings=SETTINGS.__dict__, blinding="no fitness, z-score, classifier, GO, Gale or Gumbel outputs were read",
        implementation_sha256=sha256(Path(resampling.__file__)), script_sha256=sha256(Path(__file__)))
    (BENCH/"decision.json").write_text(json.dumps(ccg._json_ready(decision), indent=2)+"\n")
    print(json.dumps(ccg._json_ready({k:decision[k] for k in ("passed", "gates", "failing_null_scenarios", "pooled_signal_fdr", "null_experiments_with_false_call",
                                                              "mean_signal_sensitivity", "seed_jaccard", "swap_null_calls", "swap_null_fraction", "reported_only")}), indent=1))


def unblind():
    prereg = check_preregistration(); decision = json.loads((BENCH/"decision.json").read_text())
    if decision["preregistration"]["sha256"] != prereg["sha256"]:raise SystemExit("decision made under a different pre-registration")
    if not decision["passed"]:raise SystemExit("phase 1 did not pass; nothing to unblind")
    raise SystemExit("unblinding: locate the zinc z-score table only now, then compare (implemented after a pass)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["blind", "unblind"]); parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    blind(args.workers) if args.phase == "blind" else unblind()


if __name__ == "__main__":main()
