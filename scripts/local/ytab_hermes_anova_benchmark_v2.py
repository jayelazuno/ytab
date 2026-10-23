#!/usr/bin/env python3
"""Pre-registered evaluation of the Hermes multi-condition tests v2 (hermes_anova_v2 and all-pairs, both with
iterative-median recentering). Same simulator, scenarios and gates as v1 (ytab_hermes_anova_benchmark.py).

dev:    development seeds (4242, 5151), 2 experiments each; not acceptance (writes development/).
blind:  fresh seed, gates in output/essentiality/multicondition/benchmark_v2/PREREGISTRATION.md; writes decision.json.
gale:   label-free real-data run on Gale 2020 pool Cg1 (parent, no-drug control, FLU); reported only.
Simulation: one library (pool) split into K conditions with replicate cultures, or two pools each split
into all conditions; each pool has its own insertion sites; counts are negative binomial around a
lognormal clone abundance (the resampling benchmark generator), with per-condition depth, dispersion
and per-site effects; BGC normalization per pool; 300 genes per experiment.
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
from ytab.essentiality import multicondition as mc
from ytab.essentiality import resampling
from ytab.normalization.methods import normalize
from ytab.normalization.ttr import sha256, write_table

OUT = ROOT/"output/essentiality/multicondition"
BENCH = OUT/"benchmark_v2"
GENES = 300
SEEDS = dict(dev=(4242, 5151), blind=71717)   # 60606 was the v1 evaluation seed
DESIGNS = {"3x2":dict(pools=1, conditions=("C", "A", "B"), reps=2), "4x2":dict(pools=1, conditions=("C", "A", "B", "D"), reps=2),
           "3x1_two_pools":dict(pools=2, conditions=("C", "A", "B"), reps=1)}
BASE = dict(m=float(np.log(30)), tau=1.5, phi=.2, design="3x2")
SCENARIOS = {
    "null_3x2":{},
    "null_3x2_unequal_depth":dict(depth={"C":(1., 1.3), "A":(.35, .3), "B":(2.5, 2.)}),
    "null_3x2_high_dispersion":dict(phi=1.),
    "null_3x2_jackpot":dict(jackpot=True),
    "null_4x2":dict(design="4x2"),
    "null_two_pools":dict(design="3x1_two_pools"),
    "null_3x2_one_noisy_condition":dict(phi=.02, phi_by={"B":.5}),
    "depletion_4x_one_condition":dict(folds={"A":.25}),
    "enrichment_4x_one_condition":dict(folds={"A":4.}),
    "opposite_effects":dict(folds={"A":.25, "B":4.}),
    "shared_depletion_4x":dict(folds={"A":.25, "B":.25}),
    "depletion_4x_two_pools":dict(design="3x1_two_pools", folds={"A":.25}),
    "complete_loss_4x2":dict(design="4x2", folds={"D":0.}),
    "null_3x2_very_noisy_condition":dict(phi=.02, phi_by={"B":1.}),   # reported only (validity boundary)
}
REPORTED_ONLY = {"null_3x2_very_noisy_condition"}
EXPOSURES = (1., 1.3, .8, 1.1)
GALE = dict(parent="SRR11561508", control=("SRR11561505", "SRR11561504"), flu=("SRR11561503", "SRR11561502"),
            gff=ROOT/"resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff")


def simulate(name, experiment, seed):
    cfg = BASE | SCENARIOS[name]; design = DESIGNS[cfg["design"]]; conds = design["conditions"]
    rng = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(sorted(SCENARIOS).index(name), experiment)))
    affected = set(rng.choice(GENES, GENES//10, replace=False).tolist()) if cfg.get("folds") else set()
    columns = [(p, c, r) for p in range(design["pools"]) for c in conds for r in range(design["reps"])]
    blocks, idx_parts = [], [[] for _ in range(GENES)]; offset = 0
    for p in range(design["pools"]):
        sites = np.clip(np.round(np.exp(rng.uniform(np.log(3), np.log(120), GENES))), 3, 120).astype(int); total = int(sites.sum())
        lam = np.exp(cfg["m"]+cfg["tau"]*rng.standard_normal(total))
        if cfg.get("jackpot"):lam *= np.where(rng.random(total) < .05, 20., 1.)
        starts = np.concatenate([[0], np.cumsum(sites)[:-1]]); effect = {c:np.ones(total) for c in conds}
        for g in affected:
            length = max(2, int(round(.5*sites[g]))); offs = np.arange(length)+sites[g]//4; offs = offs[offs < sites[g]]
            for c, f in cfg["folds"].items():effect[c][starts[g]+offs] = f
        for g in range(GENES):idx_parts[g].append(offset+starts[g]+np.arange(sites[g]))
        block = np.zeros((total, len(columns)), dtype=np.int64)
        for j, (pp, c, r) in enumerate(columns):
            if pp != p:continue
            e = cfg["depth"][c][r] if "depth" in cfg else EXPOSURES[j % len(EXPOSURES)]
            phi = cfg.get("phi_by", {}).get(c, cfg["phi"]); size = 1/phi; mean = e*lam*effect[c]
            block[:, j] = rng.negative_binomial(size, size/(size+mean))
        blocks.append(block); offset += total
    raw = np.vstack(blocks)
    genes = [dict(id=f"g{g:04d}", name=f"g{g:04d}", description="", chrom="chr1", strand="+", blocks=[(1, 2)]) for g in range(GENES)]
    idx = [np.concatenate(x) for x in idx_parts]
    folds = {c:cfg.get("folds", {}).get(c, 1.) for c in conds}
    truth = {g:(folds if g in affected else {c:1. for c in conds}) for g in range(GENES)}
    return raw, [c for _, c, _ in columns], [f"pool{p}" for p, _, _ in columns], genes, idx, truth, conds


def run_experiment(item):
    name, experiment, seed = item
    raw, conds, libs, genes, idx, truth, names = simulate(name, experiment, seed)
    values = normalize(sparse.csc_matrix(raw), "betageom", libs)[0].toarray()
    rows, diag = mc.anova(values, conds, genes, idx, libs if len(set(libs)) > 1 else None, mc.Settings(seed=seed))
    out = dict(anova=dict(unaffected=0, p05=0, calls=0, false=0, affected=0, hits=0))
    responsive = {g for g, f in truth.items() if len(set(f.values())) > 1}
    for g, r in enumerate(rows):
        a = out["anova"]
        if not r["in_family"]:a["affected" if g in responsive else "unaffected_untested"] = a.get("affected" if g in responsive else "unaffected_untested", 0)+1; continue
        called = r["call"] == "varies"; a["calls"] += called
        if g in responsive:a["affected"] += 1; a["hits"] += called
        else:a["unaffected"] += 1; a["p05"] += r["p_value"] <= .05; a["false"] += called
    pos = np.arange(len(values), dtype=np.int64)
    prow, _ = mc.all_pairs(values, conds, genes, idx, pos, libs if len(set(libs)) > 1 else None)
    p = dict(unaffected=0, p05=0, calls=0, false=0, affected=0, hits=0); any_pair = set()
    for r in prow:
        g = int(r["ORF"][1:]); fa, fb = truth[g][r["control"]], truth[g][r["treatment"]]
        if not r.get("in_family"):continue
        direction = "depleted" if fb < fa else "enriched" if fb > fa else None; called = r["call_all_pairs"] in ("depleted", "enriched")
        p["calls"] += called
        if direction:p["affected"] += 1; p["hits"] += called and r["call_all_pairs"] == direction; p["false"] += called and r["call_all_pairs"] != direction
        else:p["unaffected"] += 1; p["p05"] += r["p_value"] <= .05; p["false"] += called
        if called:any_pair.add(g)
    p["genes_any_pair_hit"] = len(any_pair & responsive); p["genes_responsive"] = len(responsive)
    out["all_pairs"] = p
    out["anova"]["permutations_total"] = int(sum(r.get("permutations", 0) for r in rows))
    out["anova"] |= dict(experiments=1, warnings=int(diag["warning"] is not None))
    return name, experiment, out, {r["ORF"]:r.get("p_value") for r in rows}


def interval(k, n):
    if n == 0:return [None, None]
    return [float(beta.ppf(.025, k, n-k+1)) if k else 0., float(beta.ppf(.975, k+1, n-k)) if k < n else 1.]


def aggregate(results):
    agg = {}
    for name, _, out, _ in results:
        for test, sc in out.items():
            a = agg.setdefault((name, test), {})
            for k, v in sc.items():a[k] = a.get(k, 0)+v
    rows = []
    for (name, test), a in sorted(agg.items()):
        rows.append(dict(scenario=name, test=test, unaffected=a["unaffected"], type1_p05=a["p05"]/a["unaffected"] if a["unaffected"] else None,
                         type1_ci=interval(a["p05"], a["unaffected"]), calls=a["calls"], fdr=a["false"]/a["calls"] if a["calls"] else 0.,
                         sensitivity=a["hits"]/a["affected"] if a["affected"] and name.split("_")[0] != "null" else None,
                         any_pair_gene_sensitivity=(a["genes_any_pair_hit"]/a["genes_responsive"]) if a.get("genes_responsive") else None,
                         noise_warnings=f"{a['warnings']}/{a['experiments']}" if "warnings" in a else None))
    return rows


def simulate_all(seed, experiments, workers):
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(run_experiment, [(s, k, seed) for s in SCENARIOS for k in range(experiments)]))


def check_preregistration():
    frozen = (BENCH/"PREREGISTRATION.sha256").read_text().split()
    if sha256(BENCH/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def blind(workers):
    prereg = check_preregistration()
    if (BENCH/"decision.json").exists():raise SystemExit("decision.json exists; the blind evaluation already ran")
    results = simulate_all(SEEDS["blind"], 8, workers); rows = aggregate(results)
    write_table(BENCH/"synthetic.tsv", list(rows[0]), [r | dict(type1_ci=";".join(f"{x:.4f}" for x in r["type1_ci"]) if r["type1_ci"][0] is not None else "") for r in rows])
    again = run_experiment(("null_3x2", 0, SEEDS["blind"]))
    first = next(r for r in results if r[0] == "null_3x2" and r[1] == 0)
    gates = {}
    for r in rows:
        if r["scenario"] in REPORTED_ONLY:continue
        noisy = r["scenario"] == "null_3x2_one_noisy_condition"
        if r["test"] == "anova" and not noisy:gates[f"type1 anova {r['scenario']}"] = dict(value=r["type1_p05"], passed=.025 <= r["type1_p05"] <= .075)
        if r["test"] == "all_pairs":gates[f"type1 all_pairs {r['scenario']}"] = dict(value=r["type1_p05"], passed=.025 <= r["type1_p05"] <= .075)
        if not r["scenario"].startswith("null") and r["test"] in ("anova", "all_pairs"):
            gates[f"fdr {r['test']} {r['scenario']}"] = dict(value=r["fdr"], passed=r["fdr"] <= .10)
    gates["determinism"] = dict(passed=first[3] == again[3])
    warn = {(n, k):out["anova"]["warnings"] for n, k, out, _ in results}
    single_pool_nulls = [n for n in SCENARIOS if n.startswith("null") and DESIGNS[(BASE | SCENARIOS[n])["design"]]["pools"] == 1]
    noisy = {"null_3x2_one_noisy_condition", "null_3x2_very_noisy_condition"}
    gates["noise warning fires in every noisy-condition experiment"] = dict(passed=all(v for (n, _), v in warn.items() if n in noisy))
    gates["noise warning silent in every other single-pool null experiment"] = dict(
        passed=not any(v for (n, _), v in warn.items() if n in single_pool_nulls and n not in noisy))
    decision = dict(preregistration=prereg, passed=all(g["passed"] for g in gates.values()), gates=gates,
                    script_sha256=sha256(Path(__file__)), module_sha256=sha256(Path(mc.__file__)))
    (BENCH/"decision.json").write_text(json.dumps(decision, indent=2)+"\n"); print(json.dumps(decision, indent=1))


def dev(workers):
    target = OUT/"development"; target.mkdir(parents=True, exist_ok=True)
    rows = aggregate(sum((simulate_all(seed, 2, workers) for seed in SEEDS["dev"]), []))
    write_table(target/"dev_v2_seeds4242_5151.tsv", list(rows[0]), [r | dict(type1_ci=str(r["type1_ci"])) for r in rows])
    for r in rows:print({k:(round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k != "type1_ci"})


def hit_counts(sample):
    path = ROOT/"output/projects/Gale_2020/create_hit_file"/sample/f"{sample}_hits.txt"
    chrom = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=(0,), dtype=str, ndmin=1)
    pos, count = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=(11, 12), dtype=np.int64, ndmin=2).T
    totals = {}
    for c, p, n in zip(chrom, pos, count):totals[(str(c), int(p))] = totals.get((str(c), int(p)), 0)+int(n)
    return totals


def gale(workers):
    target = OUT/"Gale_2020_Cg1_v2"; target.mkdir(parents=True, exist_ok=True)
    samples = [GALE["parent"], *GALE["control"], *GALE["flu"]]; conds = ["parent", "control", "control", "FLU", "FLU"]
    counts = [hit_counts(s) for s in samples]; sites = sorted(set().union(*counts))
    raw = np.array([[c.get(s, 0) for c in counts] for s in sites], dtype=np.int64)
    values = normalize(sparse.csc_matrix(raw), "betageom", ["Cg1"]*len(samples))[0].toarray()
    genes = ccg.annotations(GALE["gff"]); settings = mc.Settings(reference="control")
    idx = ccg.gene_sites(genes, sites, settings); pos = np.array([p for _, p in sites], dtype=np.int64)
    rows, diag = mc.anova(values, conds, genes, idx, settings=settings)
    write_table(target/"anova_genes.tsv", list(dict.fromkeys(k for r in rows for k in r)), rows)
    prow, pdiag = mc.all_pairs(values, conds, genes, idx, pos)
    write_table(target/"pairs_genes.tsv", list(dict.fromkeys(k for r in prow for k in r)), prow)
    import pandas as pd
    v3 = pd.read_csv(ROOT/"output/essentiality/resampling/Gale_2020_bgc_v3/genes.tsv", sep="\t")
    v3calls = {(o, c) for o, c in zip(v3.ORF, v3.call) if c in ("depleted", "enriched")}
    flu = {(r["ORF"], r["call"]) for r in prow if r["pair"] == "FLU vs control" and r["call"] in ("depleted", "enriched")}
    varies = {r["ORF"] for r in rows if r["call"] == "varies"}
    pairs = [k for k in pdiag if k != "recenter_factors"]
    by_pair = {pair:{c:sum(r["pair"] == pair and r["call_all_pairs"] == c for r in prow) for c in ("depleted", "enriched")} for pair in pairs}
    summary = dict(samples=dict(zip(samples, conds)), sites=len(sites), anova=diag, all_pairs_calls=by_pair,
                   varies_with_any_pair_call=len(varies & {r["ORF"] for r in prow if r["call_all_pairs"] in ("depleted", "enriched")}),
                   varies_top_condition={c:sum(r.get("top_condition") == c for r in rows if r["call"] == "varies") for c in ("parent", "FLU")},
                   flu_vs_control_per_pair_q_vs_published_v3=dict(this=len(flu), v3=len(v3calls), shared=len(flu & v3calls),
                                                                   jaccard=len(flu & v3calls)/len(flu | v3calls)),
                   all_pairs_recenter_factors=pdiag["recenter_factors"],
                   pair_diagnostics={k:(pdiag[k] if isinstance(pdiag[k], str) else {x:pdiag[k][x] for x in ("replicate_noise", "family_size", "calls")}) for k in pairs})
    (target/"summary.json").write_text(json.dumps(ccg._json_ready(summary), indent=2)+"\n"); print(json.dumps(ccg._json_ready(summary), indent=1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("phase", choices=["dev", "blind", "gale"]); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(); dict(dev=dev, blind=blind, gale=gale)[a.phase](a.workers)


if __name__ == "__main__":main()
