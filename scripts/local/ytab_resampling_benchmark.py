#!/usr/bin/env python3
"""Pre-registered, blinded selection of the Hermes resampling configuration, then unblinded benchmark.

blind:   synthetic calibration/sensitivity + label-free Gale diagnostics for the four candidates;
         applies PREREGISTRATION.md's gates and rule; writes selection.json. Reads no publication
         labels and no control-conditioned Gumbel output.
unblind: requires selection.json (and an unchanged pre-registration); compares the selected
         configuration with the control-conditioned Gumbel and with published Gale scores.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy import sparse
from scipy.stats import beta

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"src"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.essentiality import resampling
from ytab.normalization.methods import normalize
from ytab.normalization.ttr import sha256, write_table

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT/"output/essentiality/resampling/benchmark"
GENES = 300; EXPERIMENTS = 8; SEED = 2026
BASE = dict(m=float(np.log(30)), tau=1.5, phi=.2, exposures=(1., 1.3, .8, 1.1))
SCENARIOS = {
    "null_equal_depth": dict(exposures=(1., 1., 1., 1.)),
    "null_shallow_treatment": dict(exposures=(1., 1.3, .35, .3)),
    "null_deep_treatment": dict(exposures=(1., .8, 2.5, 2.)),
    "null_low_read_controls": dict(m=float(np.log(3))),
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
CANDIDATES = {"exact_any":dict(inference="exact", opportunities="any-support"), "exact_control":dict(inference="exact", opportunities="control-support"),
              "mc_any":dict(inference="monte-carlo", opportunities="any-support"), "mc_control":dict(inference="monte-carlo", opportunities="control-support")}
GALE = dict(controls=["Gale_2020::SRR11561505", "Gale_2020::SRR11561504"], treatments=["Gale_2020::SRR11561503", "Gale_2020::SRR11561502"],
            normalizations={"ttr":"output/normalization/gale_flu_ttr_validation_v1/methods/ttr", "totreads":"output/normalization/gale_flu_totreads_sites_v1/totreads",
                            "nzmean":"output/normalization/gale_flu_nzmean_sites_v1/nzmean"},
            gff="resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff",
            reference="resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.fna")


def check_preregistration():
    text = (BENCH/"PREREGISTRATION.sha256").read_text().split()
    if sha256(BENCH/"PREREGISTRATION.md") != text[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=text[0], frozen_utc=text[1])


def interval(k, n):
    if not n:return None, None
    return (float(beta.ppf(.025, k, n-k+1)) if k else 0., float(beta.ppf(.975, k+1, n-k)) if k < n else 1.)


def make_data(name, experiment):
    """Raw counts, TTR-normalized counts, genes and truth (+1 enriched, -1 depleted, 0 unaffected)."""
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
    e = np.asarray(config["exposures"])
    if config.get("jackpot"):
        lam = np.exp(config["m"]+config["tau"]*rng.standard_normal(sites.sum()))*np.where(rng.random(sites.sum()) < .05, 20., 1.)
        size = 1/config["phi"]; draw = lambda mean:rng.negative_binomial(size, size/(size+mean))
        _, _, genes, idx, pos = ccg.synthetic_experiment(rng, sites, e, treatments=2)
        raw = np.column_stack([draw(x*lam) for x in e])
    else:
        C, T, genes, idx, pos = ccg.synthetic_experiment(rng, sites, e, config["m"], config["tau"], config["phi"], effects, treatments=2,
                                                         treatment_phi=config.get("treatment_phi"), treatment_linear=config.get("treatment_linear"))
        raw = np.column_stack([C, T])
    scaled, factors = normalize(sparse.csc_matrix(raw), "ttr", ["synthetic"]*4)
    return raw, scaled.toarray(), np.array([r["factor"] for r in factors]), genes, idx, pos, truth


def score(rows, truth):
    out = dict(unaffected=0, p05=0, false=0, calls=0, affected=0, hits=0)
    for row, t in zip(rows, truth):
        call = row.get("call")
        if not row.get("in_family"):
            out["affected"] += t != 0; continue
        called = call in {"depleted", "enriched", "depleted run"}
        direction = 1 if call == "enriched" else -1 if call in {"depleted", "depleted run"} else 0
        out["calls"] += called
        if t == 0:
            out["unaffected"] += 1; out["p05"] += row["p_family" if "p_family" in row else "p_value"] <= .05; out["false"] += called
        else:
            out["affected"] += 1; out["hits"] += called and direction == t; out["false"] += called and direction != t
    return out


def synthetic_task(task):
    name, experiment, phase = task
    raw, scaled, factors, genes, idx, pos, truth = make_data(name, experiment)
    results = {}
    if phase == "blind":
        for cid, spec in CANDIDATES.items():
            start = time.time(); rows, _ = resampling.analyze(scaled[:, :2], scaled[:, 2:], genes, idx, pos, resampling.Settings(**spec))
            results[cid] = score(rows, truth) | dict(seconds=time.time()-start)
    else:
        selected = json.loads((BENCH/"selection.json").read_text())["selected"]
        rows, _ = resampling.analyze(scaled[:, :2], scaled[:, 2:], genes, idx, pos, resampling.Settings(**CANDIDATES[selected]))
        results["resampling"] = score(rows, truth); r_calls = {r["ORF"] for r in rows if r["call"] in {"depleted", "enriched"}}
        inverse = 1/factors; exposures = inverse/inverse[:2].mean()
        g_rows, _, _ = ccg.analyze(raw[:, :2], raw[:, 2:], exposures, genes, idx, pos, ccg.Settings())
        results["ccg"] = score(g_rows, truth); g_calls = {r["ORF"] for r in g_rows if r["call"] == "depleted run"}
        hit = {genes[g]["id"] for g in np.flatnonzero(truth != 0)}
        results["union_hits"] = len((r_calls | g_calls) & hit); results["affected_total"] = int((truth != 0).sum())
    return name, experiment, results


def aggregate(per_experiment):
    keys = ("unaffected", "p05", "false", "calls", "affected", "hits")
    total = {k:sum(r[k] for r in per_experiment) for k in keys}
    total["experiments_with_false_call"] = sum(r["false"] > 0 for r in per_experiment)
    total["p05_rate"] = total["p05"]/max(total["unaffected"], 1); total["p05_ci"] = interval(total["p05"], total["unaffected"])
    total["sensitivity"] = total["hits"]/total["affected"] if total["affected"] else None
    total["fdr"] = total["false"]/total["calls"] if total["calls"] else None
    return total


def gale_arrays(normalization, scratch):
    design = ccg.load_design(ROOT/GALE["normalizations"][normalization], normalization, ROOT/GALE["reference"], scratch, GALE["controls"], GALE["treatments"])
    genes = ccg.annotations(ROOT/GALE["gff"]); idx = ccg.gene_sites(genes, design["sites"], resampling.Settings())
    positions = np.array([p for _, p in design["sites"]], dtype=np.int64)
    return design, genes, idx, positions


def gale_task(task):
    cid, normalization, seed, swap = task
    design, genes, idx, pos = gale_arrays(normalization, BENCH/"_unused_output")
    s = design["scaled"]
    a, b = (s[:, [0, 2]], s[:, [1, 3]]) if swap else (s[:, :2], s[:, 2:])
    start = time.time(); rows, diagnostics = resampling.analyze(a, b, genes, idx, pos, resampling.Settings(**CANDIDATES[cid], seed=seed))
    return cid, normalization, seed, swap, rows, diagnostics, time.time()-start


def jaccard(a, b):return len(a & b)/len(a | b) if a | b else 1.


def blind(workers):
    prereg = check_preregistration()
    if (BENCH/"selection.json").exists():raise SystemExit("selection.json exists; phase 1 already ran (delete deliberately to rerun)")
    tasks = [(n, k, "blind") for n in SCENARIOS for k in range(EXPERIMENTS)]
    with ProcessPoolExecutor(workers) as pool:synthetic = list(pool.map(synthetic_task, tasks))
    metrics = {cid:{n:aggregate([r[2][cid] for r in synthetic if r[0] == n]) for n in SCENARIOS} for cid in CANDIDATES}
    gale_tasks = [(cid, norm, seed, False) for cid in CANDIDATES for norm in GALE["normalizations"] for seed in (2026, 2027)]
    gale_tasks += [(cid, "ttr", 2026, True) for cid in CANDIDATES]
    with ProcessPoolExecutor(workers) as pool:gale = list(pool.map(gale_task, gale_tasks))
    (BENCH/"gale_blinded").mkdir(parents=True, exist_ok=True)
    calls, family, seconds = {}, {}, {}
    for cid, norm, seed, swap, rows, diagnostics, sec in gale:
        tag = f"{cid}_{norm}_seed{seed}" + ("_swapnull" if swap else "")
        calls[tag] = {r["ORF"] for r in rows if r["call"] in {"depleted", "enriched"}}; family[tag] = diagnostics["family_size"]; seconds[tag] = sec
        write_table(BENCH/"gale_blinded"/f"{tag}.tsv", ["ORF", "Name", "n_sites", "mean_control", "mean_treatment", "log2FC", "delta_mean", "p_value", "q_bh", "call", "status"],
                    [{k:r.get(k) for k in ("ORF", "Name", "n_sites", "mean_control", "mean_treatment", "log2FC", "delta_mean", "p_value", "q_bh", "call", "status")} for r in rows])
    summary = {}
    nulls = [n for n in SCENARIOS if n.startswith("null_")]; signals = [n for n in SCENARIOS if not n.startswith("null_")]
    for cid in CANDIDATES:
        m = metrics[cid]
        gate1a = all(m[n]["p05_ci"][0] <= .05 for n in nulls)
        null_experiments = [r[2][cid] for r in synthetic if r[0] in nulls]
        fwer = sum(r["false"] > 0 for r in null_experiments)/len(null_experiments)
        false = sum(m[n]["false"] for n in signals); calls_signal = sum(m[n]["calls"] for n in signals)
        pooled_fdr = false/calls_signal if calls_signal else 0.
        seed_j = {norm:jaccard(calls[f"{cid}_{norm}_seed2026"], calls[f"{cid}_{norm}_seed2027"]) for norm in GALE["normalizations"]}
        swap_fraction = len(calls[f"{cid}_ttr_seed2026_swapnull"])/family[f"{cid}_ttr_seed2026_swapnull"]
        gates = dict(g1a_null_p05_not_above_nominal=gate1a, g1b_null_experiments_with_false_call=fwer <= .10, g1c_pooled_signal_fdr=pooled_fdr <= .10,
                     g2_seed_jaccard=min(seed_j.values()) >= .99, g3_swap_null=swap_fraction <= .005)
        summary[cid] = dict(gates=gates, passed=all(gates.values()), null_experiment_false_fraction=fwer, pooled_signal_fdr=pooled_fdr,
            seed_jaccard=seed_j, swap_null_calls=len(calls[f"{cid}_ttr_seed2026_swapnull"]), swap_null_fraction=swap_fraction,
            mean_signal_sensitivity=float(np.mean([m[n]["sensitivity"] for n in signals])),
            reported_only=dict(gale_calls={norm:len(calls[f"{cid}_{norm}_seed2026"]) for norm in GALE["normalizations"]},
                normalization_jaccard={f"{x}-{y}":jaccard(calls[f"{cid}_{x}_seed2026"], calls[f"{cid}_{y}_seed2026"])
                                       for i, x in enumerate(GALE["normalizations"]) for y in list(GALE["normalizations"])[i+1:]},
                gale_seconds={k.split(cid+"_")[1]:v for k, v in seconds.items() if k.startswith(cid+"_")}))
    passing = [c for c in CANDIDATES if summary[c]["passed"]]
    selected = None
    if passing:
        best = max(summary[c]["mean_signal_sensitivity"] for c in passing)
        near = [c for c in passing if best-summary[c]["mean_signal_sensitivity"] <= .01]
        selected = sorted(near, key=lambda c:(CANDIDATES[c]["inference"] != "exact", CANDIDATES[c]["opportunities"] != "any-support"))[0]
    rows = []
    for cid in CANDIDATES:
        for n in SCENARIOS:
            m = metrics[cid][n]
            rows.append(dict(candidate=cid, scenario=n, unaffected_p_le_05=m["p05_rate"], ci_low=m["p05_ci"][0], ci_high=m["p05_ci"][1],
                             false_calls=m["false"], calls=m["calls"], experiments_with_false_call=m["experiments_with_false_call"],
                             sensitivity=m["sensitivity"], fdr=m["fdr"]))
    write_table(BENCH/"synthetic_candidates.tsv", list(rows[0]), rows)
    (BENCH/"selection.json").write_text(json.dumps(ccg._json_ready(dict(status="success", phase="blind", preregistration=prereg,
        selected=selected, selected_settings=CANDIDATES.get(selected), passing=passing, candidates=summary,
        fixed_settings=asdict(resampling.Settings()), blinding="no publication labels or control-conditioned Gumbel outputs were read",
        implementation_sha256=sha256(Path(resampling.__file__)), script_sha256=sha256(Path(__file__)))), indent=2)+"\n")
    print(json.dumps(ccg._json_ready(dict(selected=selected, passing=passing, **{c:dict(passed=summary[c]["passed"], gates=summary[c]["gates"],
          sensitivity=round(summary[c]["mean_signal_sensitivity"], 3), swap=summary[c]["swap_null_calls"], seed_j=summary[c]["seed_jaccard"]) for c in CANDIDATES})), indent=1))


def unblind(workers):
    prereg = check_preregistration()
    selection = json.loads((BENCH/"selection.json").read_text())
    if selection["preregistration"]["sha256"] != prereg["sha256"]:raise SystemExit("selection was made under a different pre-registration")
    if not selection["selected"]:raise SystemExit("phase 1 selected no candidate; nothing to unblind")
    tasks = [(n, k, "unblind") for n in SCENARIOS for k in range(EXPERIMENTS)]
    with ProcessPoolExecutor(workers) as pool:synthetic = list(pool.map(synthetic_task, tasks))
    rows = []
    for n in SCENARIOS:
        part = [r[2] for r in synthetic if r[0] == n]
        r_m, g_m = aggregate([p["resampling"] for p in part]), aggregate([p["ccg"] for p in part])
        affected = sum(p["affected_total"] for p in part)
        rows.append(dict(scenario=n, resampling_p05=r_m["p05_rate"], ccg_p05=g_m["p05_rate"], resampling_false=r_m["false"], ccg_false=g_m["false"],
                         resampling_calls=r_m["calls"], ccg_calls=g_m["calls"], resampling_sensitivity=r_m["sensitivity"], ccg_sensitivity=g_m["sensitivity"],
                         union_sensitivity=sum(p["union_hits"] for p in part)/affected if affected else None))
    write_table(BENCH/"unblinded_synthetic_vs_ccg.tsv", list(rows[0]), rows)
    import pandas as pd
    sel = selection["selected"]
    res = pd.read_csv(BENCH/"gale_blinded"/f"{sel}_ttr_seed2026.tsv", sep="\t").set_index("ORF")
    gumbel = pd.read_csv(ROOT/"output/essentiality/control_conditioned_gumbel/Gale_2020/runs/exact_ttr_seed2026/genes.tsv", sep="\t").set_index("ORF")
    pub = pd.read_csv(ROOT/"output/essentiality/control_conditioned_gumbel/Gale_2020/publication_comparison.tsv", sep="\t").set_index("ORF")
    # Gale swap null for the Gumbel (same balanced design) for a like-for-like real-data false-positive check.
    design, genes, idx, pos = gale_arrays("ttr", BENCH/"_unused_output")
    raw, inverse = design["counts"], 1/design["factors"]
    sw = [0, 2, 1, 3]; e = inverse[sw]/inverse[sw][:2].mean()
    g_swap, _, g_swap_diag = ccg.analyze(raw[:, [0, 2]], raw[:, [1, 3]], e, genes, idx, pos, ccg.Settings())
    table = res[["Name", "n_sites", "log2FC", "p_value", "q_bh", "call"]].join(gumbel[["call", "q_bh", "r"]].rename(columns={"call":"ccg_call", "q_bh":"ccg_q"}), how="outer")
    table = table.join(pub[["gene_name", "publication_class", "published_z", "published_log2fc"]], how="left")
    table.to_csv(BENCH/"unblinded_gale_comparison.tsv", sep="\t")
    r_dep = set(res.index[res.call == "depleted"]); r_enr = set(res.index[res.call == "enriched"]); g = set(gumbel.index[gumbel.call == "depleted run"])
    strict = set(pub.index[pub.publication_class.astype(str).str.startswith("A")]); z2 = set(pub.index[pub.published_z <= -2]); z3 = set(pub.index[pub.published_z <= -3])
    zp2 = set(pub.index[pub.published_z >= 2])
    report = dict(selected=sel, resampling_depleted=len(r_dep), resampling_enriched=len(r_enr), ccg_depleted_run=len(g),
        overlap_depleted=len(r_dep & g), resampling_only=len(r_dep-g), ccg_only=len(g-r_dep), jaccard=jaccard(r_dep, g),
        strict_published_called=dict(resampling=len(r_dep & strict), ccg=len(g & strict), union=len((r_dep | g) & strict), testable=len(strict)),
        depleted_with_gale_z_le_m2=dict(resampling=len(r_dep & z2), ccg=len(g & z2)), depleted_with_gale_z_le_m3=dict(resampling=len(r_dep & z3), ccg=len(g & z3)),
        enriched_with_gale_z_ge_2=len(r_enr & zp2),
        spearman_log2fc_vs_published=float(table.log2FC.corr(table.published_log2fc, method="spearman")),
        spearman_neglogp_vs_published_z=float((-np.log10(table.p_value.clip(1e-300))).corr(-table.published_z, method="spearman")),
        swap_null=dict(resampling=selection["candidates"][sel]["swap_null_calls"], ccg=g_swap_diag["calls"]["depleted run"]))
    (BENCH/"unblinded_report.json").write_text(json.dumps(ccg._json_ready(report), indent=2)+"\n")
    print(json.dumps(ccg._json_ready(report), indent=1)); print(pd.DataFrame(rows).round(3).to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["blind", "unblind"])
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    blind(args.workers) if args.phase == "blind" else unblind(args.workers)


if __name__ == "__main__":main()
