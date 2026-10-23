#!/usr/bin/env python3
"""Pre-registered evaluation of Hermes genetic interactions (hermes_gi_v1).

dev:    development seeds (5150, 5151), 10 experiments per scenario; not acceptance (writes development/).
blind:  fresh seed, gates in output/essentiality/gi/benchmark/PREREGISTRATION.md; writes decision.json.
zn:     label-free real data, Zn screen (BGC): yH298 vs yH299 GI, and a balanced pool-swap null; reported only.

Simulator: strains A (reference) and B (query), each its own library/libraries with its own insertion sites;
lognormal clone abundance (m = log 30, tau = 1.5), negative binomial (phi 0.2); 300 genes; affected genes
change over a block of half their sites. Designs:
  split  one pool per strain, 2 control + 2 treated cultures (paired in order)
  pools  two pools per strain, each split into 1 control + 1 treated culture
BGC normalization per library. Scenario effects act on 10% of genes (a second 10% where stated).
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
from ytab.essentiality import genetic_interactions as gi
from ytab.normalization.methods import normalize
from ytab.normalization.ttr import sha256, write_table

OUT = ROOT/"output/essentiality/gi"; BENCH = OUT/"benchmark"
GENES = 300; SEEDS = dict(dev=(5150, 5151), blind=93939)
BASE = dict(m=float(np.log(30)), tau=1.5, phi=.2)
# per affected gene class: A treated fold, B treated fold, B baseline fold (both conditions); "second" = a second 10% of genes
SCENARIOS = {
    "null_no_effects":dict(),
    "null_shared_response":dict(first=(.25, .25, 1.)),
    "null_baseline_difference":dict(first=(1., 1., .25)),
    "null_shared_response_unequal_depth":dict(first=(.25, .25, 1.), depth_B=.3),
    "null_shared_response_noisy_B":dict(first=(.25, .25, 1.), phi=.02, phi_B=.5),
    "gi_aggravating":dict(first=(1., .25, 1.)),
    "gi_alleviating":dict(first=(.25, 1., 1.)),
    "gi_suppressive":dict(first=(1., 4., 1.)),
    "gi_mixed":dict(first=(.25, .25, 1.), second=(1., .25, 1.)),
    "pools_no_shared_clones":dict(first=(1., .25, 1.), fresh_clones=True),     # reported only (pools design)
}
TRUE_CLASS = {"gi_aggravating":"Aggravating", "gi_alleviating":"Alleviating", "gi_suppressive":"Suppressive", "gi_mixed":"Aggravating",
              "pools_no_shared_clones":"Aggravating"}
REPORTED_ONLY = {"pools_no_shared_clones"}
ZN = dict(normalized=ROOT/"output/normalization/zn_bgc_sites_v1", project="Zn_toxicity_screen",
          gff=ROOT/"resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff",
          fasta=ROOT/"resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.fna")


def simulate(name, design, k, seed):
    cfg = BASE | SCENARIOS[name]
    rng = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(sorted(SCENARIOS).index(name), ["split", "pools"].index(design), k)))
    cls = rng.permutation(np.repeat([0, 1, 2], [GENES-2*(GENES//10), GENES//10, GENES//10]))
    blocks, strains, conds, libs, idx, offset = [], [], [], [], [[] for _ in range(GENES)], 0
    for strain in "AB":
        fold_t = np.ones(GENES); fold_b = np.ones(GENES)
        for key, c in (("first", 1), ("second", 2)):
            if key in cfg:
                a_t, b_t, b_base = cfg[key]; fold_t[cls == c] = a_t if strain == "A" else b_t
                if strain == "B":fold_b[cls == c] = b_base
        phi = cfg.get("phi_B", cfg["phi"]) if strain == "B" else cfg["phi"]; size = 1/phi
        depth = cfg.get("depth_B", 1.) if strain == "B" else 1.
        for pool in ((1,) if design == "split" else (1, 2)):
            sites = np.clip(np.round(np.exp(rng.uniform(np.log(3), np.log(120), GENES))), 3, 120).astype(int)
            g = np.repeat(np.arange(GENES), sites); offs = np.concatenate([np.arange(n) for n in sites])
            block = (offs >= sites[g]//4) & (offs < sites[g]//4+np.maximum(2, sites[g]//2))
            lam = np.exp(cfg["m"]+cfg["tau"]*rng.standard_normal(len(g)))
            ft = np.where(block, fold_t[g], 1.); fb = np.where(block, fold_b[g], 1.)
            draw = lambda mu:rng.negative_binomial(size, size/(size+mu))
            treated_lam = (np.exp(cfg["m"]+cfg["tau"]*rng.standard_normal(len(g))) if cfg.get("fresh_clones") else lam)
            if design == "split":cols = [draw(depth*lam*fb), draw(depth*1.3*lam*fb), draw(depth*.8*treated_lam*fb*ft), draw(depth*1.1*treated_lam*fb*ft)]; cc = ["control", "control", "treated", "treated"]
            else:cols = [draw(depth*lam*fb), draw(depth*(.9 if pool == 1 else 1.2)*treated_lam*fb*ft)]; cc = ["control", "treated"]
            blocks.append(np.column_stack(cols)); strains += [strain]*len(cols); conds += cc; libs += [f"{strain}{pool}"]*len(cols)
            start = 0
            for k_, n in enumerate(sites):idx[k_].append(offset+np.arange(start, start+n)); start += n
            offset += len(g)
    raw = np.zeros((offset, len(strains)), dtype=np.int64); r = c = 0
    for b in blocks:raw[r:r+len(b), c:c+b.shape[1]] = b; r += len(b); c += b.shape[1]
    gi_genes = np.zeros(GENES, bool)
    for key, c_ in (("first", 1), ("second", 2)):
        if key in cfg and cfg[key][0] != cfg[key][1]:gi_genes |= cls == c_
    return raw, strains, conds, libs, [np.concatenate(x) for x in idx], gi_genes


def task(item):
    name, design, k, seed = item
    raw, strains, conds, libs, idx, truth = simulate(name, design, k, seed)
    values = normalize(sparse.csc_matrix(raw), "betageom", libs)[0].toarray()
    genes = [dict(id=f"g{j:04d}", name=f"g{j:04d}") for j in range(GENES)]
    rows, diag = gi.analyze(values, strains, conds, libs, genes, idx, raw=raw)
    out = dict(null=0, p05=0, calls=0, false=0, gi=0, hits=0, class_right=0, experiments=1,
               clone_warnings=int(diag["warning"] is not None))
    for t, r in zip(truth, rows):
        if r["status"] != "tested":
            out["gi"] += int(t); continue
        called = r["call"] in gi.CLASSES; out["calls"] += called
        if t:
            out["gi"] += 1; out["hits"] += called; out["class_right"] += called and r["call"] == TRUE_CLASS.get(name)
        else:out["null"] += 1; out["p05"] += r["p_value"] <= .05; out["false"] += called
    return name, design, k, out, [r.get("p_value") for r in rows]


def interval(k, n):
    return [float(beta.ppf(.025, k, n-k+1)) if k else 0., float(beta.ppf(.975, k+1, n-k)) if k < n else 1.]


def items(seed, experiments):
    return [(s, d, k, seed) for s in SCENARIOS for d in ("split", "pools") if not (s in REPORTED_ONLY and d == "split") for k in range(experiments)]


def aggregate(results):
    agg = {}
    for name, design, _, out, _ in results:
        a = agg.setdefault((name, design), {})
        for key, v in out.items():a[key] = a.get(key, 0)+v
    return [dict(scenario=n, design=d, null_genes=a["null"], type1_p05=a["p05"]/a["null"], type1_ci=interval(a["p05"], a["null"]), calls=a["calls"],
                 fdr=a["false"]/a["calls"] if a["calls"] else 0., sensitivity=a["hits"]/a["gi"] if a["gi"] else None,
                 class_accuracy=a["class_right"]/a["hits"] if a["hits"] else None, clone_warnings=f"{a['clone_warnings']}/{a['experiments']}")
            for (n, d), a in sorted(agg.items())]


def run_all(seed_list, experiments, workers):
    with ProcessPoolExecutor(workers) as pool:return [r for seed in seed_list for r in pool.map(task, items(seed, experiments))]


def dev(workers):
    target = OUT/"development"; target.mkdir(parents=True, exist_ok=True)
    rows = aggregate(run_all(SEEDS["dev"], 10, workers))
    write_table(target/"dev_final_seeds5150_5151_10exp.tsv", list(rows[0]), [r | dict(type1_ci=str(r["type1_ci"])) for r in rows])
    for r in rows:print({k:(round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k != "type1_ci"})


def check_preregistration():
    frozen = (BENCH/"PREREGISTRATION.sha256").read_text().split()
    if sha256(BENCH/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def blind(workers):
    prereg = check_preregistration()
    if (BENCH/"decision.json").exists():raise SystemExit("decision.json exists; the blind evaluation already ran")
    results = run_all((SEEDS["blind"],), 8, workers); rows = aggregate(results)
    write_table(BENCH/"synthetic.tsv", list(rows[0]), [r | dict(type1_ci=";".join(f"{x:.4f}" for x in r["type1_ci"])) for r in rows])
    gates = {}
    for r in rows:
        if r["scenario"] in REPORTED_ONLY:continue
        gates[f"type1 {r['scenario']} {r['design']}"] = dict(value=r["type1_p05"], passed=.025 <= r["type1_p05"] <= .075)
        if r["scenario"].startswith("gi_"):gates[f"fdr {r['scenario']} {r['design']}"] = dict(value=r["fdr"], passed=r["fdr"] <= .10)
    again = task(("gi_mixed", "split", 0, SEEDS["blind"])); first = next(x for x in results if x[0] == "gi_mixed" and x[1] == "split" and x[2] == 0)
    gates["determinism"] = dict(passed=again[4] == first[4])
    decision = dict(preregistration=prereg, passed=all(g["passed"] for g in gates.values()), gates=gates,
                    script_sha256=sha256(Path(__file__)), module_sha256=sha256(Path(gi.__file__)))
    (BENCH/"decision.json").write_text(json.dumps(decision, indent=2)+"\n"); print(json.dumps(decision, indent=1))


def zn(workers):
    target = OUT/"Zn_toxicity_screen"; target.mkdir(parents=True, exist_ok=True)
    pools = {"yH298":("yH298-parent-pool1", "yH298-parent-pool2"), "yH299":("yH299-parent-pool3", "yH299-parent-pool4")}
    obs = lambda p, kind:f"{ZN['project']}::{p}-{'mock' if kind == 'control' else '1_5mM-Zn-treated'}"
    real = [(obs(p, k), s, k) for s, ps in pools.items() for p in ps for k in ("control", "treated")]
    swap_groups = {"G1":("yH298-parent-pool1", "yH299-parent-pool3"), "G2":("yH298-parent-pool2", "yH299-parent-pool4")}
    swap = [(obs(p, k), g, k) for g, ps in swap_groups.items() for p in ps for k in ("control", "treated")]
    summary = {}
    for label, samples, ref, query in (("yH298_vs_yH299", real, "yH298", "yH299"), ("pool_swap_null", swap, "G1", "G2")):
        m = gi.run(ZN["normalized"], ZN["gff"], ZN["fasta"], target/label, samples, ref, query, normalization="betageom")
        d = m["diagnostics"]
        rows = list(__import__("csv").DictReader((target/label/"genes.tsv").open(), delimiter="\t"))
        p = np.array([float(r["p_value"]) for r in rows if r["status"] == "tested"])
        summary[label] = dict(tested=d["tested"], calls=d["calls"], p05_fraction=float(np.mean(p <= .05)), site_variance_scale=d["site_variance_scale"],
                              clone_sharing=d["clone_sharing"], warning=d["warning"])
    (target/"summary.json").write_text(json.dumps(ccg._json_ready(summary), indent=2)+"\n"); print(json.dumps(ccg._json_ready(summary), indent=1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("phase", choices=["dev", "blind", "zn"]); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(); dict(dev=dev, blind=blind, zn=zn)[a.phase](a.workers)


if __name__ == "__main__":main()
