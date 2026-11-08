#!/usr/bin/env python3
"""Library saturation and essentiality in four species: HMM, species-trained random forest, hybrid.

Implements output/essentiality/saturation/PREREGISTRATION.md (frozen; both phases refuse to run if it changed).
  thin-fit:  thin each library's insertion sites, build SummaryTable gene features, fit the HMM (reads no labels)
  evaluate:  labels, out-of-fold random forest and hybrid, AUCs, the saturation rule, hybrid vs RF tests
"""
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src")); sys.path.append(str(ROOT/"scripts/local"))  # after src: scripts/local/ytab.py must not shadow the package
from ytab.essentiality import model_menu as mm
from ytab.normalization.midlc import midlc
from ytab.normalization.ttr import sha256
from ytab.pipeline import summary_table_runner as st

BENCH = ROOT/"output/essentiality/saturation"
QS = [1, 2, 4, 8, 16, 32]  # keep fraction 1/q
SEEDS = [1, 2]
BOOT_SEED, BOOTSTRAPS = 52734, 2000
FEW_LABELS, DRAWS = 50, 5
SPECIES = {
    "C. glabrata": dict(project="Gale_2020", libraries=["SRR11561508", "SRR11561507", "SRR11561506"]),
    "C. albicans": dict(project="Segal_2018_albicans", libraries=["SRR7824843", "SRR7824841", "SRR7824838"]),
    "S. cerevisiae": dict(project="Michel_2017_SATAY", libraries=["ERR1533147", "ERR1533148"]),
    "S. pombe": dict(project="Pombe_Hermes", libraries=[f"ERR2641{n}" for n in range(199, 213)]),
}
SLUG = {"C. glabrata":"glabrata", "C. albicans":"albicans", "S. cerevisiae":"cerevisiae", "S. pombe":"pombe"}


def check_preregistration():
    frozen = (BENCH/"PREREGISTRATION.sha256").read_text().split()
    if sha256(BENCH/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def tasks():
    return [(sp, q, seed) for sp in SPECIES for q in QS for seed in ([0] if q == 1 else SEEDS)]


def run_dir(sp, q, seed):
    return BENCH/"runs"/SLUG[sp]/f"q{q}_s{seed}"


def sh(cmd, log):
    with open(log, "w") as handle:
        r = subprocess.run([str(c) for c in cmd], stdout=handle, stderr=subprocess.STDOUT, cwd=ROOT)
    if r.returncode:raise RuntimeError(f"failed ({log})")


def keep(lib, chrom, pos, q, seed):
    return q == 1 or int(hashlib.md5(f"{lib}:{chrom}:{pos}:{seed}".encode()).hexdigest()[:8], 16)/0xFFFFFFFF < 1/q


def thin_fit(task):
    sp, q, seed = task; cfg = SPECIES[sp]; out = run_dir(sp, q, seed)
    if (out/"run.json").exists():return json.loads((out/"run.json").read_text())
    (out/"hits").mkdir(parents=True, exist_ok=True)
    project_cfg = st.load_project_for_summary_table(ROOT/"output/projects"/cfg["project"]/"config/project.yaml")
    reference = st._reference_files(project_cfg); feature_args = st.feature_db_args(project_cfg, reference)
    fasta, gff = feature_args[feature_args.index("--reference-fasta")+1], feature_args[feature_args.index("--feature-gff")+1]
    files, combined, saturation = [], {}, []
    for lib in cfg["libraries"]:
        src = ROOT/"output/projects"/cfg["project"]/"create_hit_file"/lib/f"{lib}_hits.txt"; dst = out/"hits"/f"{lib}_hits.txt"
        counts = []
        with open(src) as f, open(dst, "w") as o:
            r = csv.reader(f, delimiter="\t"); w = csv.writer(o, delimiter="\t", lineterminator="\n"); h = next(r); w.writerow(h)
            pos, cnt = h.index("Hit position"), h.index("Hit count")
            for row in r:
                if keep(lib, row[0], row[pos], q, seed):
                    w.writerow(row); counts.append(int(row[cnt]))
                    key = (row[0], row[1], row[pos]); base = combined.setdefault(key, row[:])
                    if base is not row:base[cnt] = str(int(base[cnt])+int(row[cnt]))
        saturation.append(float(sum(counts)/midlc(np.array(counts))[0]))
        files.append(dst)
    (out/"combined").mkdir(exist_ok=True)
    with open(out/"combined/combined_hits.txt", "w") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(h); w.writerows(combined.values())
    sh([sys.executable, ROOT/"src/ytab/summary/SummaryTable.py", "--input-dir", out/"combined", "--hit-glob", "combined_hits.txt",
        "--output-dir", out/"summary", "--feature-db", st.feature_db_kind(project_cfg, reference), *feature_args[2:], "--threads", "2", "--overwrite"], out/"summary.log")
    table = out/"summary/combined.feature_table.RDF_1.csv"
    f = pd.read_csv(table, skiprows=1); f = f[f.coding_length > 0]
    sh([sys.executable, ROOT/"scripts/local/ytab_hermes_hmm.py", out/"hmm", "--fasta", fasta, "--gff", gff, *[x for p in files for x in ("--hits", p)],
        "--sequence-covariate", "learned"], out/"hmm.log")
    fit = json.loads((out/"hmm/joint/fit.json").read_text())
    genome_kb = sum(len(l.strip()) for l in open(fasta) if not l.startswith(">"))/1000
    record = dict(species=sp, q=q, seed=seed, keep_fraction=1/q, unique_sites=len({(k[0], k[2]) for k in combined}),
                  sites_per_coding_kb=float((f.hits/f.coding_length*1000).median()), sites_per_genome_kb=len({(k[0], k[2]) for k in combined})/genome_kb,
                  saturation_ratio_median=float(np.median(saturation)), hmm_converged=bool(fit["converged"]), hmm_iterations=fit["iterations"],
                  masked_fraction=fit["masked_fraction"], passed=bool(fit["converged"]) and fit["masked_fraction"] < .5, feature_table=str(table.relative_to(ROOT)))
    (out/"run.json").write_text(json.dumps(record, indent=1)); return record


def thin_fit_all(workers):
    prereg = check_preregistration()
    with ProcessPoolExecutor(max_workers=workers) as pool:records = list(pool.map(thin_fit, tasks()))
    (BENCH/"runs.json").write_text(json.dumps(dict(preregistration=prereg, script_sha256=sha256(Path(__file__)), runs=records), indent=1))
    for r in records:print(r["species"], r["q"], r["seed"], round(r["sites_per_coding_kb"], 2), r["unique_sites"], "passed" if r["passed"] else "FAILED")


# ---------------------------------------------------------------- evaluation (reads labels)

def labels(sp):
    import ytab_four_species_benchmark as fs
    if sp == "C. glabrata":
        lab = pd.read_csv(ROOT/"resources/species/glabrata/essentiality_labels/labels.tsv", sep="\t", dtype=str)
        lab = lab[(lab.tier == "silver") & lab.label.isin(["essential", "nonessential"])]
        return {g:int(l == "essential") for g, l in zip(lab.cglab_id, lab.label)}
    slug = SLUG[sp]; lab = fs.labels(slug); mask = fs.masked(slug)
    return {g:v for g, v in lab.items() if g not in mask}


def auc(score, label):
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(label, score))


def genes_table(record):
    out = BENCH/"runs"/SLUG[record["species"]]/f"q{record['q']}_s{record['seed']}"
    hmm = pd.read_csv(out/"hmm/joint/genes.tsv", sep="\t")
    feats = mm.summary_features(pd.read_csv(ROOT/record["feature_table"], skiprows=1))
    d = hmm[["ORF", "chromosome", "status"]+mm.HMM_FEATURES].merge(feats, on="ORF", how="inner")
    d["ORF"] = d.ORF.str.removeprefix("gene:")
    return d[(d.status == "called") & d[mm.RF_FEATURES+mm.HMM_FEATURES].notna().all(axis=1)].reset_index(drop=True)


def out_of_fold(d, columns):
    X, y = d[columns].to_numpy(float), d.essential.to_numpy(bool)
    folds = mm.chromosome_folds(d.chromosome.tolist()); oof = np.zeros(len(d))
    for k in np.unique(folds):
        tr, te = folds != k, folds == k
        oof[te] = mm.learner(mm.SEED).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    return oof


def evaluate(include_capped=False):
    """include_capped: exploratory sensitivity arm (not pre-registered) that also scores HMM fits stopped at the iteration
    cap with a monotone likelihood; written to sensitivity_capped/."""
    prereg = check_preregistration(); runs = json.loads((BENCH/"runs.json").read_text())["runs"]
    out = BENCH/"sensitivity_capped" if include_capped else BENCH; out.mkdir(exist_ok=True)
    rows, few, tests = [], [], []
    rng = np.random.default_rng(BOOT_SEED)
    for record in runs:
        if not (record["passed"] or include_capped and record["masked_fraction"] < .5):continue
        sp = record["species"]; lab = labels(sp)
        d = genes_table(record); d = d[d.ORF.isin(lab)].reset_index(drop=True); d["essential"] = d.ORF.map(lab).astype(bool)
        y = d.essential.to_numpy(int)
        rf, hy = out_of_fold(d, mm.RF_FEATURES), out_of_fold(d, mm.RF_FEATURES+mm.HMM_FEATURES)
        base = {k:record[k] for k in ("species", "q", "seed", "keep_fraction", "sites_per_coding_kb", "sites_per_genome_kb", "unique_sites", "saturation_ratio_median")}
        base |= dict(genes=len(d), essential=int(y.sum()), hmm_converged=record["hmm_converged"])
        rows += [dict(base, model="HMM (label-free)", auc=auc(d.P_ES_run, y)), dict(base, model="Species RF", auc=auc(rf, y)),
                 dict(base, model="Hybrid", auc=auc(hy, y))]
        if (record["q"], record["seed"]) in ((1, 0), (16, 1)):
            deltas = []
            for _ in range(BOOTSTRAPS):
                i = rng.integers(0, len(y), len(y))
                if 0 < y[i].sum() < len(y):deltas.append(auc(hy[i], y[i])-auc(rf[i], y[i]))
            deltas = np.array(deltas)
            tests.append(dict(species=sp, q=record["q"], delta=auc(hy, y)-auc(rf, y), ci_low=float(np.percentile(deltas, 2.5)),
                              ci_high=float(np.percentile(deltas, 97.5)), p_one_sided=float((deltas <= 0).mean())))
        for draw in range(DRAWS):
            r = np.random.default_rng(draw); ess, non = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
            k = max(2, round(FEW_LABELS*len(ess)/len(y)))
            train = np.r_[r.choice(ess, k, replace=False), r.choice(non, FEW_LABELS-k, replace=False)]
            test = np.setdiff1d(np.arange(len(y)), train)
            for model, cols in (("Species RF", mm.RF_FEATURES), ("Hybrid", mm.RF_FEATURES+mm.HMM_FEATURES)):
                clf = mm.learner(mm.SEED).fit(d.loc[train, cols].to_numpy(float), y[train].astype(bool))
                few.append(dict(base, model=model, draw=draw, labels=FEW_LABELS, auc=auc(clf.predict_proba(d.loc[test, cols].to_numpy(float))[:, 1], y[test])))
            few.append(dict(base, model="HMM (label-free)", draw=draw, labels=0, auc=auc(d.P_ES_run.to_numpy()[test], y[test])))
        print(sp, record["q"], record["seed"], {r["model"]:round(r["auc"], 3) for r in rows[-3:]}, flush=True)
    # Holm over the 8 pre-registered tests (m stays 8 when a test could not be run)
    order = np.argsort([t["p_one_sided"] for t in tests]); m = 2*len(SLUG); running = 0
    for rank, i in enumerate(order):
        running = max(running, min(1, (m-rank)*tests[i]["p_one_sided"])); tests[i]["p_holm"] = running; tests[i]["hybrid_better"] = running < .05
    curve = pd.DataFrame(rows); curve.to_csv(out/"auc_by_density.tsv", sep="\t", index=False)
    pd.DataFrame(few).to_csv(out/"few_labels.tsv", sep="\t", index=False); pd.DataFrame(tests).to_csv(out/"hybrid_vs_rf.tsv", sep="\t", index=False)
    mean = curve.groupby(["species", "model", "q"], as_index=False)[["auc", "sites_per_coding_kb", "sites_per_genome_kb", "saturation_ratio_median"]].mean()
    adequate = {}
    for (sp, model), g in mean.groupby(["species", "model"]):
        full = g[g.q == 1].auc.iloc[0]; ok = g[g.auc >= full-.02]
        adequate[f"{sp} | {model}"] = dict(full_auc=float(full), adequate_q=int(ok.q.max()), adequate_sites_per_coding_kb=float(ok.sites_per_coding_kb.min()),
                                           lowest_q_evaluated=int(g.q.max()))
    rec = {model:max(v["adequate_sites_per_coding_kb"] for k, v in adequate.items() if k.endswith(model)) for model in ("HMM (label-free)", "Hybrid")}
    result = dict(preregistration=prereg, script_sha256=sha256(Path(__file__)), include_capped=include_capped, adequate=adequate, recommendation_sites_per_coding_kb=rec, tests=tests)
    (out/"result.json").write_text(json.dumps(result, indent=1)); print(json.dumps({k:v for k, v in result.items() if k != "preregistration"}, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("phase", choices=["thin-fit", "evaluate"]); p.add_argument("--workers", type=int, default=6)
    p.add_argument("--include-capped", action="store_true", help="exploratory: also score fits stopped at the iteration cap")
    a = p.parse_args(); thin_fit_all(a.workers) if a.phase == "thin-fit" else evaluate(a.include_capped)


if __name__ == "__main__":main()
