#!/usr/bin/env python3
"""Method-comparison report for the Cunningham lab, step 1: every method x top-site treatment x normalization on the
two datasets with published z-score calls. Untracked study (output/comparison/ is ignored).

Datasets: Gale 2020 (BG14 fluconazole) and Nickels 2025 (BG2 and pdr1_KO, echinocandins).
Top-site treatment (normalize-project --site-filter): none, winsorize, drop-top; each a full normalization run.
Under each: z-score with every normalization (unnormalized included), with the two-control noise model and with the
noise fitted to each comparison (zscorecmp), resampling with the normalizations it accepts
(ttr, totreads, nzmean, betageom, aBGC), control-conditioned Gumbel with ttr, totreads, nzmean (conditions with fewer
than two treated libraries are declined by the method, not forced).

Output: output/comparison/report/<dataset>/<filter>/{norm,<method>_<normalization>}/ and status.tsv.
  python scripts/local/ytab_method_comparison_report_runs.py [--jobs 3]
"""
from __future__ import annotations
import argparse, csv, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[2]; PY = sys.executable; OUT = ROOT/"output/comparison/report"
DATASETS = {"Gale_2020": None, "Nickels_2025_glabrata": {"BG2", "pdr1_KO"}}
FILTERS = ["none", "winsorize", "drop-top"]
ZSCORE = ["nonorm", "totreads", "nzmean", "ttr", "betageom", "aBGC", "quantile", "zinfnb", "midlc"]
RESAMPLING = ["ttr", "totreads", "nzmean", "betageom", "aBGC"]
GUMBEL = ["ttr", "totreads", "nzmean"]


def comparisons(dataset):
    rows = [r for r in csv.DictReader(open(ROOT/"output/projects"/dataset/"config/comparison_design.csv"))
            if str(r.get("include", "")).strip().lower() in ("1", "true", "yes", "y")]
    keep = DATASETS[dataset]
    return [r["comparison_id"] for r in rows if keep is None or r["background"] in keep]


def log(line):
    with open(OUT/"status.tsv", "a") as h:h.write(line+"\n")
    print(line, flush=True)


def call(cmd, logfile):
    t = time.time(); logfile.parent.mkdir(parents=True, exist_ok=True)
    with open(logfile, "w") as h:code = subprocess.call([str(c) for c in cmd], stdout=h, stderr=subprocess.STDOUT, cwd=ROOT)
    return code, round(time.time()-t)


def normalization(dataset, how):
    base = OUT/dataset/how; norm = base/"norm"
    if not (norm/"manifest.json").exists():
        if norm.exists():subprocess.call(["rm", "-rf", str(norm)])
        code, secs = call([PY, "scripts/local/ytab_normalize_project.py", "--project-config", ROOT/"output/projects"/dataset/"config/project.yaml",
                           "--site-filter", how, "--output", norm], base/"norm.log")
        log(f"{dataset}\t{how}\tnormalization\tall\t{code}\t{secs}")
    return norm


def jobs(dataset, how, norm):
    cfg = ROOT/"output/projects"/dataset/"config/project.yaml"; gff = yaml.safe_load(open(cfg))["reference"]["gff"]
    gff = gff if Path(gff).is_absolute() else str(ROOT/gff); comps = ",".join(comparisons(dataset)); base = OUT/dataset/how; out = []
    for method, noise, ctrls in (("zscore", "control", "first-two"), ("zscorecmp", "comparison", "first-two"), ("zscoreall", "control", "all")):
        # zscorecmp: noise from each comparison itself; zscoreall: the control curve fitted on every pair of the background's controls
        for n in ZSCORE:
            out.append((dataset, how, method, n, [PY, "scripts/local/ytab_fitness_zscore.py", "--project-id", dataset, "--analysis-id", f"report_{how}_{n}",
                        "--comparison-design", ROOT/"output/projects"/dataset/"config/comparison_design.csv", "--comparisons", comps, "--noise", noise, "--noise-controls", ctrls,
                        "--output-dir", base/f"{method}_{n}", "--normalized", norm, "--normalization", n, "--gff", gff], base/f"{method}_{n}"/"zscore_manifest.json"))
    for m, norms in (("resampling", RESAMPLING), ("gumbel", GUMBEL)):
        for n in norms:
            out.append((dataset, how, m, n, [PY, "scripts/local/ytab_treated_vs_control.py", "--project-config", cfg, "--method", m, "--normalization-run", norm,
                        "--normalization", n, "--comparisons", comps, "--output", base/f"{m}_{n}"], base/f"{m}_{n}"/"manifest.json"))
    return out


def run(job):
    dataset, how, method, n, cmd, done = job
    if done.exists():return
    code, secs = call(cmd, OUT/dataset/how/f"{method}_{n}.log")
    if method.startswith("zscore"):   # per-pair tables are large and not used by the report: keep them, compressed
        for name in ("top_hit_evidence.by_pair.csv", "treated_vs_parent.by_pool.log2fc_z.csv", "gene_values.long.csv"):
            f = done.parent/"tables"/name
            if f.exists():subprocess.call(["gzip", "-f", str(f)])
    log(f"{dataset}\t{how}\t{method}\t{n}\t{code}\t{secs}")


if __name__ == "__main__":
    a = argparse.ArgumentParser(); a.add_argument("--jobs", type=int, default=3); a = a.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(3) as pool:   # the six normalization runs first
        norms = dict(zip([(d, h) for d in DATASETS for h in FILTERS], pool.map(lambda x:normalization(*x), [(d, h) for d in DATASETS for h in FILTERS])))
    all_jobs = [j for (d, h), norm in norms.items() if (norm/"manifest.json").exists() for j in jobs(d, h, norm)]
    order = {"zscore":0, "zscorecmp":0, "zscoreall":0, "resampling":1, "gumbel":2}; all_jobs.sort(key=lambda j:(order[j[2]], j[0] != "Gale_2020"))   # fast first
    log(f"#\t{len(all_jobs)} runs")
    with ThreadPoolExecutor(a.jobs) as pool:list(pool.map(run, all_jobs))
    log("#\tdone")
