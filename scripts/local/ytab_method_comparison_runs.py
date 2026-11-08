#!/usr/bin/env python3
"""Method-comparison study, step 1: run the three treated-vs-control methods with each method's default normalization
(z-score BGC, resampling BGC, Gumbel TTR) on every eligible dataset, into output/comparison/<dataset>/runs/<method>.

Eligible: included comparisons whose background has >= 2 controls (the z-score's control-control noise model and the
Gumbel's control count model both need two). Backgrounds with one control get resampling only. Project folders are not
written; the study is untracked (/output/comparison/ in .gitignore).

  python scripts/local/ytab_method_comparison_runs.py [--datasets Gale_2020,...] [--jobs 2]
"""
from __future__ import annotations
import argparse, csv, glob, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[2]; OUT = ROOT/"output/comparison"; PY = sys.executable
SKIP = {"H2O2_screen_v1", "Zn_toxicity_screen", "local_glabrata_smoke_v1"}   # unpublished / test; retired before testing


def included(d): return [r for r in csv.DictReader(open(d)) if str(r.get("include", "")).strip().lower() in ("1", "true", "yes", "y")]


def plan(project: Path) -> dict | None:
    design = project/"config/comparison_design.csv"
    if not design.exists():return None
    rows = included(design)
    if not rows:return None
    controls = {}
    for r in rows:controls.setdefault(r["background"], set()).add(r["parent_sample"])
    runs = sorted(glob.glob(str(project/"normalization/*/experiment/methods")), key=os.path.getmtime, reverse=True)
    norm = next((Path(x).parents[1] for x in runs if Path(x, "betageom").is_dir() and Path(x, "ttr").is_dir()), None)
    if norm is None:return None
    cfg = yaml.safe_load(open(project/"config/project.yaml")); gff = cfg["reference"]["gff"]
    two = [r["comparison_id"] for r in rows if len(controls[r["background"]]) >= 2]
    return dict(dataset=project.name, norm=norm, gff=gff if os.path.isabs(gff) else str(ROOT/gff), design=design,
                config=project/"config/project.yaml", all=[r["comparison_id"] for r in rows], two=two,
                controls={k: len(v) for k, v in controls.items()})


def run(cmd, log):
    t = time.time()
    with open(log, "w") as h:code = subprocess.call(cmd, stdout=h, stderr=subprocess.STDOUT, cwd=ROOT)
    return code, round(time.time()-t)


def dataset_jobs(p):
    base = OUT/p["dataset"]/"runs"; base.mkdir(parents=True, exist_ok=True)
    jobs = [("resampling", [PY, "scripts/local/ytab_treated_vs_control.py", "--project-config", str(p["config"]), "--method", "resampling",
             "--normalization-run", str(p["norm"]), "--normalization", "betageom", "--comparisons", ",".join(p["all"]), "--output", str(base/"resampling")])]
    if p["two"]:
        jobs += [("zscore", [PY, "scripts/local/ytab_fitness_zscore.py", "--project-id", p["dataset"], "--analysis-id", "comparison_zscore_BGC",
                  "--comparison-design", str(p["design"]), "--comparisons", ",".join(p["two"]), "--output-dir", str(base/"zscore"),
                  "--normalized", str(p["norm"]), "--normalization", "betageom", "--gff", p["gff"]]),
                 ("gumbel", [PY, "scripts/local/ytab_treated_vs_control.py", "--project-config", str(p["config"]), "--method", "gumbel",
                  "--normalization-run", str(p["norm"]), "--normalization", "ttr", "--comparisons", ",".join(p["two"]), "--output", str(base/"gumbel")])]
    return [(p["dataset"], m, c, base/f"{m}.log") for m, c in jobs]


def main():
    a = argparse.ArgumentParser(); a.add_argument("--datasets"); a.add_argument("--jobs", type=int, default=2); a = a.parse_args()
    projects = [Path(x) for x in sorted(glob.glob(str(ROOT/"output/projects/*"))) if Path(x).name not in SKIP]
    if a.datasets:projects = [p for p in projects if p.name in a.datasets.split(",")]
    plans = [p for p in (plan(x) for x in projects) if p]
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/"plan.json").write_text(json.dumps([{k: str(v) if isinstance(v, Path) else v for k, v in p.items()} for p in plans], indent=1))
    jobs = [j for p in plans for j in dataset_jobs(p) if not (OUT/j[0]/"runs"/j[1]/("zscore_manifest.json" if j[1] == "zscore" else "manifest.json")).exists()]
    jobs.sort(key=lambda j: (j[0] == "Nickels_2025_glabrata", j[1] == "gumbel"))   # small and fast first
    status = OUT/"runs_status.tsv"
    def go(j):
        code, secs = run(j[2], j[3])
        with open(status, "a") as h:h.write(f"{j[0]}\t{j[1]}\t{code}\t{secs}\n")
        print(f"{j[0]} {j[1]} exit={code} {secs}s", flush=True)
    print(f"{len(plans)} datasets, {len(jobs)} runs", flush=True)
    with ThreadPoolExecutor(a.jobs) as pool:list(pool.map(go, jobs))


if __name__ == "__main__":main()
