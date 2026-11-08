#!/usr/bin/env python3
"""Claims check, step 1: YTAB's standard condition comparisons on every public dataset with a comparison design.

Per project: one normalization run (aBGC, the condition-comparison default, and TTR for the Gumbel), then on all of
its included comparisons resampling (aBGC), the local-abundance z-score (aBGC) and the control-conditioned Gumbel (TTR;
it declines conditions with fewer than two treated libraries). Untracked study; output/claims_check/ is not tracked.
Output: output/claims_check/<project>/{norm,resampling_aBGC,zscore_aBGC,gumbel_ttr}/ and status.tsv.
  python scripts/local/ytab_claims_check_runs.py [--jobs 2] [--projects A,B]
"""
from __future__ import annotations
import argparse, csv, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[2]; PY = sys.executable; OUT = ROOT/"output/claims_check"
OWN = {"H2O2_screen_v1", "Zn_toxicity_screen", "local_glabrata_smoke_v1"}   # the lab's own screens: not published claims


def projects():
    out = []
    for cfg in sorted((ROOT/"output/projects").glob("*/config/comparison_design.csv")):
        name = cfg.parents[1].name
        rows = [r for r in csv.DictReader(open(cfg)) if str(r.get("include", "")).strip().lower() in ("1", "true", "yes", "y")]
        if name not in OWN and rows:out.append(name)
    return out


def call(cmd, log):
    log.parent.mkdir(parents=True, exist_ok=True); t = time.time()
    with open(log, "w") as h:code = subprocess.call([str(c) for c in cmd], stdout=h, stderr=subprocess.STDOUT, cwd=ROOT)
    return code, round(time.time()-t)


def record(project, step, code, secs):
    with open(OUT/"status.tsv", "a") as h:h.write(f"{project}\t{step}\t{code}\t{secs}\n")
    print(project, step, code, secs, flush=True)


def run(project):
    cfg = ROOT/"output/projects"/project/"config/project.yaml"; base = OUT/project; norm = base/"norm"
    if not (norm/"manifest.json").exists():
        if norm.exists():subprocess.call(["rm", "-rf", str(norm)])
        record(project, "normalization", *call([PY, "scripts/local/ytab_normalize_project.py", "--project-config", cfg, "--method", "aBGC", "--method", "ttr",
                                                "--output", norm], base/"norm.log"))
        if not (norm/"manifest.json").exists():return
    gff = yaml.safe_load(open(cfg))["reference"]["gff"]; gff = gff if Path(gff).is_absolute() else str(ROOT/gff)
    out = base/"zscore_aBGC"
    if not (out/"zscore_manifest.json").exists():
        record(project, "zscore", *call([PY, "scripts/local/ytab_fitness_zscore.py", "--project-id", project, "--analysis-id", "claims_check",
                                         "--comparison-design", ROOT/"output/projects"/project/"config/comparison_design.csv", "--output-dir", out,
                                         "--normalized", norm, "--normalization", "aBGC", "--gff", gff], base/"zscore.log"))
        for name in ("top_hit_evidence.by_pair.csv", "treated_vs_parent.by_pool.log2fc_z.csv", "gene_values.long.csv"):
            if (out/"tables"/name).exists():subprocess.call(["gzip", "-f", str(out/"tables"/name)])
    for method, n in (("resampling", "aBGC"), ("gumbel", "ttr")):
        out = base/f"{method}_{n}"
        if not (out/"manifest.json").exists():
            if out.exists():subprocess.call(["rm", "-rf", str(out)])
            record(project, method, *call([PY, "scripts/local/ytab_treated_vs_control.py", "--project-config", cfg, "--method", method,
                                           "--normalization-run", norm, "--normalization", n, "--output", out], base/f"{method}.log"))


if __name__ == "__main__":
    a = argparse.ArgumentParser(); a.add_argument("--jobs", type=int, default=2); a.add_argument("--projects")
    a = a.parse_args(); OUT.mkdir(parents=True, exist_ok=True)
    todo = a.projects.split(",") if a.projects else projects()
    print("projects:", ", ".join(todo), flush=True)
    with ThreadPoolExecutor(a.jobs) as pool:list(pool.map(run, todo))
    print("done")
