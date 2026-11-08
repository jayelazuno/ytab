#!/usr/bin/env python3
"""Method-comparison report, multi-condition arm: each drug x background dose series of Nickels 2025 tested jointly
(`ytab multicondition-project`: permutation ANOVA over control + doses, and all-pairs resampling), BGC, on the report's
normalization runs with all sites and with the top site dropped. Untracked study (output/comparison/ is ignored).

Conditions: every no-drug control library of the background is one "control" condition (the reference); each dose is
its own condition. Output: output/comparison/report/<dataset>/<filter>/multicondition_aBGC/<background>__<drug>/.
  python scripts/local/ytab_method_comparison_multicondition.py [--jobs 2]
"""
from __future__ import annotations
import argparse, csv, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]; PY = sys.executable; OUT = ROOT/"output/comparison/report"
DATASET, BACKGROUNDS, FILTERS, NORMALIZATION = "Nickels_2025_glabrata", ("BG2", "pdr1_KO"), ("none", "drop-top"), "betageom"


def series():
    """(background, drug) -> (control samples, [(sample, dose condition)])."""
    rows = [r for r in csv.DictReader(open(ROOT/"output/projects"/DATASET/"config/comparison_design.csv"))
            if r["background"] in BACKGROUNDS and str(r["include"]).strip().lower() in ("1", "true", "yes", "y")]
    out = {}
    for r in rows:
        drug = re.sub(r"_[0-9.]+ng$", "", r["treated_condition"])
        controls, treated = out.setdefault((r["background"], drug), ([], []))
        if r["parent_sample"] not in controls:controls.append(r["parent_sample"])
        if (r["treated_sample"], r["treated_condition"]) not in treated:treated.append((r["treated_sample"], r["treated_condition"]))
    return out


def run(job):
    how, (bg, drug), (controls, treated) = job
    norm = OUT/DATASET/how/"norm"; out = OUT/DATASET/how/f"multicondition_{NORMALIZATION}"/f"{bg}__{drug}"
    if (out/"manifest.json").exists():return
    if out.exists():subprocess.call(["rm", "-rf", str(out)])
    samples = [f"{DATASET}::{s}=control" for s in controls]+[f"{DATASET}::{s}={c}" for s, c in treated]
    cmd = [PY, "scripts/local/ytab_multicondition_project.py", "--project-config", ROOT/"output/projects"/DATASET/"config/project.yaml",
           "--normalized", norm/"experiment/methods", "--normalization", NORMALIZATION, "--reference-condition", "control",
           *[x for s in samples for x in ("--sample", s)], "--output", out]
    t = time.time(); out.parent.mkdir(parents=True, exist_ok=True)
    with open(out.parent/f"{bg}__{drug}.log", "w") as h:code = subprocess.call([str(c) for c in cmd], stdout=h, stderr=subprocess.STDOUT, cwd=ROOT)
    with open(OUT/"status.tsv", "a") as h:h.write(f"{DATASET}\t{how}\tmulticondition\t{bg}__{drug}\t{code}\t{round(time.time()-t)}\n")


if __name__ == "__main__":
    a = argparse.ArgumentParser(); a.add_argument("--jobs", type=int, default=2); a = a.parse_args()
    jobs = [(how, key, value) for how in FILTERS for key, value in sorted(series().items())]
    with ThreadPoolExecutor(a.jobs) as pool:list(pool.map(run, jobs))
    print("done")
