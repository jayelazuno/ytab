#!/usr/bin/env python3
"""Method-comparison study, normalization sweep on Gale 2020 (BG14, fluconazole): every method x every normalization it
accepts, on the real comparison and on the swap null (control 1 + FLU 1 vs control 2 + FLU 2: drug effects cancel, so
calls there are false calls). Output: output/comparison/Gale_2020/sweep/<design>/<method>_<normalization>/.

z-score: every method of the normalization run (including nonorm = unnormalized) and CPM from the SummaryTables;
resampling: ttr, totreads, nzmean, betageom, aBGC; Gumbel: ttr, totreads, nzmean (what each accepts; nothing forced).
"""
from __future__ import annotations
import csv, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[2]; PY = sys.executable
PROJECT = ROOT/"output/projects/Gale_2020"; NORM = PROJECT/"normalization/norm_20261022_all_methods"
OUT = ROOT/"output/comparison/Gale_2020/sweep"
C1, C2, F1, F2 = "SRR11561505", "SRR11561504", "SRR11561503", "SRR11561502"   # control 1, control 2, FLU 1, FLU 2
ZSCORE = ["nonorm", "totreads", "nzmean", "ttr", "betageom", "aBGC", "quantile", "zinfnb", "midlc", "cpm"]
RESAMPLING = ["ttr", "totreads", "nzmean", "betageom", "aBGC"]
GUMBEL = ["ttr", "totreads", "nzmean"]


def designs() -> dict[str, Path]:
    """real: the published comparison; swap_null: each side one control + one FLU culture."""
    cfg = yaml.safe_load(open(PROJECT/"config/project.yaml")); out = {}
    for name, rows in {"real": [(C1, F1, "Cg1_FLU1_vs_control1"), (C2, F2, "Cg1_FLU2_vs_control2")],
                       "swap_null": [(C1, C2, "swap_ctrl1_vs_ctrl2"), (F1, F2, "swap_flu1_vs_flu2")]}.items():
        base = OUT/name/"project"; (base/"config").mkdir(parents=True, exist_ok=True)
        with open(base/"config/comparison_design.csv", "w", newline="") as h:
            w = csv.writer(h); w.writerow(["comparison_id", "pair_id", "parent_sample", "treated_sample", "background", "pool", "parent_condition", "treated_condition", "include", "notes"])
            for p, t, cid in rows:w.writerow([cid, cid, p, t, "BG14", "Cg1", "parent", "treated", "true", f"method-comparison sweep: {name}"])
        y = dict(cfg, output_project_dir=str(base)); (base/"config/project.yaml").write_text(yaml.safe_dump(y, sort_keys=False))
        out[name] = base
    return out


def jobs(design: str, base: Path):
    cfg, comp = base/"config/project.yaml", base/"config/comparison_design.csv"; gff = yaml.safe_load(open(cfg))["reference"]["gff"]
    gff = gff if Path(gff).is_absolute() else str(ROOT/gff); out = OUT/design; js = []
    for n in ZSCORE:
        cmd = [PY, "scripts/local/ytab_fitness_zscore.py", "--project-id", "Gale_2020", "--analysis-id", f"sweep_{design}_{n}", "--comparison-design", str(comp), "--output-dir", str(out/f"zscore_{n}")]
        cmd += ([a for s in (C1, C2, F1, F2) for a in ("--sample-table", f"{s}={PROJECT}/summary/{s}/{s}.feature_table.RDF_1.csv")] if n == "cpm"
                else ["--normalized", str(NORM), "--normalization", n, "--gff", gff])
        js.append((design, "zscore", n, cmd))
    for m, norms in (("resampling", RESAMPLING), ("gumbel", GUMBEL)):
        for n in norms:
            js.append((design, m, n, [PY, "scripts/local/ytab_treated_vs_control.py", "--project-config", str(cfg), "--method", m, "--normalization-run", str(NORM),
                                      "--normalization", n, "--output", str(out/f"{m}_{n}")]))
    return js


def run(j):
    design, m, n, cmd = j; done = OUT/design/f"{m}_{n}"/("zscore_manifest.json" if m == "zscore" else "manifest.json")
    if done.exists():return
    t = time.time()
    with open(OUT/design/f"{m}_{n}.log", "w") as h:code = subprocess.call(cmd, stdout=h, stderr=subprocess.STDOUT, cwd=ROOT)
    with open(OUT/"status.tsv", "a") as h:h.write(f"{design}\t{m}\t{n}\t{code}\t{round(time.time()-t)}\n")
    print(f"{design} {m} {n} exit={code} {round(time.time()-t)}s", flush=True)


if __name__ == "__main__":
    all_jobs = [j for d, b in designs().items() for j in jobs(d, b)]
    all_jobs.sort(key=lambda j: j[1] == "gumbel")   # fast first
    print(f"{len(all_jobs)} runs", flush=True)
    with ThreadPoolExecutor(3) as pool:list(pool.map(run, all_jobs))
