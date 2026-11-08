#!/usr/bin/env python3
"""Pre-registered: condition-comparison recovery against library density (Hermes: Gale fluconazole; Ac/Ds: SATAY rapamycin).

  python scripts/local/ytab_condition_density.py freeze
  python scripts/local/ytab_condition_density.py run [--workers 3]   # thin, normalize (BGC), resampling; reads no positives
  python scripts/local/ytab_condition_density.py evaluate            # recovery, adequate density; writes result.json, levels.tsv
Add `--design pool` for the corrected design (R23): every library grown from one pool loses the same sites.

Designs: output/conditions_density/PREREGISTRATION.md (R22, each library thinned independently) and
output/conditions_density_pool/PREREGISTRATION.md (R23, thinned by pool).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.normalization.ttr import sha256

DESIGN = os.environ.get("YTAB_DENSITY_DESIGN", "independent")   # --design; read here so process-pool workers see it
OUT = ROOT/("output/conditions_density_pool" if DESIGN == "pool" else "output/conditions_density")
SCRIPTS = ROOT/"scripts/local"
QS, SEEDS = [1, 2, 4, 8, 16], [1, 2]
DATASETS = {
    "hermes_gale_fluconazole": dict(project="Gale_2020", transposon="Hermes",
                                    samples=["SRR11561505", "SRR11561504", "SRR11561503", "SRR11561502"],
                                    comparisons="Cg1_FLU1_vs_control1,Cg1_FLU2_vs_control2", direction="depleted"),
    "acds_satay_rapamycin": dict(project="Michel_2017_SATAY", transposon="Ac/Ds", samples=["ERR1533148", "ERR1533149"],
                                 comparisons="BY4741_ade2_poolWT2_treated_vs_parent", direction="enriched"),
}
RAPAMYCIN = ["FPR1", "RRD1", "TIP41", "GLN3", "SAP190", "PSP2", "CCS1", "ESL2", "PPH21", "PPM1"]


def freeze():
    stamp = OUT/"PREREGISTRATION.sha256"
    if stamp.exists():raise SystemExit("already frozen")
    stamp.write_text(f"{sha256(OUT/'PREREGISTRATION.md')} {dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}\n")
    print(stamp.read_text().strip())


def check():
    frozen = (OUT/"PREREGISTRATION.sha256").read_text().split()
    if sha256(OUT/"PREREGISTRATION.md") != frozen[0]:raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def keep(lib, chrom, pos, q, seed):   # the saturation benchmark's rule (ytab_saturation_benchmark.keep); lib is the pool in the pool design
    return q == 1 or int(hashlib.md5(f"{lib}:{chrom}:{pos}:{seed}".encode()).hexdigest()[:8], 16)/0xFFFFFFFF < 1/q


def levels():
    return [(name, q, seed) for name in DATASETS for q in QS for seed in ([0] if q == 1 else SEEDS)]


def level_dir(name, q, seed):
    return OUT/"runs"/name/f"q{q}_s{seed}"


def density(hit_files, gff):
    """Per library: median insertion sites per coding kb over genes (sites within CDS, either strand)."""
    genes = [g for g in ccg.annotations(Path(gff)) if g["blocks"]]
    out = []
    for path in hit_files:
        sites = {}
        with open(path) as handle:
            r = csv.reader(handle, delimiter="\t"); h = next(r); pos = h.index("Hit position")
            for row in r:sites.setdefault(row[0], set()).add(int(row[pos]))
        arrays = {c:np.array(sorted(s)) for c, s in sites.items()}; per_kb = []
        for g in genes:
            a = arrays.get(g["chrom"]); length = sum(b-s+1 for s, b in g["blocks"])
            n = sum(int(np.searchsorted(a, b, "right")-np.searchsorted(a, s, "left")) for s, b in g["blocks"]) if a is not None else 0
            per_kb.append(n/length*1000)
        out.append(float(np.median(per_kb)))
    return out


def run_level(task):
    name, q, seed = task; d = DATASETS[name]; out = level_dir(name, q, seed)
    if (out/"done.json").exists():return json.loads((out/"done.json").read_text())
    src_cfg = ROOT/"output/projects"/d["project"]/"config"; cfg = yaml.safe_load((src_cfg/"project.yaml").read_text())
    (out/"config").mkdir(parents=True, exist_ok=True); files = []
    with open(src_cfg/"sample_sheet.csv") as handle:pools = {r["sample"]:(r.get("pool") or r["sample"]) for r in csv.DictReader(handle)}
    for sample in d["samples"]:
        key = pools[sample] if DESIGN == "pool" else sample
        src = ROOT/"output/projects"/d["project"]/"create_hit_file"/sample/f"{sample}_hits.txt"; dst = out/"create_hit_file"/sample/f"{sample}_hits.txt"
        dst.parent.mkdir(parents=True, exist_ok=True)
        with open(src) as f, open(dst, "w") as o:
            r = csv.reader(f, delimiter="\t"); w = csv.writer(o, delimiter="\t", lineterminator="\n"); h = next(r); w.writerow(h); pos = h.index("Hit position")
            for row in r:
                if keep(key, row[0], row[pos], q, seed):w.writerow(row)
        files.append(dst)
    with open(src_cfg/"sample_sheet.csv") as handle:rows = list(csv.DictReader(handle))
    fields = list(rows[0]) + ([] if "hit_file" in rows[0] else ["hit_file"])
    with open(out/"config/sample_sheet.csv", "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields); w.writeheader()
        for r in rows:
            r = dict(r); r["include"] = "True" if r["sample"] in d["samples"] else "False"
            r["hit_file"] = str((out/"create_hit_file"/r["sample"]/f"{r['sample']}_hits.txt").relative_to(ROOT)) if r["sample"] in d["samples"] else ""
            w.writerow(r)
    shutil.copy2(src_cfg/"comparison_design.csv", out/"config/comparison_design.csv")
    rel = out.relative_to(ROOT)
    cfg |= dict(output_project_dir=str(rel), sample_sheet=str(rel/"config/sample_sheet.csv"), repo_root="../"*(len(rel.parts)+1))
    cfg["repo_root"] = cfg["repo_root"].rstrip("/")
    if "comparison_design" in cfg:cfg["comparison_design"] = str(rel/"config/comparison_design.csv")
    (out/"config/project.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    log = out/"run.log"
    with open(log, "w") as handle:
        for cmd in ([SCRIPTS/"ytab_normalize_project.py", "--project-config", out/"config/project.yaml", "--output", out/"normalization",
                     "--method", "betageom", "--normalize-together", *[x for s in d["samples"] for x in ("--sample", s)]],
                    [SCRIPTS/"ytab_treated_vs_control.py", "--project-config", out/"config/project.yaml", "--method", "resampling",
                     "--normalization-run", out/"normalization", "--normalization", "betageom", "--comparisons", d["comparisons"],
                     "--output", out/"resampling"]):
            r = subprocess.run([sys.executable, *map(str, cmd)], stdout=handle, stderr=subprocess.STDOUT, cwd=ROOT,
                               env={**os.environ, "PYTHONPATH":str(ROOT/"src")})
            if r.returncode:raise RuntimeError(f"{name} q{q} s{seed} failed; see {log}")
    rec = dict(dataset=name, q=q, seed=seed, sites_per_coding_kb=density(files, ROOT/cfg["reference"]["gff"]))
    (out/"done.json").write_text(json.dumps(rec)); return rec


def run(workers):
    check()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for rec in pool.map(run_level, levels()):print(rec, flush=True)


def evaluate():
    prereg = check(); rows = []
    strict = {r["ORF"] for r in csv.DictReader(open(ROOT/"output/essentiality/resampling/v3_gale_prioritization/strict_published_genes.tsv"), delimiter="\t")}
    full_calls = {}
    for name, q, seed in levels():
        d = DATASETS[name]; out = level_dir(name, q, seed); rec = json.loads((out/"done.json").read_text())
        genes = list(csv.DictReader(open(out/"resampling/genes.tsv"), delimiter="\t"))
        hit = (lambda g: g["ORF"] in strict) if name.startswith("hermes") else (lambda g: g["Name"] in RAPAMYCIN)
        positives = [g for g in genes if hit(g)]
        calls = {(g["ORF"], g["call"]) for g in genes if g["call"] in ("depleted", "enriched")}
        if q == 1:full_calls[name] = calls
        stat = {g["ORF"]:-float(g["log2FC"]) if d["direction"] == "depleted" else float(g["log2FC"]) for g in genes if g["log2FC"] not in ("", "nan")}
        ranked = sorted(stat, key=stat.get, reverse=True)
        rows.append(dict(dataset=name, transposon=d["transposon"], q=q, seed=seed, sites_per_coding_kb=float(np.mean(rec["sites_per_coding_kb"])),
                         recovered=sum(g["call"] == d["direction"] for g in positives), positives=len(strict) if name.startswith("hermes") else len(RAPAMYCIN),
                         depleted=sum(c == "depleted" for _, c in calls), enriched=sum(c == "enriched" for _, c in calls),
                         tested_fraction=sum(g["status"] == "tested" for g in genes)/len(genes),
                         full_calls_kept=len(calls & full_calls[name])/max(1, len(full_calls[name])), new_calls=len(calls - full_calls[name]),
                         median_rank_of_positives=float(np.median([ranked.index(g["ORF"])+1 if g["ORF"] in stat else len(genes) for g in positives]))))
    result = dict(preregistration=prereg, datasets={})
    for name in DATASETS:
        r = [x for x in rows if x["dataset"] == name]
        by_q = {q:dict(recovered=float(np.mean([x["recovered"] for x in r if x["q"] == q])),
                       sites_per_coding_kb=float(np.mean([x["sites_per_coding_kb"] for x in r if x["q"] == q]))) for q in QS}
        target = .9*by_q[1]["recovered"]; adequate = 1
        for q in QS:
            if by_q[q]["recovered"] >= target:adequate = q
            else:break
        result["datasets"][name] = dict(transposon=DATASETS[name]["transposon"], full_recovered=by_q[1]["recovered"], by_q=by_q, adequate_q=adequate,
                                        adequate_sites_per_coding_kb=by_q[adequate]["sites_per_coding_kb"],
                                        floor_reached=adequate == QS[-1])
    with open(OUT/"levels.tsv", "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    (OUT/"result.json").write_text(json.dumps(result, indent=1)); print(json.dumps(result, indent=1))


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("phase", choices=["freeze", "run", "evaluate"]); p.add_argument("--workers", type=int, default=3)
    p.add_argument("--design", choices=["independent", "pool"], default="independent")
    a = p.parse_args()
    if a.design != DESIGN:os.environ["YTAB_DENSITY_DESIGN"] = a.design; os.execv(sys.executable, [sys.executable, *sys.argv])
    {"freeze":freeze, "run":lambda: run(a.workers), "evaluate":evaluate}[a.phase]()
