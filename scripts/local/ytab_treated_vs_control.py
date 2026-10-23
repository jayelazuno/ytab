#!/usr/bin/env python3
"""Treated-vs-control gene calls by resampling or the control-conditioned Gumbel (used by the app's Treated vs Control tab).

Reads a Site Normalization run (<project>/normalization/<run>) and the comparison design. The selected comparisons are
grouped by background and treated condition; each group is one test with the group's controls and treated libraries
(resampling permutes within each library when the group spans several pools). Genome FASTA and GFF come from the project.

  resampling   site-restricted permutation test, depleted and enriched; any normalization (BGC with trimmed-sum
               recentering by default). Validation status: docs/manuscript (pre-registered benchmark).
  gumbel       control-conditioned Gumbel, depletion runs only; scalar normalizations (TTR default, total reads,
               nonzero mean); needs >= 2 controls per group.

Output (ytab.fitness.contract.validate_tvc): manifest.json, genes.tsv (all groups) and <group>/sites.tsv.gz.

  ytab treated-vs-control --project-config output/projects/P/config/project.yaml --method resampling \\
      --normalization-run output/projects/P/normalization/RUN --normalization betageom --output OUT
"""
import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT/"scripts/local"
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality.control_conditioned_gumbel import NONLINEAR_NORMALIZATIONS, NORMALIZATIONS  # noqa: E402
from ytab.fitness import contract  # noqa: E402

SCALAR = set(NORMALIZATIONS)


def groups(design: Path, comparisons: list[str] | None) -> dict[str, dict]:
    rows = [r for r in csv.DictReader(open(design)) if r.get("include", "").strip().lower() in {"1", "true", "yes", "y"}]
    if comparisons:rows = [r for r in rows if r["comparison_id"] in comparisons]
    if not rows:raise ValueError("No selected comparisons were found in the design")
    out = {}
    for r in rows:
        key = re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{r['background']}__{r['treated_condition']}")
        g = out.setdefault(key, dict(background=r["background"], condition=r["treated_condition"], comparisons=[], controls=[], treatments=[]))
        g["comparisons"].append(r["comparison_id"])
        for k, s in (("controls", r["parent_sample"]), ("treatments", r["treated_sample"])):
            if s not in g[k]:g[k].append(s)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", required=True, type=Path)
    p.add_argument("--method", required=True, choices=["resampling", "gumbel"])
    p.add_argument("--normalization-run", required=True, type=Path); p.add_argument("--normalization")
    p.add_argument("--comparisons", help="Comma-separated comparison IDs (default: all included)")
    p.add_argument("--fdr", type=float, default=.05); p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    out = a.output.resolve()
    if out.exists() and any(out.iterdir()):p.exit(1, "ERROR: Output directory is nonempty; choose a new one\n")
    normalization = a.normalization or ("betageom" if a.method == "resampling" else "ttr")
    allowed = SCALAR | set(NONLINEAR_NORMALIZATIONS) if a.method == "resampling" else set(NORMALIZATIONS)
    if normalization not in allowed:
        p.exit(1, f"ERROR: {a.method} cannot use {normalization}; choose one of {', '.join(sorted(allowed))}\n")
    cfg = yaml.safe_load(a.project_config.read_text())
    resolve = lambda x:Path(x) if Path(x).is_absolute() else ROOT/x
    fasta, gff = resolve(cfg["reference"]["fasta"]), resolve(cfg["reference"]["gff"])
    project = resolve(cfg.get("output_project_dir") or f"output/projects/{cfg['project_id']}")
    methods = a.normalization_run.resolve()/"experiment/methods"
    if not (methods/normalization/"normalized_counts.npz").is_file():
        p.exit(1, f"ERROR: Normalization run {a.normalization_run.name} has no {normalization} counts; add it in Site Normalization\n")
    obs = {r["sample"]:r for r in csv.DictReader(open(methods/"observations.tsv"), delimiter="\t")}
    try:selected = groups(project/"config/comparison_design.csv", [x for x in (a.comparisons or "").split(",") if x])
    except (ValueError, KeyError, OSError) as error:p.exit(1, f"ERROR: {error}\n")
    missing = sorted({s for g in selected.values() for s in g["controls"]+g["treatments"] if s not in obs})
    if missing:p.exit(1, f"ERROR: Not in normalization run {a.normalization_run.name}: {', '.join(missing)}; rerun Site Normalization\n")
    script = SCRIPTS/("ytab_resampling.py" if a.method == "resampling" else "ytab_control_conditioned_gumbel.py")
    results, combined = {}, []
    for i, (key, g) in enumerate(selected.items(), 1):
        print(f"{i}/{len(selected)} {a.method}: {g['background']}, {g['condition']} ({len(g['controls'])} controls, {len(g['treatments'])} treated)", file=sys.stderr, flush=True)
        command = [sys.executable, script, "--normalized", methods, "--normalization", normalization, "--gff", gff, "--reference", fasta,
                   "--output", out/key, "--fdr", a.fdr,
                   *[x for s in g["controls"] for x in ("--control", obs[s]["observation_id"])],
                   *[x for s in g["treatments"] for x in ("--treatment", obs[s]["observation_id"])]]
        if a.method == "resampling":
            libraries = {obs[s]["library_id"] for s in g["controls"]+g["treatments"]}
            if len(libraries) > 1:command.append("--stratify-libraries")
            if normalization in NONLINEAR_NORMALIZATIONS:command += ["--recenter", "trimmed-sum"]
        run = subprocess.run([str(x) for x in command], cwd=ROOT, capture_output=True, text=True)
        sys.stderr.write(run.stderr)
        if run.returncode:
            results[key] = dict(g, status="failed", message=(run.stderr.strip().splitlines() or ["failed"])[-1].removeprefix("ERROR: "))
            continue
        manifest = json.loads((out/key/"run_manifest.json").read_text())
        diagnostics = manifest.get("diagnostics", {})
        warning = (diagnostics.get("replicate_noise") or {}).get("warning")
        results[key] = dict(g, status="success", calls=diagnostics.get("calls"), warning=warning, stratified=a.method == "resampling" and "--stratify-libraries" in command)
        for row in csv.DictReader(open(out/key/"genes.tsv"), delimiter="\t"):
            if a.method == "resampling":   # plot abundance: mean normalized count per site over both conditions
                means = [float(row[k]) for k in ("mean_control", "mean_treatment") if row.get(k) not in (None, "")]
                row["mean_abundance"] = sum(means)/2 if len(means) == 2 else ""
            combined.append(dict(group=key, background=g["background"], condition=g["condition"], **row))
    if combined:
        with open(out/"genes.tsv", "w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for r in combined for k in r)), delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(combined)
    ok = [k for k, v in results.items() if v["status"] == "success"]
    status = "success" if len(ok) == len(results) else "partial" if ok else "failed"
    labels = dict(resampling="Resampling", gumbel="Control-conditioned Gumbel (depletion only)")
    manifest = dict(status=status, contract=contract.VERSION, kind="treated_vs_control", method=a.method, method_label=labels[a.method], project=cfg.get("project_id"),
                    normalization=normalization, normalization_run=str(a.normalization_run.resolve()), fdr=a.fdr, groups=results,
                    files=dict(genes="genes.tsv"))
    out.mkdir(parents=True, exist_ok=True)
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2, default=str)+"\n")
    problems = contract.validate_tvc(out)
    if problems:print("ERROR: the result does not meet the output contract: " + "; ".join(problems), file=sys.stderr); return 1
    print(json.dumps({k:manifest[k] for k in ("status", "method", "normalization")} | {k:v.get("calls") or v.get("message") for k, v in results.items()}, indent=1, default=str))
    if status == "failed":
        print("ERROR: "+"; ".join(f"{k}: {v['message']}" for k, v in results.items()), file=sys.stderr); return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
