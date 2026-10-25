#!/usr/bin/env python3
"""Multi-condition tests and their plots in one step (used by the app's Condition Comparisons tab).

Runs `ytab multicondition` (permutation ANOVA and all-pairs resampling) on a normalization run, then
`ytab multicondition-plots` (response correlation and heatmap). The genome FASTA and GFF come from the project.
A plotting failure (e.g. fewer than two non-reference conditions for a correlation) is reported, not fatal.

  ytab multicondition-project --project-config output/projects/P/config/project.yaml \\
      --normalized output/projects/P/normalization/RUN/experiment/methods --normalization betageom \\
      --sample P::ctrl1=control --sample P::flu1=FLU ... --reference-condition control --output OUT

--adopt TESTS [--adopt-plots PLOTS] turns the outputs of an earlier `ytab multicondition` (and
`ytab multicondition-plots`) run into a run of this command (folders linked, not copied; manifest as for a new run),
so results made before this command existed show in the app without being recomputed.
"""
import argparse
from collections import Counter
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT/"scripts/local"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", required=True, type=Path)
    p.add_argument("--normalized"); p.add_argument("--normalization", default="betageom")
    p.add_argument("--sample", action="append", help="OBSERVATION_ID=CONDITION; repeat")
    p.add_argument("--reference-condition")
    p.add_argument("--adopt", type=Path, help="Folder of an earlier `ytab multicondition` run (run_manifest.json, anova/pairs tables)")
    p.add_argument("--adopt-plots", type=Path, help="With --adopt: its `ytab multicondition-plots` folder")
    p.add_argument("--test", choices=["both", "anova", "all-pairs"], default="both")
    p.add_argument("--stratify-libraries", action="store_true")
    p.add_argument("--independent-libraries", action="store_true", help="Conditions are separate, independently made libraries (e.g. one per strain)")
    p.add_argument("--fdr", type=float, default=.05)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    out = a.output.resolve()
    if out.exists() and any(out.iterdir()):p.exit(1, "ERROR: Output directory is nonempty; choose a new one\n")
    cfg = yaml.safe_load(a.project_config.read_text())
    if a.adopt:
        out.mkdir(parents=True, exist_ok=True)
        link = lambda target, name:(out/name).symlink_to(os.path.relpath(target.resolve(), out))
        link(a.adopt, "tests")
        if a.adopt_plots:link(a.adopt_plots, "plots")
        write_manifest(out, cfg, a.normalization, a.normalized, plots_ok=(out/"plots/plots_manifest.json").exists(), plots_message=None, adopted=str(a.adopt.resolve()))
        return
    if not (a.normalized and a.sample and a.reference_condition):p.error("--normalized, --sample and --reference-condition are required")
    resolve = lambda x:Path(x) if Path(x).is_absolute() else ROOT/x
    fasta, gff = resolve(cfg["reference"]["fasta"]), resolve(cfg["reference"]["gff"])
    command = [sys.executable, SCRIPTS/"ytab_multicondition.py", "--normalized", a.normalized, "--normalization", a.normalization,
               "--gff", gff, "--reference", fasta, "--output", out/"tests", "--reference-condition", a.reference_condition,
               "--test", a.test, "--fdr", a.fdr, *(["--stratify-libraries"] if a.stratify_libraries else []),
               *(["--independent-libraries"] if a.independent_libraries else []),
               *[x for s in a.sample for x in ("--sample", s)]]
    print("1/2 Running the multi-condition tests", file=sys.stderr, flush=True)
    test = subprocess.run([str(x) for x in command], cwd=ROOT, capture_output=True, text=True)
    sys.stderr.write(test.stderr)
    if test.returncode:p.exit(1, (test.stderr.strip().splitlines() or ["ERROR: multicondition failed"])[-1]+"\n")
    print("2/2 Drawing response correlation and heatmap", file=sys.stderr, flush=True)
    plots = subprocess.run([sys.executable, str(SCRIPTS/"ytab_multicondition_plots.py"), "--run", str(out/"tests"), "--output", str(out/"plots")],
                           cwd=ROOT, capture_output=True, text=True)
    sys.stderr.write(plots.stderr)
    write_manifest(out, cfg, a.normalization, a.normalized, plots.returncode == 0,
                   None if plots.returncode == 0 else (plots.stderr.strip().splitlines() or [""])[-1])


def rows(path):
    with open(path, newline="") as handle:return list(csv.DictReader(handle, delimiter="\t"))


def write_manifest(out, cfg, normalization, normalized, plots_ok, plots_message, adopted=None):
    run = json.loads((out/"tests/run_manifest.json").read_text())
    plot_manifest = json.loads((out/"plots/plots_manifest.json").read_text()) if (out/"plots/plots_manifest.json").exists() else None
    anova = (run["diagnostics"].get("anova") or {})
    pairs = {k:v["calls"] for k, v in (run["diagnostics"].get("all_pairs") or {}).items() if isinstance(v, dict) and "calls" in v}
    samples = run["samples"]
    if not anova.get("calls") and (out/"tests/anova_genes.tsv").exists():   # older runs did not record the counts: count the table
        anova = dict(anova, calls=dict(Counter(r["call"] if r["status"] == "tested" else "not tested" for r in rows(out/"tests/anova_genes.tsv"))))
    if not pairs and (out/"tests/pairs_genes.tsv").exists():
        for r in rows(out/"tests/pairs_genes.tsv"):
            c = pairs.setdefault(r["pair"], {}); k = r["call_all_pairs"] if r["status"] == "tested" else "not tested"; c[k] = c.get(k, 0)+1
    conditions = run.get("conditions") or list(dict.fromkeys(s["condition"] for s in samples))   # older runs: from the samples
    recorded = run.get("normalization")   # a method name, or (ytab multicondition) a record of the method and its input
    if isinstance(recorded, dict):normalized = normalized or recorded.get("input"); recorded = recorded.get("method")
    manifest = dict(status="success", project=cfg.get("project_id"), normalization=recorded or normalization, normalized=normalized,
                    conditions=conditions, reference_condition=run["reference_condition"], stratified=bool(run.get("stratified", False)),
                    samples=samples, tests=run.get("tests") or [k for k in ("anova", "all-pairs") if (out/"tests"/f"{'anova' if k == 'anova' else 'pairs'}_genes.tsv").exists()],
                    anova_calls=anova.get("calls"), anova_warning=anova.get("warning"),
                    all_pairs_calls=pairs, recommendation=run.get("recommendation"), adopted_from=adopted,
                    plots=dict(status="success" if plots_ok else "failed", message=plots_message,
                               notes=(plot_manifest or {}).get("notes", []), responsive_genes=(plot_manifest or {}).get("responsive_genes")),
                    files=dict(anova="tests/anova_genes.tsv", pairs="tests/pairs_genes.tsv", correlation="plots/response_correlation.tsv",
                               heatmap="plots/heatmap_matrix.tsv"))
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2, default=str)+"\n")
    print(json.dumps({k:manifest[k] for k in ("status", "conditions", "anova_calls", "all_pairs_calls")}, indent=1, default=str))


if __name__ == "__main__":main()
