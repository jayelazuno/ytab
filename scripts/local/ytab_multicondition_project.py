#!/usr/bin/env python3
"""Multi-condition tests and their plots in one step (used by the app's Condition Comparisons tab).

Runs `ytab multicondition` (permutation ANOVA and all-pairs resampling) on a normalization run, then
`ytab multicondition-plots` (response correlation and heatmap). The genome FASTA and GFF come from the project.
A plotting failure (e.g. fewer than two non-reference conditions for a correlation) is reported, not fatal.

  ytab multicondition-project --project-config output/projects/P/config/project.yaml \\
      --normalized output/projects/P/normalization/RUN/experiment/methods --normalization betageom \\
      --sample P::ctrl1=control --sample P::flu1=FLU ... --reference-condition control --output OUT
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT/"scripts/local"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", required=True, type=Path)
    p.add_argument("--normalized", required=True); p.add_argument("--normalization", default="betageom")
    p.add_argument("--sample", action="append", required=True, help="OBSERVATION_ID=CONDITION; repeat")
    p.add_argument("--reference-condition", required=True)
    p.add_argument("--test", choices=["both", "anova", "all-pairs"], default="both")
    p.add_argument("--stratify-libraries", action="store_true")
    p.add_argument("--independent-libraries", action="store_true", help="Conditions are separate, independently made libraries (e.g. one per strain)")
    p.add_argument("--fdr", type=float, default=.05)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    out = a.output.resolve()
    if out.exists() and any(out.iterdir()):p.exit(1, "ERROR: Output directory is nonempty; choose a new one\n")
    cfg = yaml.safe_load(a.project_config.read_text())
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
    run = json.loads((out/"tests/run_manifest.json").read_text())
    plot_manifest = json.loads((out/"plots/plots_manifest.json").read_text()) if (out/"plots/plots_manifest.json").exists() else None
    anova = (run["diagnostics"].get("anova") or {})
    pairs = {k:v["calls"] for k, v in (run["diagnostics"].get("all_pairs") or {}).items() if isinstance(v, dict) and "calls" in v}
    manifest = dict(status="success", project=cfg.get("project_id"), normalization=a.normalization, normalized=a.normalized,
                    conditions=run["conditions"], reference_condition=run["reference_condition"], stratified=run["stratified"],
                    samples=run["samples"], tests=run["tests"], anova_calls=anova.get("calls"), anova_warning=anova.get("warning"),
                    all_pairs_calls=pairs, recommendation=run.get("recommendation"),
                    plots=dict(status="success" if plots.returncode == 0 else "failed",
                               message=None if plots.returncode == 0 else (plots.stderr.strip().splitlines() or [""])[-1],
                               notes=(plot_manifest or {}).get("notes", []), responsive_genes=(plot_manifest or {}).get("responsive_genes")),
                    files=dict(anova="tests/anova_genes.tsv", pairs="tests/pairs_genes.tsv", correlation="plots/response_correlation.tsv",
                               heatmap="plots/heatmap_matrix.tsv"))
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2, default=str)+"\n")
    print(json.dumps({k:manifest[k] for k in ("status", "conditions", "anova_calls", "all_pairs_calls")}, indent=1, default=str))


if __name__ == "__main__":main()
