#!/usr/bin/env python3
"""Essentiality for one project in one step (used by the app's Essentiality tab): Hermes window HMM, then the model table.

  1. HMM (label-free, any yeast; `ytab hermes-hmm`) on independent parent libraries of one strain, or reuse an
     earlier HMM fit (--hmm-run);
  2. model table (`ytab_essentiality_model_menu.py`): HMM calls, plus the current random forest's scores when its
     prediction table is given (--rf), plus the hybrid (random forest on gene features + HMM posteriors) when
     curated labels are also given (--labels);
  3. plot summaries: window posteriors and insertion counts along each chromosome (plots/windows.tsv.gz).

  ytab essentiality --project-config output/projects/Gale_2020/config/project.yaml \\
      --libraries SRR11561508 SRR11561507 SRR11561506 --output output/projects/Gale_2020/essentiality/run1
  ... --labels resources/species/glabrata/essentiality_labels/labels.tsv            (adds the hybrid)

With --labels and no --rf, the hybrid's gene features are built first (`ytab hybrid-inputs`) from the same parent
libraries as the HMM, normalized with --normalization (default betageom; any method, or midlc). --rf instead uses an
existing classifier prediction table.
"""
import argparse
import gzip
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import yaml

from ytab_hybrid_inputs import build as build_hybrid_inputs

ROOT_SRC = Path(__file__).resolve().parents[2]/"src"
if str(ROOT_SRC) not in sys.path:sys.path.insert(0, str(ROOT_SRC))
from ytab.essentiality import contract  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT/"scripts/local"


def step(arguments):
    result = subprocess.run([sys.executable]+[str(a) for a in arguments], cwd=ROOT, capture_output=True, text=True)
    sys.stderr.write(result.stderr)
    if result.returncode:raise ValueError((result.stderr.strip().splitlines() or ["step failed"])[-1].removeprefix("ERROR: "))
    return result.stdout


def export_windows(hmm_dir, path):
    """Per window: the HMM state posteriors, the most probable state (unmappable windows: "unmappable"), insertions."""
    fit = json.loads((hmm_dir/"fit.json").read_text()); window = fit["settings"]["window"]
    z = np.load(hmm_dir/"windows.npz")
    chromosomes = sorted({k.rsplit("_", 1)[0] for k in z.files})
    with gzip.open(path, "wt") as out:
        out.write("chromosome\tstart\tP_ES\tP_GD\tP_NE\tP_GA\tstate\tinsertions\tmappable\n")
        for c in chromosomes:
            gamma, counts, keep = z[f"{c}_gamma"], z[f"{c}_Y"].sum(1), z[f"{c}_keep"]
            for i in range(len(gamma)):
                state = contract.STATES[int(np.argmax(gamma[i]))] if keep[i] else "unmappable"
                out.write(f"{c}\t{i*window+1}\t" + "\t".join(f"{v:.4g}" for v in gamma[i]) + f"\t{state}\t{int(counts[i])}\t{int(keep[i])}\n")
    return len(chromosomes)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", required=True, type=Path)
    p.add_argument("--libraries", nargs="+", help="Independent parent libraries of one strain (sample IDs)")
    p.add_argument("--hmm-run", type=Path, help="Reuse an earlier HMM fit (a folder with genes.tsv, fit.json, windows.npz)")
    p.add_argument("--features", type=Path, help="Gene features for the hybrid: a combined SummaryTable feature table (default: built from the HMM's parents with --normalization)")
    p.add_argument("--rf", type=Path, help="Legacy classifier predictions (essentiality_predictions.<target>.csv): adds its score for comparison")
    p.add_argument("--with-classifier", action="store_true", help="Also score the built gene features with the legacy classifier (comparison only)")
    p.add_argument("--labels", type=Path, help="Curated labels (gene_id or cglab_id, label, tier), e.g. resources/species/<species>/essentiality_labels/labels.tsv: adds the hybrid")
    p.add_argument("--normalization", default="betageom", help="Normalization for the hybrid's gene features when --features is not given (method id, e.g. betageom, ttr or midlc)")
    p.add_argument("--normalization-run", type=Path, help="Site Normalization run for --normalization (default: newest with the parents)")
    p.add_argument("--train-tiers", default="silver")
    p.add_argument("--window", type=int, default=100); p.add_argument("--min-reads", type=int, default=1)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    out = a.output.resolve()
    try:
        if out.exists() and any(out.iterdir()):raise ValueError("Output directory is nonempty; choose a new one")
        if not (a.hmm_run or a.libraries):raise ValueError("Give --libraries (to fit the HMM) or --hmm-run (to reuse a fit)")
        project = a.project_config.resolve().parents[1].name
        cfg = yaml.safe_load(a.project_config.read_text())
        out.mkdir(parents=True)
        if a.hmm_run:
            hmm_dir = a.hmm_run.resolve()
            if not (hmm_dir/"genes.tsv").exists():raise ValueError(f"No HMM result in {hmm_dir}")
        else:
            print("1/3 Fitting the Hermes window HMM (several minutes)", file=sys.stderr, flush=True)
            step([SCRIPTS/"ytab_hermes_hmm.py", out/"hmm", "--project", project, "--libraries", *a.libraries,
                  "--window", a.window, "--min-reads", a.min_reads])
            hmm_dir = out/"hmm/joint"
        hybrid_inputs = None
        if a.labels and not a.features and not (a.rf and not a.with_classifier):
            parents = json.loads((hmm_dir/"fit.json").read_text())["samples"]
            print(f"Building the hybrid's gene features ({a.normalization}, parents {', '.join(parents)})", file=sys.stderr, flush=True)
            hybrid_inputs = build_hybrid_inputs(a.project_config, a.normalization, parents, a.normalization_run, with_classifier=a.with_classifier)
            a.features = Path(hybrid_inputs["features"])
            if a.with_classifier:a.rf = Path(hybrid_inputs["predictions"])
        print("2/3 Building the model table", file=sys.stderr, flush=True)
        step([SCRIPTS/"ytab_essentiality_model_menu.py", "--hmm", hmm_dir/"genes.tsv", "--out", out/"models", "--train-tiers", a.train_tiers,
              *(["--features", a.features] if a.features else []), *(["--rf", a.rf] if a.rf else []),
              *(["--labels", a.labels] if a.labels else [])])
        print("3/3 Writing plot summaries", file=sys.stderr, flush=True)
        (out/"plots").mkdir(); chromosomes = export_windows(hmm_dir, out/"plots/windows.tsv.gz")
    except (ValueError, KeyError, OSError) as error:
        p.exit(1, f"ERROR: {error}\n")
    fit = json.loads((hmm_dir/"fit.json").read_text()); models = json.loads((out/"models/summary.json").read_text())
    manifest = dict(status="success", contract=contract.VERSION, calls=models.get("hmm_calls"), project=cfg.get("project_id", project), species=cfg.get("species") or cfg.get("reference", {}).get("species"),
                    libraries=fit["samples"], hmm_run=str(hmm_dir), hmm_reused=bool(a.hmm_run), models=models["models"],
                    hybrid="hybrid" in models["models"], features=str(a.features) if a.features else None, rf=str(a.rf) if a.rf else None,
                    labels=str(a.labels) if a.labels else None,
                    hybrid_normalization=(hybrid_inputs or {}).get("source") if hybrid_inputs else ({"method":"given --features"} if a.features else {"method":"given --rf"} if a.rf else None),
                    hmm=dict(gene_calls=fit["gene_calls"], masked_fraction=fit["masked_fraction"], converged=fit["converged"],
                             iterations=fit["iterations"], notes=fit.get("notes", []), window=fit["settings"]["window"]),
                    agreement=models.get("agreement"), hybrid_training=models.get("hybrid_training"), chromosomes=chromosomes,
                    files=dict(genes="models/genes.tsv", windows="plots/windows.tsv.gz", guide="models/MODELS.md"))
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2, default=str)+"\n")
    problems = contract.validate(out)
    if problems:p.exit(1, "ERROR: the run does not meet the output contract: " + "; ".join(problems) + "\n")
    print(json.dumps({k:manifest[k] for k in ("status", "project", "libraries", "models", "hmm")}, indent=1, default=str))


if __name__ == "__main__":main()
