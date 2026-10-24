#!/usr/bin/env python3
"""Essentiality for one project in one step (used by the app's Essentiality tab): Hermes window HMM, then the model table.

  1. HMM (label-free, any yeast; `ytab hermes-hmm`) on independent parent libraries of one strain, or reuse an
     earlier HMM fit (--hmm-run);
  2. model table (`ytab_essentiality_model_menu.py`): HMM calls, plus the current random forest's scores when its
     prediction table is given (--rf), plus the hybrid (random forest on gene features + HMM posteriors) when
     curated labels are also given (--labels);
  3. plot summaries: window posteriors and insertion counts along each chromosome (plots/windows.tsv.gz); the
     insertion-target preference (plots/target_model.tsv); with --labels, ROC curves of each model against the curated
     labels (plots/roc.tsv). --refresh-plots RUN writes the last two for an earlier run.

  ytab essentiality --project-config output/projects/Gale_2020/config/project.yaml \\
      --libraries SRR11561508 SRR11561507 SRR11561506 --output output/projects/Gale_2020/essentiality/run1
  ... --labels resources/species/glabrata/essentiality_labels/labels.tsv            (adds the hybrid)

With --labels and no --rf, the hybrid's gene features are built first (`ytab hybrid-inputs`) from the same parent
libraries as the HMM, normalized with --normalization (default betageom; any method, or midlc). --rf instead uses an
existing classifier prediction table.
"""
import argparse
import csv
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


def export_target_model(hmm_dir, path, project, libraries, min_reads):
    """The insertion-target preference (position weight matrix, information per offset): the model the HMM fitted with
    when it learned one (sequence covariate "learned"); otherwise learned the same way from the same libraries' distinct
    sites as a diagnostic only (the fit used the Hermes motif). Returns True / False (used in the fit), None if unavailable."""
    import ytab_hermes_hmm as hmm
    if (hmm_dir/"target_model.tsv").exists():
        rows = list(csv.DictReader(open(hmm_dir/"target_model.tsv"), delimiter="\t")); used = True
    else:
        try:
            fasta, _ = hmm.project_reference(project)
            model = hmm.learn_target(libraries, [None]*len(libraries), hmm.hh.Settings(min_reads=min_reads), project, fasta)
        except (ValueError, OSError, KeyError):return None
        hmm.write_target_model(path, model); rows = list(csv.DictReader(open(path), delimiter="\t")); used = False
    write_rows(path, [dict(r, used_in_fit=used) for r in rows]); return used


def roc_curve(score, label):
    """False and true positive rates at each distinct score (ties form one step) and the area under the curve."""
    order = np.argsort(-score, kind="mergesort"); s, y = score[order], label[order]
    last = np.r_[np.flatnonzero(np.diff(s)), len(s)-1]
    tp, fp = np.cumsum(y)[last], np.cumsum(1-y)[last]
    fpr, tpr = np.r_[0, fp/fp[-1]], np.r_[0, tp/tp[-1]]
    return fpr, tpr, float(np.sum(np.diff(fpr)*(tpr[1:]+tpr[:-1])/2))


def export_roc(genes_path, path):
    """ROC of each model's score against the curated labels (reference_label), for each label tier and all together.
    Only genes the model was not fitted to: the HMM uses no labels, and the hybrid's training genes carry out-of-fold
    scores. Returns {model: {tier: AUC}}, or {} without labels."""
    genes = list(csv.DictReader(open(genes_path), delimiter="\t"))
    labelled = [g for g in genes if g.get("reference_label") in ("essential", "nonessential")]
    models = {"HMM (label-free)":"hmm_score", "Hybrid":"hybrid_score", "Random forest":"rf_existing_score"}
    rows, aucs = [], {}
    for model, column in models.items():
        for tier in ["all"]+sorted({g.get("reference_tier") or "" for g in labelled} - {""}):
            d = [g for g in labelled if (tier == "all" or g.get("reference_tier") == tier) and (g.get(column) or "") not in ("", "nan")]
            y = np.array([g["reference_label"] == "essential" for g in d], dtype=float)
            if len(d) < 10 or y.min() == y.max():continue
            fpr, tpr, auc = roc_curve(np.array([float(g[column]) for g in d]), y)
            keep = np.unique(np.r_[0, np.round(np.linspace(0, len(fpr)-1, min(len(fpr), 400))).astype(int), len(fpr)-1])
            rows += [dict(model=model, tier=tier, genes=len(d), essential=int(y.sum()), auc=round(auc, 4), fpr=round(float(fpr[i]), 5), tpr=round(float(tpr[i]), 5))
                     for i in keep]
            aucs.setdefault(model, {})[tier] = round(auc, 4)
    if rows:write_rows(path, rows)
    return aucs


def write_rows(path, rows):
    with open(path, "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


def plot_summaries(out, hmm_dir, project, min_reads):
    fit = json.loads((hmm_dir/"fit.json").read_text())
    used = export_target_model(hmm_dir, out/"plots/target_model.tsv", project, fit["samples"], min_reads)
    return dict(target_model=None if used is None else dict(file="plots/target_model.tsv", used_in_fit=used),
                evaluation=export_roc(out/"models/genes.tsv", out/"plots/roc.tsv") or None)


def auto_sequence_covariate(cfg, libraries):
    """learned (v3) when any library's sample-sheet transposon is not Hermes; motif (v2) otherwise, as before."""
    sheet = ROOT/cfg["sample_sheet"] if cfg.get("sample_sheet") else None
    if not sheet or not sheet.exists():return "motif"
    with open(sheet, newline="") as handle:
        transposons = {r["sample"]:(r.get("transposon") or "").strip().lower() for r in csv.DictReader(handle)}
    return "learned" if any(transposons.get(lib, "") not in ("", "hermes") for lib in libraries) else "motif"


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
    p.add_argument("--sequence-covariate", choices=("auto", "motif", "learned"), default="auto",
                   help="HMM sequence covariate: motif (Hermes nTnnnnAn, v2), learned (target preference from the libraries' own sites, v3), "
                        "auto (learned when a library's sample-sheet transposon is not Hermes, else motif)")
    p.add_argument("--output", type=Path)
    p.add_argument("--refresh-plots", metavar="RUN", type=Path, help="Write the target-preference and ROC summaries of an earlier run and exit")
    a = p.parse_args()
    if a.refresh_plots:
        run = a.refresh_plots.resolve(); m = json.loads((run/"manifest.json").read_text())
        project = a.project_config.resolve().parents[1].name
        extra = plot_summaries(run, Path(m["hmm_run"]), project, int(m.get("hmm", {}).get("settings", {}).get("min_reads", 1)))
        print(json.dumps(dict(status="success", refreshed=str(run), **extra), default=str)); return
    if not a.output:p.error("--output is required")
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
            covariate = a.sequence_covariate if a.sequence_covariate != "auto" else auto_sequence_covariate(cfg, a.libraries)
            print(f"   sequence covariate: {covariate}", file=sys.stderr, flush=True)
            step([SCRIPTS/"ytab_hermes_hmm.py", out/"hmm", "--project", project, "--libraries", *a.libraries,
                  "--window", a.window, "--min-reads", a.min_reads, "--sequence-covariate", covariate])
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
        extra = plot_summaries(out, hmm_dir, project, a.min_reads)
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
                    target_model=extra["target_model"], evaluation=extra["evaluation"],
                    files=dict(genes="models/genes.tsv", windows="plots/windows.tsv.gz", guide="models/MODELS.md", roc="plots/roc.tsv",
                               target_model="plots/target_model.tsv"))
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2, default=str)+"\n")
    problems = contract.validate(out)
    if problems:p.exit(1, "ERROR: the run does not meet the output contract: " + "; ".join(problems) + "\n")
    print(json.dumps({k:manifest[k] for k in ("status", "project", "libraries", "models", "hmm")}, indent=1, default=str))


if __name__ == "__main__":main()
