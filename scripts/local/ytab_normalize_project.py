#!/usr/bin/env python3
"""Normalize one project's insertion-site counts with chosen methods, in one step (used by the app's Normalization tab).

Runs the existing, separately validated stages in order and adds plot-ready summaries:
  1. parent site catalogue (`ytab ttr`): every included parent/control sample; library = project::background::pool;
  2. experiment matrix and methods (`ytab normalize-experiment`): samples from config/comparison_design.csv, each
     method applied to the same raw matrix, scaled within each library; unnormalized counts (nonorm) are always
     included as the baseline;
  3. summaries for plotting: per sample and method, the change from its own raw counts (median and IQR of
     log2(normalized / raw) at its insertion sites, sites kept), positive-count distributions and matched-pair diagnostics;
     and each method's native before/after view (native()): TRANSIT-style per-sample statistics, the statistic a scaling
     method equalizes, quantiles against a fitted geometric and ranked counts (BGC, aBGC, unnormalized), and saturation
     curves with the thinning target (MidLC).
Without --method every available method is run, so any of them can be shown. --refresh-summaries RUN rewrites the plot
summaries of an earlier run (no renormalization).

The normalized output for downstream tests (multicondition, genetic-interactions, resampling) is OUTPUT/experiment/methods.

Example:
  ytab normalize-project --project-config output/projects/Zn_toxicity_screen/config/project.yaml \\
      --method betageom --method ttr --output output/projects/Zn_toxicity_screen/normalization/run1
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
from scipy.stats import kurtosis, skew
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.normalization.experiment import evaluate
from ytab.normalization.methods import METHODS
from ytab.normalization.midlc import expected_occupied
from ytab.normalization.ttr import run as build_catalogue, write_table
from ytab.pipeline.project_config import load_project_config
from ytab.screen_coverage.library import is_parent_sample

# Evidence for the suggested method, not an automatic choice: in the pre-registered Hermes benchmarks BGC (then TTR)
# with recentering kept false discoveries lowest (resampling v3: FDR 7.0% vs 12.4% with TTR) and is the normalization
# the multi-condition and genetic-interaction tests were validated with.
SUGGESTED = "betageom"
PERCENTILES = np.arange(1, 100)
QQ_PROBS = np.unique(np.concatenate([np.linspace(.01, .99, 99), 1-np.logspace(-2, -4, 21)]))
SCALING = {"ttr":"density × trimmed nonzero mean", "totreads":"total reads", "nzmean":"nonzero mean", "zinfnb":"fitted NB mean (ZINB)"}
SHAPE = ("nonorm", "betageom", "aBGC")


def included(row):
    return str(row.get("include", "true")).strip().lower() in {"1", "true", "yes"}


def field(row, name):
    return (row.get(name) or row.get("guessed_"+name) or "").strip()


def catalogue_manifest(config_path, output):
    cfg = load_project_config(config_path)
    root = Path(cfg["_repo_root"]); resolve = lambda p:Path(p) if Path(p).is_absolute() else root/p
    with resolve(cfg["sample_sheet"]).open() as handle:rows = [r for r in csv.DictReader(handle) if included(r)]
    parents = [r for r in rows if is_parent_sample(r)]
    if not parents:raise ValueError("No included parent/control samples in the sample sheet")
    observations = [dict(sample=r["sample"], library_id="::".join((cfg["project_id"], field(r, "background"), field(r, "pool")))) for r in parents]
    manifest = dict(version=1, description=f"Parent catalogue for {cfg['project_id']} (app Normalization tab)",
                    sources=[dict(project_config=str(Path(config_path).resolve()), dataset=cfg["project_id"], observations=observations)])
    output.write_text(yaml.safe_dump(manifest, sort_keys=False))
    return cfg, len(parents)


def summaries(experiment, methods):
    raw = sparse.load_npz(experiment/"matrix/raw_counts.npz").tocsc()
    with (experiment/"matrix/observations.tsv").open() as handle:meta = list(csv.DictReader(handle, delimiter="\t"))
    scale, distribution = [], []
    for method in methods:
        normalized = sparse.load_npz(experiment/"methods"/method/"normalized_counts.npz").tocsc()
        with (experiment/"methods"/method/"factors.tsv").open() as handle:factors = list(csv.DictReader(handle, delimiter="\t"))
        for j, m in enumerate(meta):
            x = raw[:, j].toarray().ravel(); y = normalized[:, j].toarray().ravel(); positive = x > 0; both = positive & (y > 0)
            ratio = np.log2(y[both]/x[both]) if both.any() else np.array([np.nan])
            q1, q3 = np.percentile(ratio, [25, 75])
            scale.append(dict(method=method, observation_id=m["observation_id"], display_name=m["display_name"], role=m["role"],
                              library_id=m["library_id"], raw_total_support=int(x.sum()), occupied_sites=int(positive.sum()),
                              factor=factors[j].get("factor", ""), effective_scale=float(np.median(y[positive]/x[positive])),
                              # change from the sample's own raw counts, at sites still occupied after normalization
                              log2_shift_vs_raw=float(np.median(ratio)), log2_spread_vs_raw=float(q3-q1),
                              sites_kept=float(both.sum()/max(positive.sum(), 1))))
            values = np.percentile(y[y > 0], PERCENTILES) if (y > 0).any() else np.full(len(PERCENTILES), np.nan)
            distribution.extend(dict(method=method, observation_id=m["observation_id"], display_name=m["display_name"], role=m["role"],
                                     percentile=int(p), value=float(v)) for p, v in zip(PERCENTILES, values))
    return scale, distribution


def _rows(path):
    with path.open() as handle:return list(csv.DictReader(handle, delimiter="\t"))


def native(experiment, methods, trim):
    """Each method's native before/after view. Returns {table name: rows}."""
    raw = sparse.load_npz(experiment/"matrix/raw_counts.npz").tocsc(); meta = _rows(experiment/"matrix/observations.tsv")
    out = dict(stats=[], scaling=[], qq=[], ranked=[], saturation=[], midlc=[])
    ids = lambda m:dict(observation_id=m["observation_id"], display_name=m["display_name"], role=m["role"])
    for method in methods:
        matrix = sparse.load_npz(experiment/"methods"/method/"normalized_counts.npz").tocsc()
        factors = _rows(experiment/"methods"/method/"factors.tsv")
        for j, m in enumerate(meta):
            y = matrix[:, j].toarray().ravel().astype(float); pos = np.sort(y[y > 0]); x = raw[:, j].toarray().ravel().astype(float)
            # TRANSIT tnseq_stats columns, over every site of the experiment's catalogue
            out["stats"].append(dict(method=method, **ids(m), density=len(pos)/len(y), mean=y.mean(), nzmean=pos.mean() if len(pos) else np.nan,
                                     nzmedian=float(np.median(pos)) if len(pos) else np.nan, max=y.max(), total=y.sum(),
                                     skew=float(skew(y)), kurtosis=float(kurtosis(y))))
            if method in SCALING:   # scaling is exact: after = before x factor
                xp = x[x > 0]; f = float(factors[j]["factor"])
                before = dict(ttr=len(xp)/len(x)*float(np.mean(np.sort(xp)[int(len(xp)*trim):len(xp)-int(len(xp)*trim)])) if len(xp) else np.nan,
                              totreads=x.sum(), nzmean=xp.mean(), zinfnb=float(factors[j].get("zinb_nb_mean") or np.nan))[method]
                out["scaling"].append(dict(method=method, **ids(m), statistic=SCALING[method], before=before, after=before*f, factor=f))
            if method in SHAPE and len(pos):   # quantiles against a geometric fitted to the positive values (TRANSIT QC), ranked values
                p = 1/pos.mean(); expected = np.log1p(-QQ_PROBS)/np.log1p(-min(p, .999999))
                observed = np.quantile(pos, QQ_PROBS)
                out["qq"].extend(dict(method=method, **ids(m), probability=float(q), observed=float(o), expected=float(e))
                                 for q, o, e in zip(QQ_PROBS, observed, expected))
                ranks = np.unique(np.round(np.logspace(0, np.log10(len(pos)), 200)).astype(int))
                out["ranked"].extend(dict(method=method, **ids(m), rank=int(r), value=float(pos[-r])) for r in ranks)
            if method == "midlc":
                counts = x[x > 0]
                for frac in np.unique(np.concatenate([np.logspace(-3, 0, 60), [1.]])):
                    out["saturation"].append(dict(**ids(m), reads=float(frac*counts.sum()), sites=expected_occupied(counts, float(frac))))
                d = factors[j]
                out["midlc"].append(dict(**ids(m), raw_reads=float(counts.sum()), raw_sites=int(len(counts)), midlc_reads=float(d["midlc_reads"]),
                                         midlc_sites=len(counts)/2, saturation_ratio=float(d["saturation_ratio"]), target_ratio=float(d["group_target_ratio"]),
                                         thinned_reads=float(d["thinned_total"]), thinned_sites=int(d["thinned_occupied_sites"]),
                                         sites_kept=int(d["thinned_occupied_sites"])/max(len(counts), 1)))
    return out


def write_native(output, methods, trim):
    for name, rows in native(output/"experiment", methods, trim).items():
        if rows:write_table(output/"plots"/f"native_{name}.tsv", list(rows[0]), rows)


def method_summary(experiment, methods):
    with (experiment/"comparison_diagnostics.tsv").open() as handle:rows = list(csv.DictReader(handle, delimiter="\t"))
    out = {}
    for method in methods:
        for kind in ("control_control", "treated_control"):
            r = [x for x in rows if x["method"] == method and x["kind"] == kind and x["median_log2_ratio_shared"] not in ("", None)]
            if r:out.setdefault(method, {})[kind] = dict(pairs=len(r),
                median_abs_shift=float(np.median([abs(float(x["median_log2_ratio_shared"])) for x in r])),
                median_iqr=float(np.median([float(x["log2_ratio_iqr_shared"]) for x in r])))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--project-config")
    available = [k for k, m in METHODS.items() if m.status == "available" and k != "nonorm"]
    parser.add_argument("--method", action="append", choices=available, help="Repeat; default every available method. Unnormalized counts are always added")
    parser.add_argument("--output", help="New output directory")
    parser.add_argument("--trim", type=float, default=.05); parser.add_argument("--target", type=float, default=100)
    parser.add_argument("--refresh-summaries", metavar="RUN", help="Rewrite the plot summaries of an earlier run and exit")
    args = parser.parse_args()
    if args.refresh_summaries:
        run = Path(args.refresh_summaries).resolve(); manifest = json.loads((run/"manifest.json").read_text())
        scale, distribution = summaries(run/"experiment", manifest["methods"])
        write_table(run/"plots/effective_scale.tsv", list(scale[0]), scale)
        write_table(run/"plots/distributions.tsv", list(distribution[0]), distribution)
        write_native(run, manifest["methods"], float(manifest.get("trim", .05)))
        print(json.dumps(dict(status="success", refreshed=str(run))))
        return
    if not (args.project_config and args.output):parser.error("--project-config and --output are required")
    args.method = args.method or available
    output = Path(args.output).resolve()
    try:
        if output.exists() and any(output.iterdir()):raise ValueError("Output directory is nonempty; choose a new one")
        methods = ["nonorm"]+list(dict.fromkeys(args.method))
        output.mkdir(parents=True, exist_ok=True)
        print("1/3 Building the parent site catalogue", flush=True)
        cfg, parents = catalogue_manifest(Path(args.project_config), output/"catalogue_sources.yaml")
        build_catalogue(output/"catalogue_sources.yaml", output/"catalogue", args.trim, args.target, plots=False)
        print("2/3 Normalizing the experiment samples: "+", ".join(methods), flush=True)
        evaluate(args.project_config, output/"catalogue", output/"experiment", methods, args.trim, args.target, plots=False)
        print("3/3 Writing plot summaries", flush=True)
        scale, distribution = summaries(output/"experiment", methods)
        (output/"plots").mkdir()
        write_table(output/"plots/effective_scale.tsv", list(scale[0]), scale)
        write_table(output/"plots/distributions.tsv", list(distribution[0]), distribution)
        write_native(output, methods, args.trim)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f"ERROR: {error}\n")
    manifest = dict(status="success", project=cfg["project_id"], methods=methods, suggested=SUGGESTED, parents_in_catalogue=parents,
                    samples=len({r["observation_id"] for r in scale}), trim=args.trim, target=args.target,
                    normalized=str(output/"experiment/methods"), catalogue=str(output/"catalogue"),
                    method_summary=method_summary(output/"experiment", methods),
                    files=dict(effective_scale="plots/effective_scale.tsv", distributions="plots/distributions.tsv",
                               comparisons="experiment/comparison_diagnostics.tsv", factors="experiment/methods/method_comparison.tsv"),
                    reading="Matched controls from the same pool should differ little after normalization: a small median shift and a narrow spread of "
                            "log2 ratios at shared sites. Units differ between methods; compare relative results, not absolute values.")
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps({k:manifest[k] for k in ("status", "project", "methods", "samples", "normalized")}, indent=2))


if __name__ == "__main__":main()
