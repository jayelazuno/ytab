#!/usr/bin/env python3
"""Method-comparison study, top-site test: does removing (or capping) each gene's largest insertion site bring our calls
closer to the Cunningham lab's published z-scores? Gale 2020 (BG14) and Nickels 2025 (BG2, pdr1_KO), the in-scope strains.

z-score arms (per sample, per gene; the CDS sites with the usual N/C-terminal trimming):
  all       every site summed (what the app does today)
  drop_top  the gene's largest site removed (Nickels 2025 Methods)
  winsor    the gene's largest site lowered to its second largest (TRANSIT -winz idea, per sample)
each on CPM (raw reads per million) and on BGC (normalized site counts), and for drop_top / winsor either
  renorm=no   scaled by the library totals before the change (normalize, then drop)
  renorm=yes  rescaled so each library's gene total is what it was (drop, then normalize)
Resampling arms: BGC, with and without --winsorize (label-invariant: largest site x sample of a gene -> second largest).

Output: output/comparison/topsite/{runs,summary}. Untracked study.
  python scripts/local/ytab_method_comparison_topsite.py [--skip-resampling]
"""
from __future__ import annotations
import argparse, gzip, json, runpy, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd, yaml
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT/"src"))
from ytab.fitness import agreement, zscore  # noqa: E402
from ytab.essentiality.control_conditioned_gumbel import annotations, gene_sites  # noqa: E402
from ytab.essentiality.library_correlation import SITE_RULES, gene_key  # noqa: E402

OUT = ROOT/"output/comparison/topsite"; CMP = ROOT/"output/comparison"
SCOPE = {"Gale_2020": None, "Nickels_2025_glabrata": {"BG2", "pdr1_KO"}}


def study(dataset):
    project = ROOT/"output/projects"/dataset; design = project/"config/comparison_design.csv"
    rows = pd.read_csv(design, dtype=str); rows = rows[rows.include.str.lower().isin(["1", "true", "yes", "y"])]
    keep = SCOPE[dataset]; rows = rows if keep is None else rows[rows.background.isin(keep)]
    cfg = yaml.safe_load(open(project/"config/project.yaml")); gff = Path(cfg["reference"]["gff"]); gff = gff if gff.is_absolute() else ROOT/gff
    resampling = json.load(open(CMP/dataset/"runs/resampling/manifest.json"))   # the same normalization run as the main study
    norm = Path(resampling["normalization_run"]); norm = norm if norm.is_absolute() else ROOT/norm
    return dict(dataset=dataset, design=design, comparisons=list(rows.comparison_id), gff=gff, norm=norm, config=project/"config/project.yaml",
                published=CMP/dataset/"published.tsv")


def site_index(norm, gff):
    with gzip.open(norm/"experiment/matrix/sites.tsv.gz", "rt") as h:sites = pd.read_csv(h, sep="\t", dtype={"chromosome": str})
    genes = annotations(Path(gff)); idx = gene_sites(genes, list(zip(sites.chromosome, sites.position_1based.astype(int))), SITE_RULES)
    return [gene_key(g["id"]) for g in genes], idx


def columns(norm, path, samples):
    obs = pd.read_csv(norm/"experiment/methods/observations.tsv", sep="\t", dtype=str); m = sparse.load_npz(path).tocsc()
    return {s: np.asarray(m[:, int(np.flatnonzero(obs["sample"].to_numpy() == s)[0])].todense()).ravel() for s in samples}


def gene_values(values, keys, idx, how):
    """gene x sample sums with each gene's sites summed whole, without the largest, or with the largest capped."""
    out = {}
    for s, v in values.items():
        col = np.zeros(len(keys))
        for k, i in enumerate(idx):
            if not len(i):continue
            x = v[i]; t = x.sum()
            if how != "all" and len(x):
                top2 = np.partition(x, -2)[-2:] if len(x) >= 2 else np.array([0., x.max()])
                t -= top2[-1] if how == "drop_top" else (top2[-1]-top2[-2])
            col[k] = t
        out[s] = col
    return pd.DataFrame(out, index=keys).groupby(level=0).sum().rename_axis("feature_id")


def zscore_arms(p):
    contrasts, sample_map = zscore.load_design(p["design"], p["comparisons"]); samples = list(dict.fromkeys(sample_map["sample"]))
    keys, idx = site_index(p["norm"], p["gff"]); annot = zscore.gene_annotation_from_gff(p["gff"])
    raw = columns(p["norm"], p["norm"]/"experiment/matrix/raw_counts.npz", samples)
    bgc = columns(p["norm"], p["norm"]/"experiment/methods/betageom/normalized_counts.npz", samples)
    raw_all = gene_values(raw, keys, idx, "all")
    for unit, vals in (("cpm", raw), ("bgc", bgc)):
        base = gene_values(vals, keys, idx, "all")
        scale = (1e6/base.sum(axis=0)) if unit == "cpm" else 1.
        for how in ("all", "drop_top", "winsor"):
            g = base if how == "all" else gene_values(vals, keys, idx, how)
            for renorm in (("-",) if how == "all" else ("no", "yes")):
                v = g*scale if renorm in ("-", "no") else g*(base.sum(axis=0)/g.sum(axis=0))*scale
                arm = f"zscore_{unit}_{how}" + ("" if renorm == "-" else f"_renorm-{renorm}"); out = OUT/"runs"/p["dataset"]/arm
                if (out/"zscore_manifest.json").exists():continue
                zscore.run(p["design"], out, normalized=dict(genes=v, raw=raw_all, annotation=annot, label=f"{unit} {how} {renorm}"),
                           comparisons=p["comparisons"], project_id=p["dataset"], analysis_id=f"topsite_{arm}")
                print(p["dataset"], arm, "done", flush=True)


def resampling_arms(p):
    real = subprocess.run
    for winsor in (False, True):
        out = OUT/"runs"/p["dataset"]/("resampling_bgc_winsor" if winsor else "resampling_bgc")
        if (out/"manifest.json").exists():continue
        def patched(cmd, *a, **k):
            cmd = [str(x) for x in cmd]
            if winsor and any(x.endswith("ytab_resampling.py") for x in cmd):cmd.append("--winsorize")
            return real(cmd, *a, **k)
        subprocess.run = patched
        sys.argv = ["ytab_treated_vs_control.py", "--project-config", str(p["config"]), "--method", "resampling", "--normalization-run", str(p["norm"]),
                    "--normalization", "betageom", "--comparisons", ",".join(p["comparisons"]), "--output", str(out)]
        try:runpy.run_path(str(ROOT/"scripts/local/ytab_treated_vs_control.py"), run_name="__main__")
        except SystemExit:pass
        finally:subprocess.run = real
        print(p["dataset"], out.name, "done", flush=True)


def evaluate():
    rows = []
    for dataset in SCOPE:
        pub = agreement.read_published(CMP/dataset/"published.tsv"); keep = SCOPE[dataset]
        for arm in sorted(d for d in (OUT/"runs"/dataset).iterdir() if d.is_dir()):
            method = "zscore" if arm.name.startswith("zscore") else "resampling"
            if not ((arm/"zscore_manifest.json").exists() or (arm/"manifest.json").exists()):continue
            ours = agreement.read_method(arm, method)
            for (b, c), x in ours.groupby(["background", "condition"]):
                if keep is not None and b not in keep:continue
                y = pub[(pub.background == b) & (pub.condition == c)] if "background" in pub else pub
                m = x.merge(y, on="ORF", suffixes=("", "_pub"))
                m = m[m.tested]
                if m.empty:continue
                r = dict(dataset=dataset, background=b, condition=c, arm=arm.name, genes=len(m))
                for d in ("depleted", "enriched"):
                    truth, hit = m.call_pub.eq(d), m.call.eq(d)
                    r[f"{d}_published"] = int(truth.sum()); r[f"{d}_called"] = int(hit.sum()); r[f"{d}_recovered"] = int((truth & hit).sum())
                    r[f"{d}_recall"] = (truth & hit).sum()/truth.sum() if truth.sum() else np.nan
                ok = np.isfinite(m.score) & np.isfinite(m.score_pub)
                r["spearman_score"] = m.score[ok].rank().corr(m.score_pub[ok].rank())
                if method == "zscore":r["pearson_z"] = np.corrcoef(m.score[ok], m.score_pub[ok])[0, 1]
                ok = np.isfinite(m.effect) & np.isfinite(m.effect_pub)
                r["spearman_log2fc"] = m.effect[ok].rank().corr(m.effect_pub[ok].rank()) if ok.sum() > 2 else np.nan
                rows.append(r)
    t = pd.DataFrame(rows); (OUT/"summary").mkdir(parents=True, exist_ok=True); t.to_csv(OUT/"summary/per_condition.tsv", sep="\t", index=False)
    agg = t.groupby(["dataset", "arm"]).agg(conditions=("condition", "size"), depleted_published=("depleted_published", "sum"),
                                            depleted_recovered=("depleted_recovered", "sum"), depleted_called=("depleted_called", "sum"),
                                            enriched_published=("enriched_published", "sum"), enriched_recovered=("enriched_recovered", "sum"),
                                            enriched_called=("enriched_called", "sum"), median_spearman_z=("spearman_score", "median"),
                                            median_pearson_z=("pearson_z", "median"), median_spearman_log2fc=("spearman_log2fc", "median")).reset_index()
    agg["depleted_recall"] = agg.depleted_recovered/agg.depleted_published; agg["enriched_recall"] = agg.enriched_recovered/agg.enriched_published
    agg.to_csv(OUT/"summary/by_arm.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):print(agg.round(3).to_string(index=False))


if __name__ == "__main__":
    a = argparse.ArgumentParser(); a.add_argument("--skip-resampling", action="store_true"); a.add_argument("--evaluate-only", action="store_true"); a = a.parse_args()
    if not a.evaluate_only:
        for d in SCOPE:
            p = study(d); print(d, p["norm"].name, len(p["comparisons"]), "comparisons", flush=True); zscore_arms(p)
            if not a.skip_resampling:resampling_arms(p)
    evaluate()
