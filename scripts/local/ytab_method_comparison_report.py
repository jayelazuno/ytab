#!/usr/bin/env python3
"""Method-comparison report for the Cunningham lab, step 2: tables and figures from the grid of
ytab_method_comparison_report_runs.py, into output/comparison/report/{figures,tables}/. The report text
(report.md) is written by hand from these; build_docx() renders it as report.docx.

Published calls: the lab's per-gene z-scores, |z| >= 3 (Gale 2020; Nickels 2025 Table S1), one table per dataset in
output/comparison/<dataset>/published.tsv. Recovery = of the published hits the method tested, the share it calls in
the same direction (the Gumbel calls depletion only). Agreement = Spearman correlation of the method's signed score with
the published z over the genes both tested, per condition, median over conditions (z for the z-score, sign x -log10 p
for resampling, -(-log10 p) for the Gumbel). Call sets of Nickels conditions are pooled as gene x condition pairs.

  python scripts/local/ytab_method_comparison_report.py [--docx]
"""
from __future__ import annotations
import argparse, importlib.util, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT/"src"))
from ytab.fitness import agreement as ag  # noqa: E402

REPORT = ROOT/"output/comparison/report"; FIG = REPORT/"figures"; TAB = REPORT/"tables"
DATASETS = {"Gale_2020": ("Gale et al. 2020", None), "Nickels_2025_glabrata": ("Nickels et al. 2025", {"BG2", "pdr1_KO"})}
FILTERS = {"none": "All sites", "winsorize": "Winsorized", "drop-top": "Top site dropped"}
METHODS = {"zscore": "z-score", "zscorecmp": "z-score (comparison noise)", "resampling": "Resampling", "gumbel": "Gumbel"}
NORMS = {"zscore": ["nonorm", "totreads", "nzmean", "ttr", "midlc", "quantile", "zinfnb", "betageom", "aBGC"],
         "resampling": ["totreads", "nzmean", "ttr", "betageom", "aBGC"], "gumbel": ["totreads", "nzmean", "ttr"]}
NORMS["zscorecmp"] = NORMS["zscore"]
NORM_LABEL = dict(nonorm="Unnormalized", totreads="Total reads", nzmean="NZ mean", ttr="TTR", midlc="MidLC", quantile="Quantile",
                  zinfnb="ZINB", betageom="BGC", aBGC="aBGC", cpm="CPM (as published)")
COLOURS = {"zscore": "#0072B2", "zscorecmp": "#56B4E9", "resampling": "#E69F00", "gumbel": "#009E73", "published": "#CC79A7"}
DEFAULT = {"zscore": "betageom", "zscorecmp": "betageom", "resampling": "betageom", "gumbel": "ttr"}
MIN_TESTED = 100   # published genes a method must test in a dataset to be shown
COLLAPSE = .25   # share of tested genes with control reads and no treated reads above which a condition is collapsed (Nickels: 0.42-0.70 vs <= 0.20)

_spec = importlib.util.spec_from_file_location("summary", ROOT/"scripts/local/ytab_method_comparison_summary.py")
summary = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(summary)
style, venn, false_calls = summary.style, summary.venn, summary.false_calls


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True); fig.tight_layout()
    for ext in ("png", "pdf"):fig.savefig(FIG/f"{name}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def published(dataset):
    p = ag.read_published(ROOT/"output/comparison"/dataset/"published.tsv"); keep = DATASETS[dataset][1]
    return p if keep is None or "background" not in p else p[p.background.isin(keep)]


def result(dataset, how, method, norm):
    path = REPORT/dataset/how/f"{method}_{norm}"
    zscore = method.startswith("zscore")   # zscorecmp: the z-score with noise from each comparison
    done = (path/"zscore_manifest.json") if zscore else (path/"manifest.json")
    if not done.exists():return None
    try:x = ag.read_method(path, "zscore" if zscore else method)
    except (OSError, ValueError, KeyError):return None
    keep = DATASETS[dataset][1]
    return x if keep is None else x[x.background.isin(keep)]


def evaluate():
    rows, merged = [], {}
    for dataset in DATASETS:
        pub = published(dataset)
        for how in FILTERS:
            for method, norms in NORMS.items():
                for norm in norms:
                    x = result(dataset, how, method, norm)
                    if x is None or x.empty:continue
                    m = x[x.tested].merge(pub, on=["ORF", "background", "condition"], suffixes=("", "_pub"))
                    merged[(dataset, how, method, norm)] = m
                    r = dict(dataset=dataset, filter=how, method=method, normalization=norm, conditions=m.groupby(["background", "condition"]).ngroups,
                             genes_tested=len(m))
                    for d in ("depleted", "enriched"):
                        if d == "enriched" and method in ag.DEPLETION_ONLY:continue
                        truth, hit = m.call_pub.eq(d), m.call.eq(d)
                        r |= {f"{d}_published":int(truth.sum()), f"{d}_recovered":int((truth & hit).sum()), f"{d}_called":int(hit.sum()),
                              f"{d}_recall":(truth & hit).sum()/truth.sum() if truth.sum() else np.nan,
                              f"{d}_called_per_condition":hit.sum()/max(r["conditions"], 1)}
                    per = [g.score[ok].rank().corr(g.score_pub[ok].rank()) for _, g in m.groupby(["background", "condition"])
                           for ok in [np.isfinite(g.score) & np.isfinite(g.score_pub)] if ok.sum() > 2]
                    r["spearman_vs_published_z"] = float(np.median(per)) if per else np.nan   # median over conditions
                    rows.append(r)
    t = pd.DataFrame(rows); TAB.mkdir(parents=True, exist_ok=True); t.to_csv(TAB/"recovery_all.tsv", sep="\t", index=False)
    return t, merged


def applicable(t):
    """Drop method x dataset cells where the method tested almost none of the published genes (the Gumbel on Nickels:
    one treated library for most doses, and too few reads per site for a run of empty sites to be detectable). Shown
    as 0% recovery they would read as a failure; the report says why they are absent."""
    tested = t.groupby(["dataset", "method"]).genes_tested.transform("max")
    return t[tested >= MIN_TESTED]


def heatmaps(t, dataset):
    """Rows: method x normalization; columns: top-site treatment. Depleted recovery (with calls per condition), enriched
    recovery, agreement with the published z."""
    d = t[t.dataset == dataset]
    if d.empty:return
    labels = [(m, n) for m in METHODS for n in NORMS[m] if ((d.method == m) & (d.normalization == n)).any()]
    panels = [("depleted_recall", "Published depleted hits recovered", True), ("enriched_recall", "Published enriched hits recovered", True),
              ("spearman_vs_published_z", "Agreement with the published z\n(Spearman, median over conditions)", False)]
    fig, axes = plt.subplots(1, 3, figsize=(16, .42*len(labels)+2.2))
    for ax, (col, title, pct) in zip(axes, panels):
        grid = np.full((len(labels), len(FILTERS)), np.nan); text = [[""]*len(FILTERS) for _ in labels]
        for i, (m, n) in enumerate(labels):
            for j, f in enumerate(FILTERS):
                r = d[(d.method == m) & (d.normalization == n) & (d["filter"] == f)]
                if r.empty or col not in r or not np.isfinite(r[col].iloc[0]):continue
                v = float(r[col].iloc[0]); grid[i, j] = v
                if pct:
                    calls = r[col.replace("recall", "called_per_condition")].iloc[0]
                    text[i][j] = f"{100*v:.0f}%  ({calls:,.0f})"
                else:text[i][j] = f"{v:.2f}"
        ax.imshow(grid, cmap="Blues", vmin=0 if pct else .3, vmax=1, aspect="auto")
        for i in range(len(labels)):
            for j in range(len(FILTERS)):
                if text[i][j]:ax.text(j, i, text[i][j], ha="center", va="center", fontsize=8.5, fontweight="bold",
                                      color="white" if grid[i, j] > (.75 if pct else .8) else "black")
        ax.set_xticks(range(len(FILTERS))); ax.set_xticklabels(FILTERS.values(), fontweight="bold")
        ax.set_yticks(range(len(labels))); ax.set_yticklabels([f"{METHODS[m]} · {NORM_LABEL[n]}" for m, n in labels] if ax is axes[0] else [], fontweight="bold")
        for k in range(1, len(labels)):
            if labels[k][0] != labels[k-1][0]:ax.axhline(k-.5, color="black", lw=1.5)
        ax.set_title(title, fontweight="bold"); ax.tick_params(length=0)
    fig.text(.01, -.01, "Recovery panels: share of the published hits called in the same direction (in brackets: the method's calls in that direction per condition). Gumbel: depletion only.",
             fontsize=9, ha="left", va="top")
    fig.suptitle(f"{DATASETS[dataset][0]}: every method × normalization × top-site treatment against the published calls", fontweight="bold")
    save(fig, f"{dataset}_1_recovery_grid")


def collapsed(dataset):
    """Conditions whose treated library lost a large share of the genome: genes with control reads but no treated reads
    (resampling, default normalization, all sites), as a share of the genes tested. Above COLLAPSE the treated library
    went through a bottleneck that no test can separate from the drug effect without replicates."""
    path = REPORT/dataset/"none"/f"resampling_{DEFAULT['resampling']}"/"genes.tsv"
    if not path.exists():return pd.DataFrame(columns=["background", "condition", "genes_tested", "genes_lost", "share_lost", "collapsed"])
    g = pd.read_csv(path, sep="\t", usecols=["background", "condition", "call", "sum_control", "sum_treatment"])
    g = g[g.call.ne("not tested")]
    t = g.assign(lost=(g.sum_control > 0) & (g.sum_treatment == 0)).groupby(["background", "condition"]).agg(genes_tested=("lost", "size"), genes_lost=("lost", "sum")).reset_index()
    t["share_lost"] = t.genes_lost/t.genes_tested; t["collapsed"] = t.share_lost > COLLAPSE
    keep = DATASETS[dataset][1]; t = t if keep is None else t[t.background.isin(keep)]
    TAB.mkdir(parents=True, exist_ok=True); t.sort_values("share_lost", ascending=False).to_csv(TAB/f"{dataset}_collapsed_conditions.tsv", sep="\t", index=False)
    return t


def venns(merged, dataset):
    """Default normalization, all sites: depleted calls of each method and the published hits (gene x condition). With
    several conditions, a second row leaves out the collapsed ones (collapsed())."""
    key = lambda y:set(y.background+"|"+y.condition+"|"+y.ORF)
    shown = [m for m in METHODS if (dataset, "none", m, DEFAULT[m]) in merged]
    if "zscore" not in shown:return
    gone = collapsed(dataset); gone = set(zip(gone.background[gone.collapsed], gone.condition[gone.collapsed]))
    z = merged[(dataset, "none", "zscore", DEFAULT["zscore"])]; gone &= set(zip(z.background, z.condition))   # compared conditions only
    rows = [("All conditions", lambda y:y)]
    if gone:rows.append((f"Without the {len(gone)} collapsed conditions", lambda y:y[[(b, c) not in gone for b, c in zip(y.background, y.condition)]]))
    fig, axes = plt.subplots(len(rows), len(shown), figsize=(5.7*len(shown), 5.6*len(rows)), squeeze=False)
    for ax_row, (label, subset) in zip(axes, rows):
        z = subset(merged[(dataset, "none", "zscore", DEFAULT["zscore"])]); pub = key(z[z.call_pub.eq("depleted")])
        for ax, m in zip(ax_row, shown):
            x = subset(merged[(dataset, "none", m, DEFAULT[m])])
            venn(ax, [key(x[x.call.eq("depleted")]), pub], [f"{METHODS[m]} ({NORM_LABEL[DEFAULT[m]]})", "Published"], [COLOURS[m], COLOURS["published"]],
                 f"{METHODS[m]}: depleted calls" + (f"\n{label}" if len(rows) > 1 else ""))
    unit = "genes" if DATASETS[dataset][1] is None else "gene × condition pairs"
    fig.suptitle(f"{DATASETS[dataset][0]}: depleted calls ({unit}) at each method's default normalization vs the published hits", fontweight="bold", y=1.01)
    save(fig, f"{dataset}_2_venn_default")


def z_scatter(merged):
    """Our z-score (total reads) against the published z, all sites vs top site dropped, both datasets."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    for row, dataset in zip(axes, DATASETS):
        for ax, how in zip(row, ("none", "drop-top")):
            x = merged.get((dataset, how, "zscore", "totreads"))
            if x is None:ax.axis("off"); continue
            ok = np.isfinite(x.score) & np.isfinite(x.score_pub); x = x[ok]
            hit = x.call_pub.isin(["depleted", "enriched"])
            ax.scatter(x.score_pub[~hit], x.score[~hit], s=4, alpha=.25, color="#9aa5ad", edgecolor="none", rasterized=True)
            ax.scatter(x.score_pub[hit], x.score[hit], s=7, alpha=.6, color=COLOURS["published"], edgecolor="none", rasterized=True, label="published hit")
            lim = np.nanpercentile(np.abs(np.r_[x.score, x.score_pub]), 99.8); ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            ax.plot([-lim, lim], [-lim, lim], ls="--", color="black", lw=1)
            rho = np.median([g.score.rank().corr(g.score_pub.rank()) for _, g in x.groupby(["background", "condition"]) if len(g) > 2])
            ax.set_title(f"{DATASETS[dataset][0]} — {FILTERS[how].lower()}\nSpearman {rho:.2f} (median over conditions)", fontweight="bold")
            ax.set_xlabel("Published z"); ax.set_ylabel("Our z (total reads)"); style(ax)
    axes[0][0].legend(frameon=False, loc="upper left")
    fig.suptitle("Our z-score against the published z: summing all sites vs dropping each gene's top site", fontweight="bold")
    save(fig, "3_z_vs_published")


def error_control():
    """R27 simulated nulls: false calls per experiment and the pooled false discovery rate, by method and normalization."""
    t = false_calls(); t = t[t.normalization.isin(NORM_LABEL) | t.normalization.eq("cpm")]   # cpm: Gale's workflow as published
    t.to_csv(TAB/"error_control.tsv", sep="\t", index=False)
    order = [(m, n) for m in METHODS for n in NORMS[m] + (["cpm"] if m == "zscore" else []) if ((t.method == m) & (t.normalization == n)).any()]
    fig, axes = plt.subplots(1, 2, figsize=(15, 4.8))
    for ax, (col, lab, log) in zip(axes, [("simulated_null_false_calls", "False calls per experiment with no true effect\n(mean of 40 simulated experiments)", True),
                                          ("simulated_pooled_fdr", "False discovery rate\n(simulated screens with planted effects)", False)]):
        v = [float(t[(t.method == m) & (t.normalization == n)][col].iloc[0]) for m, n in order]
        x = np.arange(len(order)); ax.bar(x, [max(a, 1e-3) if log else a for a in v], color=[COLOURS[m] for m, _ in order], edgecolor="black", linewidth=.8*72/25.4)
        for i, a in enumerate(v):ax.text(i, max(a, 1e-3)*(1.15 if log else 1)+(0 if log else .01), f"{a:.2g}" if a < 10 else f"{a:.0f}", ha="center", va="bottom", fontsize=7.5, fontweight="bold", rotation=90 if log else 0)
        if log:ax.set_yscale("log"); ax.set_ylim(1e-3, 400)
        else:ax.axhline(.05, ls="--", color="grey", lw=1); ax.text(len(order)-.5, .06, "5% target", ha="right", fontsize=8)
        ax.set_xticks(x); ax.set_xticklabels([f"{NORM_LABEL.get(n, 'CPM')}" for _, n in order], rotation=45, ha="right"); ax.set_ylabel(lab); style(ax)
    shown = list(dict.fromkeys(m for m, _ in order))
    handles = [plt.Rectangle((0, 0), 1, 1, color=COLOURS[m], ec="black") for m in shown]
    axes[0].legend(handles, [METHODS[m] for m in shown], frameon=False, loc="upper right")
    fig.suptitle("Error control: how many calls are false when nothing (or something known) changed — pre-registered simulations", fontweight="bold")
    save(fig, "4_error_control")
    return t


def best_table(t):
    """Per dataset x method: the all-sites default normalization (BGC), aBGC on all sites and the best top-site treatment / normalization."""
    rows = []
    for (dataset, method), d in t.groupby(["dataset", "method"]):
        base = d[(d["filter"] == "none") & (d.normalization == DEFAULT[method])]
        both = d[["depleted_recall", "enriched_recall"]].mean(axis=1) if "enriched_recall" in d else d.depleted_recall   # both directions count
        best = d.assign(both=both).sort_values(["both", "spearman_vs_published_z"], ascending=False).head(1)
        abgc = d[(d["filter"] == "none") & (d.normalization == "aBGC")]
        for label, r in (("default", base), ("aBGC", abgc), ("best recovery", best)):
            if r.empty:continue
            r = r.iloc[0]
            rows.append(dict(dataset=DATASETS[dataset][0], method=METHODS[method], which=label, filter=FILTERS[r["filter"]], normalization=NORM_LABEL[r.normalization],
                             depleted=f"{r.depleted_recovered:.0f}/{r.depleted_published:.0f} ({100*r.depleted_recall:.0f}%)",
                             enriched="" if method == "gumbel" else f"{r.enriched_recovered:.0f}/{r.enriched_published:.0f} ({100*r.enriched_recall:.0f}%)",
                             depleted_calls_per_condition=round(r.depleted_called_per_condition), spearman=round(r.spearman_vs_published_z, 2)))
    b = pd.DataFrame(rows); b.to_csv(TAB/"best.tsv", sep="\t", index=False); return b


def build_docx():
    """report.md -> report.docx (pandoc, figures embedded, styles from reference.docx: Arial, bold headings) for review
    in Word; the lab copy is exported to PDF from Word."""
    md, docx, ref = REPORT/"report.md", REPORT/"report.docx", REPORT/"reference.docx"
    if not ref.exists():
        subprocess.run(["pandoc", "-o", str(ref), "--print-default-data-file", "reference.docx"], check=True)
    subprocess.run(["pandoc", str(md), "-o", str(docx), "--reference-doc", str(ref), "--resource-path", str(REPORT)], check=True, cwd=REPORT)
    return docx


if __name__ == "__main__":
    a = argparse.ArgumentParser(); a.add_argument("--docx", action="store_true", help="Only render report.md to report.docx"); a = a.parse_args()
    if a.docx:print(build_docx()); raise SystemExit
    t, merged = evaluate(); t = applicable(t)
    merged = {k:v for k, v in merged.items() if ((t.dataset == k[0]) & (t.method == k[2])).any()}
    for d in DATASETS:heatmaps(t, d); venns(merged, d)
    z_scatter(merged); error_control()
    with pd.option_context("display.width", 250, "display.max_columns", 30):print(best_table(t).to_string(index=False))
