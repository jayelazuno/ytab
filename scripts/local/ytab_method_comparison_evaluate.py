#!/usr/bin/env python3
"""Method-comparison study, step 3: agreement of the z-score, resampling and Gumbel calls on every dataset of
output/comparison (runs from ytab_method_comparison_runs.py), and against published calls where they exist
(ytab_method_comparison_published.py). Uses ytab.fitness.agreement, the same tables the app's Compare methods view reads.

Writes, under output/comparison/:
  <dataset>/agreement/{genes,pairs,sets}.tsv     per dataset (with "published" as a fourth method when available)
  summary/pairs_all.tsv                            every group x method pair of every dataset
  summary/published_recovery.tsv                   each method against the published calls: recall, precision, Spearman vs published z
  summary/*.png|pdf                                plots; RESULT.md                       the numbers in words
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]; OUT = ROOT/"output/comparison"; SUM = OUT/"summary"
sys.path.insert(0, str(ROOT/"src"))
from ytab.fitness import agreement  # noqa: E402
sys.path.append(str(ROOT/"scripts/local"))
import ytab_method_comparison_summary as summ  # noqa: E402

NAMES = {"zscore": "z-score", "resampling": "Resampling", "gumbel": "Gumbel", "published": "Published"}
COLOURS = {"zscore": "#0072B2", "resampling": "#E69F00", "gumbel": "#009E73", "published": "#CC79A7"}
# Cunningham-lab scope (user, 2026-11-05): only strains whose design meets the z-score of Gale 2020 (>= 2 controls), so
# every method is compared on data it is meant for. Single-control strains are not forced into the comparison.
CUNNINGHAM = {"Gale_2020": None, "Nickels_2025_glabrata": {"BG2", "pdr1_KO"}, "Nickels_2026_glabrata_tolerance": None}
EXCLUDED = {"Nickels_2024_glabrata": "1 control per strain", "Nickels_2025_glabrata_CBS138": "1 control", "Nickels_2026_glabrata_FKS2": "1 control"}


def in_scope(dataset: str, background: str) -> bool:
    keep = CUNNINGHAM.get(dataset); return keep is None or background in keep


def style(ax):
    for s in ("top", "right"):ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):ax.spines[s].set_linewidth(.8*72/25.4)
    ax.tick_params(width=.6*72/25.4); [t.set_fontweight("bold") for t in ax.get_xticklabels()+ax.get_yticklabels()]
    ax.xaxis.label.set_fontweight("bold"); ax.yaxis.label.set_fontweight("bold")


def save(fig, name, folder=None, exts=("png", "pdf")):
    folder = folder or SUM; folder.mkdir(parents=True, exist_ok=True); fig.tight_layout()
    for ext in exts:fig.savefig(folder/f"{name}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def finished(run: Path, method: str) -> bool:
    m = run/("zscore_manifest.json" if method == "zscore" else "manifest.json")
    if not m.exists():return False
    j = json.loads(m.read_text()); return method == "zscore" or any(g.get("status") == "success" for g in j.get("groups", {}).values())


def per_dataset():
    rows, recovery = [], []
    for d in sorted(p for p in OUT.iterdir() if (p/"runs").is_dir() and p.name not in EXCLUDED):
        inputs = {m: d/"runs"/m for m in ("zscore", "resampling", "gumbel") if finished(d/"runs"/m, m)}
        pub = d/"published.tsv"
        if len(inputs) + pub.exists() < 2:continue
        agreement.run(inputs, d/"agreement", pub if pub.exists() else None, sources={"dataset": d.name})
        p = pd.read_csv(d/"agreement/pairs.tsv", sep="\t"); p.insert(0, "dataset", d.name)
        rows.append(p[[in_scope(d.name, b) for b in p.background.astype(str)]])
        if pub.exists():
            g = pd.read_csv(d/"agreement/genes.tsv", sep="\t", low_memory=False)
            g = g[g.published_tested.astype(str).eq("True") & np.array([in_scope(d.name, b) for b in g.background.astype(str)], dtype=bool)]
            for m in inputs:
                for grp, x in g[g[f"{m}_tested"].astype(str).eq("True")].groupby("group"):
                    for direction in ("depleted", "enriched"):
                        if direction == "enriched" and m in agreement.DEPLETION_ONLY:continue
                        truth = x.published_call.eq(direction); hit = x[f"{m}_call"].eq(direction)
                        recovery.append(dict(dataset=d.name, group=grp, method=m, direction=direction, published=int(truth.sum()), called=int(hit.sum()),
                                             both=int((truth & hit).sum()), recall=(truth & hit).sum()/truth.sum() if truth.sum() else np.nan,
                                             precision=(truth & hit).sum()/hit.sum() if hit.sum() else np.nan,
                                             spearman_vs_published_z=x[f"{m}_score"].rank().corr(x.published_score.rank())))
    pairs = pd.concat(rows, ignore_index=True); rec = pd.DataFrame(recovery)
    SUM.mkdir(parents=True, exist_ok=True); pairs.to_csv(SUM/"pairs_all.tsv", sep="\t", index=False); rec.to_csv(SUM/"published_recovery.tsv", sep="\t", index=False)
    return pairs, rec


def strip(ax, data, groups, value, colour_of, ylabel):
    for i, k in enumerate(groups):
        v = data[k][value].dropna().to_numpy(); x = i + np.random.default_rng(1).uniform(-.18, .18, len(v))
        ax.scatter(x, v, s=14, alpha=.55, color=colour_of(k), edgecolor="none")
        if len(v):ax.hlines(np.median(v), i-.3, i+.3, color="black", lw=2)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels([NAMES.get(k, k) if isinstance(k, str) else " vs\n".join(NAMES[m] for m in k) for k in groups])
    ax.set_ylabel(ylabel); style(ax)


def plots(pairs, rec):
    methods = [p for p in pairs if True]
    within = pairs[(pairs.method_a != "published") & (pairs.method_b != "published")]
    keys = [k for k in [("zscore", "resampling"), ("zscore", "gumbel"), ("resampling", "gumbel")] if ((within.method_a == k[0]) & (within.method_b == k[1])).any()]
    data = {k: within[(within.method_a == k[0]) & (within.method_b == k[1])] for k in keys}
    mix = lambda k: "#555555"
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    strip(axes[0], data, keys, "depleted_a_recovered_by_b", mix, "Of the first method's depleted calls,\nshare the second also calls"); axes[0].set_ylim(0, 1)
    strip(axes[1], data, keys, "depleted_b_recovered_by_a", mix, "Of the second method's depleted calls,\nshare the first also calls"); axes[1].set_ylim(0, 1)
    strip(axes[2], data, keys, "spearman_score", mix, "Spearman, signed significance"); axes[2].set_ylim(-.2, 1)
    fig.suptitle("Method vs method, every background x condition of every dataset (bar: median; Gumbel pairs depleted only)", fontweight="bold"); save(fig, "A_method_agreement")
    if len(rec):
        ms = [m for m in ("zscore", "resampling", "gumbel") if m in set(rec.method)]
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
        for ax, (val, lab) in zip(axes, [("recall", "Recall of published hits"), ("precision", "Precision vs published hits"), ("spearman_vs_published_z", "Spearman vs published z")]):
            d = {m: rec[(rec.method == m) & (rec.direction == "depleted")] for m in ms}; strip(ax, d, ms, val, COLOURS.get, lab); ax.set_ylim(0 if val != "spearman_vs_published_z" else -.2, 1)
        fig.suptitle("Against the Cunningham lab's published z-score calls (|z| ≥ 3), depleted, per condition (bar: median)", fontweight="bold"); save(fig, "C_vs_published")


def upset(dataset, group, direction="depleted", folder=None, name=None, exts=("png", "pdf")):
    s = pd.read_csv(OUT/dataset/"agreement/sets.tsv", sep="\t"); s = s[(s.group == group) & (s.direction == direction)].sort_values("genes", ascending=False)
    if s.empty:return
    methods = [m for m in ("zscore", "resampling", "gumbel", "published") if s.methods.str.contains(m).any()]
    fig, (top, bot) = plt.subplots(2, 1, figsize=(max(6, .55*len(s)+2), 5.5), gridspec_kw=dict(height_ratios=[3, 1.2]), sharex=True)
    x = np.arange(len(s)); members = [m.split("&") for m in s.methods]
    colour = [COLOURS[m[0]] if len(m) == 1 else ("#087f85" if len(m) == len(methods) else "#9aa5ad") for m in members]
    top.bar(x, s.genes, color=colour, edgecolor="black", linewidth=.8*72/25.4)
    for i, v in enumerate(s.genes):top.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=8, fontweight="bold")
    top.set_ylabel("Genes"); style(top); top.set_title(f"{dataset} — {group.replace('__', ', ')}: {direction} calls", fontweight="bold")
    for r, m in enumerate(methods):
        on = [m in mm for mm in members]; bot.scatter(x, [r]*len(x), s=60, c=["black" if o else "#dddddd" for o in on])
    for i, mm in enumerate(members):
        rows = [methods.index(m) for m in mm]
        if len(rows) > 1:bot.vlines(i, min(rows), max(rows), color="black", lw=2)
    bot.set_yticks(range(len(methods))); bot.set_yticklabels([NAMES[m] for m in methods], fontweight="bold"); bot.set_xticks([]); bot.invert_yaxis()
    for sp in bot.spines.values():sp.set_visible(False)
    save(fig, name or f"D_upset_{dataset}_{group}", folder, exts)


def report(pairs, rec):
    w = pairs[(pairs.method_a != "published") & (pairs.method_b != "published")]
    lines = ["# Method comparison: z-score vs resampling vs control-conditioned Gumbel", "",
             "Each method with its default normalization (z-score BGC, resampling BGC, Gumbel TTR), same comparisons, run by",
             "scripts/local/ytab_method_comparison_runs.py; agreement by ytab.fitness.agreement (the app's Compare methods tables).",
             "Published calls: Cunningham-lab z-scores, |z| >= 3 (Gale 2020; Nickels & Cunningham 2025 Table S1, which includes the",
             "Nickels 2024 CBS138 screens). Untracked; for our own checks.", "",
             "Scope: Cunningham-lab strains whose design meets the z-score of Gale 2020 (>= 2 no-drug controls): Gale 2020 BG14,",
             "Nickels 2025 BG2 and pdr1_KO, Nickels 2026 tolerance BG2. Not compared (single control, not forced): " + "; ".join(f"{k} ({v})" for k, v in EXCLUDED.items())
             + ", and Nickels 2025 PDR1_GOF. Other datasets: method vs method only, each method where its requirements are met.", "",
             ]
    for title, cohort in (("Cunningham lab, four strains", w[w.dataset.isin(list(CUNNINGHAM))]), ("Other datasets", w[~w.dataset.isin(list(CUNNINGHAM))])):
        lines += [f"## Method vs method: {title} (median over background x condition)", "",
                  "| A vs B | groups | Spearman effect | Spearman significance | depleted calls A / B | A's depleted also by B | B's depleted also by A | enriched calls A / B |", "|---|---|---|---|---|---|---|---|"]
        for (a, b), x in cohort.groupby(["method_a", "method_b"]):
            enr = "depletion only" if x.enriched_n_a.isna().all() else f"{x.enriched_n_a.median():.0f} / {x.enriched_n_b.median():.0f}"
            lines.append(f"| {NAMES[a]} vs {NAMES[b]} | {len(x)} | {x.spearman_effect.median():.2f} | {x.spearman_score.median():.2f} | {x.depleted_n_a.median():.0f} / {x.depleted_n_b.median():.0f} | "
                         f"{x.depleted_a_recovered_by_b.median():.0%} | {x.depleted_b_recovered_by_a.median():.0%} | {enr} |")
        lines.append("")
    if len(rec):
        lines += ["", "## Against published calls (median over conditions)", "", "| dataset | method | direction | conditions | recall | precision | Spearman vs published z |", "|---|---|---|---|---|---|---|"]
        for (d, m, dr), x in rec.groupby(["dataset", "method", "direction"]):
            lines.append(f"| {d} | {NAMES[m]} | {dr} | {len(x)} | {x.recall.median():.2f} | {x.precision.median():.2f} | {x.spearman_vs_published_z.median():.2f} |")
    status = OUT/"runs_status.tsv"
    if status.exists():
        s = pd.read_csv(status, sep="\t", names=["dataset", "method", "exit", "seconds"]); bad = s[s.exit != 0]
        lines += ["", "## Runs", "", f"{len(s)} runs, {len(bad)} not completed: " + ("; ".join(f"{r.dataset} {r.method}" for r in bad.itertuples()) or "none") + ".",
                  "The Gumbel needs controls and treatments from one insertion library (it compares runs of sites); datasets whose libraries",
                  "differ, or backgrounds with a single control, get resampling only."]
    (OUT/"RESULT.md").write_text("\n".join(lines)+"\n")


def safe(x): return "".join(c if c.isalnum() or c in "-_." else "_" for c in x)


def dataset_plots(d: Path):
    """output/comparison/<dataset>/plots: an overview of every condition, and per background x condition the agreement
    matrix, effect and significance scatters of every method pair, UpSet of the depleted and enriched calls, and each
    method against the published z where there are published calls."""
    a = d/"agreement"
    if not (a/"pairs.tsv").exists():return 0
    pairs = pd.read_csv(a/"pairs.tsv", sep="\t"); genes = pd.read_csv(a/"genes.tsv", sep="\t", low_memory=False)
    keep = np.array([in_scope(d.name, b) for b in pairs.background.astype(str)], dtype=bool); pairs = pairs[keep]
    genes = genes[np.array([in_scope(d.name, b) for b in genes.background.astype(str)], dtype=bool)]
    methods = [m for m in ("zscore", "resampling", "gumbel", "published") if f"{m}_tested" in genes]
    out = d/"plots"; n = 0
    # overview: of each method's depleted calls, the share each other method also calls, per condition
    both = pd.concat([pairs.assign(pair=pairs.method_a.map(NAMES)+" in "+pairs.method_b.map(NAMES), value=pairs.depleted_a_recovered_by_b),
                      pairs.assign(pair=pairs.method_b.map(NAMES)+" in "+pairs.method_a.map(NAMES), value=pairs.depleted_b_recovered_by_a)])
    h = both.pivot_table(index="group", columns="pair", values="value")
    if len(h):
        fig, ax = plt.subplots(figsize=(1.6*h.shape[1]+3, .32*len(h)+1.8)); im = ax.imshow(h.to_numpy(dtype=float), vmin=0, vmax=1, cmap="GnBu", aspect="auto")
        for (i, j), v in np.ndenumerate(h.to_numpy(dtype=float)):ax.text(j, i, "" if np.isnan(v) else f"{v:.2f}", ha="center", va="center", fontsize=8, fontweight="bold")
        ax.set_xticks(range(h.shape[1])); ax.set_xticklabels(h.columns, rotation=25, ha="right", fontweight="bold")
        ax.set_yticks(range(len(h))); ax.set_yticklabels([g.replace("__", ", ") for g in h.index], fontweight="bold", fontsize=8)
        fig.colorbar(im, ax=ax, label="share of the first method's depleted calls\nalso called by the second"); ax.set_title(f"{d.name}: depleted-call recovery per condition", fontweight="bold")
        save(fig, "overview_recovery_depleted", out); n += 1
    for grp, g in genes.groupby("group"):
        gdir = out/safe(grp); p = pairs[pairs.group == grp]; title = f"{d.name}: {grp.replace('__', ', ')}"
        # agreement matrix (depleted)
        ms = [m for m in methods if g[f"{m}_tested"].astype(str).eq("True").any()]
        if len(ms) >= 2:
            M = np.full((len(ms), len(ms)), np.nan); L = [["" for _ in ms] for _ in ms]
            for i, m1 in enumerate(ms):
                for j, m2 in enumerate(ms):
                    if i == j:L[i][j] = f"{NAMES[m1]}"; continue
                    r = p[(p.method_a == m1) & (p.method_b == m2)]; flip = False
                    if not len(r):r = p[(p.method_a == m2) & (p.method_b == m1)]; flip = True
                    if len(r):
                        r = r.iloc[0]; share = r.depleted_b_recovered_by_a if flip else r.depleted_a_recovered_by_b
                        mine, both_n = (r.depleted_n_b if flip else r.depleted_n_a), r.depleted_n_both
                        M[i, j] = share; L[i][j] = f"{share:.0%}\n{both_n:.0f} of {mine:.0f}"
            fig, ax = plt.subplots(figsize=(1.5*len(ms)+2.5, 1.3*len(ms)+1.2)); im = ax.imshow(M, vmin=0, vmax=1, cmap="GnBu")
            for i in range(len(ms)):
                for j in range(len(ms)):ax.text(j, i, L[i][j], ha="center", va="center", fontsize=8, fontweight="bold")
            ax.set_xticks(range(len(ms))); ax.set_xticklabels([NAMES[m] for m in ms], fontweight="bold"); ax.set_yticks(range(len(ms))); ax.set_yticklabels([NAMES[m] for m in ms], fontweight="bold")
            ax.set_ylabel("Calls of", fontweight="bold"); ax.set_xlabel("Also called by", fontweight="bold")
            fig.colorbar(im, ax=ax, label="row's depleted calls also called by column"); ax.set_title(title, fontweight="bold", fontsize=10); save(fig, "1_recovery_matrix", gdir, ("png",)); n += 1
        # scatters of every pair, effect and significance
        prs = [(m1, m2) for i, m1 in enumerate(ms) for m2 in ms[i+1:] if not (m2 == "published" and g["published_effect"].isna().all() and False)]
        for what, label in (("effect", "effect"), ("score", "signed significance")):
            use = [(m1, m2) for m1, m2 in prs if g[f"{m1}_{what}"].notna().any() and g[f"{m2}_{what}"].notna().any()]
            if not use:continue
            cols = min(3, len(use)); rows_n = -(-len(use)//cols)
            fig, axes = plt.subplots(rows_n, cols, figsize=(4.4*cols, 4*rows_n), squeeze=False)
            for ax in axes.flat[len(use):]:ax.set_visible(False)
            for ax, (m1, m2) in zip(axes.flat, use):
                x = g[g[f"{m1}_tested"].astype(str).eq("True") & g[f"{m2}_tested"].astype(str).eq("True")]
                x = x[np.isfinite(x[f"{m1}_{what}"]) & np.isfinite(x[f"{m2}_{what}"])]
                called = ["depleted"] if "gumbel" in (m1, m2) else ["depleted", "enriched"]
                c1 = x[f"{m1}_call"].isin(called); c2 = x[f"{m2}_call"].isin(called)
                for mask, col, lab in ((~c1 & ~c2, "#d9d9d9", "neither"), (c1 & ~c2, COLOURS[m1], f"{NAMES[m1]} only"), (~c1 & c2, COLOURS[m2], f"{NAMES[m2]} only"), (c1 & c2, "#222222", "both")):
                    ax.scatter(x.loc[mask, f"{m1}_{what}"], x.loc[mask, f"{m2}_{what}"], s=7, alpha=.6, color=col, edgecolor="none", label=f"{lab} ({mask.sum()})")
                rho = x[f"{m1}_{what}"].rank().corr(x[f"{m2}_{what}"].rank())
                ax.axhline(0, color="#bbbbbb", lw=.6); ax.axvline(0, color="#bbbbbb", lw=.6)
                ax.set_xlabel(NAMES[m1]); ax.set_ylabel(NAMES[m2]); ax.set_title(f"ρ = {rho:.2f} ({len(x):,} genes)", fontsize=9, fontweight="bold"); ax.legend(fontsize=7, frameon=False); style(ax)
            fig.suptitle(f"{title}: {label}", fontweight="bold", fontsize=10); save(fig, f"2_scatter_{what}", gdir, ("png",)); n += 1
        for direction in ("depleted", "enriched"):
            if (a/"sets.tsv").exists():upset(d.name, grp, direction, gdir, f"3_upset_{direction}", ("png",)); n += 1
            use = [m for m in ms if not (direction == "enriched" and m in agreement.DEPLETION_ONLY)]
            use = (["zscore", "resampling", "published"] if {"zscore", "resampling", "published"} <= set(use) else use)[:3]
            if len(use) >= 2:
                t = g[np.logical_and.reduce([g[f"{m}_tested"].astype(str).eq("True") for m in use])]
                fig, ax = plt.subplots(figsize=(6, 5.4))
                summ.venn(ax, [set(t.ORF[t[f"{m}_call"] == direction]) for m in use], [NAMES[m] for m in use], [COLOURS[m] for m in use], f"{title}\n{direction} calls")
                save(fig, f"3_venn_{direction}", gdir, ("png",)); n += 1
    return n


def main():
    pairs, rec = per_dataset(); plots(pairs, rec); report(pairs, rec)
    for d in sorted(p for p in OUT.iterdir() if (p/"agreement").is_dir() and p.name not in EXCLUDED):
        print(f"{d.name}: {dataset_plots(d)} plots in {d/'plots'}")
    for d, grp in [("Gale_2020", "BG14__treated"), ("Nickels_2025_glabrata", "BG2__micafungin_8ng"), ("Nickels_2026_glabrata_tolerance", None)]:
        if (OUT/d/"agreement/sets.tsv").exists():
            g = grp or pd.read_csv(OUT/d/"agreement/sets.tsv", sep="\t").group.iloc[0]; upset(d, g)
    print((OUT/"RESULT.md").read_text())


if __name__ == "__main__":main()
