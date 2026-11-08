#!/usr/bin/env python3
"""Method-comparison study, summary of the Gale 2020 normalization sweep in three plots and a report
(output/comparison/Gale_2020/summary/):

  1_recovery_by_normalization   recall of Gale's published hits per method x normalization, with the number of calls
  2_false_calls_by_normalization  false calls: pre-registered R27 simulated nulls (mean per experiment; pooled FDR over
                                  planted-effect screens) and the real-data swap null
  3_venn_default                  calls at each method's default normalization vs published (UpSet versions alongside)

Published hits (Gale et al. 2020, per-gene z and log2FC): strict = z <= -3 and log2FC <= -1 (19 depleted; z >= 3 and
log2FC >= 1, 185 enriched), the paper's significance with a 2-fold change; |z| >= 3 alone gives 37 / 284.
Recovery = of the published hits tested by the method, the share it calls in the same direction. No Jaccard (capped by
unequal list sizes). The Gumbel calls depletion only.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT/"src"))
from ytab.fitness import agreement as ag  # noqa: E402

SWEEP = ROOT/"output/comparison/Gale_2020/sweep"; OUT = ROOT/"output/comparison/Gale_2020/summary"
R27 = ROOT/"output/treated_vs_control/benchmark_r27_2026-10-31/metrics.tsv"
PUB = ROOT/"output/essentiality/control_conditioned_gumbel/Gale_2020/publication_comparison.tsv"
NORMS = ["nonorm", "cpm", "totreads", "nzmean", "ttr", "midlc", "quantile", "zinfnb", "betageom", "aBGC"]
NORM_LABEL = dict(nonorm="None", cpm="CPM", totreads="Total reads", nzmean="NZ mean", ttr="TTR", midlc="MidLC", quantile="Quantile",
                  zinfnb="ZINB", betageom="BGC", aBGC="aBGC")
METHODS = {"zscore": "z-score", "resampling": "Resampling", "gumbel": "Gumbel"}
COLOURS = {"zscore": "#0072B2", "resampling": "#E69F00", "gumbel": "#009E73", "published": "#CC79A7"}
DEFAULT = {"zscore": "betageom", "resampling": "betageom", "gumbel": "ttr"}


def style(ax):
    for s in ("top", "right"):ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):ax.spines[s].set_linewidth(.8*72/25.4)
    ax.tick_params(width=.6*72/25.4)
    for t in ax.get_xticklabels()+ax.get_yticklabels():t.set_fontweight("bold")
    ax.xaxis.label.set_fontweight("bold"); ax.yaxis.label.set_fontweight("bold")


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True); fig.tight_layout()
    for ext in ("png", "pdf"):fig.savefig(OUT/f"{name}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def published():
    p = pd.read_csv(PUB, sep="\t"); p = p[np.isfinite(p.published_z)]
    z, l = p.published_z, p.published_log2fc
    return {("depleted", "strict"): set(p.ORF[(z <= -3) & (l <= -1)]), ("depleted", "|z|>=3"): set(p.ORF[z <= -3]),
            ("enriched", "strict"): set(p.ORF[(z >= 3) & (l >= 1)]), ("enriched", "|z|>=3"): set(p.ORF[z >= 3])}, p.set_index("ORF").published_z


def runs():
    out = {}
    for d in sorted((SWEEP/"real").glob("*_*")):
        if not d.is_dir() or d.name == "project":continue
        m, n = d.name.split("_", 1)
        marker = d/("zscore_manifest.json" if m == "zscore" else "manifest.json")
        if marker.exists():
            try:out[(m, n)] = ag.read_method(d, m)
            except Exception as e:print(f"skip {d.name}: {e}")
    return out


def recovery_table(rs, pub, pz):
    rows = []
    for (m, n), x in rs.items():
        x = x[x.tested].set_index("ORF")
        rank = x.score.abs().rank(ascending=False) if m == "zscore" else (x.score.abs().rank(ascending=False))
        for (direction, kind), hits in pub.items():
            if direction == "enriched" and m == "gumbel":continue
            tested = hits & set(x.index); called = set(x.index[x.call == direction])
            signed = x.score if direction == "enriched" else -x.score
            r = signed.rank(ascending=False)
            rows.append(dict(method=m, normalization=n, direction=direction, published=kind, published_tested=len(tested),
                             recovered=len(tested & called), recall=len(tested & called)/len(tested) if tested else np.nan, calls=len(called),
                             median_rank_of_published=float(r.reindex(list(tested)).median()) if tested else np.nan,
                             published_in_top50=int((r.reindex(list(tested)) <= 50).sum()) if tested else 0,
                             spearman_vs_published_z=float(x.score.rank().corr(pz.reindex(x.index).rank()))))
    return pd.DataFrame(rows)


def plot_recovery(t):
    fig, axes = plt.subplots(2, 2, figsize=(13, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[1.3, 1]))
    xs = {n: i for i, n in enumerate(NORMS)}; off = {"zscore": -.22, "resampling": 0, "gumbel": .22}
    for j, direction in enumerate(("depleted", "enriched")):
        top, bot = axes[0, j], axes[1, j]
        for m in METHODS:
            if direction == "enriched" and m == "gumbel":continue
            for kind, filled in (("|z|>=3", False), ("strict", True)):   # strict drawn last and smaller, so both show at the same recall
                d = t[(t.method == m) & (t.direction == direction) & (t.published == kind)]
                d = d[d.normalization.isin(NORMS)]; x = [xs[n]+off[m] for n in d.normalization]
                top.scatter(x, d.recall, s=30 if filled else 90, marker="o", facecolor=COLOURS[m] if filled else "white", edgecolor=COLOURS[m], linewidth=1.6,
                            label=f"{METHODS[m]}, {'strict' if filled else '|z| ≥ 3'}" if j == 0 else None, zorder=3)
            d = t[(t.method == m) & (t.direction == direction) & (t.published == "strict") & t.normalization.isin(NORMS)]
            bot.bar([xs[n]+off[m] for n in d.normalization], d.calls, width=.21, color=COLOURS[m], edgecolor="black", linewidth=.8*72/25.4/2)
        n_strict = int(t[(t.direction == direction) & (t.published == "strict")].published_tested.max())
        n_loose = int(t[(t.direction == direction) & (t.published == "|z|>=3")].published_tested.max())
        top.set_ylim(-.03, 1.05); top.set_ylabel("Recall of published hits"); top.set_title(f"{direction.capitalize()}: published {n_strict} strict / {n_loose} |z| ≥ 3",
                                                                                               fontweight="bold")
        bot.set_ylabel(f"Genes called {direction}"); style(top); style(bot)
    for ax in axes[1]:ax.set_xticks(range(len(NORMS))); ax.set_xticklabels([NORM_LABEL[n] for n in NORMS], rotation=35, ha="right")
    axes[0, 0].legend(fontsize=8, frameon=False, ncol=2, loc="lower left")
    fig.suptitle("Gale 2020 fluconazole: recovery of the published hits by normalization (bars: genes each method calls)", fontweight="bold")
    save(fig, "1_recovery_by_normalization")


def false_calls():
    m = pd.read_csv(R27, sep="\t", low_memory=False, keep_default_na=False, na_values={"auc": [""], "fdp": [""]})   # the scenario "null" is not NA
    for c in ("calls", "false_calls"):m[c] = pd.to_numeric(m[c], errors="coerce")
    arm = m.arm.str.extract(r"^(zscore|resampling|gumbel)_(.+)$"); m["method"], m["normalization"] = arm[0], arm[1].replace({"bgc": "betageom"})
    m = m[m.method.notna() & ~m.normalization.str.endswith("pairs", na=False)]
    null = m[m.scenario == "null"].groupby(["method", "normalization"]).calls.mean().rename("simulated_null_false_calls")
    eff = m[m.scenario != "null"].groupby(["method", "normalization"])[["false_calls", "calls"]].sum()
    fdr = (eff.false_calls/eff.calls).rename("simulated_pooled_fdr")
    t = pd.concat([null, fdr], axis=1).reset_index()
    swap = SWEEP/"swap_null_arms.tsv"
    if swap.exists():
        s = pd.read_csv(swap, sep="\t").assign(normalization=lambda x: x.normalization.replace({"bgc": "betageom"}))
        t = t.merge(s[["method", "normalization", "calls"]].rename(columns=dict(calls="swap_null_calls")), how="outer", on=["method", "normalization"])
    return t


def plot_false_calls(t):
    fig, axes = plt.subplots(1, 3 if "swap_null_calls" in t else 2, figsize=(15, 4.4), sharex=True)
    xs = {n: i for i, n in enumerate(NORMS)}; off = {"zscore": -.25, "resampling": 0, "gumbel": .25}
    panels = [("simulated_null_false_calls", "False calls per null experiment\n(R27 simulations, mean)"), ("simulated_pooled_fdr", "False discovery rate\n(R27 planted-effect screens)")]
    if "swap_null_calls" in t:panels.append(("swap_null_calls", "Calls in the real swap null\n(control+FLU vs control+FLU)"))
    notes = {"swap_null_calls": "z-score: 0, but not a calibration test — its noise\nmodel is fitted on a control and a FLU culture here"}
    for ax, (col, lab) in zip(axes, panels):
        for m in METHODS:
            d = t[(t.method == m) & t.normalization.isin(NORMS) & t[col].notna()]
            ax.bar([xs[n]+off[m] for n in d.normalization], d[col], width=.24, color=COLOURS[m], edgecolor="black", linewidth=.6, label=METHODS[m])
            top = t[col].max()
            for n, v in zip(d.normalization, d[col]):   # values too small to see as bars
                if v < .04*top:ax.text(xs[n]+off[m], top*.02, f"{v:.2g}", rotation=90, ha="center", va="bottom", fontsize=7, color=COLOURS[m], fontweight="bold")
        if col == "simulated_pooled_fdr":ax.axhline(.1, ls="--", color="grey", lw=1); ax.text(len(NORMS)-.5, .105, "10%", ha="right", fontsize=8)
        if col in notes:ax.text(.03, .97, notes[col], transform=ax.transAxes, ha="left", va="top", fontsize=8, color=COLOURS["zscore"])
        ax.set_ylabel(lab); ax.set_xticks(range(len(NORMS))); ax.set_xticklabels([NORM_LABEL[n] for n in NORMS], rotation=35, ha="right"); style(ax)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Are the extra calls justified? False calls by method and normalization (z-score calls are not error-controlled)", fontweight="bold")
    save(fig, "2_false_calls_by_normalization")


def venn(ax, sets, names, colours, title):
    k = len(sets); cen = [(-.6, 0), (.6, 0)] if k == 2 else [(-.65, .45), (.65, .45), (0, -.6)]
    for (x, y), c in zip(cen, colours):ax.add_patch(Circle((x, y), 1.15, facecolor=c, alpha=.25, edgecolor="black", lw=1.2))
    def n(inside):
        s = set.intersection(*[sets[i] for i in inside]); out = set.union(*[sets[i] for i in range(k) if i not in inside]) if len(inside) < k else set()
        return len(s - out)
    spots = ({(0,): (-1.15, 0), (1,): (1.15, 0), (0, 1): (0, 0)} if k == 2 else
             {(0,): (-1.2, .85), (1,): (1.2, .85), (2,): (0, -1.25), (0, 1): (0, .95), (0, 2): (-.75, -.25), (1, 2): (.75, -.25), (0, 1, 2): (0, .1)})
    for inside, (x, y) in spots.items():ax.text(x, y, f"{n(inside):,}", ha="center", va="center", fontsize=12, fontweight="bold")
    for (x, y), name, s in zip(cen, names, sets):ax.text(x*1.9, y*1.9+(.55 if y >= 0 else -.6), f"{name}\n({len(s):,})", ha="center", va="center", fontsize=10, fontweight="bold")
    ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.4, 2.4); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontweight="bold")


def upset(ax, sets, names, title):
    k = len(sets); members = {}
    for g in set.union(*sets):members.setdefault(tuple(i for i in range(k) if g in sets[i]), []).append(g)
    combos = sorted(members, key=lambda c: -len(members[c])); x = np.arange(len(combos))
    ax.bar(x, [len(members[c]) for c in combos], color="#9aa5ad", edgecolor="black", linewidth=.8)
    for i, c in enumerate(combos):
        ax.text(i, len(members[c]), str(len(members[c])), ha="center", va="bottom", fontsize=8, fontweight="bold")
        for r in range(k):ax.scatter(i, -(r+1)*max(len(v) for v in members.values())*.12, s=40, color="black" if r in c else "#dddddd", clip_on=False)
    ax.set_xticks([]); ax.set_title(title, fontweight="bold"); style(ax)
    ax.set_yticks(list(ax.get_yticks()[ax.get_yticks() >= 0]) + [-(r+1)*max(len(v) for v in members.values())*.12 for r in range(k)])
    ax.set_yticklabels([f"{int(v)}" for v in ax.get_yticks()[: -k]] + names)


def plot_venn(rs, pub):
    calls = {m: (rs[(m, DEFAULT[m])] if (m, DEFAULT[m]) in rs else None) for m in METHODS}
    def called(m, d):
        x = calls[m]; return set(x.ORF[(x.call == d) & x.tested]) if x is not None else set()
    panels = [("depleted", ["zscore", "resampling"], "strict"), ("enriched", ["zscore", "resampling"], "strict"), ("depleted", ["gumbel"], "strict")]
    for kind, fn in (("venn", venn), ("upset", upset)):
        fig, axes = plt.subplots(1, 3, figsize=(17, 5.6))
        for ax, (d, ms, p) in zip(axes, panels):
            sets = [called(m, d) for m in ms]+[pub[(d, p)]]; names = [f"{METHODS[m]} ({NORM_LABEL[DEFAULT[m]]})" for m in ms]+["Published (strict)"]
            title = f"{d.capitalize()} calls"
            if kind == "venn":fn(ax, sets, names, [COLOURS[m] for m in ms]+[COLOURS["published"]], title)
            else:fn(ax, sets, names, title)
        fig.suptitle("Gale 2020 fluconazole: calls at each method's default normalization vs the published hits", fontweight="bold", y=1.04)
        save(fig, "3_venn_default" if kind == "venn" else "3_upset_default")


def main():
    pub, pz = published(); rs = runs()
    if not rs:raise SystemExit("No sweep results yet")
    t = recovery_table(rs, pub, pz); f = false_calls()
    OUT.mkdir(parents=True, exist_ok=True); t.to_csv(OUT/"recovery_by_normalization.tsv", sep="\t", index=False); f.to_csv(OUT/"false_calls_by_normalization.tsv", sep="\t", index=False)
    plot_recovery(t); plot_false_calls(f); plot_venn(rs, pub)
    s = t[t.published == "strict"].pivot_table(index=["method", "normalization"], columns="direction", values=["recall", "calls", "published_in_top50"])
    print(s.round(2).to_string()); print(f.round(3).to_string())


if __name__ == "__main__":main()
