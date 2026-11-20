#!/usr/bin/env python3
"""Comet plot from a treated-vs-control run: write the plot-ready table, then draw it.

  python ytab_comet.py RUN_DIR --gene-sets SETS.gmt --orthologs TABLE.txt -o OUT_DIR

RUN_DIR is a treated_vs_parent run directory. The table is written first and the figure is drawn from it, so
the figure can be redrawn, restyled or replaced without recomputing anything.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from ytab.fitness.comet import Settings, comet_table, ortholog_names, UNASSIGNED

# Colours are assigned in the GMT's declared order. Grey is reserved for genes in no set, which must stay
# visually behind the overlays however many sets are declared.
PALETTE = ["#d4a017", "#2e8b3d", "#e07a28", "#2f6fb5", "#8e44ad", "#c0392b", "#16a085", "#7f8c8d"]
BODY = "#1a1a1a"
GREY = "#9aa0a6"


def load_run(run_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Read one run's results and the comparisons that went into them.

    A run directory holds one results table. Its hash-suffixed subdirectories are byte-identical copies of that
    same table (verified on Gale 2020: all five match the top level and each other), not per-comparison splits,
    so reading them would multiply every gene by the number of copies. One run is therefore one panel, and a
    multi-panel comet is built from several run directories.
    """
    results = run_dir / "treated_vs_parent_results.csv"
    if not results.exists():
        raise SystemExit(f"No treated_vs_parent_results.csv in {run_dir}")
    d = pd.read_csv(results)
    summary_path = run_dir / "treated_vs_parent_comparison_summary.csv"
    summary = pd.read_csv(summary_path) if summary_path.exists() else None
    return d, summary


def panel_label(run_dir: Path, summary: pd.DataFrame | None) -> str:
    """Name a panel after the background it compares, falling back to the run directory."""
    if summary is not None and "background" in summary.columns:
        backgrounds = [b for b in dict.fromkeys(summary["background"].astype(str)) if b and b != "nan"]
        if len(backgrounds) == 1:
            return backgrounds[0]
        if backgrounds:
            return ", ".join(backgrounds)
    return run_dir.name


def _place_labels(ax, rows) -> None:
    """Draw each label at the first candidate offset that does not overlap one already placed.

    Without this, genes that sit close together on the log axes get labels written on top of each other, which
    reads as a single nonsense name (ERG6 over ERG3 reads as "ERG63"). The candidates are tried in a fixed
    order, so the placement is deterministic: the same table always produces the same figure.
    """
    placed = []
    # Short offsets first: at this distance the label sits against its own point and needs no leader line.
    # Leader lines are what garble a dense figure -- the boxes can be collision-free while a line from one
    # label crosses the text of another, which reads as a corrupted gene name (ERG6 struck through by ERG3's
    # line). A line is therefore drawn only for the long fallback offsets, where the label has travelled far
    # enough that its owner is genuinely ambiguous.
    candidates = [(7, 3), (7, -8), (-7, 4), (-7, -9), (0, 9), (0, -11),
                  (22, 7), (-22, 7), (22, -11), (-22, -11)]
    LEADER_FROM = 6                                # index at which offsets become long enough to need a line
    renderer = ax.figure.canvas.get_renderer()
    for _, r in rows.iterrows():
        for i, (dx, dy) in enumerate(candidates):
            kw = dict(arrowprops=dict(arrowstyle="-", lw=.4, color="#777", shrinkA=0, shrinkB=2)) \
                 if i >= LEADER_FROM else {}
            ann = ax.annotate(r["label_name"], xy=(r["x"], r["y"]), xytext=(dx, dy),
                              textcoords="offset points", fontsize=5.8, style="italic", color="#111",
                              ha="left" if dx >= 0 else "right", va="center", **kw)
            box = ann.get_window_extent(renderer=renderer).expanded(1.12, 1.30)
            if not any(box.overlaps(b) for b in placed):
                placed.append(box)
                break
            ann.remove()
        else:                                   # nowhere free: drop the label rather than print a collision
            continue


def draw(table: pd.DataFrame, out: Path, *, title: str = "", panel_size: float = 2.9) -> None:
    """One panel per comparison, shared log-log axes, gene sets over a dark unassigned body."""
    panels = list(dict.fromkeys(table["comparison_id"]))
    sets = [s for s in dict.fromkeys(table["gene_set"]) if s != UNASSIGNED]
    colours = {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(sets)}

    lo = 1.0
    hi = float(max(table["x"].max(), table["y"].max())) * 1.7
    fig, axes = plt.subplots(len(panels), 1, figsize=(panel_size + 0.9, panel_size * len(panels)),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes[:, 0]

    for ax, panel in zip(axes, panels):
        d = table[table["comparison_id"] == panel]
        base, over = d[d["gene_set"] == UNASSIGNED], d[d["gene_set"] != UNASSIGNED]
        ax.plot([lo, hi], [lo, hi], color=GREY, lw=.7, ls="--", zorder=1)
        ax.scatter(base["x"], base["y"], s=1.9, c=BODY, linewidths=0, alpha=.5, zorder=2)
        for name in sets:
            g = over[over["gene_set"] == name]
            if not g.empty:
                ax.scatter(g["x"], g["y"], s=11, c=colours[name], edgecolors="white",
                           linewidths=.35, zorder=3)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")
        label = str(d["background"].iloc[0]) if str(d["background"].iloc[0]) else panel
        ax.text(.04, .95, label, transform=ax.transAxes, fontsize=9, fontweight="bold",
                style="italic", va="top")
        # The house style: black axis lines at 0.8 mm, ticks at 0.6 mm, bold black axis labels.
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("black"); ax.spines[side].set_linewidth(2.27)
        ax.tick_params(colors="black", width=1.70, labelsize=7)
        ax.grid(True, which="major", color="#ebebeb", lw=.45, zorder=0)

    axes[-1].set_xlabel("Control (reads/gene)", fontweight="bold", color="black", fontsize=8.5)
    axes[len(axes)//2].set_ylabel("Treated (reads/gene)", fontweight="bold", color="black", fontsize=8.5)

    handles = [Line2D([], [], marker="o", ls="", ms=3.4, mfc=BODY, mec="none", label="all genes")]
    handles += [Line2D([], [], marker="o", ls="", ms=4.2, mfc=colours[s], mec="white", mew=.3, label=s)
                for s in sets]
    # Outside the axes: an in-plot legend covers the top-right of the comet, which is where the
    # best-supported genes sit.
    fig.legend(handles=handles, frameon=False, fontsize=6.2, loc="center left",
               bbox_to_anchor=(.86, .5), handletextpad=.35, labelspacing=.5)
    if title:
        fig.suptitle(title, fontsize=9, fontweight="bold")

    # Labels are measured in display coordinates, so every transform that moves an axis must already have
    # happened: tight_layout first, then a draw to build the renderer, then placement. Placing before
    # tight_layout silently defeats the collision test, because the boxes it compared are no longer where
    # the text ends up.
    fig.tight_layout()
    fig.canvas.draw()
    for ax, panel in zip(axes, panels):
        d = table[table["comparison_id"] == panel]
        _place_labels(ax, d[d["label"]].sort_values("off_diagonal", key=abs, ascending=False))

    for ext in ("png", "pdf"):
        fig.savefig(out.with_suffix("." + ext), dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run_dirs", type=Path, nargs="+", help="one panel per run directory")
    p.add_argument("-o", "--output", type=Path,
                   help="default: <first run directory>/comet, where the app looks for it")
    p.add_argument("--gene-sets", type=Path)
    p.add_argument("--orthologs", type=Path)
    p.add_argument("--id-column", default="cglab_id")
    p.add_argument("--name-column", default="scer_gene_name")
    p.add_argument("--set-order", nargs="*", default=[])
    p.add_argument("--label-sets", nargs="*", default=[])
    p.add_argument("--label-count", type=int, default=12)
    p.add_argument("--label-min-reads", type=float, default=10.0)
    p.add_argument("--direction", choices=("both", "depleted", "enriched"), default="both",
                   help="which hits may be labelled; every set member is drawn either way")
    p.add_argument("--label-genes", nargs="*", default=[])
    p.add_argument("--min-reads", type=float, default=0.0)
    p.add_argument("--labels", nargs="*", default=[], help="panel labels; default: each run's background")
    p.add_argument("--title", default="")
    a = p.parse_args(argv)

    # The app finds a run's comet table at <run>/comet/comet_table.csv, so that is where one lands by
    # default. An explicit -o is for figures built from several runs, which belong to no single run.
    if a.output is None:
        a.output = a.run_dirs[0] / "comet"
    if a.labels and len(a.labels) != len(a.run_dirs):
        raise SystemExit(f"--labels has {len(a.labels)} entries for {len(a.run_dirs)} run directories")
    orth = ortholog_names(a.orthologs, a.id_column, a.name_column) if a.orthologs else {}
    settings = Settings(gene_sets=a.gene_sets, set_order=tuple(a.set_order),
                        label_sets=tuple(a.label_sets), label_count=a.label_count,
                        label_genes=tuple(a.label_genes), min_reads=a.min_reads,
                        label_min_reads=a.label_min_reads, direction=a.direction)

    panels = []
    for i, run_dir in enumerate(a.run_dirs):
        results, summary = load_run(run_dir)
        label = a.labels[i] if a.labels else panel_label(run_dir, summary)
        results["comparison_id"] = label
        t = comet_table(results, settings, orthologs=orth, comparisons=None)
        t["background"] = label
        panels.append(t)
    table = pd.concat(panels, ignore_index=True)

    a.output.mkdir(parents=True, exist_ok=True)
    table.to_csv(a.output / "comet_table.csv", index=False)
    (a.output / "comet_run_metadata.json").write_text(json.dumps(dict(
        run_dirs=[str(r) for r in a.run_dirs], genes=int(table.feature_id.nunique()),
        comparisons=list(dict.fromkeys(table.comparison_id)),
        gene_sets=str(a.gene_sets) if a.gene_sets else None,
        orthologs=str(a.orthologs) if a.orthologs else None,
        set_counts={k: int(v) for k, v in table.gene_set.value_counts().items()},
        settings={k: (list(v) if isinstance(v, tuple) else str(v) if isinstance(v, Path) else v)
                  for k, v in settings.__dict__.items()}), indent=2) + "\n")
    draw(table, a.output / "comet_plot", title=a.title)
    print(f"{a.output}/comet_table.csv  ({len(table)} rows, "
          f"{table.feature_id.nunique()} genes, {table.comparison_id.nunique()} comparisons)")
    for k, v in table.gene_set.value_counts().items():
        print(f"  {k:16s} {v}")
    print(f"{a.output}/comet_plot.png / .pdf")


if __name__ == "__main__":
    raise SystemExit(main())
