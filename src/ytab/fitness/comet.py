"""Comet plots: a gene's reads under treatment against its reads in the control, one point per gene.

The name is descriptive. Genes that do not respond sit on the diagonal, and because read totals span orders of
magnitude the bulk of them forms a dense elongated body along it; the genes that do respond spray away from that
body, so the figure reads as a comet with a tail. It is the view that shows, in one panel per comparison, both how
much evidence a gene carries (how far along the diagonal it sits) and how much it moved (how far off the diagonal),
which a volcano plot separates and a bare hit list discards.

Three things the plot needs that a results table does not already carry:

  gene sets      the overlay colours. Sets are declared in GMT files (the MSigDB interchange format: one set per
                 line, name, description, then members), so a lab can supply its own without touching code. Members
                 are matched against the gene's own identifier and name AND against its ortholog's name in a
                 reference species, because in a non-model species most genes carry no name of their own: on
                 C. glabrata 4,812 of 5,234 genes reach an S. cerevisiae ortholog while only ~746 have a name in
                 the reference annotation, so sets declared on ortholog names cover the genome and sets declared on
                 native names cover a seventh of it. Ortholog resolution is therefore how the sets become usable,
                 not a convenience.

  set priority   a gene may be in several sets. Rather than draw it twice, the first set in the declared order wins
                 and the assignment is recorded, so the figure cannot silently depend on dictionary ordering.

  labels         which points carry a name. Chosen by distance from the diagonal within the highlighted sets, or
                 given explicitly, never by a threshold on the statistic: the plot is read for where genes sit, and
                 a label rule tied to a call would hide a gene that moved a long way without reaching significance.

Nothing here computes a statistic. The table this writes is derived entirely from an existing treated-vs-control
run, and the plot is drawn from the table, so the figure and the numbers cannot disagree.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import numpy as np
import pandas as pd

PSEUDOCOUNT = 1.0          # so a gene with no reads in one condition still plots on a log axis
UNASSIGNED = "other"       # the set name used for genes in none of the declared sets


@dataclass(frozen=True)
class Settings:
    """How the comet table is built. Every field is recorded in the run metadata."""
    gene_sets: Path | None = None            # GMT file; None -> every gene is UNASSIGNED
    set_order: tuple[str, ...] = ()          # priority, first match wins; () -> GMT file order
    label_sets: tuple[str, ...] = ()         # sets eligible for labels; () -> all declared sets
    label_count: int = 12                    # labels per gene set, by distance from the diagonal
    direction: str = "both"                  # both | depleted | enriched -- which hits may be labelled
    label_min_reads: float = 10.0            # a gene must reach this in EITHER condition to be labelled
    label_genes: tuple[str, ...] = ()        # explicit names/ids; these are always labelled
    min_reads: float = 0.0                   # drop genes below this in BOTH conditions
    pseudocount: float = PSEUDOCOUNT


def read_gmt(path: Path) -> dict[str, dict]:
    """Parse a GMT file: name, description, then members, tab separated. Blank lines and # comments ignored."""
    sets: dict[str, dict] = {}
    for line in Path(path).read_text().splitlines():
        line = line.rstrip("\r")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 3:
            raise ValueError(f"GMT line needs name, description and at least one member: {line[:60]!r}")
        name, description, members = parts[0], parts[1], [p for p in parts[2:] if p]
        if not name:
            raise ValueError("GMT set name may not be empty")
        if name in sets:
            raise ValueError(f"GMT declares set {name!r} twice")
        # Matching is case-insensitive because gene names are written inconsistently across annotations
        # (Erg11, ERG11), while the identifiers they match against are not.
        sets[name] = dict(description=description, members={m.upper() for m in members})
    if not sets:
        raise ValueError(f"No gene sets found in {path}")
    return sets


def ortholog_names(table: Path, id_column: str, name_column: str) -> dict[str, str]:
    """Map a gene identifier to its ortholog's gene name, from an orthology table."""
    o = pd.read_csv(table, sep="\t", dtype=str, keep_default_na=False)
    for column in (id_column, name_column):
        if column not in o.columns:
            raise ValueError(f"Orthology table has no column {column!r}: {table}")
    out: dict[str, str] = {}
    for gene, name in zip(o[id_column], o[name_column]):
        gene, name = gene.strip(), name.strip()
        # "NA" and "protein_coding" are placeholders these tables use for "no name", not names.
        if gene and name and name.upper() not in {"NA", "PROTEIN_CODING"} and gene not in out:
            out[gene] = name
    return out


def assign_sets(frame: pd.DataFrame, sets: dict[str, dict], order: tuple[str, ...]) -> pd.Series:
    """Label each gene with the first declared set it belongs to, matching on any of its names."""
    names = order or tuple(sets)
    unknown = [n for n in names if n not in sets]
    if unknown:
        raise ValueError(f"set_order names sets that are not in the GMT: {unknown}")
    # Each gene's candidate keys: its id, its own name, and its ortholog's name.
    keys = [
        {str(v).upper() for v in row if v and str(v) != "nan"}
        for row in zip(frame["feature_id"], frame["display_name"], frame["ortholog_name"])
    ]
    assigned = []
    for candidates in keys:
        for name in names:
            if candidates & sets[name]["members"]:
                assigned.append(name)
                break
        else:
            assigned.append(UNASSIGNED)
    return pd.Series(assigned, index=frame.index, dtype="object")


def choose_labels(frame: pd.DataFrame, settings: Settings) -> pd.Series:
    """Mark the genes whose names are drawn, by distance from the diagonal within the eligible sets."""
    explicit = {g.upper() for g in settings.label_genes}
    label = pd.Series(False, index=frame.index)
    if explicit:
        for column in ("feature_id", "display_name", "ortholog_name"):
            label |= frame[column].astype(str).str.upper().isin(explicit)
    eligible = set(settings.label_sets) or (set(frame["gene_set"]) - {UNASSIGNED})
    if settings.label_count > 0 and eligible:
        for _, block in frame.groupby("comparison_id", sort=False):
            pool = block[block["gene_set"].isin(eligible)]
            # A gene sitting at one or two reads in both conditions has a large off-diagonal distance for no
            # reason but counting noise, and those genes pile up at the origin where their labels overlap each
            # other. Requiring real support in at least one condition removes the pile-up and the noise together.
            pool = pool[(pool["control_reads"] >= settings.label_min_reads) |
                        (pool["treated_reads"] >= settings.label_min_reads)]
            # Which side of the diagonal may be named, matching the app's "Calls shown". Every set member stays
            # drawn whichever direction is chosen: the overlay exists to show where the whole family sits.
            if settings.direction == "depleted":
                pool = pool[pool["off_diagonal"] < 0]
            elif settings.direction == "enriched":
                pool = pool[pool["off_diagonal"] > 0]
            # Top hits **per set**: with one pool a large family crowds out every small one, so a 90-gene
            # ribosome set would take every label and leave none for the 29-gene ergosterol set beside it.
            for name, members in pool.groupby("gene_set", sort=False):
                top = members["off_diagonal"].abs().nlargest(settings.label_count).index
                label.loc[top] = True
    return label


def comet_table(results: pd.DataFrame, settings: Settings, *,
                orthologs: dict[str, str] | None = None,
                comparisons: pd.DataFrame | None = None) -> pd.DataFrame:
    """Build the plot-ready table from a treated-vs-control results frame.

    One row per gene per comparison, carrying everything the figure draws and nothing it does not.
    """
    required = {"feature_id", "control_mean", "treated_mean"}
    missing = required - set(results.columns)
    if missing:
        raise ValueError(f"Results table missing columns {sorted(missing)}")

    d = results.copy()
    d["display_name"] = (d["common_name"] if "common_name" in d else d["feature_id"]).astype(str)
    # A display name that is just the locus id is not a name; the plot should fall back to the ortholog.
    d.loc[d["display_name"].str.strip() == "", "display_name"] = d["feature_id"]
    d["ortholog_name"] = d["feature_id"].map(orthologs or {}).fillna("")

    if "comparison_id" not in d.columns:
        d["comparison_id"] = "all"
    if comparisons is not None and "background" in comparisons.columns:
        lookup = dict(zip(comparisons["comparison_id"], comparisons["background"]))
        d["background"] = d["comparison_id"].map(lookup).fillna("")
    elif "background" not in d.columns:
        d["background"] = ""

    c = settings.pseudocount
    d["control_reads"] = pd.to_numeric(d["control_mean"], errors="coerce").fillna(0.0)
    d["treated_reads"] = pd.to_numeric(d["treated_mean"], errors="coerce").fillna(0.0)
    if settings.min_reads > 0:
        keep = (d["control_reads"] >= settings.min_reads) | (d["treated_reads"] >= settings.min_reads)
        d = d[keep].copy()
    d["x"] = d["control_reads"] + c
    d["y"] = d["treated_reads"] + c
    # The plot's own quantity: vertical distance from the diagonal on log axes. Recomputed here from the
    # plotted coordinates rather than taken from the run's log2FC, so a point can never sit off the diagonal
    # while its label says it did not move.
    d["off_diagonal"] = np.log2(d["y"]) - np.log2(d["x"])

    sets = read_gmt(settings.gene_sets) if settings.gene_sets else {}
    d["gene_set"] = assign_sets(d, sets, settings.set_order) if sets else UNASSIGNED
    d["call"] = d["final_call"].astype(str) if "final_call" in d.columns else ""
    d["label"] = choose_labels(d, settings)
    d["label_name"] = np.where(d["label"],
                               np.where(d["ortholog_name"].astype(bool), d["ortholog_name"], d["display_name"]),
                               "")

    columns = ["comparison_id", "background", "feature_id", "display_name", "ortholog_name",
               "control_reads", "treated_reads", "x", "y", "off_diagonal", "call",
               "gene_set", "label", "label_name"]
    return d[columns].sort_values(["comparison_id", "gene_set", "feature_id"]).reset_index(drop=True)
