"""Top-site treatment of a library's raw counts, applied before any normalization.

A gene's largest insertion site can carry most of its reads without many mutants behind it: a PCR jackpot, an
insertion that does not disrupt the gene, a selection marker or the donor locus. Within each gene's coding sequence:
  none       counts unchanged (default)
  winsorize  the largest site is lowered to the gene's second-largest count (a gene with one site is unchanged);
             adapted from TRANSIT's -winz option, per library
  drop-top   the largest site is removed (Nickels et al. 2025); a gene with one site loses it
Each library is treated on its own, so treated and control libraries get the same rule. Counts stay integers.
"""
from __future__ import annotations

FILTERS = ("none", "winsorize", "drop-top")


def gene_groups(genes: list[dict], sites: list[tuple]) -> list:
    """Site indices of each gene's coding sequence (whole CDS, last 3 bases excluded)."""
    from ytab.essentiality.control_conditioned_gumbel import gene_sites
    from ytab.essentiality.library_correlation import SITE_RULES
    return [idx for idx in gene_sites(genes, sites, SITE_RULES) if len(idx)]


def apply(counts: dict, genes: list[dict], how: str) -> tuple[dict, dict]:
    """counts: (chromosome, position) -> reads of one library. Returns (treated counts, summary)."""
    if how not in FILTERS:raise ValueError(f"Unknown site filter {how!r}; choose one of {', '.join(FILTERS)}")
    summary = dict(site_filter=how, genes_changed=0, reads_removed=0)
    if how == "none" or not counts:return dict(counts), summary
    keys = list(counts); out = dict(counts)
    for idx in gene_groups(genes, keys):
        values = sorted(((counts[keys[i]], keys[i]) for i in idx), reverse=True)
        top, key = values[0]
        if how == "winsorize":
            if len(values) < 2 or values[1][0] == top:continue
            new = values[1][0]
        else:
            new = 0
        if top <= 0:continue
        summary["genes_changed"] += 1; summary["reads_removed"] += int(top-new)
        if new:out[key] = new
        else:out.pop(key, None)
    return out, summary
