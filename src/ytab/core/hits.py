"""Reading YTAB hit files: the per-coordinate insertion counts a mapping run writes."""

import csv
from pathlib import Path

import numpy as np


def read_hits(path, min_reads=1):
    """{chromosome: (sorted positions, reads)} from a YTAB *_hits.txt or a table with chromosome/position/reads columns.
    Counts at one position (both strands) are summed."""
    path = Path(path); text = path.read_text().splitlines()
    delimiter = "\t" if "\t" in text[0] else ","
    rows = csv.DictReader(text, delimiter=delimiter); fields = {f.lower().strip():f for f in rows.fieldnames}
    pick = lambda *names:next((fields[n] for n in names if n in fields), None)
    chrom, pos, reads = pick("chromosome", "chrom", "chr"), pick("hit position", "position", "pos"), pick("hit count", "reads", "count", "counts")
    if not (chrom and pos and reads):raise ValueError(f"{path}: need chromosome, position and read-count columns (YTAB hit file or chromosome/position/reads)")
    totals = {}
    for r in rows:
        key = (r[chrom], int(float(r[pos]))); totals[key] = totals.get(key, 0)+int(float(r[reads]))
    out = {}
    for (c, p), n in totals.items():
        if n >= min_reads:out.setdefault(c, []).append((p, n))
    return {c:(np.array([p for p, _ in v]), np.array([n for _, n in v], dtype=float)) for c, v in ((c, sorted(v)) for c, v in out.items())}
