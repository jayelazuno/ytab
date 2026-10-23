"""Hybrid essentiality inputs from any normalization method.

The hybrid's gene features come from parent hit files that are combined, summarized (SummaryTable) and scored by
the classifier. MidLC writes such hit files itself (sample_normalization/<target tag>/). For any other method, this
module normalizes the parents with it (from a Site Normalization run's parent catalogue) and writes the same files
under sample_normalization/<method>__<run>__<parents>/, so the combine, combined-summary and classifier steps treat
every method alike.

SummaryTable needs whole-number counts. A normalized count is rounded to the nearest whole number, and a site with
insertions in the raw hit file keeps at least 1, so site-based features (hits, freedom and neighbourhood indices,
longest free interval) match the data and read-based features carry the method's scaling.

MidLC is different: it thins reads to a common saturation, so its counts are already whole numbers and a site thinned
to zero is meant to disappear. Its parents are normalized together (one group), and each site's thinned count is
shared between its strand rows by largest remainder, so site totals are exact and zero sites are left out.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[3]
SOURCE_FILE = "normalization_source.json"
ROUNDING = "nearest whole number; sites with raw insertions kept at >= 1"
THINNING = "thinned counts kept as whole numbers; strands share a site's count by largest remainder; thinned-out sites dropped"
THINNING_METHODS = {"midlc"}


def method_tag(method: str, run_dir: Path, samples: list[str]) -> str:
    """<method>__<normalization run>__<parent set>: different parents give different inputs, so they get their own tag."""
    parents = hashlib.sha1(",".join(sorted(samples)).encode()).hexdigest()[:6]
    return f"{method}__{Path(run_dir).name}__{parents}"


def catalogue_runs(project_root: Path, observations: list[str]) -> list[Path]:
    """Site normalization runs whose parent catalogue holds every observation (project::sample), newest first."""
    root = Path(project_root)/"normalization"; found = []
    for run in (sorted(root.glob("*"), key=lambda p:p.stat().st_mtime, reverse=True) if root.is_dir() else []):
        table = run/"catalogue/observations.tsv"
        if (run/"catalogue/raw_counts.npz").is_file() and table.is_file():
            present = {row["observation_id"] for row in csv.DictReader(open(table), delimiter="\t")}
            if set(observations) <= present:found.append(run)
    return found


def round_keep_observed(values: np.ndarray) -> np.ndarray:
    return np.maximum(np.floor(np.asarray(values, dtype=float)+.5), 1).astype(np.int64)


def share_by_largest_remainder(index: np.ndarray, row_raw: np.ndarray, site_counts: np.ndarray) -> np.ndarray:
    """Whole-number shares of each site's count among its rows, proportional to the rows' raw reads; sums are exact."""
    share = site_counts[index]*row_raw/np.bincount(index, weights=row_raw, minlength=len(site_counts))[index]
    values = np.floor(share).astype(np.int64)
    short = (site_counts - np.bincount(index, weights=values, minlength=len(site_counts))).astype(np.int64)
    order = np.lexsort((-(share-values), index))
    rank = np.empty(len(index), dtype=np.int64); seen = {}
    for i in order:
        rank[i] = seen.get(index[i], 0); seen[index[i]] = rank[i]+1
    return values + (rank < short[index])


def _column(matrix: sparse.spmatrix, j: int) -> np.ndarray:
    return np.asarray(matrix.tocsc()[:, j].todense()).ravel()


def _observation_columns(path: Path) -> dict[str, int]:
    return {row["sample"]:i for i, row in enumerate(csv.DictReader(open(path), delimiter="\t"))}


def normalize_parents(project_root: Path, run_dir: Path, method: str, project_id: str, samples: list[str], raw_hits: dict[str, Path],
                      force: bool = False, python: str = sys.executable) -> dict:
    """Normalize the parents with the method (from the run's catalogue) and write <sample>_normalized_hits.txt.

    A sample's raw hit file, summed over strands, must match the catalogue's raw count at every site; otherwise the
    catalogue was built from different hit files and nothing is written."""
    run_dir = Path(run_dir); tag = method_tag(method, run_dir, samples)
    target_dir = Path(project_root)/"sample_normalization"/tag; source_path = target_dir/SOURCE_FILE
    if source_path.is_file() and not force:return json.loads(source_path.read_text())
    catalogue = run_dir/"catalogue"; engine = target_dir/"engine"
    if engine.exists():shutil.rmtree(engine)
    command = [python, str(ROOT/"scripts/local/ytab_normalize.py"), "--catalogue", str(catalogue), "--output", str(engine),
               "--method", method, "--group-by", "cohort" if method in THINNING_METHODS else "library_id", *[x for s in samples for x in ("--observation", f"{project_id}::{s}")]]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:raise ValueError(f"{method} normalization failed: "+((result.stderr or result.stdout).strip().splitlines() or ["no output"])[-1])
    with gzip.open(catalogue/"sites.tsv.gz", "rt") as handle:
        sites = {(row["chromosome"], int(row["position_1based"])):int(row["site_index"]) for row in csv.DictReader(handle, delimiter="\t")}
    raw = sparse.load_npz(catalogue/"raw_counts.npz"); normalized = sparse.load_npz(engine/method/"normalized_counts.npz")
    raw_columns, norm_columns = _observation_columns(catalogue/"observations.tsv"), _observation_columns(engine/"observations.tsv")
    totals = {}
    for sample in samples:
        raw_col, norm_col = _column(raw, raw_columns[sample]), _column(normalized, norm_columns[sample])
        with open(raw_hits[sample]) as handle:
            reader = csv.reader(handle, delimiter="\t"); header = next(reader); rows = list(reader)
        pos, count = header.index("Hit position"), header.index("Hit count")
        index = np.array([sites.get((r[0], int(r[pos])), -1) for r in rows]); row_raw = np.array([int(r[count]) for r in rows])
        site_raw = np.bincount(index[index >= 0], weights=row_raw[index >= 0], minlength=len(raw_col))
        if (index < 0).any() or not np.array_equal(site_raw[index], raw_col[index]):
            raise ValueError(f"{sample}: its hit file does not match the catalogue of {run_dir.name}; rerun Site Normalization")
        # A site's W and C strands are separate rows; the site's normalized count is shared in proportion to their raw reads.
        row_norm = norm_col[index]*row_raw/raw_col[index]
        thinning = method in THINNING_METHODS
        values = share_by_largest_remainder(index, row_raw, norm_col.astype(np.int64)) if thinning else round_keep_observed(row_norm)
        out = target_dir/sample/f"{sample}_normalized_hits.txt"; out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(header)
            for r, v in zip(rows, values):
                if v or not thinning:writer.writerow([*r[:count], int(v), *r[count+1:]])
        totals[sample] = dict(rows=len(rows), sites=int(len(set(index.tolist()))), raw_reads=int(row_raw.sum()),
                              normalized_reads=round(float(row_norm.sum()), 3), written_reads=int(values.sum()),
                              raised_to_one=0 if thinning else int((row_norm < .5).sum()), rows_dropped=int((values == 0).sum()) if thinning else 0)
    source = dict(method=method, normalization_run=str(run_dir), tag=tag, samples=list(samples), rounding=THINNING if method in THINNING_METHODS else ROUNDING,
                  command=command, created=datetime.now(timezone.utc).isoformat(), totals=totals)
    source_path.write_text(json.dumps(source, indent=2)+"\n")
    return source


def read_source(project_root: Path, tag: str) -> dict | None:
    path = Path(project_root)/"sample_normalization"/tag/SOURCE_FILE
    return json.loads(path.read_text()) if path.is_file() else None
