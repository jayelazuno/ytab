"""Principal components of libraries from a gene x library matrix (QC and normalization plots).

Values are log2(1 + x), the genes restricted to those with a value in at least one library and then to the most variable
(500 by default, as DESeq2's plotPCA), each gene centred (not scaled), and the libraries projected on the components.
Writes nothing; callers write the rows.
"""
from __future__ import annotations

import numpy as np

N_COMPONENTS = 4
TOP_GENES = 500


def pca(values: np.ndarray, top: int = TOP_GENES, components: int = N_COMPONENTS) -> tuple[np.ndarray, np.ndarray, int] | None:
    """values: genes x libraries (non-negative). Returns (scores libraries x k, fraction of variance per component,
    genes used), or None when fewer than 3 libraries or 3 usable genes."""
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or x.shape[1] < 3:return None
    x = np.log2(1+x[np.isfinite(x).all(axis=1) & (x > 0).any(axis=1)])
    if len(x) < 3:return None
    x = x[np.argsort(x.var(axis=1))[::-1][:top]]
    centred = (x-x.mean(axis=1, keepdims=True)).T   # libraries x genes
    u, s, _ = np.linalg.svd(centred, full_matrices=False)
    total = float((s**2).sum()); k = min(components, len(s))
    if total <= 0:return None
    return u[:, :k]*s[:k], (s[:k]**2)/total, len(x)


def rows(values: np.ndarray, libraries: list[dict], group: dict | None = None, top: int = TOP_GENES) -> tuple[list[dict], list[dict]]:
    """(score rows: one per library with its metadata and PC1..PCk, variance rows: pc, variance_explained, genes).
    `group` (e.g. {"method": "ttr"}) is added to every row."""
    result = pca(values, top)
    if result is None:return [], []
    scores, variance, genes = result; group = group or {}
    score_rows = [group | lib | {f"PC{k+1}":round(float(v), 5) for k, v in enumerate(scores[i])} for i, lib in enumerate(libraries)]
    variance_rows = [group | dict(pc=f"PC{k+1}", variance_explained=round(float(v), 5), genes=genes) for k, v in enumerate(variance)]
    return score_rows, variance_rows
