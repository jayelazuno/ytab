"""MidLC common-saturation normalization: thin libraries to the same depth relative to their own complexity.

MidLC_j is the read depth at which library j is expected to show half of its occupied sites. Under binomial thinning
with retention p a site with c reads stays occupied with probability 1-(1-p)^c, so MidLC_j = p_half * total_j where
sum_s [1-(1-p_half)^c_s] = occupied_j / 2 (solved exactly by bisection, no subsampling grid). The saturation ratio
total_j / MidLC_j says how far past half-saturation a library was sequenced. Every library in a normalization group is
thinned to the group's smallest ratio, reads_target_j = min_ratio * MidLC_j, by binomial thinning with
p_j = reads_target_j / total_j (<= 1): libraries are only ever downsampled, the least saturated one is unchanged, and
the counts stay whole numbers. Follows the depth correction of Gale et al. 2020 with a data-derived common target.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse

SEED = 1


def expected_occupied(counts: np.ndarray, p: float) -> float:
    return float(np.sum(-np.expm1(counts*np.log1p(-p)))) if p < 1 else float(len(counts))


def midlc(counts: np.ndarray) -> tuple[float, float]:
    """(MidLC in reads, retention fraction at half-saturation) for one library's positive site counts."""
    counts = np.asarray(counts, dtype=float); counts = counts[counts > 0]
    half = len(counts)/2; lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = (lo+hi)/2
        lo, hi = (mid, hi) if expected_occupied(counts, mid) < half else (lo, mid)
    return hi*counts.sum(), hi


def thin(matrix: sparse.csc_matrix, groups: list, seed: int = SEED) -> tuple[sparse.csc_matrix, list[dict]]:
    matrix = sparse.csc_matrix(matrix, copy=True); n = matrix.shape[1]
    total = np.asarray(matrix.sum(axis=0)).ravel()
    est = [midlc(matrix.data[matrix.indptr[j]:matrix.indptr[j+1]]) for j in range(n)]
    ratio = np.array([total[j]/est[j][0] for j in range(n)])
    rng = np.random.default_rng(seed); out = matrix.astype(np.int64); details = []
    for j, group in enumerate(groups):
        members = [i for i, g in enumerate(groups) if g == group]
        target_ratio = ratio[members].min(); p = min(1.0, target_ratio/ratio[j])
        block = slice(out.indptr[j], out.indptr[j+1])
        if p < 1:out.data[block] = rng.binomial(out.data[block], p)
        details.append(dict(midlc_reads=float(est[j][0]), saturation_ratio=float(ratio[j]), group_target_ratio=float(target_ratio),
                            reads_target=float(target_ratio*est[j][0]), thinning_p=float(p), thinned_total=int(out.data[block].sum()),
                            thinned_occupied_sites=int(np.count_nonzero(out.data[block]))))
    out.eliminate_zeros()
    return out, details
