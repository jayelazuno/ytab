"""Deterministic, zero-preserving positive-count quantile adaptation for Hermes.

Not a literal port of TRANSIT QuantileNorm's dense-rank lookup. Each explicitly
selected library group defines a mean positive empirical quantile function.
Zeros remain undetected; tied counts share the mean target over their ranks.
"""
import numpy as np
from scipy import sparse


def normalize_positive_quantiles(matrix, groups):
    matrix = sparse.csc_matrix(matrix, dtype=float, copy=True)
    matrix.eliminate_zeros()
    if matrix.shape[1] != len(groups) or np.any(~np.isfinite(matrix.data)) or np.any(matrix.data < 0):
        raise ValueError("Finite nonnegative counts and one explicit group per sample required")
    output = matrix.copy()
    for group in dict.fromkeys(groups):
        members = [i for i, value in enumerate(groups) if value == group]
        values = [matrix.data[matrix.indptr[i]:matrix.indptr[i+1]] for i in members]
        if any(len(x) == 0 for x in values):
            raise ValueError("Empty samples cannot define positive-count quantiles")
        if len(members) == 1:
            continue
        sorted_values = [np.sort(x) for x in values]
        for sample, original in zip(members, values):
            order = np.argsort(original, kind="stable")
            ordered = original[order]
            # Quantile knots: endpoints for n>1; singleton is a constant function.
            q = np.linspace(0, 1, len(ordered)) if len(ordered) > 1 else np.array([.5])
            target = np.zeros(len(q), dtype=float)
            for reference in sorted_values:
                target += (np.interp(q, np.linspace(0, 1, len(reference)), reference)
                    if len(reference) > 1 else reference[0]) / len(members)
            # Average the reference values across every tied input rank block.
            starts = np.r_[0, np.flatnonzero(np.diff(ordered)) + 1]
            stops = np.r_[starts[1:], len(ordered)]
            target = np.repeat(np.add.reduceat(target, starts) / (stops-starts), stops-starts)
            transformed = np.empty(len(target))
            transformed[order] = target
            output.data[output.indptr[sample]:output.indptr[sample+1]] = transformed
    return output.tocsr()
