"""Small numeric helpers shared across YTAB's analyses: Benjamini-Hochberg correction and JSON-safe conversion of
numpy scalars for manifests."""

import numpy as np


def benjamini_hochberg(p):
    p = np.asarray(p, dtype=float); order = np.argsort(p); ranked = p[order]*len(p)/np.arange(1, len(p)+1)
    adjusted = np.empty_like(p); adjusted[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return adjusted


def json_ready(value):
    if isinstance(value, dict):return {str(k):json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):return [json_ready(v) for v in value]
    if isinstance(value, (np.bool_,)):return bool(value)
    if isinstance(value, np.integer):return int(value)
    if isinstance(value, (float, np.floating)):return float(value) if np.isfinite(value) else None
    return value
