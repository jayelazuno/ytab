import numpy as np
from scipy import sparse

from ytab.normalization.methods import METHODS, normalize
from ytab.normalization.midlc import midlc
from ytab.pipeline.hybrid_inputs import share_by_largest_remainder


def _libraries(seed=0):
    rng = np.random.default_rng(seed)
    lam = rng.gamma(.5, 20, 20000)
    return np.c_[rng.poisson(lam), rng.poisson(lam*4), rng.poisson(lam*2)]


def test_midlc_is_half_saturation_depth():
    counts = _libraries()[:, 0]; counts = counts[counts > 0]
    reads, p = midlc(counts)
    assert np.isclose(reads, p*counts.sum())
    kept = [np.count_nonzero(np.random.default_rng(s).binomial(counts, p)) for s in range(5)]
    assert abs(np.mean(kept) - len(counts)/2) < .01*len(counts)


def test_common_saturation_downsamples_only_to_least_saturated():
    assert METHODS["midlc"].status == "available"
    m = sparse.csc_matrix(_libraries())
    out, rows = normalize(m, "midlc", ["g", "g", "g"])
    ratios = [r["saturation_ratio"] for r in rows]
    least = int(np.argmin(ratios))
    assert rows[least]["thinning_p"] == 1 and rows[least]["thinned_total"] == rows[least]["raw_total_support"]
    assert all(r["thinning_p"] <= 1 and r["thinned_total"] <= r["raw_total_support"] for r in rows)
    assert np.all(out.data == np.floor(out.data)) and np.all(out.toarray() <= m.toarray())
    # each library is thinned to target ratio x its own MidLC (binomial noise only)
    for r in rows:
        assert np.isclose(r["reads_target"], min(ratios)*r["midlc_reads"])
        assert abs(r["thinned_total"]-r["reads_target"]) < 5*np.sqrt(r["reads_target"])


def test_groups_are_thinned_separately():
    m = sparse.csc_matrix(_libraries())
    _, rows = normalize(m, "midlc", ["a", "b", "b"])
    assert rows[0]["thinning_p"] == 1


def test_strand_shares_sum_to_site_counts():
    index = np.array([0, 0, 1, 2, 2]); row_raw = np.array([3, 1, 5, 1, 1]); site = np.array([3, 0, 1])
    values = share_by_largest_remainder(index, row_raw, site)
    assert np.array_equal(np.bincount(index, weights=values, minlength=3), site)
    assert values[2] == 0 and values.min() >= 0
