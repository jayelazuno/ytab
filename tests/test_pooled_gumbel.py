import numpy as np
import pytest
from scipy.stats import chi2

from ytab.essentiality.control_conditioned_gumbel import synthetic_experiment
from ytab.essentiality.pooled_gumbel import Settings, analyze, partial_conjunction, required_pools


def pools(seed, effects=None, genes=80, sites=30, P=4, **model):
    rng = np.random.default_rng(seed); mocks, treated = [], []
    for j in range(P):
        C, T, gene_list, idx, pos = synthetic_experiment(rng, np.full(genes, sites), np.array([1., 1.1]), effects=effects, treatments=1, **model)
        mocks.append(C[:, 0]); treated.append(T)
    return np.column_stack(mocks), np.column_stack(treated), np.column_stack([np.ones(P), np.full(P, 1.1)]), gene_list, idx, pos


def test_partial_conjunction_and_majority_rule():
    assert [required_pools(n, Settings()) for n in (2, 3, 4, 5, 6)] == [2, 2, 3, 3, 4]
    assert required_pools(4, Settings(min_pools=4)) == 4 and required_pools(2, Settings(min_pools=3)) == 2
    p = [1e-6, 1e-5, .3, .8]
    assert partial_conjunction(p, 1) == pytest.approx(chi2.sf(-2*np.log(p).sum(), 8))
    assert partial_conjunction(p, 3) == pytest.approx(chi2.sf(-2*np.log([.3, .8]).sum(), 4))
    assert partial_conjunction(p, 3) > .1 > partial_conjunction(p, 2)   # two strong pools do not establish three
    with pytest.raises(ValueError):Settings(direction="up").validate()


def test_null_is_calibrated_and_replicated_depletion_detected():
    M, T, E, genes, idx, pos = pools(1)
    rows, diagnostics = analyze(M, T, E, genes, idx, pos)
    assert diagnostics["calls"]["depleted run (replicated)"] == 0
    for fit in diagnostics["pool_fits"].values():
        assert abs(fit["flagged_sites"]-fit["expected_flagged"]) < 5*np.sqrt(fit["expected_flagged"])
    effects = np.ones(80*30)
    for g in range(0, 80, 8):effects[g*30+5:g*30+25] = 0.
    M, T, E, genes, idx, pos = pools(2, effects=effects)
    rows = analyze(M, T, E, genes, idx, pos)[0]
    called = {r["ORF"] for r in rows if r["call"] == "depleted run (replicated)"}; truth = {genes[g]["id"] for g in range(0, 80, 8)}
    assert len(called & truth) >= 8 and len(called-truth) <= 1


def test_single_pool_artifact_is_not_called_but_reported():
    M, T, E, genes, idx, pos = pools(3)
    T = T.copy(); T[idx[5][5:25], 0] = 0   # complete loss in one pool only (e.g. a clone-specific artifact)
    row = {r["ORF"]:r for r in analyze(M, T, E, genes, idx, pos)[0]}[genes[5]["id"]]
    assert row["pools_p_le_05"] == 1 and row["call"] == "not significant" and row["required_pools"] == 3


def test_enrichment_direction():
    effects = np.ones(80*30)
    for g in range(0, 80, 8):effects[g*30+5:g*30+25] = 6.
    M, T, E, genes, idx, pos = pools(4, effects=effects)
    rows = analyze(M, T, E, genes, idx, pos, Settings(direction="enriched"))[0]
    called = {r["ORF"] for r in rows if r["call"] == "enriched run (replicated)"}
    assert len(called & {genes[g]["id"] for g in range(0, 80, 8)}) >= 7
