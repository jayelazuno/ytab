import numpy as np
from scipy import stats

from ytab.essentiality import pooled_gene_ql as pq


def experiment(seed, P=4, genes=600, sites=20, phi=.3, effect=None, effect_pools=None):
    """Site-level NB counts: gene g's sites share a pool-specific representation; treated means x effect[g]."""
    rng = np.random.default_rng(seed); n = genes*sites
    gene_level = np.repeat(np.exp(rng.normal(np.log(30), 1., genes)), sites)
    mocks, treated = [], []
    for p in range(P):
        pool = gene_level*np.repeat(np.exp(rng.normal(0, .4, genes)), sites)*rng.gamma(2, .5, n)
        factor = np.ones(n) if effect is None or (effect_pools is not None and p not in effect_pools) else np.repeat(effect, sites)
        draw = lambda mean:rng.negative_binomial(1/phi, 1/(1+phi*mean))
        mocks.append(draw(pool)); treated.append(draw(pool*factor))
    layout = [dict(id=f"g{g}", name=f"g{g}", description="", chrom="c", strand="+") for g in range(genes)]
    return np.column_stack(mocks), np.column_stack(treated), layout, [np.arange(g*sites, (g+1)*sites) for g in range(genes)]


def test_glm_matches_closed_form_for_poisson_limit():
    y = np.array([[10., 20., 5., 15.]]); X = pq.design_matrix(2)
    beta, mu, _ = pq.fit_glm(y, np.zeros(4), X, 1e-8)
    # saturated pool intercepts with a shared log fold change: ratio of treated to mock totals at phi -> 0
    assert np.isclose(np.exp(beta[0, -1]), (5+15)/(10+20), rtol=1e-4)
    assert np.allclose(mu.sum(), y.sum())


def test_squeeze_variances_recovers_prior_df():
    rng = np.random.default_rng(1); df = 3; s2 = stats.chi2.rvs(df, size=20000, random_state=rng)/df*stats.invgamma.rvs(5, size=20000, random_state=rng)*4
    _, d0, _ = pq.squeeze_variances(s2, df, np.zeros(len(s2)))
    assert 6 < d0 < 14   # inverse-gamma(5) prior on variance = scaled inverse chi-square with 10 df


def test_null_is_calibrated_for_both_count_types():
    for counts in ("reads", "clones"):
        p = []
        for seed in range(3):
            m, t, genes, idx = experiment(seed)
            rows, _ = pq.analyze(m, t, np.ones((4, 2)), genes, idx, pq.Settings(counts=counts, min_rpm=50.))
            p += [r["p_value"] for r in rows if r["in_family"]]
        p = np.array(p); assert .02 < np.mean(p <= .05) < .075, (counts, np.mean(p <= .05))


def test_replicated_depletion_is_called_and_single_pool_artifact_is_not():
    effect = np.ones(600); effect[:60] = .1
    m, t, genes, idx = experiment(10, effect=effect)
    rows, _ = pq.analyze(m, t, np.ones((4, 2)), genes, idx, pq.Settings(counts="reads"))
    called = {r["ORF"] for r in rows if r["call"] == "depleted (replicated)"}
    assert len(called & {f"g{g}" for g in range(60)}) >= 40 and len(called - {f"g{g}" for g in range(60)}) <= 6   # BH: FDR controlled in expectation
    m, t, genes, idx = experiment(11, effect=np.where(np.arange(600) < 60, .02, 1.), effect_pools={0})
    rows, _ = pq.analyze(m, t, np.ones((4, 2)), genes, idx, pq.Settings(counts="reads"))
    assert sum(r["call"] == "depleted (replicated)" for r in rows[:60]) <= 3


def test_enrichment_direction_and_filter():
    effect = np.ones(600); effect[:60] = 8.
    m, t, genes, idx = experiment(12, effect=effect)
    rows, _ = pq.analyze(m, t, np.ones((4, 2)), genes, idx, pq.Settings(counts="reads", direction="enriched"))
    assert sum(r["call"] == "enriched (replicated)" for r in rows[:60]) >= 40
    m[:20] = 0   # gene 0 has no mock reads in any pool
    rows, _ = pq.analyze(m, t, np.ones((4, 2)), genes, idx, pq.Settings(counts="reads"))
    assert not rows[0]["in_family"] and rows[0]["call"] == "not tested"
