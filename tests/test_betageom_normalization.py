import json
import runpy
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from scipy.special import beta as beta_function
from scipy.stats import trim_mean

from ytab.normalization.betageom import beta_geometric_cdf, correct, fit_parameters
from ytab.normalization.methods import normalize, run_choices
from ytab.normalization.ttr import ttr_parameters


def test_cdf_matches_beta_geometric_pmf_and_simulation():
    alpha,beta = 2.,5.
    x = np.arange(30)
    pmf = np.array([0.] + [beta_function(alpha+1,beta+k-1)/beta_function(alpha,beta) for k in x[1:]])
    np.testing.assert_allclose(beta_geometric_cdf(x,alpha,beta), np.cumsum(pmf), atol=1e-14)
    rng = np.random.default_rng(8)
    draws = rng.geometric(rng.beta(alpha,beta,size=200000))
    np.testing.assert_allclose(beta_geometric_cdf([1,3,10,20],alpha,beta), [(draws<=k).mean() for k in [1,3,10,20]],atol=.004)


def test_fit_estimators_match_inspected_transit_source():
    positive = np.array([1,2,3,5,20,100,1000])
    n = 30
    actual = fit_parameters(positive,n)
    ordered = np.sort(np.r_[np.zeros(n-len(positive)),positive])
    rho = max(1/trim_mean(ordered+1,.001),.0001)
    q = np.arange(n)/float(n)
    A = np.sum(np.log(1-q)**2)/np.sum(ordered*np.log(1-q))
    kappa = max((2*np.exp(A)-1)/(np.exp(A)+rho-1),10)
    assert actual["rho"] == pytest.approx(rho)
    assert actual["A"] == pytest.approx(A)
    assert actual["kappa"] == pytest.approx(kappa)
    # The upstream family really is a mixture of scipy's positive-support geom.
    path=Path(__file__).resolve().parents[1]/"docs/codex/transit/src/pytransit/norm_tools.py"
    if path.exists():
        source=runpy.run_path(str(path))
        assert source["BetaGeomNorm"].name == "betageom"
        state=np.random.get_state()
        try:
            np.random.seed(3)
            _,fitted=source["BetaGeomNorm"].normalize(np.vstack([ordered,ordered*2]),bgsamples=1000,doTTR=False)
        finally:
            np.random.set_state(state)
        np.testing.assert_allclose(fitted[0],[actual["rho"],actual["kappa"]],rtol=1e-12)


def test_zeros_raw_counts_ties_and_post_ttr():
    raw = sparse.csr_matrix([[0,1],[1,3],[2,3],[10,20],[1000,600],[0,0]])
    before = raw.copy()
    unscaled,detail = correct(raw,["g","g"],post_ttr=False)
    scaled,params = correct(raw,["g","g"],post_ttr=True)
    assert (raw != before).nnz == 0
    assert np.array_equal(raw.toarray()==0,scaled.toarray()==0)
    assert unscaled.toarray()[1,1] == unscaled.toarray()[2,1]
    for j in range(2):
        positive = unscaled.toarray()[:,j]
        positive = positive[positive>0]
        expected = ttr_parameters(positive,raw.shape[0])["factor"]
        np.testing.assert_allclose(scaled.toarray()[:,j],unscaled.toarray()[:,j]*expected)
        assert params[j]["bgc_post_ttr_factor"] == pytest.approx(expected)
        assert detail[j]["bgc_post_ttr_factor"] == 1
    result, rows = normalize(raw,"betageom",["g","g"])
    np.testing.assert_allclose(result.toarray(),scaled.toarray())
    assert rows[0]["factor"] == "" and rows[0]["bgc_cdf"] == "analytic"


def test_simulation_seed_order_and_analytic_agreement():
    raw = np.array([[0,0],[1,2],[3,5],[10,15],[100,80],[1000,1200]])
    a,detail = correct(raw,["g","g"],cdf="simulation",seed=17,draws=100000,post_ttr=False)
    b,_ = correct(raw,["g","g"],cdf="simulation",seed=17,draws=100000,post_ttr=False)
    c,_ = correct(raw[:,::-1],["g","g"],cdf="simulation",seed=17,draws=100000,post_ttr=False)
    np.testing.assert_array_equal(a.toarray(),b.toarray())
    np.testing.assert_array_equal(a.toarray(),c.toarray()[:,::-1])
    exact,_ = correct(raw,["g","g"],post_ttr=False)
    # Monte Carlo percentile noise matters most at the upper tail; compare scale-relative errors.
    np.testing.assert_allclose(a.toarray(),exact.toarray(),rtol=.08,atol=3)
    assert detail[0]["bgc_draws"] == 100000


def test_groups_and_catalogue_zero_sensitivity_are_explicit():
    raw = np.array([[0,0,0],[1,2,100],[3,5,200],[10,15,800],[100,80,1000]])
    a,details = correct(raw,["a","a","b"],post_ttr=False)
    b,_ = correct(raw[:,:2],["a","a"],post_ttr=False)
    np.testing.assert_allclose(a.toarray()[:,:2],b.toarray())
    assert details[2]["bgc_target_geometric_mean"] != details[0]["bgc_target_geometric_mean"]
    before = fit_parameters([1,3,10,100],5)
    after = fit_parameters([1,3,10,100],50)
    assert before["rho"] != after["rho"]  # Zero-only catalogue changes affect the fit.
    for kwargs in (dict(seed=-1),dict(draws=0),dict(cdf="invalid")):
        with pytest.raises(ValueError):correct(raw,["a"]*3,**kwargs)
    with pytest.raises(ValueError):fit_parameters([],10)
    with pytest.raises(ValueError):correct([[0,1]],["a","a"])


def test_shared_runner_records_parameters_with_other_methods(tmp_path):
    from test_normalization_choices import make_catalogue
    catalogue = tmp_path/"catalogue"
    make_catalogue(catalogue)
    out = tmp_path/"results"
    result = run_choices(catalogue,out,["nonorm","ttr","betageom"],"library_id",bgc_options=dict(post_ttr=False))
    assert result["bgc_parameters"]["post_ttr"] is False
    assert result["bgc_parameters"]["fit_universe"] == "full catalogue including zeros"
    assert result["betageom_implementation_sha256"]
    assert (out/"betageom/normalized_counts.npz").exists()
    assert "bgc_rho" in (out/"method_comparison.tsv").read_text()


def test_audit_cli_checks_real_output_and_seed_reproducibility(tmp_path):
    from test_normalization_choices import make_catalogue
    catalogue = tmp_path/"catalogue"
    make_catalogue(catalogue)
    out = tmp_path/"results"
    run_choices(catalogue,out,["betageom"],"library_id")
    script = Path(__file__).resolve().parents[1]/"scripts/local/ytab_bgc_audit.py"
    audit = runpy.run_path(str(script))["audit"]
    report = audit(out,tmp_path/"audit",seed=9,draws=1000)
    assert report["simulation_seed_reproduction"] and report["zero_pattern_preserved"]
    assert (tmp_path/"audit/diagnostics.tsv").exists()
    with pytest.raises(ValueError,match="new output"):
        audit(out,tmp_path/"audit",draws=1000)
