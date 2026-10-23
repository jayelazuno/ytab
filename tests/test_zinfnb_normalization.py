import numpy as np
import pytest
from scipy.stats import nbinom
from ytab.normalization.zinfnb import fit, profile
from ytab.normalization.methods import normalize, run_choices


def test_proper_likelihood_matches_direct_mixture():
    values=np.array([1,2,7]); frequencies=np.array([4,3,1]); zeros=12
    loss,pi=profile(values,frequencies,zeros,np.log([4.,2.]))
    p=2/6
    expected=-(zeros*np.log(pi+(1-pi)*nbinom.pmf(0,2,p))+
               np.dot(frequencies,np.log1p(-pi)+nbinom.logpmf(values,2,p)))/20
    assert loss==pytest.approx(expected)
    mass=pi+(1-pi)*nbinom.pmf(np.arange(1000),2,p).sum()
    assert mass==pytest.approx(1.)


def test_synthetic_recovery_determinism_zeros_and_raw_preservation():
    rng=np.random.default_rng(81)
    x=rng.negative_binomial(2,.2,40000)
    x[rng.random(len(x))<.35]=0
    raw=x[:,None].copy()
    result,rows=normalize(raw,"zinfnb",["pool"])
    row=rows[0]
    assert row["zinb_pi"]==pytest.approx(.35,abs=.025)
    assert row["zinb_nb_mean"]==pytest.approx(8,rel=.06)
    assert row["zinb_size"]==pytest.approx(2,rel=.1)
    assert row["zinb_cdf_max_discrepancy"]<.015
    np.testing.assert_array_equal(raw[:,0],x)
    np.testing.assert_array_equal(result.toarray()==0,raw==0)
    np.testing.assert_allclose(result.toarray(),raw*row["factor"])
    assert fit(x[x>0][::-1],len(x))==fit(x[x>0],len(x))


def test_catalogue_sensitivity_group_independence_and_validation():
    x=np.array([1,1,2,3,5,8,13,20]*10)
    a=fit(x,100); b=fit(x,200)
    assert b["zinb_pi"]>a["zinb_pi"]
    _,rows=normalize(np.column_stack([x,x]),"zinfnb",["a","b"])
    assert rows[0]["factor"]==rows[1]["factor"]
    with pytest.raises(ValueError):fit([],10)
    with pytest.raises(ValueError):fit([1.5],2)
    with pytest.raises(ValueError,match="empty samples"):normalize([[0]],"zinfnb",["g"])


def test_failed_optimization_is_not_silent(monkeypatch):
    from types import SimpleNamespace
    from ytab.normalization import zinfnb
    monkeypatch.setattr(zinfnb,"minimize",lambda *a,**k:SimpleNamespace(success=False))
    with pytest.raises(ValueError,match="no fallback"):fit([1,2],4)


def test_nb_boundary_mean_and_parameter_order_invariance():
    # Counts deliberately have fewer zeros than a highly dispersed NB expects.
    x=np.array([1,2,3,4,5,6,7,8]*100)
    params=fit(x,len(x))
    assert params["zinb_pi"]==0
    assert params["zinb_boundary"]
    assert params["zinb_nb_mean"]==pytest.approx(x.mean(),rel=1e-3)
    raw=np.column_stack([x,x*2])
    a,rows=normalize(raw,"zinfnb",["g","g"])
    b,reverse=normalize(raw[:,::-1],"zinfnb",["g","g"])
    np.testing.assert_allclose(a.toarray(),b.toarray()[:,::-1],rtol=1e-12)
    assert rows[0]["factor"]==reverse[1]["factor"]


def test_profiled_zero_weight_is_likelihood_optimum():
    from scipy.optimize import minimize_scalar
    x=np.array([1,2,3,10]); h=np.array([30,20,10,5]); zeros=100
    loss,pi=profile(x,h,zeros,np.log([5.,1.5]))
    p=1.5/6.5
    def objective(w):
        return -(zeros*np.log(w+(1-w)*nbinom.pmf(0,1.5,p))+
                 np.dot(h,np.log1p(-w)+nbinom.logpmf(x,1.5,p)))/165
    optimum=minimize_scalar(objective,bounds=(0,1-1e-12),method="bounded")
    assert pi==pytest.approx(optimum.x,abs=1e-5)
    assert loss==pytest.approx(optimum.fun,abs=1e-10)


def test_runner_provenance_and_gene_aggregation(tmp_path,monkeypatch):
    from test_normalization_genes import fixture
    from ytab.normalization.experiment import evaluate
    from ytab.normalization.gene_validation import evaluate as genes
    project,_=fixture(tmp_path,monkeypatch)
    validation=tmp_path/"assay"
    evaluate(project/"config/project.yaml",tmp_path/"catalogue",validation,["nonorm","ttr"],plots=False)
    site=tmp_path/"zinb"
    manifest=run_choices(validation/"matrix",site,["zinfnb"],"library_id")
    assert manifest["zinfnb_implementation_sha256"]
    result=genes(validation,tmp_path/"genes",site_normalization=site,plots=False)
    assert result["raw_aggregation_exact"]
    assert result["scaling_before_after_aggregation_agrees"]
    assert "zinfnb" in result["methods"]
