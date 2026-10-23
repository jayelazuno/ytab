import csv
import runpy
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import geom,trim_mean

from ytab.normalization.adaptive_bgc import DEFAULT_FRACTIONS,correct,discrepancy,fit_candidate,mapped_counts
from ytab.normalization.methods import normalize,run_choices


def counts():
    return np.array([[0,0],[1,2],[2,4],[4,8],[8,16],[16,32],[32,64],[64,128],[128,256]])


def test_positive_fit_matches_transit_expressions_and_grid():
    assert DEFAULT_FRACTIONS==tuple(i/100 for i in range(0,31,2))
    z=counts()[1:,0]
    actual=fit_candidate(z,.1)
    q=np.arange(len(z))/len(z)
    A=np.sum(np.log(1-q)**2)/np.sum(z*np.log(1-q))
    rho=1/trim_mean(z,.1)
    kappa=(2*np.exp(A)-1)/(np.exp(A)+rho-1)
    np.testing.assert_allclose([actual["A"],actual["rho"],actual["kappa"]],[A,rho,kappa],rtol=1e-12)
    with pytest.raises(ValueError):fit_candidate([1,1,1],0)


def test_score_tail_conventions_are_exact_and_explicit():
    y=np.array([1.,2,4,10])
    theoretical=geom.ppf(np.arange(1,5)/4,1/y.mean())
    theoretical[-1]=max(theoretical[:-1])+100
    expected=max((y-theoretical)**2/theoretical)
    assert discrepancy(y,"endpoint")==pytest.approx(expected)
    midpoint=geom.ppf((np.arange(4)+.5)/4,1/y.mean())
    assert discrepancy(y,"midpoint")==pytest.approx(max((y-midpoint)**2/midpoint))
    assert discrepancy(y,"endpoint")!=discrepancy(y,"midpoint")


def test_selection_candidate_reporting_zeros_and_monotonicity():
    raw=counts()
    before=raw.copy()
    result,details=correct(raw,["g","g"],post_ttr=False)
    np.testing.assert_array_equal(raw,before)
    assert np.array_equal(raw==0,result.toarray()==0)
    for j,detail in enumerate(details):
        candidates=detail["_abgc_candidates"]
        assert len(candidates)==16
        valid=[r for r in candidates if r["status"]=="valid"]
        winner=min(valid,key=lambda r:(r["score"],r["trim_fraction"]))
        assert winner["selected"] and detail["abgc_selected_trim"]==winner["trim_fraction"]
        assert sum(r["selected"] for r in candidates)==1
        assert np.all(np.diff(result.toarray()[1:,j])>=0)


def test_extra_zeros_do_not_change_positive_fit_or_mapping():
    raw=counts()
    padded=np.vstack([raw,np.zeros_like(raw)])
    a,da=correct(raw,["g","g"],post_ttr=False)
    b,db=correct(padded,["g","g"],post_ttr=False)
    np.testing.assert_array_equal(a.toarray(),b.toarray()[:len(raw)])
    assert [r["abgc_selected_trim"] for r in da]==[r["abgc_selected_trim"] for r in db]
    a,_=correct(raw,["g","g"])
    b,_=correct(padded,["g","g"])
    np.testing.assert_allclose(b.toarray()[:len(raw)],a.toarray()*2)


def test_simulation_reproducible_and_sample_order_independent():
    raw=counts()
    options=dict(cdf="simulation",seed=7,draws=3000,fractions=[0,.1,.2],post_ttr=False)
    a,da=correct(raw,["g","g"],**options)
    b,db=correct(raw[:,::-1],["g","g"],**options)
    np.testing.assert_array_equal(a.toarray(),b.toarray()[:,::-1])
    assert da[0]["abgc_selected_trim"]==db[1]["abgc_selected_trim"]
    assert da[0]["abgc_draws"]==3000


def test_invalid_candidates_and_inputs_fail_clearly():
    for raw in ([[0],[1]],[[1],[1],[1]]):
        with pytest.raises(ValueError,match="occupied coordinates|no valid"):
            correct(raw,["g"])
    for kwargs in (dict(tail="bad"),dict(seed=-1),dict(draws=0),dict(fractions=[.5]),dict(fractions=[0,0])):
        with pytest.raises(ValueError):correct(counts(),["g","g"],**kwargs)
    _,details=correct([[1],[1],[1],[2]],["g"])
    assert details[0]["abgc_invalid_candidates"]>0 and details[0]["abgc_valid_candidates"]>0
    assert all(r["reason"] for r in details[0]["_abgc_candidates"] if r["status"]=="invalid")


def test_shared_runner_exports_all_candidates_and_parameters(tmp_path):
    from test_normalization_choices import make_catalogue
    catalogue=tmp_path/"catalogue"
    make_catalogue(catalogue)
    output=tmp_path/"adaptive"
    result=run_choices(catalogue,output,["nonorm","aBGC"],"library_id")
    assert result["abgc_parameters"]["fit_universe"]=="positive counts only"
    with (output/"aBGC/candidates.tsv").open() as handle:
        rows=list(csv.DictReader(handle,delimiter="\t"))
    assert len(rows)==32
    assert sum(r["selected"]=="True" for r in rows)==2
    assert "_abgc_candidates" not in (output/"aBGC/factors.tsv").read_text()
    normalized,params=normalize(counts(),"aBGC",["g","g"])
    assert params[0]["factor"]=="" and params[0]["abgc_tail"]=="endpoint"


def test_analytic_one_read_boundary_and_mixture_quantile_bound():
    unique=np.array([1.,2.,3.,10.,100.,1000.])
    for rho,kappa in [(.02728847340291065,55.72826355514959),(.0012792628544724266,1113.6796070267583),(.2,2.)]:
        fit=dict(rho=rho,alpha=kappa*rho,beta=kappa*(1-rho))
        values,probability=mapped_counts(unique,fit,"analytic",0,1000,"fixture",.1,"candidate")
        assert values[0]==1 and probability[0]==rho
        assert np.all(values<=unique) and np.all(np.diff(values)>=0)


def test_audit_exports_tail_simulation_and_zero_padding_checks(tmp_path):
    from test_normalization_choices import make_catalogue
    catalogue=tmp_path/"catalogue"
    make_catalogue(catalogue)
    site=tmp_path/"adaptive"
    run_choices(catalogue,site,["aBGC"],"library_id")
    script=Path(__file__).resolve().parents[1]/"scripts/local/ytab_abgc_audit.py"
    audit=runpy.run_path(str(script))["audit"]
    report=audit(site,tmp_path/"audit",draws=1000)
    assert report["positive_fit_zero_padding_invariant"] and report["simulation_seed_reproduction"]
    assert (tmp_path/"audit/diagnostics.tsv").exists()


def test_earlier_tail_name_is_an_alias():
    y=np.array([1.,2.,3.,5.,8.])
    assert discrepancy(y,"transit")==discrepancy(y,"endpoint")
