import numpy as np
import pandas as pd
from scipy import sparse

from ytab.benchmarks import treated_vs_control as tvc


def _calibration(rng):
    lam = rng.gamma(.8, 30, 3000)
    raw = np.column_stack([rng.poisson(lam*e) for e in (1., .9, 1.1, .8)])
    sites = [("c", i*10) for i in range(len(lam))]
    return tvc.calibrate(sparse.csc_matrix(raw), sites, [0, 1], [2, 3])


def test_calibration_and_simulation_truth():
    rng = np.random.default_rng(0); cal = _calibration(rng)
    assert len(cal.exposures) == 4 and cal.abundance.min() > 0
    genes = [list(range(i*10, i*10+10)) for i in range(len(cal.abundance)//10)]
    raw, truth = tvc.simulate(cal, genes, "mixed", np.random.default_rng(1))
    assert raw.shape == (len(cal.abundance), 4) and raw.dtype == np.int64
    assert (truth != 0).sum() == round(.05*len(genes)) and set(np.unique(truth)) <= {-1, 0, 1}
    null_raw, null_truth = tvc.simulate(cal, genes, "null", np.random.default_rng(1))
    assert not null_truth.any()


def test_score_counts_wrong_direction_as_false():
    truth = pd.Series([-1, -1, 1, 0, 0, 0], index=list("abcdef"))
    result = pd.DataFrame(dict(score=[-5., -3., -2., 0., .1, 4.], call=[-1, 0, -1, 0, 0, 1], p=[.001, .2, .01, .9, .5, .02]), index=truth.index)
    s = tvc.score(result, truth)
    assert s["calls"] == 3 and s["correct_calls"] == 1 and s["false_calls"] == 2
    assert s["unaffected_tested"] == 3 and s["unaffected_p05"] == 1


def test_paired_tests():
    diff = np.full(30, .02)
    assert tvc.sign_flip(diff, flips=2000) < .01
    lo, hi = tvc.bootstrap_ci(diff); assert np.isclose(lo, .02) and np.isclose(hi, .02)
    adj = tvc.holm({"a": .01, "b": .04, "c": .03})
    assert adj == {"a": .03, "c": .06, "b": .06}


def test_false_discovery_proportion_and_asymmetric_scenario():
    truth = pd.Series([-1, 1, 0, 0], index=list("abcd"))
    none = pd.DataFrame(dict(score=0., call=0, p=np.nan), index=truth.index)
    assert tvc.score(none, truth)["fdp"] == 0
    some = none.assign(call=[-1, -1, 1, 0])   # one correct, one wrong direction, one on an unaffected gene
    assert np.isclose(tvc.score(some, truth)["fdp"], 2/3)
    rng = np.random.default_rng(0); cal = _calibration(rng)
    genes = [list(range(i*10, i*10+10)) for i in range(len(cal.abundance)//10)]
    _, t = tvc.simulate(cal, genes, "widespread_asymmetric", np.random.default_rng(2))
    assert abs((t == 1).sum()/(t != 0).sum() - .8) < .05 and (t != 0).sum() == round(.2*len(genes))
