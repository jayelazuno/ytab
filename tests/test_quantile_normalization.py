import numpy as np
from scipy import sparse

from ytab.normalization.methods import normalize


def test_quantile_equal_lengths_ties_and_zeros():
    raw = np.array([[0,0],[1,2],[1,4],[5,8]])
    result, rows = normalize(raw, "quantile", ["pool","pool"])
    # Mean sorted positives: [1.5,2.5,6.5]. First sample's tied 1s -> mean 2.
    np.testing.assert_allclose(result.toarray(), [[0,0],[2,1.5],[2,2.5],[6.5,6.5]])
    assert all(r["factor"] == "" for r in rows)  # Nonlinear, not a scalar factor.
    assert np.array_equal(result.toarray()==0, raw==0)


def test_quantile_unequal_density_and_singletons():
    raw = np.array([[0,4],[2,8],[6,12],[0,0]])
    result, _ = normalize(raw, "quantile", ["pool","pool"])
    np.testing.assert_allclose(result.toarray(), [[0,3],[3,6],[9,9],[0,0]])
    one, _ = normalize([[2,0],[0,6]], "quantile", ["pool","pool"])
    np.testing.assert_allclose(one.toarray(), [[4,0],[0,4]])


def test_quantile_strata_reproducibility_and_single_sample_identity():
    raw = sparse.csr_matrix([[0,0,100],[1,2,200],[1,4,0],[5,8,300]])
    before = raw.copy()
    a,_ = normalize(raw,"quantile",["a","a","b"])
    b,_ = normalize(raw,"quantile",["a","a","b"])
    np.testing.assert_array_equal(a.toarray(),b.toarray())
    np.testing.assert_array_equal(a.toarray()[:,2], raw.toarray()[:,2])
    assert (raw != before).nnz == 0
    permutation = [2,0,1]
    c,_ = normalize(raw[:,permutation], "quantile", ["b","a","a"])
    np.testing.assert_allclose(c.toarray()[:,np.argsort(permutation)], a.toarray())


def test_quantile_monotonicity_equal_input_ties_and_scale():
    raw = np.array([[0,1],[1,5],[2,5],[2,8],[10,50]])
    a,_ = normalize(raw,"quantile",["a","a"])
    b,_ = normalize(raw*3,"quantile",["a","a"])
    np.testing.assert_allclose(b.toarray(),a.toarray()*3)
    values = a.toarray()
    for i in range(2):
        assert np.all(np.diff(values[:,i])>=0)
    assert values[2,0] == values[3,0] and values[1,1] == values[2,1]


def test_nonlinear_transform_must_precede_gene_aggregation():
    raw = sparse.csr_matrix([[1,100],[3,2],[10,4]])
    mapping = sparse.csr_matrix([[1,1,0],[0,0,1]])
    sites,_ = normalize(raw,"quantile",["pool","pool"])
    genes_wrong,_ = normalize(mapping @ raw,"quantile",["pool","pool"])
    assert not np.allclose((mapping @ sites).toarray(), genes_wrong.toarray())
