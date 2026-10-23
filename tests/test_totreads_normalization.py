import csv
import gzip

import numpy as np
import pytest
from scipy import sparse

from ytab.normalization.methods import normalize,run_choices
from ytab.normalization import gene_validation as gv


def test_hand_calculation_equalizes_totals_not_positive_means():
    raw=np.array([[2,10],[4,20],[0,30],[0,0]])
    result,params=normalize(raw,"totreads",["g","g"])
    np.testing.assert_allclose([r["factor"] for r in params],[33/6,33/60])
    np.testing.assert_allclose(np.asarray(result.sum(axis=0)).ravel(),[33,33])
    positive_means=[result.toarray()[:,j][raw[:,j]>0].mean() for j in range(2)]
    assert positive_means[0]!=positive_means[1]
    assert np.array_equal(raw==0,result.toarray()==0)


def test_zero_padding_invariance_outliers_and_raw_preservation():
    raw=np.array([[1,2],[3,6],[0,0]])
    before=raw.copy()
    a,pa=normalize(raw,"totreads",["g","g"])
    b,pb=normalize(np.vstack([raw,np.zeros((50,2))]),"totreads",["g","g"])
    np.testing.assert_array_equal(raw,before)
    np.testing.assert_allclose(a.toarray(),b.toarray()[:len(raw)])
    assert [r["factor"] for r in pa]==[r["factor"] for r in pb]
    outlier=raw.copy();outlier[1,0]=3000
    _,pc=normalize(outlier,"totreads",["g","g"])
    assert not np.isclose(pa[0]["factor"],pc[0]["factor"])


def test_group_targets_single_sample_and_empty_rejection():
    raw=np.array([[1,2,100,200,7],[3,6,300,600,0]])
    result,params=normalize(raw,"totreads",["a","a","b","b","c"])
    np.testing.assert_allclose(np.asarray(result.sum(axis=0)).ravel(),[6,6,600,600,7])
    assert params[-1]["factor"]==1
    with pytest.raises(ValueError,match="empty samples"):
        normalize([[0,1]],"totreads",["g","g"])


def test_site_total_not_gene_cpm_denominator(tmp_path,monkeypatch):
    from test_normalization_genes import fixture
    from ytab.normalization.experiment import evaluate as evaluate_sites
    project,_=fixture(tmp_path,monkeypatch)
    hit=project/"create_hit_file/c1/c1_hits.txt"
    with hit.open("a") as handle:handle.write("chr1\t4\t300\n")
    validation=tmp_path/"assay"
    evaluate_sites(project/"config/project.yaml",tmp_path/"catalogue",validation,["nonorm","ttr"],plots=False)
    site=tmp_path/"totreads"
    run_choices(validation/"matrix",site,["totreads"],"library_id")
    output=tmp_path/"genes"
    result=gv.evaluate(validation,output,site_normalization=site,plots=False)
    assert result["methods"]==["nonorm","ttr","cpm","totreads"]
    with gzip.open(output/"totreads_gene_counts.tsv.gz","rt") as handle:
        records=list(csv.DictReader(handle,delimiter="\t"))
    raw_sites=sparse.load_npz(validation/"matrix/raw_counts.npz")
    _,params=normalize(raw_sites,"totreads",["g"]*4)
    raw_genes=np.array([[10,20,5,10],[20,40,5,10],[0,0,0,0]])
    _,gene_params=normalize(raw_genes,"totreads",["g"]*4)
    assert not np.isclose(params[0]["factor"],gene_params[0]["factor"])
    assert float(records[0]["c1"])==pytest.approx(10*params[0]["factor"])
    # Noncoding support influences scaling but is not added to raw gene sums.
    with gzip.open(output/"nonorm_gene_counts.tsv.gz","rt") as handle:
        raw_records=list(csv.DictReader(handle,delimiter="\t"))
    assert float(raw_records[0]["c1"])==10
