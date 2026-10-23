import csv
import gzip

import numpy as np
from scipy import sparse

from ytab.normalization.methods import normalize,run_choices
from ytab.normalization import gene_validation as gv


def test_hand_calculation_equalizes_positive_mean_not_density_or_total():
    raw=np.array([[2,10],[4,20],[0,30],[0,0]])
    result,params=normalize(raw,"nzmean",["g","g"])
    # Positive means 3 and 20; equally weighted group mean 11.5.
    np.testing.assert_allclose([r["factor"] for r in params],[11.5/3,11.5/20])
    actual=result.toarray()
    np.testing.assert_allclose([actual[:,j][raw[:,j]>0].mean() for j in range(2)],[11.5,11.5])
    assert not np.isclose(actual[:,0].sum(),actual[:,1].sum())
    assert np.array_equal(raw==0,actual==0)


def test_zero_padding_does_not_change_factors_or_counts_and_outliers_do():
    raw=np.array([[1,2],[3,6],[0,0]])
    a,pa=normalize(raw,"nzmean",["g","g"])
    b,pb=normalize(np.vstack([raw,np.zeros((50,2))]),"nzmean",["g","g"])
    np.testing.assert_allclose(a.toarray(),b.toarray()[:len(raw)])
    assert [r["factor"] for r in pa]==[r["factor"] for r in pb]
    outlier=raw.copy(); outlier[1,0]=3000
    _,pc=normalize(outlier,"nzmean",["g","g"])
    assert not np.isclose(pa[0]["factor"],pc[0]["factor"])


def test_pool_groups_and_single_sample_identity():
    raw=np.array([[1,2,100,200,7],[3,6,300,600,0]])
    actual,params=normalize(raw,"nzmean",["a","a","b","b","c"])
    np.testing.assert_allclose([r["factor"] for r in params],[1.5,.75,1.5,.75,1])
    np.testing.assert_array_equal(actual.toarray()[:,4],raw[:,4])


def test_nzmean_gene_validation_reuses_site_factors_not_gene_nonzero_means(tmp_path,monkeypatch):
    from test_normalization_genes import fixture
    _,validation=fixture(tmp_path,monkeypatch)
    site=tmp_path/"nzmean"
    run_choices(validation/"matrix",site,["nzmean"],"library_id")
    output=tmp_path/"genes"
    result=gv.evaluate(validation,output,site_normalization=site,plots=False)
    assert result["methods"]==["nonorm","ttr","cpm","nzmean"]
    with gzip.open(output/"nzmean_gene_counts.tsv.gz","rt") as handle:
        rows=list(csv.DictReader(handle,delimiter="\t"))
    expected=sparse.load_npz(site/"nzmean/normalized_counts.npz").toarray()
    np.testing.assert_allclose([[float(r[k]) for k in ["c1","c2","t1","t2"]] for r in rows],expected)


def test_noncoding_site_support_affects_factor_but_not_raw_gene_sums(tmp_path,monkeypatch):
    from test_normalization_genes import fixture
    from ytab.normalization.experiment import evaluate as evaluate_sites
    project,_=fixture(tmp_path,monkeypatch)
    hit=project/"create_hit_file/c1/c1_hits.txt"
    with hit.open("a") as handle:handle.write("chr1\t4\t300\n")
    validation=tmp_path/"assay_with_noncoding"
    evaluate_sites(project/"config/project.yaml",tmp_path/"catalogue",validation,["nonorm","ttr"],plots=False)
    choices=tmp_path/"nzmean"
    run_choices(validation/"matrix",choices,["nzmean"],"library_id")
    output=tmp_path/"genes"
    gv.evaluate(validation,output,site_normalization=choices,plots=False)
    with gzip.open(output/"nzmean_gene_counts.tsv.gz","rt") as handle:
        records=list(csv.DictReader(handle,delimiter="\t"))
    raw_sites=sparse.load_npz(validation/"matrix/raw_counts.npz")
    _,site_params=normalize(raw_sites,"nzmean",["pool"]*4)
    raw_genes=np.array([[10,20,5,10],[20,40,5,10],[0,0,0,0]])
    _,gene_params=normalize(raw_genes,"nzmean",["pool"]*4)
    assert not np.isclose(site_params[0]["factor"],gene_params[0]["factor"])
    assert float(records[0]["c1"])==10*site_params[0]["factor"]
