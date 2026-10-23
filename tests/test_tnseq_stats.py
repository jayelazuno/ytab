import csv
import json
import numpy as np
import pytest
from scipy import stats
from ytab.qc.tnseq_stats import statistics, run, TNSEQ_STATS_COLUMNS


def test_transit_definitions_and_median_precision():
    x=np.array([1,2,4,9])
    result=statistics(x,10)
    assert result["density"]==.4
    assert result["mean_ct"]==1.6
    assert result["NZmean"]==4
    assert result["NZmedian"]==3
    assert result["skewness"]==pytest.approx(stats.skew(x))
    assert result["kurtosis"]==pytest.approx(stats.kurtosis(x))
    assert statistics([1,2],3)["NZmedian"]==1.5


def test_pti_matches_upstream_and_handles_ties_empty_short():
    x=np.arange(1,601)**2
    y=sorted(x,reverse=True)
    expected=np.median([np.log2((y[m]-y[2*m])/(y[2*m]-y[4*m])) for m in range(10,100)])
    assert statistics(x,700)["pickands_tail_index"]==pytest.approx(expected)
    for x,u in (([],1000),([1]*500,500),([1,2],2)):
        result=statistics(x,u)
        assert result["pickands_tail_index"] is None
        assert result["pti_valid_terms"]==0
    assert statistics([1]*500,500)["skewness"] is None
    with pytest.raises(ValueError):statistics([-1],10)


def test_excel_format_reuse_and_mismatch_guard(tmp_path):
    from test_normalization_choices import make_catalogue
    catalogue=tmp_path/"catalogue"; make_catalogue(catalogue)
    metadata=catalogue/"observations.tsv"
    with metadata.open() as handle:
        records=list(csv.DictReader(handle,delimiter="\t"))
    for record in records:record["sample"]=record["observation_id"]
    with metadata.open("w") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(records[0]),delimiter="\t")
        writer.writeheader();writer.writerows(records)
    qc=tmp_path/"qc"; qc.mkdir()
    (qc/"library_diagnostics.summary.csv").write_text(
        "sample,total_reads,unique_sites,jackpot_frac_reads,midlc_est\ncontrol,3,2,0.2,100\n")
    (qc/"control").mkdir()
    curve=qc/"control/control.midlc.csv"; curve.write_text("Reads Sampled,Unique Sites Trial1\n1,1\n")
    result=run(catalogue,tmp_path/"out",{"study":qc})
    assert result["count_state"]=="raw; no normalization"
    with (tmp_path/"out/tnseq_stats.tsv").open() as handle:
        reader=csv.DictReader(handle,delimiter="\t"); rows=list(reader)
        assert reader.fieldnames[:10]==TNSEQ_STATS_COLUMNS
    assert rows[0]["hermes_jackpot_frac_reads"]=="0.2"
    assert rows[0]["hermes_midlc_csv_file"]==str(curve)
    assert rows[1]["hermes_qc_status"]=="sample unavailable"
    assert str(curve) in result["source_hashes"]
    with pytest.raises(ValueError,match="empty output"):run(catalogue,tmp_path/"out")
    (qc/"library_diagnostics.summary.csv").write_text("sample,total_reads,unique_sites\ncontrol,99,2\n")
    run(catalogue,tmp_path/"bad",{"study":qc})
    with (tmp_path/"bad/tnseq_stats.tsv").open() as handle:
        rows=list(csv.DictReader(handle,delimiter="\t"))
    assert rows[0]["hermes_qc_status"]=="count mismatch; saved metrics withheld"
    assert "hermes_jackpot_frac_reads" not in rows[0]
