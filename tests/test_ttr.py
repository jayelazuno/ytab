import csv
import gzip
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np
import pytest
from scipy import sparse
import yaml

from ytab.normalization.ttr import read_counts, ttr_parameters, expected_sites_at_depth, run


def test_ttr_zeros_scaling_and_trim():
    counts=np.array([1,2,4,8,16,32,64,128,256,1024]*2)
    p=ttr_parameters(counts,100)
    assert p["density"]==.2 and p["trimmed_sites_each_tail"]==1
    assert p["trimmed_nonzero_mean"]==np.mean(np.sort(counts)[1:-1])
    assert p["normalized_density_times_trimmed_mean"]==pytest.approx(100)
    larger=ttr_parameters(counts,200)
    assert larger["factor"]==pytest.approx(p["factor"]*2)
    a,b=ttr_parameters(counts,100),ttr_parameters(counts*3,100)
    assert a["factor"]/b["factor"]==pytest.approx(3)
    for invalid in ([],[0],[1,-1],[np.nan]):
        with pytest.raises(ValueError):ttr_parameters(invalid,100)
    with pytest.raises(ValueError):ttr_parameters([1],1,trim=.5)


def test_numerical_parity_with_local_transit():
    source=Path(__file__).resolve().parents[1]/"docs/codex/transit/src/pytransit/norm_tools.py"
    if not source.exists():pytest.skip("Optional upstream TRANSIT checkout unavailable")
    transit=runpy.run_path(str(source))["TTRNorm"]
    data=np.array([[0,1,2,3,10,100],[1,0,0,30,20,1000]],dtype=float)
    normalized,factors=transit.normalize(data)
    for i,row in enumerate(data):
        result=ttr_parameters(row[row>0],len(row))
        assert result["factor"]==pytest.approx(factors[i,0])
        np.testing.assert_allclose(row*result["factor"],normalized[i])


def test_exact_depth_rarefaction_against_enumeration():
    labels=["a","a","a","b","b","c"]
    for depth in range(7):
        expected=np.mean([len({labels[i] for i in subset}) for subset in itertools.combinations(range(6),depth)])
        assert expected_sites_at_depth([3,2,1],depth)==pytest.approx(expected)


def fixture_manifest(tmp_path):
    reference=tmp_path/"ref.fa"
    reference.write_text(">chr1\n"+"A"*100+"\n")
    sources=[]
    for name,counts in (("study-a",[(1,2),(2,3)]),("study-b",[(2,4),(3,5)])):
        project=tmp_path/name
        config=project/"config"
        config.mkdir(parents=True)
        sheet=config/"sample_sheet.csv"
        sheet.write_text("sample,include,library_role,pool\nparent,True,parent,1\ntreated,True,treated,1\n")
        path=project/"create_hit_file/parent/parent_hits.txt"
        path.parent.mkdir(parents=True)
        path.write_text("Chromosome\tSource\tHit position\tHit count\n"+"".join(f"chr1\tW\t{pos}\t{n}\n" for pos,n in counts))
        cfg=dict(repo_root=str(tmp_path),project_id=name,output_project_dir=name,sample_sheet=str(sheet),reference=dict(fasta=str(reference)))
        (config/"project.yaml").write_text(yaml.safe_dump(cfg))
        sources.append(dict(project_config=str(config/"project.yaml"),dataset=name,
          observations=[dict(sample="parent",library_id=name+"-pool")]))
    manifest=tmp_path/"manifest.yaml"
    manifest.write_text(yaml.safe_dump(dict(version=1,repo_root=str(tmp_path),sources=sources)))
    return manifest


def test_pipeline_matrix_provenance_overlap_and_input_preservation(tmp_path):
    manifest=fixture_manifest(tmp_path)
    output=tmp_path/"results"
    result=run(manifest,output,plots=False)
    assert result["universe_sites"]==3 and result["independent_libraries"]==2
    matrix=sparse.load_npz(output/"raw_counts.npz").toarray()
    np.testing.assert_array_equal(matrix,[[2,0],[3,4],[0,5]])
    normalized=sparse.load_npz(output/"ttr_normalized_counts.npz").toarray()
    assert np.array_equal(matrix==0,normalized==0)
    with gzip.open(output/"raw_counts.tsv.gz","rt") as handle:
        rows=list(csv.reader(handle,delimiter="\t"))
    assert rows[1][2:]==["2","0"] and rows[3][2:]==["0","5"]
    with (output/"dataset_overlap.tsv").open() as handle:
        overlap=list(csv.DictReader(handle,delimiter="\t"))[0]
    assert overlap["shared_sites"]=="1" and overlap["left_only"]==overlap["right_only"]=="1"
    assert json.loads((output/"run_manifest.json").read_text())["reference_sha256"]
    with pytest.raises(ValueError,match="nonempty"):run(manifest,output,plots=False)


def test_repeated_observations_are_not_independent_libraries(tmp_path):
    manifest=fixture_manifest(tmp_path)
    data=yaml.safe_load(manifest.read_text())
    data["sources"][1]["dataset"]=data["sources"][0]["dataset"]
    data["sources"][1]["observations"][0]["library_id"]=data["sources"][0]["observations"][0]["library_id"]
    manifest.write_text(yaml.safe_dump(data))
    result=run(manifest,tmp_path/"repeated",plots=False)
    assert result["observations"]==2 and result["independent_libraries"]==1


def test_strands_sum_never_cancel_and_treatment_rejected(tmp_path):
    path=tmp_path/"hits.txt"
    path.write_text("Chromosome\tSource\tHit position\tHit count\nchr1\tW\t5\t3\nchr1\tC\t5\t4\n")
    counts,strands=read_counts(path)
    assert counts[("chr1",5)]==7 and strands[("chr1",5)]=={"W","C"}
    path.write_text("Chromosome\tHit position\tHit count\nchr1\t5\t-1\n")
    with pytest.raises(ValueError):read_counts(path)
    manifest=fixture_manifest(tmp_path)
    data=yaml.safe_load(manifest.read_text())
    data["sources"][0]["observations"][0]["sample"]="treated"
    manifest.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError,match="untreated"):run(manifest,tmp_path/"treated-results",plots=False)


def test_reference_mismatch_fails_without_writing_results(tmp_path):
    manifest=fixture_manifest(tmp_path)
    data=yaml.safe_load(manifest.read_text())
    config=Path(data["sources"][1]["project_config"])
    cfg=yaml.safe_load(config.read_text())
    other=tmp_path/"different_reference.fa"
    other.write_text(">chr1\n"+"C"*100+"\n")
    cfg["reference"]["fasta"]=str(other)
    config.write_text(yaml.safe_dump(cfg))
    output=tmp_path/"incompatible"
    with pytest.raises(ValueError,match="reference FASTA checksum"):run(manifest,output,plots=False)
    assert not output.exists()


def test_duplicate_observation_rejected(tmp_path):
    manifest=fixture_manifest(tmp_path)
    data=yaml.safe_load(manifest.read_text())
    data["sources"].append(data["sources"][0])
    manifest.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError,match="unique"):run(manifest,tmp_path/"duplicates",plots=False)
