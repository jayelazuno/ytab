import csv
import json
from pathlib import Path
import runpy
import subprocess
import sys

import numpy as np
import pytest
from scipy import sparse

from ytab.normalization.methods import METHODS, method_catalogue, normalize, run_choices


def test_registry_describes_status_assumptions_and_scope():
    cards=method_catalogue()
    assert len({c["id"] for c in cards})==len(cards)
    assert {c["id"] for c in cards if c["status"]=="available"}=={"nonorm","ttr","totreads","nzmean","quantile","betageom","aBGC","zinfnb","midlc"}
    for card in cards:
        for field in ("description","formula","may_suit","caution","input_level","diagnostics","source"):
            assert card[field]
    assert METHODS["cpm"].status=="existing_workflow" and METHODS["midlc"].status=="available"


def test_all_available_methods_against_transit():
    path=Path(__file__).resolve().parents[1]/"docs/codex/transit/src/pytransit/norm_tools.py"
    if not path.exists():pytest.skip("Optional local TRANSIT checkout unavailable")
    upstream=runpy.run_path(str(path))["methods"]
    matrix=np.array([[0,2,1],[2,4,0],[10,8,3],[100,60,11]],dtype=float)
    for method,source in (("ttr","TTR"),("totreads","totreads"),("nzmean","nzmean"),("nonorm","nonorm")):
        result,rows=normalize(matrix,method,["group"]*3)
        expected,factors=upstream[source].normalize(matrix.T.copy())
        np.testing.assert_allclose(result.toarray(),expected.T)
        if method!="nonorm":np.testing.assert_allclose([r["factor"] for r in rows],factors.ravel())
        assert np.array_equal(result.toarray()==0,matrix==0)


def test_grouping_preserves_pool_strata_and_raw_inputs():
    matrix=sparse.csr_matrix([[10,20,100,300],[0,0,0,0]])
    before=matrix.copy()
    result,rows=normalize(matrix,"totreads",["pool-a","pool-a","pool-b","pool-b"])
    np.testing.assert_allclose(result.toarray()[0],[15,15,200,200])
    assert (before!=matrix).nnz==0
    assert all(r["group_samples"]==2 for r in rows)
    with pytest.raises(ValueError,match="not selectable"):normalize(matrix,"emphist",["g"]*4)
    with pytest.raises(ValueError,match="not selectable"):normalize(matrix,"cpm",["g"]*4)
    with pytest.raises(ValueError,match="empty samples"):normalize([[1,0]],"ttr",["g","g"])
    with pytest.raises(ValueError,match="raw integer"):normalize([[1.5,2]],"totreads",["g","g"])


def make_catalogue(root):
    root.mkdir()
    raw=sparse.csr_matrix([[1,2],[2,4],[0,3]])
    sparse.save_npz(root/"raw_counts.npz",raw)
    with (root/"observations.tsv").open("w") as handle:
        writer=csv.DictWriter(handle,fieldnames=["observation_id","display_name","library_id","project","dataset"],delimiter="\t")
        writer.writeheader()
        writer.writerows([dict(observation_id=x,display_name=x,library_id="pool",project="study",dataset="data") for x in ["control","treatment"]])
    (root/"sites.tsv.gz").write_bytes(b"fixture coordinate identity")
    (root/"run_manifest.json").write_text(json.dumps(dict(status="success",universe_sites=3,reference_sha256="fixture-reference")))
    return raw


def test_method_outputs_are_separate_and_each_starts_from_raw(tmp_path):
    catalogue=tmp_path/"catalogue"
    raw=make_catalogue(catalogue)
    output=tmp_path/"choices"
    result=run_choices(catalogue,output,["ttr","totreads","nzmean","nonorm"],"library_id")
    for method in result["methods"]:
        actual=sparse.load_npz(output/method/"normalized_counts.npz")
        expected,_=normalize(raw,method,["pool","pool"])
        np.testing.assert_allclose(actual.toarray(),expected.toarray())
        assert (output/method/"method.json").exists()
    assert (sparse.load_npz(catalogue/"raw_counts.npz")!=raw).nnz==0
    assert result["software_versions"] and len(result["source_hashes"] )==4
    from hashlib import sha256
    assert sha256((output/"implementation_snapshot.py").read_bytes()).hexdigest()==result["implementation_sha256"]
    with pytest.raises(ValueError,match="nonempty"):run_choices(catalogue,output,["ttr"],"library_id")
    with pytest.raises(ValueError,match="not selectable"):run_choices(catalogue,tmp_path/"planned",["ttr","emphist"],"library_id")
    assert not (tmp_path/"planned").exists()


def test_selection_order_and_missing_method_fail_clearly(tmp_path):
    catalogue=tmp_path/"catalogue"
    raw=make_catalogue(catalogue)
    output=tmp_path/"selected"
    result=run_choices(catalogue,output,["nonorm"],"library_id",["treatment"])
    assert result["observations"]==["treatment"]
    np.testing.assert_array_equal(sparse.load_npz(output/"nonorm/normalized_counts.npz").toarray(),raw.toarray()[:,1:])
    script=Path(__file__).resolve().parents[1]/"scripts/local/ytab_normalize.py"
    call=subprocess.run([sys.executable,str(script),"--catalogue",str(catalogue),"--output",str(tmp_path/"missing")],capture_output=True,text=True)
    assert call.returncode!=0 and "explicitly" in call.stderr
    describe=subprocess.run([sys.executable,str(script),"--describe","ttr"],capture_output=True,text=True,check=True)
    assert json.loads(describe.stdout)["status"]=="available"
