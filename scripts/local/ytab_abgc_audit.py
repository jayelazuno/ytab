#!/usr/bin/env python3
"""Audit adaptive BGC endpoint and simulation sensitivity on an existing site run."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
from ytab.normalization.adaptive_bgc import correct
from ytab.normalization.ttr import sha256,write_table


def audit(site_run,output,seed=0,draws=200000):
    site_run,output=Path(site_run).resolve(),Path(output).resolve()
    if output.exists() and any(output.iterdir()):raise ValueError("Choose a new output directory")
    source=json.loads((site_run/"run_manifest.json").read_text())
    params=source["abgc_parameters"]
    if source.get("status")!="success" or params["cdf"]!="analytic":raise ValueError("Audit requires a completed analytic adaptive BGC run")
    catalogue=Path(source["catalogue"])
    for name,digest in source["source_hashes"].items():
        if sha256(catalogue/name)!=digest:raise ValueError("Source inputs changed")
    with (site_run/"observations.tsv").open() as handle:
        metadata=list(csv.DictReader(handle,delimiter="\t"))
    if source["group_by"]!="library_id" or source["observations"]!=[r["observation_id"] for r in metadata]:
        raise ValueError("Audit requires full ordered library_id observations")
    groups=[r["library_id"] for r in metadata]
    raw=sparse.load_npz(catalogue/"raw_counts.npz").tocsc()
    actual=sparse.load_npz(site_run/"aBGC/normalized_counts.npz").tocsc()
    options={k:params[k] for k in ("cdf","seed","draws","post_ttr","tail","fractions","trim","target")}
    exact,selected=correct(raw,groups,**options)
    exact=exact.tocsc()
    np.testing.assert_allclose(exact.data,actual.data,rtol=1e-12)
    np.testing.assert_array_equal(raw.indptr,actual.indptr)
    np.testing.assert_array_equal(raw.indices,actual.indices)
    alternative="midpoint" if options["tail"] in {"endpoint","transit"} else "endpoint"
    other,other_details=correct(raw,groups,**(options|dict(tail=alternative)))
    sim,sim_details=correct(raw,groups,**(options|dict(cdf="simulation",seed=seed,draws=draws)))
    repeat,_=correct(raw,groups,**(options|dict(cdf="simulation",seed=seed,draws=draws)))
    sim,repeat=sim.tocsc(),repeat.tocsc()
    np.testing.assert_array_equal(sim.data,repeat.data)
    extra=max(1,int(np.ceil(raw.shape[0]*.1)))
    padded=sparse.vstack([raw,sparse.csc_matrix((extra,raw.shape[1]))]).tocsc()
    pad,pad_details=correct(padded,groups,**options)
    pad=pad[:raw.shape[0]].tocsc()
    np.testing.assert_allclose(pad.data,actual.data*((raw.shape[0]+extra)/raw.shape[0] if options["post_ttr"] else 1),rtol=1e-12)
    results=[]
    for j,meta in enumerate(metadata):
        x=raw.data[raw.indptr[j]:raw.indptr[j+1]]
        y=actual.data[actual.indptr[j]:actual.indptr[j+1]]
        z=sim.data[sim.indptr[j]:sim.indptr[j+1]]
        if np.any(np.diff(y[np.argsort(x,kind="stable")]) < -1e-10):raise ValueError("Nonmonotone transform")
        assert selected[j]["abgc_selected_trim"]==pad_details[j]["abgc_selected_trim"]
        results.append(dict(observation_id=meta["observation_id"],selected_trim=selected[j]["abgc_selected_trim"],
            selected_score=selected[j]["abgc_score"],valid_candidates=selected[j]["abgc_valid_candidates"],
            invalid_candidates=selected[j]["abgc_invalid_candidates"],alternative_tail=alternative,
            alternative_tail_selected_trim=other_details[j]["abgc_selected_trim"],
            simulation_selected_trim=sim_details[j]["abgc_selected_trim"],
            simulation_median_relative_count_difference=float(np.median(np.abs(z-y)/y)),
            simulation_99pct_relative_count_difference=float(np.quantile(np.abs(z-y)/y,.99)),
            positive_fit_unchanged_with_extra_zero_rows=True))
    output.mkdir(parents=True,exist_ok=True)
    write_table(output/"diagnostics.tsv",list(results[0]),results)
    report=dict(status="success",site_run=str(site_run),seed=seed,draws=draws,exact_reproduction=True,
        simulation_seed_reproduction=True,zero_pattern_preserved=True,positive_fit_zero_padding_invariant=True,
        source_hashes_unchanged=True,source_manifest_sha256=sha256(site_run/"run_manifest.json"),
        implementation_sha256=sha256(__file__),
        caution="Selection scores and tail/simulation sensitivity are not p-values or evidence of biological correctness; endpoint rule is explicit.")
    (output/"run_manifest.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-run",required=True)
    parser.add_argument("--output",required=True)
    parser.add_argument("--seed",type=int,default=0)
    parser.add_argument("--draws",type=int,default=200000)
    args=parser.parse_args()
    try:result=audit(args.site_run,args.output,args.seed,args.draws)
    except (ValueError,KeyError,OSError,AssertionError) as error:parser.exit(1,f"ERROR: {error}\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
