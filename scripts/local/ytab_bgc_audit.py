#!/usr/bin/env python3
"""Audit completed analytic BGC site results: simulation, fit and catalogue sensitivity."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
from ytab.normalization.betageom import beta_geometric_cdf,correct,fit_parameters
from ytab.normalization.ttr import sha256,write_table


def audit(site_run, output, seed=0, draws=200000):
    site_run,output=Path(site_run).resolve(),Path(output).resolve()
    if output.exists() and any(output.iterdir()):raise ValueError("Choose a new output directory")
    manifest=json.loads((site_run/"run_manifest.json").read_text())
    params=manifest["bgc_parameters"]
    if manifest.get("status")!="success" or params["cdf"]!="analytic":raise ValueError("Audit requires a completed analytic BGC result")
    catalogue=Path(manifest["catalogue"])
    for name,digest in manifest["source_hashes"].items():
        if sha256(catalogue/name)!=digest:raise ValueError("Source inputs changed")
    with (site_run/"observations.tsv").open() as handle:
        metadata=list(csv.DictReader(handle,delimiter="\t"))
    groups=[r["library_id"] for r in metadata]
    raw=sparse.load_npz(catalogue/"raw_counts.npz").tocsc()
    actual=sparse.load_npz(site_run/"betageom/normalized_counts.npz").tocsc()
    if manifest["group_by"]!="library_id" or manifest["observations"]!=[r["observation_id"] for r in metadata]:
        raise ValueError("Audit requires full ordered library_id observations")
    exact,details=correct(raw,groups,post_ttr=params["post_ttr"],trim=params["trim"],target=params["target"])
    exact=exact.tocsc()
    if not np.array_equal(exact.indptr,actual.indptr) or not np.array_equal(exact.indices,actual.indices):raise ValueError("Site patterns differ")
    np.testing.assert_allclose(exact.data,actual.data,rtol=1e-12)
    simulation,_=correct(raw,groups,cdf="simulation",seed=seed,draws=draws,post_ttr=params["post_ttr"],trim=params["trim"],target=params["target"])
    repeated,_=correct(raw,groups,cdf="simulation",seed=seed,draws=draws,post_ttr=params["post_ttr"],trim=params["trim"],target=params["target"])
    simulation,repeated=simulation.tocsc(),repeated.tocsc()
    np.testing.assert_array_equal(simulation.data,repeated.data)
    np.testing.assert_array_equal(raw.indptr,actual.indptr)
    np.testing.assert_array_equal(raw.indices,actual.indices)
    results=[]
    for j,meta in enumerate(metadata):
        x=raw.data[raw.indptr[j]:raw.indptr[j+1]]
        y=actual.data[actual.indptr[j]:actual.indptr[j+1]]
        z=simulation.data[simulation.indptr[j]:simulation.indptr[j+1]]
        order=np.argsort(x,kind="stable")
        if np.any(np.diff(y[order]) < -1e-10):raise ValueError("Nonmonotone mapping")
        fit=fit_parameters(x,raw.shape[0])
        expanded=fit_parameters(x,raw.shape[0]+int(np.ceil(raw.shape[0]*.1)))
        unique,frequency=np.unique(x,return_counts=True)
        empirical=np.cumsum(frequency)/len(x)
        cdf=beta_geometric_cdf(unique,fit["alpha"],fit["beta"])
        results.append(dict(observation_id=meta["observation_id"],positive_sites=len(x),
            descriptive_max_abs_positive_cdf_difference=float(np.max(np.abs(empirical-cdf))),
            simulation_median_relative_count_difference=float(np.median(np.abs(z-y)/y)),
            simulation_99pct_relative_count_difference=float(np.quantile(np.abs(z-y)/y,.99)),
            rho=fit["rho"],rho_with_10pct_extra_zero_rows=expanded["rho"],
            kappa=fit["kappa"],kappa_with_10pct_extra_zero_rows=expanded["kappa"],
            raw_top_1pct_fraction=details[j]["bgc_raw_top_1pct_support_fraction"],
            corrected_top_1pct_fraction=details[j]["bgc_corrected_top_1pct_support_fraction"]))
    output.mkdir(parents=True,exist_ok=True)
    write_table(output/"diagnostics.tsv",list(results[0]),results)
    report=dict(status="success",site_run=str(site_run),seed=seed,draws=draws,
        exact_reproduction=True,simulation_seed_reproduction=True,zero_pattern_preserved=True,
        source_hashes_unchanged=True,positive_transform_monotone=True,
        source_manifest_sha256=sha256(site_run/"run_manifest.json"),implementation_sha256=sha256(__file__),
        caution="CDF difference is descriptive, not a goodness-of-fit p-value. Added zero rows are a fit-sensitivity diagnostic, not a different biological catalogue or a recommendation.")
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
