#!/usr/bin/env python3
"""Hermes site-restricted resampling, adapted from TRANSIT's resampling test (DeJesus et al. 2015) to empirical Hermes sites.

Input is a YTAB `normalize` output: ttr (default), totreads, nzmean, or the beta-geometric corrections
betageom / aBGC (BGC; recommended with --recenter trimmed-sum when replicate noise differs between
conditions). No normalization is run here. A warning is printed when replicate noise is unequal.
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
from ytab.essentiality.control_conditioned_gumbel import NONLINEAR_NORMALIZATIONS, NORMALIZATIONS
from ytab.essentiality.resampling import Settings,run

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--normalized",required=True,help="YTAB normalize output directory, or one method folder inside it (e.g. .../ttr)")
    parser.add_argument("--normalization",choices=NORMALIZATIONS+NONLINEAR_NORMALIZATIONS,help="Method folder to use when --normalized is a run directory (default ttr)")
    for name in ("gff","reference","output"):
        parser.add_argument("--"+name,required=True)
    parser.add_argument("--control",action="append",required=True,help="Control observation ID; repeat for each replicate")
    parser.add_argument("--treatment",action="append",required=True,help="Treatment observation ID; repeat for each replicate")
    defaults=Settings()
    parser.add_argument("--inference",choices=["exact","monte-carlo"],default=defaults.inference)
    parser.add_argument("--opportunities",choices=["any-support","control-support"],default=defaults.opportunities)
    for name in ("samples","grid_bins","min_sites","seed"):
        parser.add_argument("--"+name.replace("_","-"),type=int,default=getattr(defaults,name))
    parser.add_argument("--alternative",choices=["two-sided","depleted","enriched"],default=defaults.alternative)
    parser.add_argument("--recenter",choices=["none","median-gene","trimmed-sum"],default=defaults.recenter,
        help="Rescale the treatment once so typical genes have no change (guards against normalization shifts under widespread effects)")
    parser.add_argument("--winsorize",action="store_true",help="Label-invariant winsorization: the largest site x sample value of each gene is replaced by the second largest")
    for name in ("pseudocount","fdr","nterm","cterm","recenter_trim","asymmetry_warning"):
        parser.add_argument("--"+name.replace("_","-"),type=float,default=getattr(defaults,name))
    parser.add_argument("--stratify-libraries",action="store_true",help="Several independent libraries (e.g. paired pools): permute labels only within each library_id")
    parser.add_argument("--include-stop",action="store_true")
    args=parser.parse_args()
    values={k:getattr(args,k) for k in Settings.__dataclass_fields__ if k!="ignore_stop"}
    values["ignore_stop"]=not args.include_stop
    try:
        result=run(args.normalized,args.gff,args.reference,args.output,args.control,args.treatment,Settings(**values),args.normalization,args.stratify_libraries)
    except (ValueError,KeyError,OSError) as error:parser.exit(1,f"ERROR: {error}\n")
    print(json.dumps(dict(status=result["status"],analysis=result["analysis"],normalization=result["normalization"]["method"],**result["diagnostics"]),indent=2))
    warning=(result["diagnostics"].get("replicate_noise") or {}).get("warning")
    if warning:print(f"WARNING: {warning}",file=sys.stderr)

if __name__=="__main__":main()
