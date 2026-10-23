#!/usr/bin/env python3
"""Explicit normalization choices and method explanations; no automatic replacement of existing workflows."""
import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
from ytab.normalization.methods import METHODS,method_catalogue,run_choices

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list-methods",action="store_true")
    parser.add_argument("--describe",choices=list(METHODS))
    parser.add_argument("--catalogue",help="Existing empirical catalogue directory containing raw_counts.npz")
    parser.add_argument("--output",help="New, separate output directory")
    parser.add_argument("--method",action="append",choices=list(METHODS),help="Repeat to compare methods; no default method")
    parser.add_argument("--group-by",choices=["library_id","project","dataset","cohort"],help="Required explicit scaling groups; library_id preserves pool stratification")
    parser.add_argument("--observation",action="append",help="Optional explicit project::sample selection; repeat as needed")
    parser.add_argument("--trim",type=float,default=.05)
    parser.add_argument("--target",type=float,default=100)
    parser.add_argument("--bgc-cdf",choices=["analytic","simulation"],default="analytic")
    parser.add_argument("--bgc-seed",type=int,default=0)
    parser.add_argument("--bgc-draws",type=int,default=200000)
    parser.add_argument("--bgc-no-ttr",action="store_true",help="Disable BGC's explicitly documented final TTR scaling")
    parser.add_argument("--abgc-cdf",choices=["analytic","simulation"],default="analytic")
    parser.add_argument("--abgc-seed",type=int,default=0)
    parser.add_argument("--abgc-draws",type=int,default=200000)
    parser.add_argument("--abgc-tail",choices=["endpoint","midpoint","transit"],default="endpoint",help="Explicit geometric scoring quantiles: endpoint (r/n; 'transit' is an accepted earlier name) or finite midpoint sensitivity")
    parser.add_argument("--abgc-fractions",help="Optional comma-separated trim fractions; default 0,0.02,...,0.30")
    parser.add_argument("--abgc-no-ttr",action="store_true")
    args=parser.parse_args()
    if args.list_methods:
        print(json.dumps(method_catalogue(),indent=2));return
    if args.describe:
        print(json.dumps(asdict(METHODS[args.describe]),indent=2));return
    if not all((args.catalogue,args.output,args.method,args.group_by)):
        parser.error("Choose --catalogue, --output, --method and --group-by explicitly, or use --list-methods/--describe")
    try:
        adaptive=dict(cdf=args.abgc_cdf,seed=args.abgc_seed,draws=args.abgc_draws,tail=args.abgc_tail,post_ttr=not args.abgc_no_ttr)
        if args.abgc_fractions is not None:adaptive["fractions"]=[float(f) for f in args.abgc_fractions.split(",")]
        result=run_choices(args.catalogue,args.output,args.method,args.group_by,args.observation,args.trim,args.target,
            dict(cdf=args.bgc_cdf,seed=args.bgc_seed,draws=args.bgc_draws,post_ttr=not args.bgc_no_ttr),adaptive)
    except (ValueError,KeyError,OSError) as error:
        parser.exit(1,f"ERROR: {error}\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":main()
