#!/usr/bin/env python3
"""Build a parent site universe and evaluate site-level TTR, without changing other analyses."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
from ytab.normalization.ttr import run

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest",required=True,help="Explicit source projects, parent observations and independent library IDs")
    parser.add_argument("--output",required=True,help="New, versioned output directory")
    parser.add_argument("--trim",type=float,default=.05,help="Fraction trimmed from each tail of positive counts (default .05)")
    parser.add_argument("--target",type=float,default=100)
    parser.add_argument("--no-plots",action="store_true")
    args=parser.parse_args()
    try:
        result=run(args.manifest,args.output,args.trim,args.target,not args.no_plots,
                   progress=lambda message:print(message,file=sys.stderr,flush=True))
    except (ValueError,KeyError,OSError) as error:
        parser.exit(1,f"ERROR: {error}\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
