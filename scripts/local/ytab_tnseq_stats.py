#!/usr/bin/env python3
"""Export per-library Tn-seq count statistics (adapted from TRANSIT's tnseq_stats) plus saved Hermes QC (Excel-readable TSV)."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
from ytab.qc.tnseq_stats import run

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue",required=True,help="Existing raw coordinate matrix and metadata directory")
    parser.add_argument("--output",required=True,help="New output directory")
    parser.add_argument("--diagnostics",action="append",default=[],metavar="PROJECT=DIR",
                        help="Repeat for existing library_diagnostics directories; exact metadata project key")
    args=parser.parse_args()
    try:
        pairs=[entry.split("=",1) for entry in args.diagnostics]
        if any(len(pair)!=2 or not all(pair) for pair in pairs) or len({p[0] for p in pairs})!=len(pairs):
            raise ValueError("Use unique --diagnostics PROJECT=DIR mappings")
        print(json.dumps(run(args.catalogue,args.output,dict(pairs)),indent=2))
    except (ValueError,KeyError,OSError) as error:
        parser.exit(1,f"ERROR: {error}\n")

if __name__=="__main__":main()
