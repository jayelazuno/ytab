#!/usr/bin/env python3
"""Identify individual Hermes mutants from local Sanger files or archives."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
from ytab.mutant_screen.identify import identify, identify_current
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference-config","--project-config",dest="project_config",required=True,type=Path,help="Configuration selecting the alignment reference; does not link screen results")
    p.add_argument("--input",nargs="+",required=True,type=Path)
    current=Path(__file__).resolve().parents[2]/"output/mutant_screen/current"
    p.add_argument("--output",type=Path,default=current,help="Default: shared current identification loaded by Shiny; specify another directory for standalone results")
    p.add_argument("--hermes",type=Path)
    for name,default in (("threads",2),("min-length",80),("min-mapq",20),("quality",20),("overlap",12),("window",200),("junction-min-overlap",20),("junction-gene-window",100),("junction-gene-min-bases",50)):
        p.add_argument("--"+name,type=int,default=default)
    p.add_argument("--min-identity",type=float,default=.9)
    p.add_argument("--error-rate",type=float,default=.15)
    p.add_argument("--fitness-result",type=Path)
    p.add_argument("--link-project-config",type=Path)
    p.add_argument("--replace-current",action="store_true")
    p.add_argument("--junction-read-end",choices=("start","end","require"),default="start")
    p.add_argument("--condition",default="");p.add_argument("--pool",default="")
    a=vars(p.parse_args());a["inputs"]=a.pop("input")
    replace=a.pop("replace_current")
    runner=identify_current if replace or a["output"].resolve()==current.resolve() else identify
    try:
        result=runner(**a);print(json.dumps(dict(status="success",output_dir=result["output_dir"],reads=len(result["results"]),reused=result.get("reused",False))));return 0
    except Exception as e:
        print(json.dumps(dict(status="failure",error=str(e))),file=sys.stderr);return 1
if __name__=="__main__":raise SystemExit(main())
