#!/usr/bin/env python3
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"));from ytab.pipeline.project_report import build_project_report
def main():
 p=argparse.ArgumentParser(description="Build a lightweight HTML summary from existing YTAB outputs.");p.add_argument("--project-config",type=Path,required=True);p.add_argument("--output",type=Path);p.add_argument("--force",action="store_true");p.add_argument("--open-browser",action="store_true");a=p.parse_args()
 try:r=build_project_report(a.project_config,a.output,a.force,a.open_browser);print(f"Status: {r['status']}\nReport: {r['report']}\nMetadata: {r['metadata']}");return 0
 except Exception as e:print(f"ERROR: {e}",file=sys.stderr);return 1
if __name__=="__main__":raise SystemExit(main())
