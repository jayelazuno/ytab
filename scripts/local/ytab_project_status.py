#!/usr/bin/env python3
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT/"src"))
from ytab.pipeline.project_status import build_project_status,write_project_status,determine_next_ready_stage
def main():
 p=argparse.ArgumentParser(description="Inspect YTAB project stage status without running analysis."); p.add_argument("--project-config",type=Path,required=True); p.add_argument("--json",action="store_true"); p.add_argument("--refresh",action="store_true"); p.add_argument("--profile",choices=("core","classifier","fitness","all"),default="all"); p.add_argument("--show-next",action="store_true"); a=p.parse_args()
 try: s=build_project_status(a.project_config); paths=write_project_status(a.project_config,s)
 except Exception as e: print(f"ERROR: {e}",file=sys.stderr); return 1
 if a.json: print(json.dumps(s,indent=2))
 else:
  print(f"Project ID: {s['project_id']}\nSpecies: {s['species']}\nProject config: {s['project_config']}"); print(f"{'stage':24} {'branch':10} {'status':12} {'done':>9}  blocked reason")
  for x in s["stages"]: print(f"{x['stage']:24} {x['branch']:10} {x['status']:12} {x['completed']:>4}/{x['expected']:<4}  {x['blocked_reason'] or x['message']}")
  if a.show_next: print("Next ready stage:",determine_next_ready_stage(s,a.profile) or "none")
  for w in s["warnings"]: print("WARNING:",w)
 print(f"JSON: {paths['json']}\nCSV: {paths['csv']}"); return 0
if __name__=="__main__": raise SystemExit(main())
