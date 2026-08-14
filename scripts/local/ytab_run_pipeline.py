#!/usr/bin/env python3
import argparse,json,shlex,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT/"src"))
from ytab.pipeline.orchestrator import build_pipeline_plan,run_pipeline_plan
from ytab.pipeline.project_status import load_project_status_context
def main():
 p=argparse.ArgumentParser(description="Run existing YTAB stages sequentially with cache-aware runners."); p.add_argument("--project-config",type=Path,required=True); p.add_argument("--profile",choices=("core","classifier","fitness","all"),default="all"); p.add_argument("--target",default="recommended"); p.add_argument("--analysis-id",default="treated_vs_parent_raw_cpm"); p.add_argument("--stages"); p.add_argument("--from-stage"); p.add_argument("--to-stage"); p.add_argument("--threads",type=int); p.add_argument("--force",action="store_true"); p.add_argument("--dry-run",action="store_true"); p.add_argument("--keep-going",action="store_true"); p.add_argument("--print-plan",action="store_true"); a=p.parse_args(); stages=a.stages.split(",") if a.stages else None
 try:
  c=load_project_status_context(a.project_config); plan=build_pipeline_plan(a.project_config,a.profile,a.target,a.analysis_id,stages,a.from_stage,a.to_stage)
  print(f"Project ID: {c['config']['project_id']}\nSpecies: {c['config']['species']}\nProfile: {a.profile}");
  if a.profile in {"classifier","all"}: print("Classifier target:",a.target)
  if a.profile in {"fitness","all"}: print("Fitness analysis ID:",a.analysis_id,"(raw-summary / CPM in R)")
  print("Plan:"); [print(f"  {i+1}. {x['stage']} [{x['current_status']}] ({x['branch']})") for i,x in enumerate(plan)]
  result=run_pipeline_plan(a.project_config,a.profile,a.target,a.analysis_id,a.threads,stages,a.from_stage,a.to_stage,a.force,a.dry_run,a.keep_going)
 except Exception as e: print(f"ERROR: {e}",file=sys.stderr); return 1
 print("Status:",result["status"]); print("Log:",result["log_file"]); [print(f"{x['stage']}: {x['status']} | {shlex.join(x['command'])}") for x in result["results"]]
 return 1 if result["status"]=="failed" else 0
if __name__=="__main__": raise SystemExit(main())
