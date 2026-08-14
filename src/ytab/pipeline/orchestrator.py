"""Sequential orchestration of existing YTAB stage CLIs."""
from __future__ import annotations
import json, os, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path
from .project_status import DEPENDENCIES, build_project_status, load_project_status_context, write_project_status

ORDER=["reference","mapfastq","create_hit_file","summary","library_diagnostics","sample_normalization",
       "summary_normalized","combined_hits","summary_combined","classifier","comparison_design","treated_vs_parent"]
PROFILES={"core":ORDER[:5],"classifier":ORDER[:10],"fitness":ORDER[:5]+ORDER[10:],"all":ORDER}
CLI={"reference":"ytab_prepare_reference.py","mapfastq":"ytab_run_mapfastq.py","create_hit_file":"ytab_run_create_hit_file.py",
 "summary":"ytab_run_summary_table.py","library_diagnostics":"ytab_run_library_diagnostics.py",
 "sample_normalization":"ytab_run_sample_normalization.py","summary_normalized":"ytab_run_summary_normalized.py",
 "combined_hits":"ytab_run_combine_hits.py","summary_combined":"ytab_run_summary_combined.py","classifier":"ytab_run_classifier.py",
 "comparison_design":"ytab_init_comparison_design.py","treated_vs_parent":"ytab_run_treated_vs_parent.py"}
def _now(): return datetime.now(timezone.utc).isoformat()
def _stamp(): return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

def _command(stage,context,target,analysis_id,threads,force,dry_run,keep_going):
    root=context["root"]; cfg=str(context["project_config"]); script=root/"scripts"/"local"/CLI[stage]
    if stage=="reference": cmd=[sys.executable,str(script),"--species",str(context["config"]["species"]),"--threads",str(threads),"--repo-root",str(root)]
    else: cmd=[sys.executable,str(script),"--project-config",cfg]
    if stage in {"mapfastq","create_hit_file","summary","library_diagnostics","sample_normalization","summary_normalized","summary_combined"}: cmd += ["--threads",str(threads)]
    if stage=="sample_normalization": cmd += ["--targets","auto","--sample-mode","parents"]
    if stage=="summary_normalized": cmd += ["--targets",target,"--sample-mode","parents"]
    if stage in {"combined_hits","summary_combined","classifier"}: cmd += ["--target",target]
    if stage=="treated_vs_parent": cmd += ["--analysis-id",analysis_id,"--input-mode","raw-summary"]
    if force: cmd.append("--overwrite" if stage=="comparison_design" else "--force")
    if dry_run and stage not in {"reference","comparison_design"}: cmd.append("--dry-run")
    if keep_going and stage in {"mapfastq","create_hit_file","summary","sample_normalization","summary_normalized","treated_vs_parent"}: cmd.append("--keep-going")
    return cmd

def build_pipeline_plan(project_config:Path,profile:str,target="recommended",analysis_id="treated_vs_parent_raw_cpm",stages=None,from_stage=None,to_stage=None):
    if profile not in PROFILES: raise ValueError(f"Unsupported profile: {profile}")
    context=load_project_status_context(project_config); selected=list(PROFILES[profile])
    if stages:
        unknown=set(stages)-set(ORDER)
        if unknown: raise ValueError("Unknown stages: "+", ".join(sorted(unknown)))
        selected=[s for s in ORDER if s in stages]
    if from_stage:
        if from_stage not in selected: raise ValueError(f"from-stage is not in plan: {from_stage}")
        selected=selected[selected.index(from_stage):]
    if to_stage:
        if to_stage not in selected: raise ValueError(f"to-stage is not in plan: {to_stage}")
        selected=selected[:selected.index(to_stage)+1]
    status=build_project_status(project_config); lookup={x["stage"]:x for x in status["stages"]}
    return [{"stage":s,"branch":lookup[s]["branch"],"current_status":lookup[s]["status"],"dependencies":DEPENDENCIES.get(s,[]),
             "target":target if s in {"summary_normalized","combined_hits","summary_combined","classifier"} else None,
             "analysis_id":analysis_id if s=="treated_vs_parent" else None} for s in selected]

def validate_pipeline_plan(plan,context):
    errors=[]
    for item in plan:
        if item["stage"]=="treated_vs_parent" and item.get("target"): errors.append("Fitness analysis must not receive a classifier normalization target.")
    return {"valid":not errors,"errors":errors,"warnings":[]}

def execute_pipeline_stage(stage,context,threads,force,dry_run,keep_going,target="recommended",analysis_id="treated_vs_parent_raw_cpm",log=None):
    command=_command(stage["stage"],context,target,analysis_id,threads,force,dry_run,keep_going)
    if dry_run: return {"stage":stage["stage"],"status":"dry-run","command":command,"returncode":0}
    if stage.get("current_status")=="success" and not force:
        return {"stage":stage["stage"],"status":"skipped","command":command,"returncode":0,
                "message":"Successful project status found; existing stage outputs reused."}
    completed=subprocess.run(command,cwd=context["root"],stdout=log or subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    return {"stage":stage["stage"],"status":"success" if completed.returncode==0 else "failed","command":command,
            "returncode":completed.returncode,"output":completed.stdout if log is None else None}

def run_pipeline_plan(project_config:Path,profile="all",target="recommended",analysis_id="treated_vs_parent_raw_cpm",threads=None,
                      stages=None,from_stage=None,to_stage=None,force=False,dry_run=False,keep_going=False):
    context=load_project_status_context(project_config); threads=int(threads or min(2,int(context["config"].get("threads",2))))
    plan=build_pipeline_plan(project_config,profile,target,analysis_id,stages,from_stage,to_stage); validation=validate_pipeline_plan(plan,context)
    if not validation["valid"]: raise ValueError("\n".join(validation["errors"]))
    orch=context["manifests"]/"orchestrator"; orch.mkdir(parents=True,exist_ok=True); logs=context["project"]/"logs"/"orchestrator"; logs.mkdir(parents=True,exist_ok=True)
    stamp=_stamp(); plan_file=orch/"pipeline_plan.json"; run_file=orch/f"pipeline_run.{stamp}.json"; current=orch/"current_job.json"; log_file=logs/f"pipeline_run.{stamp}.log"
    plan_file.write_text(json.dumps({"project_id":context["config"]["project_id"],"profile":profile,"target":target,"analysis_id":analysis_id,"plan":plan},indent=2)+"\n")
    started=time.monotonic(); results=[]; warnings=[]
    if "sample_normalization" in [x["stage"] for x in plan]:
        status=build_project_status(project_config); diag=next(x for x in status["stages"] if x["stage"]=="library_diagnostics")
        if diag["status"]!="success": warnings.append("Normalization requested before raw LibraryDiagnostics completed.")
    with log_file.open("w",encoding="utf-8") as log:
      for item in plan:
        command=_command(item["stage"],context,target,analysis_id,threads,force,dry_run,keep_going)
        job={"project_id":context["config"]["project_id"],"profile":profile,"current_stage":item["stage"],"command":command,
             "pid":os.getpid(),"status":"running","start_time":_now(),"last_update":_now(),"cancel_requested":False,"message":""}
        current.write_text(json.dumps(job,indent=2)+"\n")
        log.write(f"\n=== {item['stage']} ===\nCommand: {' '.join(command)}\n"); log.flush()
        result=execute_pipeline_stage(item,context,threads,force,dry_run,keep_going,target,analysis_id,log); results.append(result)
        write_project_status(project_config,build_project_status(project_config))
        if result["status"]=="failed" and not keep_going: break
    final="failed" if any(r["status"]=="failed" for r in results) else ("dry-run" if dry_run else "success")
    job.update(status=final,current_stage=None,last_update=_now(),message="Pipeline run completed." if final!="failed" else "One or more stages failed.")
    current.write_text(json.dumps(job,indent=2)+"\n")
    run={"project_id":context["config"]["project_id"],"species":context["config"]["species"],"profile":profile,"target":target,
         "analysis_id":analysis_id,"threads":threads,"status":final,"plan_file":str(plan_file),"log_file":str(log_file),"results":results,
         "warnings":warnings,"start_time":job["start_time"],"end_time":_now(),"elapsed_seconds":round(time.monotonic()-started,3)}
    run_file.write_text(json.dumps(run,indent=2)+"\n"); run["run_metadata_file"]=str(run_file); return run
