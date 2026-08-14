"""Read-only, manifest-first project status inspection."""
from __future__ import annotations
import csv, json
from datetime import datetime, timezone
from pathlib import Path
from .mapfastq_runner import load_project_for_mapping, get_included_samples

STAGES=("project_config","reference","mapfastq","create_hit_file","summary","library_diagnostics",
        "sample_normalization","summary_normalized","combined_hits","summary_combined","classifier",
        "comparison_design","treated_vs_parent")
BRANCH={"project_config":"setup","reference":"setup","mapfastq":"core","create_hit_file":"core",
        "summary":"core","library_diagnostics":"qc","sample_normalization":"classifier",
        "summary_normalized":"classifier","combined_hits":"classifier","summary_combined":"classifier",
        "classifier":"classifier","comparison_design":"fitness","treated_vs_parent":"fitness"}
DEPENDENCIES={"reference":["project_config"],"mapfastq":["reference"],"create_hit_file":["mapfastq"],
 "summary":["create_hit_file"],"library_diagnostics":["create_hit_file"],"sample_normalization":["create_hit_file"],
 "summary_normalized":["sample_normalization"],"combined_hits":["summary_normalized"],
 "summary_combined":["combined_hits"],"classifier":["summary_combined"],
 "comparison_design":["summary"],"treated_vs_parent":["summary","comparison_design"]}
SAMPLE_STATUS={"mapfastq":"mapfastq/mapfastq_status.csv","create_hit_file":"create_hit_file/create_hit_file_status.csv",
               "summary":"summary/summary_status.csv","library_diagnostics":"library_diagnostics/library_diagnostics_status.csv"}
TARGET_STATUS={"sample_normalization":"sample_normalization/sample_normalization_status.csv",
 "summary_normalized":"summary_normalized/summary_normalized_status.csv","combined_hits":"combined_hits/combined_hits_status.csv",
 "summary_combined":"summary_combined/summary_combined_status.csv","classifier":"classifier/classifier_status.csv",
 "treated_vs_parent":"treated_vs_parent/treated_vs_parent_status.csv"}
COMPLETE_STATUSES={"success","skipped","cached","imported_success","imported_or_not_required"}

def _now(): return datetime.now(timezone.utc).isoformat()
def _bool(v): return v if isinstance(v,bool) else str(v).lower() in {"true","1","yes","y"}
def _rows(path):
    if not path.is_file(): return []
    with path.open(newline="",encoding="utf-8") as h: return list(csv.DictReader(h))
def _item(stage,status="not_started",**kw):
    base={"stage":stage,"branch":BRANCH[stage],"status":status,"completed":0,"expected":0,"successful":0,
          "failed":0,"skipped":0,"blocked_reason":"","primary_output":"","manifest":"","last_updated":"","message":""}
    base.update(kw); return base

def _is_complete_status(status):
    return str(status or "").lower() in COMPLETE_STATUSES

def _start_stage(context):
    return str((context.get("config") or {}).get("start_stage") or "fastq").strip().lower()

def _sample_hit_files_complete(context):
    complete=0
    for sample in context["samples"]:
        name=str(sample.get("sample") or "")
        hits=context["project"]/"create_hit_file"/name/f"{name}_hits.txt"
        if hits.is_file() and hits.stat().st_size>0:
            complete+=1
    return complete

def load_project_status_context(project_config: Path)->dict:
    config=load_project_for_mapping(project_config); root=Path(config["_repo_root"])
    project=(root/str(config["output_project_dir"])).resolve() if not Path(str(config["output_project_dir"])).is_absolute() else Path(config["output_project_dir"])
    export=(root/str(config["output_export_dir"])).resolve() if not Path(str(config["output_export_dir"])).is_absolute() else Path(config["output_export_dir"])
    return {"config":config,"root":root,"project":project,"export":export,"manifests":project/"manifests",
            "samples":get_included_samples(config),"project_config":Path(config["_project_config"])}

def inspect_reference_status(context):
    ref=context["config"].get("reference") or {}; fasta=Path(ref.get("fasta") or ""); prefix=Path(ref.get("bowtie2_index_prefix") or "")
    root=context["root"]; fasta=fasta if fasta.is_absolute() else root/fasta
    prefix=prefix if prefix.is_absolute() else root/prefix
    index=any(all(Path(str(prefix)+s+e).is_file() for s in (".1",".2",".3",".4",".rev.1",".rev.2")) for e in (".bt2",".bt2l"))
    ok=fasta.is_file() and index; manifest=context["project"]/"config"/"reference_resolved.json"
    return _item("reference","success" if ok else ("blocked" if not fasta.is_file() else "ready"),completed=int(ok),expected=1,
      successful=int(ok),blocked_reason="Reference FASTA is missing." if not fasta.is_file() else "",primary_output=str(fasta),manifest=str(manifest),
      last_updated=datetime.fromtimestamp(manifest.stat().st_mtime,timezone.utc).isoformat() if manifest.is_file() else "")

def _summarize(stage,rows,expected,primary,manifest):
    statuses=[str(r.get("status","")).lower() for r in rows]
    completed=[_is_complete_status(s) and "dry run" not in str(r.get("message","")).lower() for s,r in zip(statuses,rows)]
    good=sum(completed); failed=sum(s=="failed" for s in statuses)
    if failed: state="failed"
    elif expected and good>=expected and statuses and all(s=="imported_or_not_required" for s in statuses): state="imported_or_not_required"
    elif expected and good>=expected and statuses and all(s=="imported_success" for s in statuses): state="imported_success"
    elif expected and good>=expected: state="success"
    elif good: state="partial"
    elif rows: state="not_started"
    else: state="not_started"
    return _item(stage,state,completed=good,expected=expected,successful=sum(s=="success" for s in statuses),failed=failed,
      skipped=sum(s in {"skipped","imported_or_not_required"} for s in statuses),primary_output=str(primary),manifest=str(manifest),
      last_updated=datetime.fromtimestamp(manifest.stat().st_mtime,timezone.utc).isoformat() if manifest.is_file() else "")

def inspect_sample_stage_status(context,stage_name):
    status=context["manifests"]/SAMPLE_STATUS[stage_name]; rows=_rows(status); expected=len(context["samples"])
    output=context["project"]/{"mapfastq":"mapfastq","create_hit_file":"create_hit_file","summary":"summary","library_diagnostics":"library_diagnostics"}[stage_name]
    if stage_name=="library_diagnostics": expected=1
    if _start_stage(context)=="create_hit_file" and stage_name=="mapfastq":
        if rows:
            return _summarize(stage_name,rows,expected,output,status)
        return _item(stage_name,"imported_or_not_required",completed=expected,expected=expected,
          successful=0,skipped=expected,primary_output=str(output),manifest=str(status),
          message="Project starts from existing CreateHitFile outputs; MapFastq is not required.")
    if _start_stage(context)=="create_hit_file" and stage_name=="create_hit_file" and not rows:
        complete=_sample_hit_files_complete(context)
        return _item(stage_name,"imported_success" if expected and complete>=expected else "blocked",
          completed=complete,expected=expected,successful=complete,primary_output=str(output),manifest=str(status),
          blocked_reason="" if complete>=expected else "Imported hit files are missing.")
    return _summarize(stage_name,rows,expected,output,status)

def inspect_target_stage_status(context,stage_name):
    status=context["manifests"]/TARGET_STATUS[stage_name]; rows=_rows(status)
    keys=[]
    for row in rows:
        key=row.get("target_tag") or row.get("target") or row.get("analysis_id") or row.get("sample")
        if key: keys.append(key)
    tags=sorted(set(keys)); grouped={tag:[str(r.get("status","")).lower() for r in rows if (r.get("target_tag") or r.get("target") or r.get("analysis_id") or r.get("sample"))==tag] for tag in tags}
    completed=sum(bool(values) and all(_is_complete_status(v) for v in values) for values in grouped.values())
    failed=sum(any(v=="failed" for v in values) for values in grouped.values()); expected=max(1,len(tags))
    state="failed" if failed else ("success" if completed>=expected else ("partial" if completed else "not_started"))
    result=_item(stage_name,state,completed=completed,expected=expected,successful=completed,failed=failed,
      skipped=sum(all(v=="skipped" for v in values) for values in grouped.values()),primary_output=str(context["project"]/{"sample_normalization":"sample_normalization","summary_normalized":"summary_normalized","combined_hits":"combined_hits","summary_combined":"summary_combined","classifier":"classifier","treated_vs_parent":"treated_vs_parent"}[stage_name]),manifest=str(status),last_updated=datetime.fromtimestamp(status.stat().st_mtime,timezone.utc).isoformat() if status.is_file() else "")
    result["message"]=("Available analysis IDs: " if stage_name=="treated_vs_parent" else "Available targets: ")+", ".join(sorted(set(keys))) if keys else ""
    return result

def inspect_project_stage_status(context,stage_name):
    if stage_name=="project_config": return _item(stage_name,"success",completed=1,expected=1,successful=1,primary_output=str(context["project_config"]),manifest=str(context["project_config"]))
    if stage_name=="comparison_design":
        path=context["project"]/"config"/"comparison_design.csv"; issues=context["project"]/"config"/"comparison_design_issues.csv"
        bad=len(_rows(issues)); state="success" if path.is_file() and not bad else ("failed" if bad else "not_started")
        return _item(stage_name,state,completed=int(state=="success"),expected=1,successful=int(state=="success"),failed=int(bool(bad)),primary_output=str(path),manifest=str(issues),message=f"{bad} unresolved issue(s)." if bad else "")
    raise ValueError(stage_name)

def build_project_status(project_config: Path)->dict:
    c=load_project_status_context(project_config); items=[]
    for stage in STAGES:
        if stage=="reference": item=inspect_reference_status(c)
        elif stage in SAMPLE_STATUS: item=inspect_sample_stage_status(c,stage)
        elif stage in TARGET_STATUS: item=inspect_target_stage_status(c,stage)
        else: item=inspect_project_stage_status(c,stage)
        items.append(item)
    lookup={x["stage"]:x for x in items}
    current=c["manifests"]/"orchestrator"/"current_job.json"
    if current.is_file():
        try:
            job=json.loads(current.read_text()); active=job.get("current_stage")
            if job.get("status")=="running" and active in lookup: lookup[active]["status"]="running"
        except (OSError,json.JSONDecodeError): pass
    for stage,deps in DEPENDENCIES.items():
        item=lookup[stage]
        if item["status"]=="not_started":
            missing=[d for d in deps if not _is_complete_status(lookup[d]["status"])]
            if missing: item["status"]="blocked"; item["blocked_reason"]="Requires: "+", ".join(missing)
            else: item["status"]="ready"
        if item["status"]=="success" and item["primary_output"] and not Path(item["primary_output"]).exists(): item["status"]="stale"; item["message"]="Manifest reports success but primary output is missing."
    warnings=[]
    if "smoke" in str(c["config"].get("project_id","")).lower(): warnings.append("This shallow smoke project validates software behavior; scientific calls are not final biological conclusions.")
    return {"project_id":c["config"].get("project_id"),"species":c["config"].get("species"),"project_config":str(c["project_config"]),"generated_at":_now(),"stages":items,"warnings":warnings}

def write_project_status(project_config: Path,status:dict)->dict:
    c=load_project_status_context(project_config); j=c["manifests"]/"project_status.json"; t=c["manifests"]/"project_status.csv"; j.parent.mkdir(parents=True,exist_ok=True)
    j.write_text(json.dumps(status,indent=2)+"\n",encoding="utf-8")
    with t.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=list(_item("project_config")),lineterminator="\n"); w.writeheader(); w.writerows(status["stages"])
    return {"json":j,"csv":t}

def determine_next_ready_stage(status,profile="all"):
    allowed={"core":{"reference","mapfastq","create_hit_file","summary","library_diagnostics"},"classifier":{"reference","mapfastq","create_hit_file","summary","library_diagnostics","sample_normalization","summary_normalized","combined_hits","summary_combined","classifier"},"fitness":{"reference","mapfastq","create_hit_file","summary","library_diagnostics","comparison_design","treated_vs_parent"},"all":set(STAGES)}[profile]
    return next((x["stage"] for x in status["stages"] if x["stage"] in allowed and x["status"]=="ready"),None)
