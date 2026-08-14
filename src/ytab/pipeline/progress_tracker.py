"""Atomic, reusable progress records for local YTAB jobs."""
from __future__ import annotations

import json, os, secrets, statistics, tempfile
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1
FINAL = {"success", "cached", "dry_run_success", "failed", "cancelled", "partial"}

def now_iso() -> str: return datetime.now(timezone.utc).isoformat()
def parse_time(value: str | None) -> datetime | None:
    if not value: return None
    try: return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except (TypeError, ValueError): return None
def make_job_id(stage: str) -> str:
    stamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{stage}_{stamp}_{secrets.token_hex(3)}"

def atomic_write_json(path: Path, data: dict[str, Any]) -> None:
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix=f".{path.name}.",suffix=".tmp",dir=path.parent)
    try:
        with os.fdopen(fd,"w",encoding="utf-8") as handle:
            json.dump(data,handle,indent=2);handle.write("\n");handle.flush();os.fsync(handle.fileno())
        os.replace(tmp,path)
    finally:
        try: Path(tmp).unlink(missing_ok=True)
        except OSError: pass

def load_progress_state(path: Path) -> dict[str, Any] | None:
    try:
        value=json.loads(Path(path).read_text(encoding="utf-8"))
        return value if isinstance(value,dict) else None
    except (OSError,json.JSONDecodeError): return None

def _item(name: str,input_file: str="") -> dict[str,Any]:
    size=0
    try: size=Path(input_file).stat().st_size if input_file else 0
    except OSError: pass
    return {"item":name,"input_file":input_file,"input_size_bytes":size,"status":"queued","started_at":None,"finished_at":None,"elapsed_seconds":0.0,"output_files":[],"message":"","error_message":""}

def create_progress_state(job_id: str,project_id: str,stage: str,items: list[str],command: list[str]|None=None,input_files: dict[str,str]|None=None,profile: str|None=None) -> dict[str,Any]:
    timestamp=now_iso();records=[_item(name,(input_files or {}).get(name,"")) for name in items]
    state={"schema_version":SCHEMA_VERSION,"job_id":job_id,"project_id":project_id,"stage":stage,"profile":profile,"status":"queued","command":command or [],"selected_items":items,"total_items":len(items),"processed_items":0,"successful_items":0,"skipped_items":0,"failed_items":0,"remaining_items":len(items),"current_item":None,"current_item_index":None,"current_item_started_at":None,"current_phase":"queued","job_started_at":timestamp,"updated_at":timestamp,"job_elapsed_seconds":0.0,"current_item_elapsed_seconds":0.0,"progress_fraction":0.0,"progress_percent":0.0,"eta_seconds":None,"estimated_completion_time":None,"eta_method":"completed items in current run","eta_confidence":"unavailable","message":"Queued","cancel_requested":False,"pid":os.getpid(),"items":records}
    return refresh_progress(state)

def estimate_remaining_time(state: dict[str,Any]) -> tuple[float|None,str]:
    completed=[x for x in state["items"] if x["status"]=="success" and (x.get("elapsed_seconds") or 0)>0]
    if not completed:return None,"unavailable"
    rates=[x["elapsed_seconds"]/x["input_size_bytes"] for x in completed if (x.get("input_size_bytes") or 0)>0]
    durations=[x["elapsed_seconds"] for x in completed]
    rate=statistics.median(rates) if rates else None;duration=statistics.median(durations)
    remaining=0.0
    for item in state["items"]:
        if item["status"] not in {"queued","running","starting"}:continue
        expected=(item.get("input_size_bytes") or 0)*rate if rate and item.get("input_size_bytes") else duration
        if item["status"] in {"running","starting"}:expected=max(expected-(state.get("current_item_elapsed_seconds") or 0),0)
        remaining+=expected
    n=len(completed);confidence="low" if n==1 else "medium" if n<=3 else "high"
    if len(rates)>=2 and statistics.mean(rates)>0 and statistics.pstdev(rates)/statistics.mean(rates)>0.75:confidence="low"
    return round(remaining,1),confidence

def refresh_progress(state: dict[str,Any]) -> dict[str,Any]:
    now=datetime.now(timezone.utc);start=parse_time(state.get("job_started_at"));current=parse_time(state.get("current_item_started_at"))
    state["updated_at"]=now.isoformat();state["job_elapsed_seconds"]=round((now-start).total_seconds(),1) if start else 0.0;state["current_item_elapsed_seconds"]=round((now-current).total_seconds(),1) if current else 0.0
    statuses=[x["status"] for x in state["items"]];state["successful_items"]=statuses.count("success");state["skipped_items"]=statuses.count("skipped");state["failed_items"]=statuses.count("failed");state["processed_items"]=sum(s in {"success","skipped","failed","cancelled"} for s in statuses);state["remaining_items"]=max(state["total_items"]-state["processed_items"],0);state["progress_fraction"]=state["processed_items"]/state["total_items"] if state["total_items"] else 1.0;state["progress_percent"]=round(100*state["progress_fraction"],1)
    eta,confidence=estimate_remaining_time(state);state["eta_seconds"]=eta;state["eta_confidence"]=confidence;state["estimated_completion_time"]=(now+timedelta(seconds=eta)).isoformat() if eta is not None else None
    return state

class ProgressTracker:
    def __init__(self,path:Path,state:dict[str,Any]):self.path=Path(path);self.state=state;self.write()
    @classmethod
    def create(cls,path:Path,job_id:str,project_id:str,stage:str,items:list[str],**kwargs):return cls(path,create_progress_state(job_id,project_id,stage,items,**kwargs))
    def write(self):refresh_progress(self.state);atomic_write_json(self.path,self.state);return self.state
    def start(self,message="Running"):self.state.update(status="running",current_phase="starting",message=message);return self.write()
    def start_item(self,item:str,index:int,phase="starting"):
        row=next(x for x in self.state["items"] if x["item"]==item);row.update(status="running",started_at=now_iso(),message="Running")
        self.state.update(status="running",current_item=item,current_item_index=index,current_item_started_at=row["started_at"],current_phase=phase,message=f"Processing {item}");return self.write()
    def phase(self,phase:str,message:str=""):self.state.update(current_phase=phase,message=message or phase);return self.write()
    def finish_item(self,item:str,status="success",output_files=None,message="",error_message=""):
        row=next(x for x in self.state["items"] if x["item"]==item);finished=now_iso();started=parse_time(row.get("started_at"));row.update(status=status,finished_at=finished,elapsed_seconds=round((parse_time(finished)-started).total_seconds(),3) if started else 0.0,output_files=output_files or [],message=message,error_message=error_message)
        self.state.update(current_item=None,current_item_started_at=None,current_item_elapsed_seconds=0,current_phase="between items",message=message or f"{item}: {status}");return self.write()
    def heartbeat(self):return self.write()
    def cancel(self,message="Cancellation requested"):
        self.state.update(status="cancelled",cancel_requested=True,current_phase="cancelled",message=message);return self.write()
    def finalize(self,status:str|None=None,message=""):
        if status is None:status="failed" if self.state["failed_items"] else "success"
        self.state.update(status=status,current_item=None,current_item_started_at=None,current_phase="complete" if status=="success" else status,message=message or status.title());return self.write()
