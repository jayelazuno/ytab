#!/usr/bin/env python3
"""Non-scientific worker used to validate persistent YTAB progress."""
from __future__ import annotations
import argparse,time
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from ytab.pipeline.progress_tracker import ProgressTracker,make_job_id
def main()->int:
 p=argparse.ArgumentParser();p.add_argument("--progress-file",type=Path,required=True);p.add_argument("--items",required=True);p.add_argument("--sleep-seconds",type=float,default=.1);p.add_argument("--fail-item");p.add_argument("--job-id");a=p.parse_args();items=[x for x in a.items.split(",") if x];tracker=ProgressTracker.create(a.progress_file,a.job_id or make_job_id("test"),"test_project","test_worker",items);tracker.start("Test worker running")
 for i,item in enumerate(items,1):
  tracker.start_item(item,i,"simulating work");time.sleep(a.sleep_seconds)
  if item==a.fail_item:tracker.finish_item(item,"failed",error_message="Simulated failure")
  else:tracker.finish_item(item,"success",message="Simulated success")
 tracker.finalize("partial" if a.fail_item else "success");return 0
if __name__=="__main__":raise SystemExit(main())
