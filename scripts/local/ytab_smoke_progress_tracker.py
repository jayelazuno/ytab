#!/usr/bin/env python3
from pathlib import Path
import tempfile,time,sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from ytab.pipeline.progress_tracker import *
def main():
 with tempfile.TemporaryDirectory(prefix="YTAB progress path with spaces ") as d:
  path=Path(d)/"job progress.json";t=ProgressTracker.create(path,"test_job","project","mapfastq",["a","b","c","d"],input_files={"a":__file__,"b":__file__,"c":__file__,"d":__file__});assert path.is_file();assert load_progress_state(path)["eta_seconds"] is None
  t.start();t.start_item("a",1);time.sleep(.02);t.finish_item("a","success");s=load_progress_state(path);assert s["eta_seconds"] is not None and s["eta_confidence"]=="low" and s["progress_fraction"]==.25
  t.start_item("b",2);t.finish_item("b","skipped");t.start_item("c",3);t.finish_item("c","failed",error_message="test");assert load_progress_state(path)["processed_items"]==3
  t.start_item("d",4);t.cancel();assert load_progress_state(path)["status"]=="cancelled";t.finalize("partial");assert load_progress_state(path)["status"]=="partial"
  path.write_text("{malformed",encoding="utf-8");assert load_progress_state(path) is None
 print("PASS");return 0
if __name__=="__main__":raise SystemExit(main())
