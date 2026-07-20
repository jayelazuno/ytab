#!/usr/bin/env python3
import argparse, shlex, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT/"src"))
from ytab.pipeline.treated_vs_parent_runner import run_treated_vs_parent_project

def main():
    p=argparse.ArgumentParser(description="Run raw-SummaryTable CPM treated-versus-parent fitness analysis.")
    p.add_argument("--project-config",type=Path,required=True); p.add_argument("--analysis-id",default="treated_vs_parent_raw_cpm")
    p.add_argument("--comparison-design",type=Path); p.add_argument("--comparisons",help="Comma-separated comparison or pair IDs")
    p.add_argument("--input-mode",choices=("auto","raw-summary"),default="auto"); p.add_argument("--annotate-classifier",action="store_true")
    p.add_argument("--classifier-target"); p.add_argument("--force",action="store_true"); p.add_argument("--dry-run",action="store_true")
    p.add_argument("--keep-going",action="store_true"); p.add_argument("--extra-r-arg",action="append",default=[]); a=p.parse_args()
    try:
        result=run_treated_vs_parent_project(a.project_config,a.analysis_id,a.comparison_design,
          a.comparisons.split(",") if a.comparisons else None,a.input_mode,a.classifier_target,a.annotate_classifier,
          a.force,a.dry_run,a.keep_going,a.extra_r_arg)
    except Exception as exc: print(f"ERROR: {exc}",file=sys.stderr); return 1
    print(f"Project ID: {result['project_id']}"); print(f"Species: {result['species']}"); print(f"Analysis ID: {result['analysis_id']}")
    print(f"Comparison design: {result['comparison_design_file']}"); print(f"Requested input mode: {a.input_mode}"); print("Resolved input mode: raw-summary")
    for sample,path in result["resolved_inputs"]["sample_files"].items(): print(f"Raw SummaryTable: {sample} = {path}")
    inspection=result["script_inspection"]; print(f"CPM normalization detected: {inspection['cpm_normalization_detected']}")
    print("CPM evidence: "+"; ".join(inspection["cpm_detection_evidence"])); print(f"Classifier annotation requested: {a.annotate_classifier}")
    print("Command: "+shlex.join(result["command_run"])); print(f"Status: {result['status']}"); print(f"Output directory: {result['output_dir']}")
    if result.get("stable_result_table"): print(f"Stable result table: {result['stable_result_table']}"); print(f"Result rows: {result.get('result_row_count','')}"); print("Scientific columns: "+", ".join(result.get("scientific_columns",[])))
    return 0 if result["status"] in {"success","skipped","dry-run"} else 1
if __name__=="__main__": raise SystemExit(main())
