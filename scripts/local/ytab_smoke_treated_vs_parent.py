#!/usr/bin/env python3
import argparse, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT/"src"))
from ytab.pipeline.treated_vs_parent_runner import run_treated_vs_parent_project
def main():
    p=argparse.ArgumentParser(description="Run one comparison through treated-versus-parent analysis.")
    p.add_argument("--project-config",type=Path,required=True); p.add_argument("--comparison-id",required=True)
    p.add_argument("--input-mode",choices=("auto","raw-summary"),default="auto"); p.add_argument("--force",action="store_true"); a=p.parse_args()
    try: result=run_treated_vs_parent_project(a.project_config,"smoke_"+a.comparison_id,comparisons=[a.comparison_id],input_mode=a.input_mode,force=a.force)
    except Exception as exc: print(f"FAIL: {exc}",file=sys.stderr); return 1
    row=result["selected_comparisons"][0]; print(f"Parent: {row['parent_sample']}\nTreated: {row['treated_sample']}")
    for sample,path in result["resolved_inputs"]["sample_files"].items(): print(f"Raw SummaryTable: {sample} = {path}")
    print(f"CPM normalization detected: {result['script_inspection']['cpm_normalization_detected']}")
    print(f"Output directory: {result['output_dir']}\nStable result table: {result.get('stable_result_table','')}\nRows: {result.get('result_row_count','')}\nManifest: {result.get('manifest_file','')}\n"+("PASS" if result["status"] in {"success","skipped"} else "FAIL"))
    return 0 if result["status"] in {"success","skipped"} else 1
if __name__=="__main__": raise SystemExit(main())
