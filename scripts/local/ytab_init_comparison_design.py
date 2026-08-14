#!/usr/bin/env python3
import argparse, csv, hashlib, json, sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from ytab.pipeline.treated_vs_parent_runner import (infer_comparison_design, load_comparison_design,
    load_project_for_treated_vs_parent, validate_comparison_design)

def main():
    parser=argparse.ArgumentParser(description="Create and validate a parent-treated comparison design.")
    parser.add_argument("--project-config", type=Path, required=True); parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--print-design", action="store_true"); args=parser.parse_args()
    config=load_project_for_treated_vs_parent(args.project_config); path=infer_comparison_design(config,args.overwrite)
    rows=load_comparison_design(path); result=validate_comparison_design(config,rows)
    root=Path(config["_repo_root"]);sample_sheet=Path(config.get("sample_sheet") or "")
    if not sample_sheet.is_absolute():sample_sheet=(root/sample_sheet).resolve()
    summary_refs=[]
    for row in rows:
        for key in ("parent_sample","treated_sample"):
            sample=row.get(key,"");directory=Path(config["output_project_dir"])/"summary"/sample
            if not directory.is_absolute():directory=(root/directory).resolve()
            summary_refs.extend(str(item) for item in sorted(directory.glob("*.feature_table*")) if item.is_file())
    sha=lambda value:hashlib.sha256(Path(value).read_bytes()).hexdigest() if Path(value).is_file() else None
    metadata={"project_id":config.get("project_id"),"generated_at":datetime.now(timezone.utc).isoformat(),"design_version":1,"valid_comparison_count":len(rows) if result["valid"] else 0,"included_comparison_count":result["included_count"],"design_issues":result["issues"],"sample_sheet_sha256":sha(sample_sheet),"raw_summary_inputs":summary_refs}
    path.with_name("comparison_design.metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(f"Comparison design: {path}")
    if args.print_design:
        for row in rows: print(f"{row['comparison_id']}: {row['parent_sample']} -> {row['treated_sample']} (background={row['background']}, pool={row['pool']}, include={row['include']})")
    for issue in result["issues"]: print(f"ISSUE {issue['issue']}: {issue['message']}")
    print("PASS" if result["valid"] else "FAIL"); return 0 if result["valid"] else 1
if __name__ == "__main__": raise SystemExit(main())
