#!/usr/bin/env python3
import argparse, csv, sys
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
    print(f"Comparison design: {path}")
    if args.print_design:
        for row in rows: print(f"{row['comparison_id']}: {row['parent_sample']} -> {row['treated_sample']} (background={row['background']}, pool={row['pool']}, include={row['include']})")
    for issue in result["issues"]: print(f"ISSUE {issue['issue']}: {issue['message']}")
    print("PASS" if result["valid"] else "FAIL"); return 0 if result["valid"] else 1
if __name__ == "__main__": raise SystemExit(main())
