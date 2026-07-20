#!/usr/bin/env python3
"""Run one SummaryTable smoke sample without downstream analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.summary_table_runner import run_summary_table_project


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        summary = run_summary_table_project(
            args.project_config, samples=[args.sample], threads=args.threads,
            force=args.force,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    result = summary["results"][0]
    manifest = Path(result["output_dir"]).parents[1] / "manifests" / "summary" / f"{args.sample}.summary_manifest.json"
    print(f"Input hits file: {result['input_hits_file']}")
    print(f"Output directory: {result['output_dir']}")
    print("Detected feature tables: " + (", ".join(result["detected_feature_tables"]) or "none"))
    print(f"Manifest path: {manifest}")
    if result["status"] in {"success", "skipped"}:
        print(f"PASS: SummaryTable {result['status']}")
        return 0
    print(f"FAIL: {result.get('error_message') or 'SummaryTable failed'}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
