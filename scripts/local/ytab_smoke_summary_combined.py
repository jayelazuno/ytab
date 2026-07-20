#!/usr/bin/env python3
"""Smoke-test SummaryTable on one combined normalized parent hit file."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.summary_combined_runner import run_summary_combined_project


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--target", default="recommended")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        summary = run_summary_combined_project(
            args.project_config, target=args.target, threads=args.threads, force=args.force,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print(f"Selected target: {summary['target']} ({summary['target_tag']})")
    print(f"Input combined hits: {summary['input_combined_hits_file']}")
    print(f"Stable combined feature table: {summary.get('stable_combined_feature_table') or 'not produced'}")
    print(f"Summary stats: {summary.get('summary_stats_file') or 'not produced'}")
    print(f"Manifest path: {summary['manifest_path']}")
    stable = summary.get("stable_combined_feature_table")
    if summary["status"] in {"success", "skipped"} and stable and Path(stable).is_file():
        print("PASS: combined normalized parent SummaryTable")
        return 0
    print(f"FAIL: {summary.get('error_message') or 'stable combined feature table was not detected'}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
