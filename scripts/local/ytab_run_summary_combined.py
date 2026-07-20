#!/usr/bin/env python3
"""Run SummaryTable on one combined normalized parent hit file."""

from __future__ import annotations

import argparse
import shlex
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
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        summary = run_summary_combined_project(
            args.project_config, target=args.target, threads=args.threads,
            force=args.force, dry_run=args.dry_run,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"Project id: {summary['project_id']}")
    print(f"Species: {summary['species']}")
    print(f"Target: {summary['target']}")
    print(f"Target tag: {summary['target_tag']}")
    print(f"Input combined hits: {summary['input_combined_hits_file']}")
    print("Command: " + shlex.join(summary["command_run"]))
    print(f"Stable combined feature table: {summary.get('stable_combined_feature_table') or 'not produced'}")
    print(f"Summary stats: {summary.get('summary_stats_file') or 'not produced'}")
    print(f"Status: {summary['status']}")
    print(f"Manifest: {summary['manifest_path']}")
    if summary.get("warnings"):
        print("Warnings: " + "; ".join(summary["warnings"]))
    if summary.get("error_message"):
        print(f"ERROR: {summary['error_message']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
