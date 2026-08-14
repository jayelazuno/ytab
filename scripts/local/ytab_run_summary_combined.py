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
from ytab.pipeline.progress_tracker import ProgressTracker, make_job_id
from ytab.pipeline.summary_combined_runner import (
    load_project_for_summary_combined,
    resolve_summary_combined_target,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--target", default="recommended")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--job-id")
    parser.add_argument("--progress-file", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    tracker = None
    try:
        config = load_project_for_summary_combined(args.project_config)
        target = resolve_summary_combined_target(config, args.target)
        job_id = args.job_id or make_job_id("summary_combined")
        progress_file = args.progress_file or args.project_config.resolve().parents[1] / "manifests" / "jobs" / f"{job_id}.progress.json"
        tracker = ProgressTracker.create(progress_file, job_id, str(config.get("project_id")),
                                         "summary_combined", [target["target_tag"]])
        tracker.state.update(dry_run=args.dry_run,
                             execution_mode="preview" if args.dry_run else "run",
                             target=target["target"], target_tag=target["target_tag"])
        tracker.start("Validating combined hit file")
        tracker.start_item(target["target_tag"], 1, "validating combined hit file")
        tracker.phase("loading reference annotation", "Loading configured reference annotation")
        summary = run_summary_combined_project(
            args.project_config, target=args.target, threads=args.threads,
            force=args.force, dry_run=args.dry_run,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        if tracker is not None:
            tracker.finish_item(tracker.state["items"][0]["item"], "failed", error_message=str(exc))
            tracker.finalize("failed", str(exc))
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
        tracker.finish_item(summary["target_tag"], "failed",
                            error_message=summary["error_message"])
        tracker.finalize("failed", summary["error_message"])
        print(f"ERROR: {summary['error_message']}", file=sys.stderr)
        return 1
    tracker.phase("selecting stable feature table", "Validating stable combined feature table")
    tracker.finish_item(summary["target_tag"],
                        "skipped" if summary["status"] == "skipped" else "success",
                        [summary.get("stable_combined_feature_table")] if summary.get("stable_combined_feature_table") else [],
                        "Preview complete; SummaryTable was not executed." if args.dry_run else
                        "Cached feature table reused." if summary["status"] == "skipped" else
                        "Combined parent feature table complete")
    final = "dry_run_success" if args.dry_run else "cached" if summary["status"] == "skipped" else "success"
    tracker.finalize(final,
                     "Preview complete" if args.dry_run else
                     "Cached combined feature table reused" if final == "cached" else
                     "Combined SummaryTable complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
