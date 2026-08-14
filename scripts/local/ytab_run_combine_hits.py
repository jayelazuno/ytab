#!/usr/bin/env python3
"""Combine normalized parent hit files for one selected normalization target."""

from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.combine_hits_runner import (
    find_normalized_parent_hits,
    get_parent_samples,
    load_project_for_combine_hits,
    resolve_combine_target,
    run_combine_hits_project,
)
from ytab.pipeline.progress_tracker import ProgressTracker, make_job_id


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--target", default="recommended")
    parser.add_argument("--samples", help="Comma-separated included parent sample names")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--job-id")
    parser.add_argument("--progress-file", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    samples = [part.strip() for part in args.samples.split(",") if part.strip()] if args.samples else None
    tracker = None
    try:
        config = load_project_for_combine_hits(args.project_config)
        selected = get_parent_samples(config, samples)
        target = resolve_combine_target(config, args.target)
        inputs = find_normalized_parent_hits(config, target["target_tag"], selected)
        job_id = args.job_id or make_job_id("combined_hits")
        progress_file = args.progress_file or args.project_config.resolve().parents[1] / "manifests" / "jobs" / f"{job_id}.progress.json"
        tracker = ProgressTracker.create(progress_file, job_id, str(config.get("project_id")),
                                         "combined_hits", [target["target_tag"]])
        tracker.state.update(dry_run=args.dry_run,
                             execution_mode="preview" if args.dry_run else "run",
                             target=target["target"], target_tag=target["target_tag"],
                             parent_samples=[str(row["sample"]) for row in selected],
                             parent_count=len(selected))
        tracker.start("Validating normalized hit files")
        tracker.start_item(target["target_tag"], 1, "validating normalized hit files")
        tracker.phase("loading parent hits", "Loading normalized parent hit files")
        summary = run_combine_hits_project(
            args.project_config, target=args.target, samples=samples,
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
    print("Selected parent samples: " + ", ".join(summary["selected_samples"]))
    print("Input normalized hit files:")
    for path in inputs:
        print(f"  {path}")
    print(f"Output combined hit file: {summary['combined_hits_file']}")
    print("Command: " + shlex.join(summary["command_run"]))
    print(f"Status: {summary['status']}")
    print(f"Manifest: {summary['manifest_path']}")
    if summary.get("total_combined_sites") is not None:
        print(f"Total combined sites: {summary['total_combined_sites']}")
        print(f"Total combined reads: {summary['total_combined_reads']}")
    if summary.get("warnings"):
        print("Warnings: " + "; ".join(summary["warnings"]))
    if summary.get("error_message"):
        tracker.finish_item(summary["target_tag"], "failed",
                            error_message=summary["error_message"])
        tracker.finalize("failed", summary["error_message"])
        print(f"ERROR: {summary['error_message']}", file=sys.stderr)
        return 1
    tracker.phase("writing combined hit file", "Validating stable combined parent hit output")
    tracker.finish_item(summary["target_tag"],
                        "skipped" if summary["status"] == "skipped" else "success",
                        [summary["combined_hits_file"]],
                        "Preview complete; parent hits were not combined." if args.dry_run else
                        "Cached combined parent library reused." if summary["status"] == "skipped" else
                        "Parent libraries combined")
    final = "dry_run_success" if args.dry_run else "cached" if summary["status"] == "skipped" else "success"
    tracker.finalize(final,
                     "Preview complete" if args.dry_run else
                     "Cached combined parent library reused" if final == "cached" else
                     "Combine Parent Libraries complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
