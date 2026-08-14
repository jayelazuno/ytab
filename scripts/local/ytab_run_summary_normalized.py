#!/usr/bin/env python3
"""Run SummaryTable on normalized hit files and evaluate normalization targets."""

from __future__ import annotations

import argparse
import json
import shlex
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.summary_normalized_runner import (
    find_normalized_hits_file,
    get_summary_normalized_samples,
    load_project_for_summary_normalized,
    resolve_normalization_targets,
    run_summary_normalized_project,
)
from ytab.pipeline.progress_tracker import ProgressTracker, make_job_id


def retention(value: str) -> float:
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be numeric") from None
    if not 0 < number <= 1:
        raise argparse.ArgumentTypeError("must be greater than 0 and at most 1")
    return number


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--targets", default="recommended")
    parser.add_argument("--sample-mode", choices=("parents", "all", "treated"), default="parents")
    parser.add_argument("--samples", help="Comma-separated included sample names; overrides sample mode")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--keep-going", action="store_true")
    parser.add_argument("--min-feature-retention", type=retention, default=0.95)
    parser.add_argument("--job-id")
    parser.add_argument("--progress-file", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    samples = [part.strip() for part in args.samples.split(",") if part.strip()] if args.samples else None
    tracker = None
    try:
        config = load_project_for_summary_normalized(args.project_config)
        selected = get_summary_normalized_samples(config, args.sample_mode, samples)
        target_infos = resolve_normalization_targets(config, args.targets)
        print(f"Project id: {config.get('project_id')}")
        print(f"Species: {config.get('species')}")
        print(f"Sample mode: {'explicit' if samples is not None else args.sample_mode}")
        print("Selected samples: " + ", ".join(str(row["sample"]) for row in selected))
        print("Selected targets: " + ", ".join(str(row["target"]) for row in target_infos))
        print("Target tags: " + ", ".join(str(row["target_tag"]) for row in target_infos))
        print("Input normalized hit files:")
        for target in target_infos:
            for sample in selected:
                print(f"  {find_normalized_hits_file(config, target['target_tag'], sample)}")
        items = [f"{target['target_tag']} × {sample['sample']}"
                 for target in target_infos for sample in selected]
        job_id = args.job_id or make_job_id("summary_normalized")
        progress_file = args.progress_file or args.project_config.resolve().parents[1] / "manifests" / "jobs" / f"{job_id}.progress.json"
        tracker = ProgressTracker.create(progress_file, job_id, str(config.get("project_id")),
                                         "summary_normalized", ["feature-level target evaluation"])
        tracker.state.update(dry_run=args.dry_run,
                             execution_mode="preview" if args.dry_run else "run",
                             target_parent_items=items, target_count=len(target_infos),
                             parent_count=len(selected))
        tracker.start("Validating normalized parent hit files")
        tracker.start_item("feature-level target evaluation", 1, "validating normalized parent hit files")
        tracker.phase("running normalized SummaryTable", "Evaluating target × parent SummaryTable inputs")
        summary = run_summary_normalized_project(
            args.project_config, targets=args.targets, sample_mode=args.sample_mode,
            samples=samples, threads=args.threads, force=args.force, dry_run=args.dry_run,
            keep_going=args.keep_going, min_feature_retention=args.min_feature_retention,
        )
    except (FileNotFoundError, OSError, ValueError, json.JSONDecodeError) as exc:
        if tracker is not None:
            tracker.finish_item("feature-level target evaluation", "failed", error_message=str(exc))
            tracker.finalize("failed", str(exc))
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"Output normalized summary directory: {summary['output_dir']}")
    for result in summary["results"]:
        print(f"\n{result['target_tag']}/{result['sample']}: {result['status']}")
        print("Command: " + shlex.join(result["command_run"]))
        print("Feature tables: " + (", ".join(result["feature_tables"]) or "none"))
        print(f"Manifest: {result['manifest_path']}")
        if result.get("error_message"):
            print(f"ERROR: {result['error_message']}", file=sys.stderr)
    if summary["evaluation"]:
        print(f"Target evaluation: {summary['evaluation']['evaluation_path']}")
        print(f"Target summary: {summary['evaluation']['target_summary_path']}")
        feature = summary["evaluation"]["feature_recommendation"]
        print(f"Feature recommendation: {feature['recommended_target']} ({feature['recommended_target_tag']})")
        print(f"Feature recommendation file: {summary['evaluation']['feature_recommendation_path']}")
    print(
        f"Summary: success={summary['success']} failed={summary['failed']} skipped={summary['skipped']}"
    )
    final = "dry_run_success" if args.dry_run else "failed" if summary["failed"] else "success"
    item_status = "failed" if summary["failed"] else "skipped" if summary["success"] == 0 else "success"
    tracker.phase("collecting retention outputs", "Loading stored site- and feature-retention tables")
    tracker.finish_item("feature-level target evaluation", item_status,
                        [path for result in summary["results"] for path in result["feature_tables"]],
                        "Preview complete; SummaryTable was not executed." if args.dry_run else
                        "Feature-level target evaluation complete")
    if not args.dry_run and item_status == "skipped":
        final = "cached"
    tracker.finalize(final, "Preview complete" if args.dry_run else
                     "Cached target evaluation reused" if final == "cached" else
                     "Normalized SummaryTable complete" if not summary["failed"] else
                     "Normalized SummaryTable failed")
    return 0 if args.dry_run or summary["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
