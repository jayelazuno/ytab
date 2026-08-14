#!/usr/bin/env python3
"""Run essentiality classification on a combined parent feature table."""

from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.classifier_runner import (
    load_project_for_classifier,
    resolve_classifier_target,
    run_classifier_project,
    save_reviewed_classifier_target,
)
from ytab.pipeline.progress_tracker import ProgressTracker, make_job_id


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--target", default="recommended")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--save-final-target", action="store_true")
    parser.add_argument("--save-final-target-only", action="store_true",
                        help="Record an existing reviewed successful result without rerunning classification")
    parser.add_argument("--job-id")
    parser.add_argument("--progress-file", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.save_final_target_only:
        if args.dry_run or args.force or args.save_final_target:
            print("ERROR: --save-final-target-only cannot be combined with execution flags.", file=sys.stderr)
            return 2
        try:
            saved = save_reviewed_classifier_target(args.project_config, args.target)
        except (FileNotFoundError, OSError, ValueError) as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        print(f"Final classifier target: {saved['target_tag']}")
        print(f"Reviewed prediction rows: {saved['prediction_rows']}")
        print(f"Final target record: {saved['final_target_files']['text']}")
        return 0
    tracker = None
    try:
        config = load_project_for_classifier(args.project_config)
        target = resolve_classifier_target(config, args.target)
        job_id = args.job_id or make_job_id("classifier")
        progress_file = args.progress_file or args.project_config.resolve().parents[1] / "manifests" / "jobs" / f"{job_id}.progress.json"
        tracker = ProgressTracker.create(progress_file, job_id, str(config.get("project_id")),
                                         "classifier", [target["target_tag"]])
        tracker.state.update(dry_run=args.dry_run,
                             execution_mode="preview" if args.dry_run else "run",
                             target=target["target"], target_tag=target["target_tag"],
                             seed=0 if args.seed is None else args.seed)
        tracker.start("Validating inputs")
        tracker.start_item(target["target_tag"], 1, "validating inputs")
        tracker.phase("loading classifier resources", "Loading combined feature table and required classifier resources")
        tracker.phase("fitting classifier", "Running the validated essentiality classifier")
        result = run_classifier_project(
            args.project_config, target=args.target, seed=args.seed, force=args.force,
            dry_run=args.dry_run, save_final_target=args.save_final_target,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        if tracker is not None:
            tracker.finish_item(tracker.state["items"][0]["item"], "failed", error_message=str(exc))
            tracker.finalize("failed", str(exc))
        print(f"ERROR: {exc}", file=sys.stderr); return 2
    print(f"Project id: {result['project_id']}")
    print(f"Species: {result['species']}")
    print(f"Requested target: {args.target}")
    print(f"Resolved target: {result['target']} ({result['target_tag']})")
    print(f"Target resolution source: {result['target_resolution_source']}")
    print(f"Input combined feature table: {result['input_combined_feature_table']}")
    print(f"Input feature count: {result['input_feature_count']}")
    print("Classifier resources:")
    for name, entries in result["classifier_resources_resolved"].items():
        print(f"  {name}: " + ", ".join(item["path"] for item in entries))
    print(f"Classifier output directory: {result['output_dir']}")
    print(f"Save final target: {'yes' if args.save_final_target else 'no'}")
    print("Command: " + shlex.join(result["command_run"]))
    print(f"Status: {result['status']}")
    print(f"Stable prediction table: {result.get('stable_prediction_table') or 'not produced'}")
    print(f"Output prediction count: {result.get('output_prediction_count') if result.get('output_prediction_count') is not None else 'not available'}")
    if result.get("class_counts"):
        print("Prediction labels: " + ", ".join(f"{key}={value}" for key, value in sorted(result["class_counts"].items())))
    print(f"Manifest: {result['manifest_path']}")
    if result.get("warnings"): print("Warnings: " + "; ".join(result["warnings"]))
    if result.get("error_message"):
        tracker.finish_item(result["target_tag"], "failed",
                            error_message=result["error_message"])
        tracker.finalize("failed", result["error_message"])
        print(f"ERROR: {result['error_message']}", file=sys.stderr); return 1
    tracker.phase("writing results", "Validating stable prediction table and manifest")
    tracker.finish_item(result["target_tag"],
                        "skipped" if result["status"] == "skipped" else "success",
                        [result["stable_prediction_table"]] if result.get("stable_prediction_table") else [],
                        "Preview complete; classifier was not executed." if args.dry_run else
                        "Cached predictions reused." if result["status"] == "skipped" else
                        "Essentiality classifier complete")
    final = "dry_run_success" if args.dry_run else "cached" if result["status"] == "skipped" else "success"
    tracker.finalize(final,
                     "Preview complete" if args.dry_run else
                     "Cached classifier predictions reused" if final == "cached" else
                     "Essentiality Classifier complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
