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
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    samples = [part.strip() for part in args.samples.split(",") if part.strip()] if args.samples else None
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
        summary = run_summary_normalized_project(
            args.project_config, targets=args.targets, sample_mode=args.sample_mode,
            samples=samples, threads=args.threads, force=args.force, dry_run=args.dry_run,
            keep_going=args.keep_going, min_feature_retention=args.min_feature_retention,
        )
    except (FileNotFoundError, OSError, ValueError, json.JSONDecodeError) as exc:
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
    return 0 if args.dry_run or summary["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
