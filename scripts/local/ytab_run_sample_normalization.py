#!/usr/bin/env python3
"""Explore local MidLC normalization targets and insertion-site retention."""

from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.normalization_runner import (
    get_normalization_samples,
    load_project_for_normalization,
    normalize_target_tag,
    parse_targets,
    run_normalization_project,
)


def positive(value: str) -> float:
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be numeric") from None
    if not number > 0 or number == float("inf"):
        raise argparse.ArgumentTypeError("must be a finite positive number")
    return number


def retention(value: str) -> float:
    number = positive(value)
    if number > 1:
        raise argparse.ArgumentTypeError("must be at most 1")
    return number


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--targets", default="auto", help="auto or comma-separated positive numeric values")
    parser.add_argument("--sample-mode", choices=("parents", "all", "treated"), default="parents")
    parser.add_argument("--samples", help="Comma-separated included sample names; overrides sample mode")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--keep-going", action="store_true")
    parser.add_argument("--min-site-retention", type=retention, default=0.95)
    parser.add_argument("--auto-min-target", type=positive)
    parser.add_argument("--auto-max-target", type=positive)
    parser.add_argument("--auto-step", type=positive, default=5.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    samples = [part.strip() for part in args.samples.split(",") if part.strip()] if args.samples else None
    try:
        parsed = parse_targets(args.targets)
        config = load_project_for_normalization(args.project_config)
        selected = get_normalization_samples(config, args.sample_mode, samples)
        print(f"Project id: {config.get('project_id')}")
        print(f"Species: {config.get('species')}")
        print(f"Sample mode: {'explicit' if samples is not None else args.sample_mode}")
        print("Selected samples: " + ", ".join(str(row["sample"]) for row in selected))
        print(f"Target mode: {'auto' if parsed == 'auto' else 'manual'}")
        if parsed == "auto":
            print(
                f"Auto settings: minimum={args.auto_min_target or 'adaptive'} "
                f"maximum={args.auto_max_target or 'adaptive'} step={args.auto_step} "
                f"minimum site retention={args.min_site_retention}"
            )
        else:
            print("Targets: " + ", ".join(str(value) for value in parsed))
            print("Target tags: " + ", ".join(normalize_target_tag(value) for value in parsed))
        print("Input hit files:")
        for sample in selected:
            name = str(sample["sample"])
            print(f"  {Path(config['output_project_dir']) / 'create_hit_file' / name / f'{name}_hits.txt'}")
        summary = run_normalization_project(
            args.project_config, targets=args.targets, sample_mode=args.sample_mode,
            samples=samples, threads=args.threads, force=args.force, dry_run=args.dry_run,
            keep_going=args.keep_going, min_site_retention=args.min_site_retention,
            auto_min_target=args.auto_min_target, auto_max_target=args.auto_max_target,
            auto_step=args.auto_step,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"Output normalization directory: {summary['output_dir']}")
    for result in summary["results"]:
        print(f"\n{result['target_tag']}: {result['status']}")
        print("Command: " + shlex.join(result["command_run"]))
        print("Normalized hit files: " + (", ".join(result["normalized_hit_files"]) or "none"))
        print(f"Manifest: {result['manifest_path']}")
        if result.get("error_message"):
            print(f"ERROR: {result['error_message']}", file=sys.stderr)
    if summary["recommendation"]:
        recommendation = summary["recommendation"]
        print(
            f"Recommended target: {recommendation['recommended_target']} "
            f"({recommendation['recommended_target_tag']})"
        )
        print(f"Recommendation: {recommendation['recommendation_path']}")
    print(f"Comparison table: {summary['comparison_path'] or 'not available (dry run or no outputs)'}")
    return 0 if args.dry_run or summary["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
