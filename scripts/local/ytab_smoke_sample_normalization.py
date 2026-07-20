#!/usr/bin/env python3
"""Smoke-test auto or manual MidLC normalization on selected parent samples."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.normalization_runner import run_normalization_project


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--targets", default="auto")
    parser.add_argument("--sample-mode", choices=("parents", "all", "treated"), default="parents")
    parser.add_argument("--samples", help="Comma-separated included sample names")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--min-site-retention", type=float, default=0.95)
    parser.add_argument("--auto-min-target", type=float)
    parser.add_argument("--auto-max-target", type=float)
    parser.add_argument("--auto-step", type=float, default=5.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    samples = [part.strip() for part in args.samples.split(",") if part.strip()] if args.samples else None
    try:
        summary = run_normalization_project(
            args.project_config, targets=args.targets, sample_mode=args.sample_mode,
            samples=samples, threads=args.threads, force=args.force,
            min_site_retention=args.min_site_retention, auto_min_target=args.auto_min_target,
            auto_max_target=args.auto_max_target, auto_step=args.auto_step,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print("Selected samples: " + ", ".join(summary["selected_samples"]))
    print("Input hit files:")
    for path in summary["input_hit_files"]:
        print(f"  {path}")
    print(f"Target mode: {summary['target_mode']}")
    print(f"Output directory: {summary['output_dir']}")
    normalized = [path for result in summary["results"] for path in result["normalized_hit_files"]]
    print("Detected normalized hit files: " + (", ".join(normalized) or "none"))
    if summary["recommendation"]:
        print(f"Recommendation file: {summary['recommendation']['recommendation_path']}")
        print(
            f"Recommended target: {summary['recommendation']['recommended_target']} "
            f"({summary['recommendation']['recommended_target_tag']})"
        )
    print("Manifest path: " + ", ".join(result["manifest_path"] for result in summary["results"]))
    if summary["failed"] == 0 and normalized:
        print("PASS: sample normalization")
        return 0
    print("FAIL: normalization failed or normalized hit files were not detected")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
