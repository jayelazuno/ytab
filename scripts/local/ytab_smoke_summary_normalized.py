#!/usr/bin/env python3
"""Smoke-test normalized SummaryTable and feature-level target evaluation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.summary_normalized_runner import run_summary_normalized_project


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--targets", default="recommended")
    parser.add_argument("--sample-mode", choices=("parents", "all", "treated"), default="parents")
    parser.add_argument("--samples", help="Comma-separated included sample names")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--min-feature-retention", type=float, default=0.95)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    samples = [part.strip() for part in args.samples.split(",") if part.strip()] if args.samples else None
    try:
        summary = run_summary_normalized_project(
            args.project_config, targets=args.targets, sample_mode=args.sample_mode,
            samples=samples, threads=args.threads, force=args.force,
            min_feature_retention=args.min_feature_retention,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print("Selected targets: " + ", ".join(row["target_tag"] for row in summary["target_infos"]))
    print("Selected samples: " + ", ".join(summary["selected_samples"]))
    print("Input normalized hit files:")
    for result in summary["results"]:
        print(f"  {result['input_normalized_hits_file']}")
    print(f"Output directory: {summary['output_dir']}")
    tables = [path for result in summary["results"] for path in result["feature_tables"]]
    print("Detected normalized feature tables: " + (", ".join(tables) or "none"))
    evaluation = summary["evaluation"]
    print(f"Target evaluation file: {evaluation['evaluation_path'] if evaluation else 'none'}")
    print(f"Feature recommendation file: {evaluation['feature_recommendation_path'] if evaluation else 'none'}")
    if summary["failed"] == 0 and len(tables) == len(summary["results"]) and evaluation:
        print("PASS: normalized SummaryTable")
        return 0
    print("FAIL: normalized SummaryTable failed or outputs were not detected")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
