#!/usr/bin/env python3
"""Run essentiality classification on a combined parent feature table."""

from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.classifier_runner import run_classifier_project


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--target", default="recommended")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--save-final-target", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        result = run_classifier_project(
            args.project_config, target=args.target, seed=args.seed, force=args.force,
            dry_run=args.dry_run, save_final_target=args.save_final_target,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
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
        print(f"ERROR: {result['error_message']}", file=sys.stderr); return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
