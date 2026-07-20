#!/usr/bin/env python3
"""Smoke-test essentiality classification for one combined parent target."""

from __future__ import annotations

import argparse
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
    parser.add_argument("--save-final-target", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        result = run_classifier_project(
            args.project_config, target=args.target, seed=args.seed, force=args.force,
            save_final_target=args.save_final_target,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr); return 1
    print(f"Selected target: {result['target']} ({result['target_tag']})")
    print(f"Combined feature table: {result['input_combined_feature_table']}")
    print("Required classifier resources:")
    for name, entries in result["classifier_resources_resolved"].items():
        print(f"  {name}: " + ", ".join(item["path"] for item in entries))
    print(f"Stable prediction table: {result.get('stable_prediction_table') or 'not produced'}")
    print(f"Output rows: {result.get('output_prediction_count') if result.get('output_prediction_count') is not None else 'not available'}")
    if result.get("class_counts"):
        print("Prediction labels: " + ", ".join(f"{key}={value}" for key, value in sorted(result["class_counts"].items())))
    print(f"Manifest path: {result['manifest_path']}")
    stable = result.get("stable_prediction_table")
    if result["status"] in {"success", "skipped"} and stable and Path(stable).is_file():
        print("PASS: essentiality classifier"); return 0
    print(f"FAIL: {result.get('error_message') or 'stable prediction table was not detected'}"); return 1


if __name__ == "__main__":
    raise SystemExit(main())
