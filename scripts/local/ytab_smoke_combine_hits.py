#!/usr/bin/env python3
"""Smoke-test combining normalized parent hit files for one target."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.combine_hits_runner import run_combine_hits_project


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--target", default="recommended")
    parser.add_argument("--samples", help="Comma-separated included parent sample names")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    samples = [part.strip() for part in args.samples.split(",") if part.strip()] if args.samples else None
    try:
        summary = run_combine_hits_project(
            args.project_config, target=args.target, samples=samples, force=args.force,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print(f"Selected target: {summary['target']} ({summary['target_tag']})")
    print("Selected parent samples: " + ", ".join(summary["selected_samples"]))
    print("Input normalized hit files:")
    for path in summary["input_normalized_hit_files"]:
        print(f"  {path}")
    print(f"Output combined hit file: {summary['combined_hits_file']}")
    print(f"Manifest path: {summary['manifest_path']}")
    if summary["status"] in {"success", "skipped"} and Path(summary["combined_hits_file"]).is_file():
        print(f"Total combined sites: {summary.get('total_combined_sites', 'unknown')}")
        print("PASS: combined normalized parent hits")
        return 0
    print(f"FAIL: {summary.get('error_message') or 'combined hit file was not detected'}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
