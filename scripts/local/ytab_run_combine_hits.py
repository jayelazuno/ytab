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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--target", default="recommended")
    parser.add_argument("--samples", help="Comma-separated included parent sample names")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    samples = [part.strip() for part in args.samples.split(",") if part.strip()] if args.samples else None
    try:
        config = load_project_for_combine_hits(args.project_config)
        selected = get_parent_samples(config, samples)
        target = resolve_combine_target(config, args.target)
        inputs = find_normalized_parent_hits(config, target["target_tag"], selected)
        summary = run_combine_hits_project(
            args.project_config, target=args.target, samples=samples,
            force=args.force, dry_run=args.dry_run,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
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
        print(f"ERROR: {summary['error_message']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
