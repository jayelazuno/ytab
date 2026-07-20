#!/usr/bin/env python3
"""Smoke-test raw LibraryDiagnostics on two selected project samples."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.library_diagnostics_runner import (
    get_included_samples,
    load_project_for_library_diagnostics,
    run_library_diagnostics_project,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--samples", help="Comma-separated included sample names; defaults to first two")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        config = load_project_for_library_diagnostics(args.project_config)
        included = get_included_samples(config)
        selected = (
            [value.strip() for value in args.samples.split(",") if value.strip()]
            if args.samples else [str(row["sample"]) for row in included[:2]]
        )
        if len(selected) < 2:
            raise ValueError("Smoke test requires at least two included samples.")
        summary = run_library_diagnostics_project(
            args.project_config, samples=selected, threads=args.threads, force=args.force,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print("Input hit files:")
    for path in summary["input_hit_files"]:
        print(f"  {path}")
    print(f"Output directory: {summary['output_dir']}")
    summaries = [path for path in summary["detected_outputs"] if path.endswith("summary.csv")]
    print("Detected summary outputs: " + (", ".join(summaries) or "none"))
    print(f"Manifest path: {summary['manifest_path']}")
    if summary["status"] in {"success", "skipped"} and summaries:
        print(f"PASS: LibraryDiagnostics {summary['status']}")
        return 0
    print(f"FAIL: {summary.get('error_message') or 'diagnostic summaries were not detected'}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
