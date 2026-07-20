#!/usr/bin/env python3
"""Run raw LibraryDiagnostics across selected local project hit files."""

from __future__ import annotations

import argparse
import shlex
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
    parser.add_argument("--samples", help="Comma-separated included sample names")
    parser.add_argument("--threads", type=int, help="Reserved runner resource setting")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selected = [value.strip() for value in args.samples.split(",") if value.strip()] if args.samples else None
    try:
        config = load_project_for_library_diagnostics(args.project_config)
        included = get_included_samples(config)
        summary = run_library_diagnostics_project(
            args.project_config, samples=selected, threads=args.threads,
            force=args.force, dry_run=args.dry_run,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"Project id: {summary['project_id']}")
    print(f"Species: {summary['species']}")
    print(f"Included samples: {len(included)}")
    print("Selected samples: " + ", ".join(summary["selected_samples"]))
    print("Input hit files:")
    for path in summary["input_hit_files"]:
        print(f"  {path}")
    print(f"Output diagnostics directory: {summary['output_dir']}")
    print(f"Export directory: {summary['export_dir']}")
    print("Command: " + shlex.join(summary["command_run"]))
    print(f"Status: {summary['status']}")
    print(f"Manifest: {summary['manifest_path']}")
    if summary.get("warnings"):
        print("Warnings: " + "; ".join(summary["warnings"]))
    if summary.get("error_message"):
        print(f"ERROR: {summary['error_message']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
