#!/usr/bin/env python3
"""Run CreateHitFile locally, one included project sample at a time."""

from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.create_hit_file_runner import (
    get_included_samples,
    load_project_for_create_hit_file,
    run_create_hit_file_project,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--samples", help="Comma-separated included sample names")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--keep-going", action="store_true")
    parser.add_argument("--job-id")
    parser.add_argument("--progress-file", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selected = [value.strip() for value in args.samples.split(",") if value.strip()] if args.samples else None
    try:
        config = load_project_for_create_hit_file(args.project_config)
        included = get_included_samples(config)
        root = Path(config["_repo_root"])
        project_dir = Path(config["output_project_dir"])
        if not project_dir.is_absolute():
            project_dir = root / project_dir
        print(f"Project id: {config.get('project_id')}")
        print(f"Species: {config.get('species')}")
        print(f"Included samples: {len(included)}")
        print("Selected samples: " + ", ".join(selected or [str(row['sample']) for row in included]))
        print(f"Input BAM directory: {project_dir / 'mapfastq'}")
        print(f"Output create_hit_file directory: {project_dir / 'create_hit_file'}")
        print(f"Threads: {args.threads if args.threads is not None else config.get('threads')}")
        print(f"Dry run: {'yes' if args.dry_run else 'no'}")
        summary = run_create_hit_file_project(
            args.project_config, samples=selected, threads=args.threads,
            force=args.force, dry_run=args.dry_run, keep_going=args.keep_going,
            job_id=args.job_id, progress_file=args.progress_file,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    for result in summary["results"]:
        print(f"\n{result['sample']}: {result['status']}")
        print("Command: " + shlex.join(result["command_run"]))
        print(f"Input BAM: {result['input_bam']}")
        print(f"Hits file: {result['hits_file']}")
        if result.get("error_message"):
            print(f"ERROR: {result['error_message']}")
    print(
        f"\nSummary: success={summary['success']} failed={summary['failed']} "
        f"skipped={summary['skipped']}"
    )
    return 0 if args.dry_run or summary["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
