#!/usr/bin/env python3
"""Run MapFastq locally, one included project sample at a time."""

from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.mapfastq_runner import (
    get_included_samples,
    load_project_for_mapping,
    run_mapfastq_project,
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
    parser.add_argument("--junction-mate", choices=("R1","R2"), help="Override the designated Hermes/junction read for selected samples (default R1)")
    parser.add_argument("--require-mate-alignment", action=argparse.BooleanOptionalAction, default=None, help="Require concordant pairs; default uses only the designated junction mate")
    parser.add_argument("--hermes-cleanup", action=argparse.BooleanOptionalAction, default=None, help="Automatically trim observed Hermes/adaptor sequence (default enabled)")
    parser.add_argument("--transposon", help="Transposon profile in resources/transposons (default: the project's, else hermes)")
    parser.add_argument("--library-design", help="Library design within the transposon profile (default: the profile's default)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selected = [value.strip() for value in args.samples.split(",") if value.strip()] if args.samples else None
    try:
        config = load_project_for_mapping(args.project_config)
        included = get_included_samples(config)
        print(f"Project id: {config.get('project_id')}")
        print(f"Species: {config.get('species')}")
        print(f"Bowtie2 index prefix: {(config.get('reference') or {}).get('bowtie2_index_prefix')}")
        print(f"Included samples: {len(included)}")
        print("Selected samples: " + ", ".join(selected or [str(row['sample']) for row in included]))
        print(f"Threads: {args.threads if args.threads is not None else config.get('threads')}")
        print(f"Dry run: {'yes' if args.dry_run else 'no'}")
        summary = run_mapfastq_project(
            args.project_config, samples=selected, threads=args.threads,
            force=args.force, dry_run=args.dry_run, keep_going=args.keep_going,
            job_id=args.job_id, progress_file=args.progress_file,
            junction_mate=args.junction_mate, require_mate_alignment=args.require_mate_alignment,
            hermes_cleanup=args.hermes_cleanup,
            transposon=args.transposon, library_design=args.library_design,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    for result in summary["results"]:
        print(f"\n{result['sample']}: {result['status']}")
        print("Command: " + shlex.join(result["command_run"]))
        print(f"Output directory: {result['output_dir']}")
        print(f"Manifest: {Path(result['output_dir']).parents[1] / 'manifests' / 'mapfastq' / (result['sample'] + '.mapfastq_manifest.json')}")
        if result.get("error_message"):
            print(f"ERROR: {result['error_message']}")
    print(
        f"\nSummary: success={summary['success']} failed={summary['failed']} "
        f"skipped={summary['skipped']}"
    )
    return 0 if args.dry_run or summary["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
