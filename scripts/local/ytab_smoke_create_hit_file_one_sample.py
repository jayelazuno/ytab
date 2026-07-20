#!/usr/bin/env python3
"""Run one CreateHitFile smoke sample without downstream analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.create_hit_file_runner import run_create_hit_file_project


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        summary = run_create_hit_file_project(
            args.project_config, samples=[args.sample], threads=args.threads,
            force=args.force,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    result = summary["results"][0]
    manifest = Path(result["output_dir"]).parents[1] / "manifests" / "create_hit_file" / f"{args.sample}.create_hit_file_manifest.json"
    print(f"Input BAM: {result['input_bam']}")
    print(f"Output directory: {result['output_dir']}")
    print(f"Hits file: {result['hits_file']}")
    print(f"Manifest path: {manifest}")
    if result["status"] in {"success", "skipped"}:
        print(f"PASS: CreateHitFile {result['status']}")
        return 0
    print(f"FAIL: {result.get('error_message') or 'CreateHitFile failed'}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
