#!/usr/bin/env python3
"""Run a one-sample local MapFastq smoke test without downstream steps."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.mapfastq_runner import run_mapfastq_project


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
        summary = run_mapfastq_project(
            args.project_config, samples=[args.sample], threads=args.threads,
            force=args.force,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    result = summary["results"][0]
    output_dir = Path(result["output_dir"])
    manifest = output_dir.parents[1] / "manifests" / "mapfastq" / f"{args.sample}.mapfastq_manifest.json"
    print(f"Output directory: {output_dir}")
    print(f"Manifest path: {manifest}")
    if result["status"] in {"success", "skipped"}:
        print(f"PASS: sample mapping {result['status']}")
        return 0
    print(f"FAIL: {result.get('error_message') or 'mapping failed'}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
