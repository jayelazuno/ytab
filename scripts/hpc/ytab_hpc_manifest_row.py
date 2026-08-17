#!/usr/bin/env python3
"""Read one row from hpc_sample_manifest.csv for SGE array jobs."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--task-id", type=int)
    group.add_argument("--sample")
    parser.add_argument("--format", choices=("shell", "json"), default="shell")
    parser.add_argument("--count", action="store_true")
    return parser.parse_args()


def safe_key(name: str) -> str:
    return "YTAB_" + re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").upper()


def main() -> int:
    args = parse_args()
    if not args.manifest.is_file():
        raise SystemExit(f"ERROR: manifest is missing: {args.manifest}")
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if args.count:
        print(len(rows))
        return 0
    if args.task_id is not None:
        index = args.task_id - 1
        if index < 0 or index >= len(rows):
            raise SystemExit(f"ERROR: SGE_TASK_ID {args.task_id} is outside manifest row range 1-{len(rows)}")
        row = rows[index]
    else:
        matches = [row for row in rows if row.get("sample") == args.sample]
        if len(matches) != 1:
            raise SystemExit(f"ERROR: expected one manifest row for sample {args.sample}; found {len(matches)}")
        row = matches[0]
    if args.format == "json":
        print(json.dumps(row, indent=2))
    else:
        for key, value in row.items():
            print(f"export {safe_key(key)}={shlex.quote(str(value or ''))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
