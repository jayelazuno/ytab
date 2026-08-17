#!/usr/bin/env python3
"""Write one sample-level HPC status row and a JSON manifest."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


FIELDS = {
    "mapfastq": [
        "sample", "status", "layout", "fastq_1", "fastq_2", "output_dir",
        "log_file", "elapsed_seconds", "message",
    ],
    "create_hit_file": [
        "sample", "status", "input_bam", "hits_file", "output_dir",
        "log_file", "elapsed_seconds", "message",
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=sorted(FIELDS), required=True)
    parser.add_argument("--project-dir", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--status", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--message", default="")
    parser.add_argument("--elapsed-seconds", default="0")
    parser.add_argument("--layout", default="")
    parser.add_argument("--fastq-1", default="")
    parser.add_argument("--fastq-2", default="")
    parser.add_argument("--input-bam", default="")
    parser.add_argument("--hits-file", default="")
    parser.add_argument("--manifest-extra-json", default="{}")
    return parser.parse_args()


def update_csv(path: Path, fieldnames: list[str], row: dict[str, str]) -> None:
    existing: dict[str, dict[str, str]] = {}
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            existing = {item.get("sample", ""): item for item in csv.DictReader(handle)}
    existing[row["sample"]] = {key: row.get(key, "") for key in fieldnames}
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for sample in sorted(existing):
            writer.writerow(existing[sample])


def main() -> int:
    args = parse_args()
    fieldnames = FIELDS[args.stage]
    status_path = args.project_dir / "manifests" / args.stage / f"{args.stage}_status.csv"
    row = {
        "sample": args.sample,
        "status": args.status,
        "layout": args.layout,
        "fastq_1": args.fastq_1,
        "fastq_2": args.fastq_2,
        "input_bam": args.input_bam,
        "hits_file": args.hits_file,
        "output_dir": args.output_dir,
        "log_file": args.log_file,
        "elapsed_seconds": args.elapsed_seconds,
        "message": args.message,
    }
    update_csv(status_path, fieldnames, row)
    try:
        extra = json.loads(args.manifest_extra_json)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"ERROR: invalid manifest JSON: {exc}") from exc
    manifest_dir = args.project_dir / "manifests" / args.stage
    manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = manifest_dir / f"{args.sample}.{args.stage}_manifest.json"
    manifest = {
        **extra,
        "sample": args.sample,
        "stage": args.stage,
        "status": args.status,
        "output_dir": args.output_dir,
        "log_file": args.log_file,
        "message": args.message,
        "last_updated": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Status: {status_path}")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
