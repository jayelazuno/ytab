#!/usr/bin/env python3
"""Smoke-check a YTAB project imported from existing CreateHitFile outputs."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ytab.pipeline.project_status import build_project_status


def rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", type=Path, default=ROOT / "output/projects/H2O2_screen_v1/config/project.yaml")
    parser.add_argument("--expected-hit-files", type=int, default=8)
    args = parser.parse_args()

    config_path = args.project_config.resolve()
    if not config_path.is_file():
        print(f"ERROR: project config missing: {config_path}", file=sys.stderr)
        return 1
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if config.get("start_stage") != "create_hit_file":
        print("ERROR: imported project must set start_stage: create_hit_file", file=sys.stderr)
        return 1
    project_dir = ROOT / config["output_project_dir"]
    sample_sheet = ROOT / config["sample_sheet"]
    samples = rows(sample_sheet)
    if len(samples) != args.expected_hit_files:
        print(f"ERROR: expected {args.expected_hit_files} imported samples, found {len(samples)}", file=sys.stderr)
        return 1
    missing = [row["sample"] for row in samples if not (ROOT / row["hit_file"]).is_file()]
    if missing:
        print("ERROR: imported hit files missing: " + ", ".join(missing), file=sys.stderr)
        return 1
    if any(row.get("fastq_1") or row.get("fastq_2") for row in samples):
        print("ERROR: imported hit project should not require FASTQ columns", file=sys.stderr)
        return 1
    map_status = {row["sample"]: row["status"] for row in rows(project_dir / "manifests/mapfastq/mapfastq_status.csv")}
    hit_status = {row["sample"]: row["status"] for row in rows(project_dir / "manifests/create_hit_file/create_hit_file_status.csv")}
    if set(map_status.values()) != {"imported_or_not_required"}:
        print("ERROR: MapFastq import status must be imported_or_not_required", file=sys.stderr)
        return 1
    if set(hit_status.values()) != {"imported_success"}:
        print("ERROR: CreateHitFile import status must be imported_success", file=sys.stderr)
        return 1
    status = build_project_status(config_path)
    stages = {row["stage"]: row["status"] for row in status["stages"]}
    if stages.get("mapfastq") != "imported_or_not_required" or stages.get("create_hit_file") != "imported_success":
        print(json.dumps(stages, indent=2), file=sys.stderr)
        print("ERROR: project status did not preserve imported stage states", file=sys.stderr)
        return 1
    print(f"Import smoke passed: {len(samples)} hit files imported; MapFastq not required.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
