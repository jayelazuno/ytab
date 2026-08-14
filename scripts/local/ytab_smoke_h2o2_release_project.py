#!/usr/bin/env python3
"""Smoke-check the H2O2 release project after downstream-from-hit execution."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ytab.pipeline.project_status import build_project_status


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", type=Path, default=ROOT / "output/projects/H2O2_screen_v1/config/project.yaml")
    args = parser.parse_args()
    config_path = args.project_config.resolve()
    if not config_path.is_file():
        print(f"ERROR: project config missing: {config_path}", file=sys.stderr)
        return 1
    project = config_path.parents[1]
    required = [
        project / "config/sample_sheet.csv",
        project / "config/reference_resolved.json",
        project / "config/comparison_design.csv",
        project / "manifests/create_hit_file/import_existing_hit_project_manifest.json",
        project / "summary/summary_stats.all_samples.csv",
        project / "manifests/project_status.json",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        print("ERROR: required release outputs missing:\n" + "\n".join(missing), file=sys.stderr)
        return 1
    if (project / "config/final_classifier_target.txt").exists():
        print("ERROR: final classifier target was saved without explicit review.", file=sys.stderr)
        return 1
    status = build_project_status(config_path)
    stages = {row["stage"]: row["status"] for row in status["stages"]}
    if stages.get("mapfastq") != "imported_or_not_required" or stages.get("create_hit_file") != "imported_success":
        print(json.dumps(stages, indent=2), file=sys.stderr)
        print("ERROR: H2O2 release project no longer records imported entry stages.", file=sys.stderr)
        return 1
    downstream = ["summary", "library_diagnostics", "sample_normalization", "summary_normalized", "combined_hits", "summary_combined", "classifier", "comparison_design", "treated_vs_parent"]
    incomplete = [stage for stage in downstream if stages.get(stage) not in {"success", "cached", "skipped"}]
    if incomplete:
        print("ERROR: downstream stages incomplete: " + ", ".join(f"{s}={stages.get(s)}" for s in incomplete), file=sys.stderr)
        return 1
    print("H2O2 release smoke passed: imported hit files and downstream release stages are complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
