#!/usr/bin/env python3
"""Record the working classifier target without marking it as final/reviewed."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ytab.pipeline.classifier_runner import find_combined_feature_table, load_project_for_classifier, resolve_classifier_target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--target", default="recommended")
    args = parser.parse_args()
    try:
        config = load_project_for_classifier(args.project_config)
        target = resolve_classifier_target(config, args.target)
        feature_table = find_combined_feature_table(config, target["target_tag"])
        configured_project = Path(str(config["output_project_dir"]))
        project = configured_project if configured_project.is_absolute() else Path(config["_repo_root"]) / configured_project
        path = project / "config" / "working_classifier_target.json"
        data = {
            "project_id": config.get("project_id"),
            "record_type": "working_classifier_target",
            "final_reviewed_target": False,
            "target": target["target"],
            "target_tag": target["target_tag"],
            "target_resolution_source": target.get("source", ""),
            "combined_feature_table": str(feature_table),
            "recorded_at": datetime.now(timezone.utc).isoformat(),
            "message": "Working target for downstream release analysis; not a saved final classifier target.",
        }
        path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"Working classifier target: {data['target_tag']}")
    print(f"Record: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
