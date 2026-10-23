#!/usr/bin/env python3
"""Read-only project inventory worker; also usable from the command line."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from ytab.screen_coverage.project_summary import project_summary

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        result = project_summary(args.project_config)
    except Exception as exc:
        result = {"error": str(exc), "project_config": str(Path(args.project_config).resolve())}
    Path(args.output).write_text(json.dumps(result, allow_nan=False))
    if "error" in result:
        sys.exit(1)
