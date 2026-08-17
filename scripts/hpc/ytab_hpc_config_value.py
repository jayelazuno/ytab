#!/usr/bin/env python3
"""Print a scalar value from a YTAB project YAML config."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--key", required=True, help="Dot-separated key, e.g. reference.fasta")
    parser.add_argument("--resolve-repo", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.project_config.is_file():
        raise SystemExit(f"ERROR: project config is missing: {args.project_config}")
    data = yaml.safe_load(args.project_config.read_text(encoding="utf-8"))
    value = data
    for part in args.key.split("."):
        if not isinstance(value, dict) or part not in value:
            raise SystemExit(f"ERROR: key not found in project config: {args.key}")
        value = value[part]
    if isinstance(value, (dict, list)):
        print(json.dumps(value))
        return 0
    text = "" if value is None else str(value)
    if args.resolve_repo and text and not Path(text).is_absolute():
        raw_root = Path(str(data.get("repo_root") or args.project_config.resolve().parents[3]))
        repo_root = raw_root if raw_root.is_absolute() else (args.project_config.resolve().parent / raw_root).resolve()
        text = str((repo_root / text).resolve())
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
