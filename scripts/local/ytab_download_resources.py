#!/usr/bin/env python3
"""Placeholder manager for future manifest-backed YTAB resource downloads."""

from __future__ import annotations

from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = REPO_ROOT / "resources" / "resource_manifest.yaml"
NO_MANIFEST_MESSAGE = (
    "No resource manifest found. Copy reference files into "
    "resources/species/<species>/reference_genome/ or add resources/resource_manifest.yaml."
)


def main() -> int:
    if not MANIFEST_PATH.is_file():
        print(NO_MANIFEST_MESSAGE)
        return 0
    try:
        with MANIFEST_PATH.open(encoding="utf-8") as handle:
            manifest = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        print(f"ERROR: could not read {MANIFEST_PATH}: {exc}")
        return 1
    if not isinstance(manifest, dict):
        print(f"ERROR: resource manifest must contain a YAML mapping: {MANIFEST_PATH}")
        return 1
    print(f"Resource manifest found: {MANIFEST_PATH}")
    print("Manifest download execution is not implemented yet; no files were downloaded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
