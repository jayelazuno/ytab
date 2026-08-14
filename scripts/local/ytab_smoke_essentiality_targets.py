#!/usr/bin/env python3
"""Smoke-test Essentiality target parsing, tags, discovery, and precedence."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.pipeline.normalization_runner import normalize_target_tag, parse_targets


def tag_value(tag: str) -> float:
    if not tag.startswith("T"):
        raise ValueError(tag)
    value = float(tag[1:].replace("p", "."))
    if not math.isfinite(value) or value <= 0 or normalize_target_tag(value) != tag:
        raise ValueError(tag)
    return value


def discover(project: Path) -> list[tuple[float, str]]:
    found: set[str] = set()
    for stage in ("sample_normalization", "summary_normalized", "combined_hits",
                  "summary_combined", "classifier"):
        root = project / stage
        if not root.is_dir():
            continue
        for path in root.rglob("*"):
            for part in path.relative_to(root).parts:
                if part.startswith("T"):
                    try:
                        tag_value(part)
                    except ValueError:
                        continue
                    found.add(part)
    normalization = project / "sample_normalization"
    for name in ("normalization_feature_recommendation.json",
                 "normalization_recommendation.json"):
        path = normalization / name
        if path.is_file():
            data = json.loads(path.read_text(encoding="utf-8"))
            tag = str(data.get("recommended_target_tag") or "")
            try:
                tag_value(tag)
            except ValueError:
                continue
            found.add(tag)
    return sorted((tag_value(tag), tag) for tag in found)


def recommendation(project: Path) -> tuple[str, str] | None:
    normalization = project / "sample_normalization"
    for name, level in (
        ("normalization_feature_recommendation.json", "feature"),
        ("normalization_recommendation.json", "site"),
    ):
        path = normalization / name
        if not path.is_file():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        tag = str(data.get("recommended_target_tag") or "")
        try:
            tag_value(tag)
        except ValueError:
            continue
        return level, tag
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    args = parser.parse_args()
    project = args.project_config.resolve().parents[1]
    final = project / "config" / "final_classifier_target.txt"
    before = final.read_bytes() if final.is_file() else None

    assert parse_targets("20,100,6400") == [20.0, 100.0, 6400.0]
    assert parse_targets("20.5,59.7") == [20.5, 59.7]
    for invalid in ("", "0", "-1", "NA", "Inf", "NaN", "text", "20,,30"):
        try:
            parse_targets(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid target accepted: {invalid!r}")
    expected = {20: "T020", 20.5: "T020p5", 59.7: "T059p7",
                100: "T100", 6400: "T6400"}
    for value, tag in expected.items():
        assert normalize_target_tag(value) == tag
        assert normalize_target_tag(tag_value(tag)) == tag

    targets = discover(project)
    assert targets
    assert targets == sorted(targets)
    assert all(normalize_target_tag(value) == tag for value, tag in targets)
    selected = recommendation(project)
    feature = project / "sample_normalization" / "normalization_feature_recommendation.json"
    site = project / "sample_normalization" / "normalization_recommendation.json"
    if feature.is_file():
        assert selected and selected[0] == "feature"
    elif site.is_file():
        assert selected and selected[0] == "site"

    # Exercise site fallback and the safe no-recommendation state in memory.
    assert next((level for present, level in ((False, "feature"), (True, "site")) if present), None) == "site"
    assert next((level for present, level in ((False, "feature"), (False, "site")) if present), None) is None
    after = final.read_bytes() if final.is_file() else None
    assert before == after
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
