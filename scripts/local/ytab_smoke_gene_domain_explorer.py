#!/usr/bin/env python3
"""Smoke-test the YTAB-native Gene & Domain Insertion Explorer."""

from __future__ import annotations

import argparse
import csv
import inspect
import re
import struct
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.domain_explorer import domain_insertion_explorer as explorer  # noqa: E402


def png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) < 24 or not header.startswith(b"\x89PNG\r\n\x1a\n"):
        raise AssertionError(f"not a PNG file: {path}")
    return struct.unpack(">II", header[16:24])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-config",
        type=Path,
        default=REPO_ROOT / "output/projects/H2O2_screen_v1/config/project.yaml",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        lookup = explorer.build_gene_lookup(args.project_config)
        records = lookup["records"]
        if not records:
            raise AssertionError("gene lookup is empty")
        eligible_records = [r for r in records if r.gene_id and r.chromosome and r.start > 0 and r.end > r.start]
        gene = min(eligible_records, key=lambda r: r.end - r.start) if eligible_records else records[0]
        candidates = explorer.query_gene(args.project_config, gene.gene_id)
        if not candidates:
            raise AssertionError(f"query did not return selected gene: {gene.gene_id}")
        tracks = explorer.list_insertion_tracks(args.project_config, "raw")
        if not tracks:
            raise AssertionError("no raw insertion tracks found")
        all_tracks = explorer.resolve_track_preset(args.project_config, "raw", "all")
        parent_tracks = explorer.resolve_track_preset(args.project_config, "raw", "parents")
        treated_tracks = explorer.resolve_track_preset(args.project_config, "raw", "treated")
        if not parent_tracks or not treated_tracks:
            raise AssertionError("parent/treated presets did not resolve")
        first_treated = next((i for i, track in enumerate(all_tracks) if track.get("role") == "treated"), None)
        last_parent = max((i for i, track in enumerate(all_tracks) if track.get("role") == "parent"), default=-1)
        if first_treated is not None and last_parent >= first_treated:
            raise AssertionError("all-track preset did not order parent tracks before treated tracks")
        selected_tracks = [track["sample"] for track in (parent_tracks[:1] + treated_tracks[:1])]
        if len(selected_tracks) < 2:
            raise AssertionError("smoke test requires at least two preset-selected tracks")
        pool1_pair = None
        try:
            pool1_pair = explorer.resolve_track_preset(args.project_config, "raw", "pool1_pair")
            if len(pool1_pair) >= 2 and pool1_pair[0].get("role") != "parent":
                raise AssertionError("pool1_pair did not place parent before treated")
        except ValueError:
            pool1_pair = None
        result = explorer.run_gene_domain_explorer(
            args.project_config,
            gene.gene_id,
            samples=[track["sample"] for track in pool1_pair] if pool1_pair else selected_tracks,
            track_source="raw",
            track_preset="custom" if pool1_pair is None else "pool1_pair",
            flank_bp=1000,
            width_px=1800,
            dpi=150,
            label_mode="full",
            show_site_counts=True,
        )
        if result.get("status") not in {"success", "cached"}:
            raise AssertionError(f"plot generation failed: {result}")
        manifest = result["manifest"]
        for key in ("figure_path", "table_path", "manifest_path"):
            path = Path(manifest[key])
            if not path.is_file():
                raise AssertionError(f"missing {key}: {path}")
            if key == "figure_path" and path.stat().st_size <= 0:
                raise AssertionError(f"empty PNG: {path}")
        width, height = png_dimensions(Path(manifest["figure_path"]))
        if width < 1200 or height < 500:
            raise AssertionError(f"unexpected PNG dimensions: {width}x{height}")
        if int(manifest.get("track_count") or 0) < 2:
            raise AssertionError("multi-track plot was not generated")
        if manifest.get("count_definition") != "insertion sites inside gene_start..gene_end only; flank insertions are displayed but not counted":
            raise AssertionError("manifest count definition is missing or incorrect")
        with Path(manifest["table_path"]).open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            columns = set(reader.fieldnames or [])
            table_rows = list(reader)
        for column in ("in_gene", "region_class", "displayed_region_start", "displayed_region_end", "gene_start", "gene_end"):
            if column not in columns:
                raise AssertionError(f"missing insertion table column: {column}")
        table_gene_counts: dict[str, int] = {}
        for row in table_rows:
            source = row.get("source_file", "").lower()
            if source.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".bai")):
                raise AssertionError(f"unexpected sequence/alignment input required: {source}")
            if str(row.get("in_gene", "")).lower() in {"true", "1", "yes"}:
                table_gene_counts[row["sample"]] = table_gene_counts.get(row["sample"], 0) + 1
        for track in manifest.get("track_summaries", []):
            expected = int(track.get("insertions_in_gene") or 0)
            observed = table_gene_counts.get(track["sample"], 0)
            if expected != observed:
                raise AssertionError(f"gene-body count mismatch for {track['sample']}: manifest={expected}, table={observed}")
        source_file = Path(inspect.getsourcefile(explorer) or "")
        if "codex" in source_file.parts:
            raise AssertionError(f"production module imported from codex path: {source_file}")
        source_text = source_file.read_text(encoding="utf-8")
        forbidden = (
            r"^\s*import\s+DomainFigures\b",
            r"^\s*from\s+DomainFigures\b",
            r"^\s*import\s+Organisms\b",
            r"^\s*from\s+Organisms\b",
            r"^\s*import\s+SortedCollection\b",
            r"^\s*from\s+SortedCollection\b",
            r"^\s*import\s+RangeSet\b",
            r"^\s*from\s+RangeSet\b",
        )
        hits = [pattern for pattern in forbidden if re.search(pattern, source_text, re.MULTILINE)]
        if hits:
            raise AssertionError(f"legacy imports detected: {hits}")
    except Exception as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
