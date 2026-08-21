#!/usr/bin/env python3
"""Create an R2 FASTQ manifest for the Zn Hermes BLAST diagnostic."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Any


R2_PATTERNS = (
    "*R2*.fastq.gz",
    "*R2*.fq.gz",
    "*_R2_*.fastq.gz",
    "*_R2_*.fq.gz",
    "*R2*.fastq",
    "*R2*.fq",
)


def normalize(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def clean_text(value: Any) -> str:
    return str(value or "").strip()


def infer_role(text: str) -> str:
    lower = text.lower()
    if re.search(r"treated|h2o2|zn|1[._-]?5\s*m?m", lower):
        return "Zn-treated"
    if re.search(r"mock|control|parent", lower):
        return "mock/control"
    return "unknown"


def infer_pool(text: str) -> str:
    match = re.search(r"pool[\s_-]*(\d+)", text, flags=re.IGNORECASE)
    return f"pool{match.group(1)}" if match else ""


def infer_background(text: str) -> str:
    match = re.search(r"\byH\d+\b", text, flags=re.IGNORECASE)
    return match.group(0) if match else ""


def read_sample_sheet(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample)
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        return [{str(k): clean_text(v) for k, v in row.items()} for row in reader]


def row_values(row: dict[str, str]) -> list[str]:
    values = []
    preferred = (
        "sample", "sample_id", "Sample", "Sample_ID", "sample_name", "SampleName",
        "name", "Name", "fastq", "fastq_2", "R2", "read2"
    )
    for key in preferred:
        if key in row and row[key]:
            values.append(row[key])
    values.extend(v for v in row.values() if v)
    return list(dict.fromkeys(values))


def match_sample_sheet(fastq: Path, rows: list[dict[str, str]]) -> tuple[str, dict[str, str] | None, str]:
    basename = fastq.name
    stem = basename
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    stem = re.sub(r"(^|[_\-.])R2([_\-.].*|$)", "", stem, flags=re.IGNORECASE)
    fastq_key = normalize(basename)
    stem_key = normalize(stem)
    best: tuple[int, dict[str, str] | None, str] = (0, None, "")
    for row in rows:
        for value in row_values(row):
            key = normalize(value)
            if not key:
                continue
            score = 0
            if key and key in fastq_key:
                score = len(key)
            elif stem_key and stem_key in key:
                score = len(stem_key)
            if score > best[0]:
                best = (score, row, value)
    if best[1] is not None:
        return best[2], best[1], "matched sampleSheetUsed.csv"
    return stem, None, "no sample-sheet match; sample_id inferred from FASTQ basename"


def find_r2_fastqs(r2_dir: Path) -> list[Path]:
    hits: dict[Path, None] = {}
    for pattern in R2_PATTERNS:
        for path in r2_dir.rglob(pattern):
            name = path.name.lower()
            if "r1" in name and "r2" not in name:
                continue
            if re.search(r"(^|[^a-z0-9])r2([^a-z0-9]|$)|read2", name):
                hits[path.resolve()] = None
            elif "r2" in name:
                hits[path.resolve()] = None
    return sorted(hits)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r2-dir", type=Path, default=Path("/nfsscratch/jayelazuno/workspace/EO46-Zn-tox/He_26110"))
    parser.add_argument("--sample-sheet", type=Path, default=Path("/nfsscratch/jayelazuno/workspace/EO46-Zn-tox/He_26110/sampleSheetUsed.csv"))
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.r2_dir.is_dir():
        raise SystemExit(f"ERROR: R2 FASTQ directory is not visible: {args.r2_dir}")
    rows = read_sample_sheet(args.sample_sheet)
    fastqs = find_r2_fastqs(args.r2_dir)
    if not fastqs:
        raise SystemExit(f"ERROR: no R2 FASTQ files found under {args.r2_dir}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "array_index", "sample_id", "r2_fastq", "fastq_basename",
        "condition_or_role", "pool", "background", "sample_sheet_match", "notes",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for index, fastq in enumerate(fastqs, start=1):
            sample_id, row, note = match_sample_sheet(fastq, rows)
            context = " ".join([sample_id, fastq.name] + (row_values(row) if row else []))
            writer.writerow({
                "array_index": index,
                "sample_id": sample_id,
                "r2_fastq": str(fastq),
                "fastq_basename": fastq.name,
                "condition_or_role": infer_role(context),
                "pool": infer_pool(context),
                "background": infer_background(context),
                "sample_sheet_match": "yes" if row else "no",
                "notes": note,
            })
    print(f"Wrote {len(fastqs)} R2 FASTQ rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
