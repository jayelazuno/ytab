#!/usr/bin/env python3
"""Combine per-sample R2 Hermes BLAST summaries."""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "nan"))
    except ValueError:
        return float("nan")


def write_delimited(path: Path, rows: list[dict[str, str]], delimiter: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-sample-dir", type=Path, required=True)
    parser.add_argument("--summary-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    paths = sorted(args.per_sample_dir.glob("*.r2_hermes_blast_summary.tsv"))
    if not paths:
        raise SystemExit(f"ERROR: no per-sample summaries found in {args.per_sample_dir}")
    rows: list[dict[str, str]] = []
    for path in paths:
        data = read_tsv(path)
        if data:
            rows.extend(data)
    if not rows:
        raise SystemExit("ERROR: per-sample summary files were empty")
    args.summary_dir.mkdir(parents=True, exist_ok=True)
    write_delimited(args.summary_dir / "r2_hermes_blast_summary.csv", rows, ",")
    write_delimited(args.summary_dir / "r2_hermes_blast_summary.tsv", rows, "\t")
    compact_rows = []
    for row in rows:
        compact_rows.append({
            "Sample": row.get("sample_id", ""),
            "Role": row.get("condition_or_role", ""),
            "Pool": row.get("pool", ""),
            "Background": row.get("background", ""),
            "Total R2 reads": row.get("total_r2_reads", ""),
            "Hermes-positive R2 reads": row.get("hermes_positive_reads", ""),
            "Percent Hermes-positive R2 reads": row.get("percent_r2_reads_hermes_positive", ""),
        })
    write_delimited(args.summary_dir / "r2_hermes_blast_compact.csv", compact_rows, ",")
    percents = [
        numeric(row, "percent_r2_reads_hermes_positive")
        for row in rows
        if row.get("status", "") == "success"
    ]
    percents = [x for x in percents if x == x]
    failed = [row for row in rows if row.get("status", "") != "success"]
    top = sorted(rows, key=lambda r: numeric(r, "percent_r2_reads_hermes_positive"), reverse=True)[:5]
    lines = [
        "R2 Hermes BLAST diagnostic summary",
        "",
        f"R2 FASTQ files processed: {len(rows)}",
        f"Successful samples: {sum(row.get('status') == 'success' for row in rows)}",
        f"Failed samples: {len(failed)}",
    ]
    if percents:
        lines.extend([
            f"Median percent Hermes-positive R2 reads: {statistics.median(percents):.4f}",
            f"Minimum percent Hermes-positive R2 reads: {min(percents):.4f}",
            f"Maximum percent Hermes-positive R2 reads: {max(percents):.4f}",
            "",
            "Samples with highest percent Hermes-positive R2 reads:",
        ])
        for row in top:
            lines.append(
                f"  {row.get('sample_id', '')}: "
                f"{numeric(row, 'percent_r2_reads_hermes_positive'):.4f}% "
                f"({row.get('hermes_positive_reads', '')}/{row.get('total_r2_reads', '')})"
            )
    if failed:
        lines.append("")
        lines.append("Samples with failed status:")
        for row in failed:
            lines.append(f"  {row.get('sample_id', '')}: {row.get('notes', '')}")
    lines.extend([
        "",
        "Interpretation guidance:",
        "If a large percentage of R2 reads contain Hermes CDS sequence, this supports the hypothesis that R2 reads require Hermes/transposon clipping before mapping.",
        "If only a small percentage of R2 reads contain Hermes CDS sequence, the low mapped-read count likely has another cause and should be investigated separately.",
        "This BLAST diagnostic does not by itself change the pipeline.",
        "Any clipping/remapping decision should be made as a separate follow-up step after reviewing the BLAST results.",
    ])
    (args.summary_dir / "r2_hermes_blast_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote combined summaries to {args.summary_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
