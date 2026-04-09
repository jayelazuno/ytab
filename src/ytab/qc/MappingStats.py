
#!/usr/bin/env python3
"""
MappingStats.py  - BAM-based mapping QC summary for YTAB

WHAT IT DOES
------------
Given a mapping run directory containing one subdirectory per sample, where each
sample directory contains a coordinate-sorted BAM file (`*.sorted.bam`), compute
a compact mapping/QC summary table directly from the BAM using `samtools view`.

For each sample, the script extracts:

1) Total alignment records
   - `total_records`
   - Counts all alignment records in the BAM

2) Primary alignment records
   - `primary_records`
   - Excludes secondary and supplementary alignments

3) Primary mapped and unmapped records
   - `primary_mapped`
   - `primary_unmapped`

4) Percent mapped
   - `percent_mapped`
   - Computed as:
       primary_mapped / primary_records * 100

5) Duplicate-marked primary records
   - `primary_duplicates`
   - `percent_duplicates`
   - Based on the BAM duplicate flag
   - If duplicates were not marked upstream, these values will remain 0

6) MAPQ-based mapping quality summary
   - `mapq_ge20`
   - `percent_mapq_ge20`
   - `avg_mapq_mapped_primary`
   - By default, MAPQ >= 20 is used as the high-confidence threshold
   - Threshold can be changed with `--mapq-threshold`

7) BAM provenance
   - `bam_path`
   - Absolute path to the BAM used for the summary row

WHY THIS EXISTS
---------------
This script provides a stable, app-friendly mapping QC table that can be used
for downstream visualization in YTAB without depending on aligner log files.
It is intended as the mapping-QC summary layer for the app.

CURRENT DEVELOPMENT USE
-----------------------
At this stage, we are starting with a single smoke-test sample only.

That means the script can be run on a run directory containing just one sample
subdirectory, and it will produce a one-row CSV. This is enough to build the
initial app tab and plotting logic. Later, the exact same script can be rerun
on the full mapping run with all samples.

INPUT EXPECTATION
-----------------
The script expects:

  <run_dir>/<sample_name>/<sample_name>.sorted.bam

More generally, each sample directory must contain exactly one:

  *.sorted.bam

Example smoke-test layout:

  results/smoketests/mapfastq/
    yH298-parent-pool1/
      yH298-parent-pool1.sorted.bam
      yH298-parent-pool1.sorted.bam.bai
      yH298-parent-pool1_log.txt

The log file is ignored by this script.

OUTPUT
------
Writes one combined CSV file with one row per sample:

  <output_csv>

Columns written:

  sample
  total_records
  primary_records
  primary_mapped
  primary_unmapped
  percent_mapped
  primary_duplicates
  percent_duplicates
  mapq_ge20
  percent_mapq_ge20
  avg_mapq_mapped_primary
  bam_path

EXAMPLE USAGE
-------------
Single smoke-test sample:

  python src/ytab/qc/MappingStats.py \
    --run-dir /Users/jayelazuno/workspace/ytab/results/smoketests/mapfastq \
    --output-csv /Users/jayelazuno/workspace/ytab/results/smoketests/mapfastq/mapping_stats.csv

Larger run with multiple sample subdirectories:

  python src/ytab/qc/MappingStats.py \
    --run-dir /Users/jayelazuno/workspace/ytab/results/runs/run_2026-03-12_mapfastq \
    --output-csv /Users/jayelazuno/workspace/ytab/results/runs/run_2026-03-12_mapfastq/mapping_stats.csv \
    --threads 8

NOTES
-----
- This script is BAM-only and does not parse Bowtie2 log files.
- `total_records` includes all alignment records in the BAM.
- `primary_records` excludes secondary and supplementary alignments.
- Duplicate metrics are only meaningful if duplicate flags were set upstream.
- The script is intended to live under `src/ytab/qc/` and serve as a reusable
  QC summary step for both local runs and future pipeline execution.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path
from typing import Dict, List


SECONDARY = 0x100
SUPPLEMENTARY = 0x800
UNMAPPED = 0x4
DUPLICATE = 0x400


def run_samtools_view(bam_path: Path, threads: int = 1):
    cmd = ["samtools", "view", "-@", str(threads), str(bam_path)]
    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )


def compute_bam_stats(
    bam_path: Path,
    mapq_threshold: int = 20,
    threads: int = 1,
) -> Dict[str, object]:
    total_records = 0
    primary_records = 0
    primary_mapped = 0
    primary_unmapped = 0
    primary_duplicates = 0
    mapq_ge20 = 0
    sum_mapq_mapped_primary = 0

    proc = run_samtools_view(bam_path, threads=threads)

    assert proc.stdout is not None
    for line in proc.stdout:
        if not line.strip():
            continue

        total_records += 1

        fields = line.rstrip("\n").split("\t")
        flag = int(fields[1])
        mapq = int(fields[4])

        is_primary = not (flag & SECONDARY) and not (flag & SUPPLEMENTARY)
        if not is_primary:
            continue

        primary_records += 1

        if flag & DUPLICATE:
            primary_duplicates += 1

        if flag & UNMAPPED:
            primary_unmapped += 1
        else:
            primary_mapped += 1
            sum_mapq_mapped_primary += mapq
            if mapq >= mapq_threshold:
                mapq_ge20 += 1

    stderr = proc.stderr.read() if proc.stderr is not None else ""
    returncode = proc.wait()
    if returncode != 0:
        raise RuntimeError(
            f"samtools view failed for {bam_path}\n"
            f"Command exit code: {returncode}\n"
            f"stderr:\n{stderr}"
        )

    percent_mapped = round((primary_mapped / primary_records) * 100, 3) if primary_records else 0.0
    percent_duplicates = round((primary_duplicates / primary_records) * 100, 3) if primary_records else 0.0
    percent_mapq_ge20 = round((mapq_ge20 / primary_mapped) * 100, 3) if primary_mapped else 0.0
    avg_mapq_mapped_primary = round((sum_mapq_mapped_primary / primary_mapped), 2) if primary_mapped else 0.0

    return {
        "total_records": total_records,
        "primary_records": primary_records,
        "primary_mapped": primary_mapped,
        "primary_unmapped": primary_unmapped,
        "percent_mapped": percent_mapped,
        "primary_duplicates": primary_duplicates,
        "percent_duplicates": percent_duplicates,
        "mapq_ge20": mapq_ge20,
        "percent_mapq_ge20": percent_mapq_ge20,
        "avg_mapq_mapped_primary": avg_mapq_mapped_primary,
    }


def find_bam(sample_dir: Path) -> Path | None:
    bams = sorted(sample_dir.glob("*.sorted.bam"))
    if not bams:
        return None
    if len(bams) > 1:
        raise RuntimeError(f"Expected one BAM in {sample_dir}, found {len(bams)}")
    return bams[0]


def collect_rows(run_dir: Path, mapq_threshold: int, threads: int) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []

    for sample_dir in sorted(run_dir.iterdir()):
        if not sample_dir.is_dir():
            continue

        bam_path = find_bam(sample_dir)
        if bam_path is None:
            continue

        sample = sample_dir.name
        stats = compute_bam_stats(
            bam_path=bam_path,
            mapq_threshold=mapq_threshold,
            threads=threads,
        )

        row = {
            "sample": sample,
            **stats,
            "bam_path": str(bam_path.resolve()),
        }
        rows.append(row)

    if not rows:
        raise RuntimeError(f"No sample BAMs found under: {run_dir}")

    return rows


def write_csv(rows: List[Dict[str, object]], output_csv: Path) -> None:
    fieldnames = [
        "sample",
        "total_records",
        "primary_records",
        "primary_mapped",
        "primary_unmapped",
        "percent_mapped",
        "primary_duplicates",
        "percent_duplicates",
        "mapq_ge20",
        "percent_mapq_ge20",
        "avg_mapq_mapped_primary",
        "bam_path",
    ]

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract mapping QC statistics from BAM files using samtools."
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Directory containing one subdirectory per sample, each with one *.sorted.bam file.",
    )
    parser.add_argument(
        "--output-csv",
        required=True,
        help="Output CSV path.",
    )
    parser.add_argument(
        "--mapq-threshold",
        type=int,
        default=20,
        help="MAPQ threshold for the mapq_ge20-style metric (default: 20).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads to pass to samtools view (default: 1).",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    output_csv = Path(args.output_csv).resolve()

    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    rows = collect_rows(
        run_dir=run_dir,
        mapq_threshold=args.mapq_threshold,
        threads=args.threads,
    )
    write_csv(rows, output_csv)
    print(f"Wrote {len(rows)} rows to {output_csv}")


if __name__ == "__main__":
    main()
