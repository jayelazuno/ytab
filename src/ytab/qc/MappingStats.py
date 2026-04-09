
#!/usr/bin/env python3
"""
MappingStats.py  - BAM-based mapping QC summary for YTAB

WHAT IT DOES
------------
Provides reusable functions to compute mapping statistics from a full
alignment file (SAM or BAM) using samtools view.

This is intended to be called during the MapFastq stage, immediately after
Bowtie2 writes the full SAM and before MAPQ/flag filtering is applied.

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


NOTES
-----
- `primary_duplicates` depends on duplicate flags already being present.
  For current YTAB/Hermes mapping outputs this will usually remain 0.
- `mapq_ge_threshold` uses the threshold supplied by MapFastq (default 20).
"""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Dict, List


SECONDARY = 0x100
SUPPLEMENTARY = 0x800
UNMAPPED = 0x4
DUPLICATE = 0x400


def _run_samtools_view(alignment_path: Path, threads: int = 1) -> subprocess.Popen:
    cmd = ["samtools", "view", "-@", str(threads), str(alignment_path)]
    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )


def compute_alignment_stats(
    alignment_path: str | Path,
    sample: str,
    mapq_threshold: int = 20,
    threads: int = 1,
) -> Dict[str, object]:
    alignment_path = Path(alignment_path).resolve()

    total_records = 0
    primary_records = 0
    primary_mapped = 0
    primary_unmapped = 0
    primary_duplicates = 0
    mapq_ge_threshold = 0
    sum_mapq_mapped_primary = 0

    proc = _run_samtools_view(alignment_path, threads=threads)

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
                mapq_ge_threshold += 1

    stderr = proc.stderr.read() if proc.stderr is not None else ""
    returncode = proc.wait()
    if returncode != 0:
        raise RuntimeError(
            f"samtools view failed for {alignment_path}\n"
            f"exit code: {returncode}\n"
            f"stderr:\n{stderr}"
        )

    percent_mapped = round((primary_mapped / primary_records) * 100, 3) if primary_records else 0.0
    percent_duplicates = round((primary_duplicates / primary_records) * 100, 3) if primary_records else 0.0
    percent_mapq_ge_threshold = round((mapq_ge_threshold / primary_mapped) * 100, 3) if primary_mapped else 0.0
    avg_mapq_mapped_primary = round((sum_mapq_mapped_primary / primary_mapped), 2) if primary_mapped else 0.0

    return {
        "sample": sample,
        "total_records": total_records,
        "primary_records": primary_records,
        "primary_mapped": primary_mapped,
        "primary_unmapped": primary_unmapped,
        "percent_mapped": percent_mapped,
        "primary_duplicates": primary_duplicates,
        "percent_duplicates": percent_duplicates,
        "mapq_ge_threshold": mapq_ge_threshold,
        "percent_mapq_ge_threshold": percent_mapq_ge_threshold,
        "avg_mapq_mapped_primary": avg_mapq_mapped_primary,
        "alignment_path": str(alignment_path),
    }


def write_mapping_stats_csv(row: Dict[str, object], output_csv: str | Path) -> None:
    output_csv = Path(output_csv).resolve()
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames: List[str] = [
        "sample",
        "total_records",
        "primary_records",
        "primary_mapped",
        "primary_unmapped",
        "percent_mapped",
        "primary_duplicates",
        "percent_duplicates",
        "mapq_ge_threshold",
        "percent_mapq_ge_threshold",
        "avg_mapq_mapped_primary",
        "alignment_path",
    ]

    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)