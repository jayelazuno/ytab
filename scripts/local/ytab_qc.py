#!/usr/bin/env python3
"""Library QC table for a project: every per-library QC metric in one flagged, commented table plus plot data.

  python scripts/local/ytab.py qc --project-config output/projects/P/config/project.yaml --output output/projects/P/qc/RUN

Writes qc_table.tsv (one row per library), qc_flags.tsv (one row per library x metric: value, flag, reason, source,
comment, normalization hint), qc_table.xlsx (coloured, with the reasons as cell notes), metrics.tsv and plots/
(count histogram, QQ against a geometric fit and ranked counts). See ytab.qc.library_qc for the rules.
"""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"src"))
from ytab.qc.library_qc import build  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    try:
        manifest = build(a.project_config, a.output)
    except (ValueError, KeyError, OSError) as error:
        p.exit(1, f"ERROR: {error}\n")
    print(json.dumps({k:manifest[k] for k in ("status", "project", "libraries", "libraries_without_hits", "flag_counts")}, indent=1))


if __name__ == "__main__":main()
