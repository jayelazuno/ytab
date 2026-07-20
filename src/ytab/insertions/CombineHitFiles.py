
#!/usr/bin/env python3

"""
CombineHitFiles.py

Combine multiple YTAB/Tn-seq hit files into one SummaryTable-compatible hit file.

Typical use:
    - After MidLC normalization, each parent pool has a normalized hit file.
    - To follow Gale-style classifier logic, combine normalized parent hit files
      first, then run SummaryTable once on the combined hit file, then run
      Classifier on the combined feature table.

The script combines by:
    chromosome column + position column

and sums:
    reads/count column

Important:
    This version preserves all original columns in the combined output.

For duplicate insertion sites:
    - The read/count column is summed.
    - The chromosome and position columns are kept.
    - All other columns are retained.
    - If a retained annotation column has one unique value across contributing
      rows, that value is kept.
    - If it has multiple unique values, the unique values are collapsed using ";".

Output:
    <sample_name>_hits.txt
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


CHROM_ALIASES = {
    "chromosome",
    "chrom",
    "chr",
    "seqname",
    "seqid",
    "contig",
}

POSITION_ALIASES = {
    "position",
    "pos",
    "site",
    "coordinate",
    "coord",
    "insertion_position",
    "insertion_site",

    # LibraryDiagnostics / normalized hit files
    "hit_position",
    "hit_pos",
    "hit_site",
}

READS_ALIASES = {
    "reads",
    "read",
    "count",
    "counts",
    "n_reads",
    "read_count",

    # LibraryDiagnostics / normalized hit files
    "hit_count",
    "hit_counts",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine multiple hit files by chromosome/position, summing reads, "
            "while preserving all original annotation columns."
        )
    )

    parser.add_argument(
        "--hits-txt",
        nargs="*",
        default=[],
        help="Input hit files to combine.",
    )

    parser.add_argument(
        "--hits-list",
        default=None,
        help=(
            "Optional text file containing one hit-file path per line. "
            "Blank lines and # comments are ignored."
        ),
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for the combined hit file.",
    )

    parser.add_argument(
        "--sample-name",
        required=True,
        help=(
            "Name for the combined pseudo-sample. Output will be "
            "<sample-name>_hits.txt."
        ),
    )

    parser.add_argument(
        "--output-file",
        default=None,
        help=(
            "Optional output filename or path for the combined hit table. "
            "Default remains <sample-name>_hits.txt under --outdir."
        ),
    )

    parser.add_argument(
        "--chrom-col",
        default=None,
        help="Optional explicit chromosome column name.",
    )

    parser.add_argument(
        "--position-col",
        default=None,
        help="Optional explicit position column name.",
    )

    parser.add_argument(
        "--reads-col",
        default=None,
        help="Optional explicit reads/count column name.",
    )

    parser.add_argument(
        "--keep-zero",
        action="store_true",
        help="Keep rows with zero reads. Default is to drop reads <= 0.",
    )

    parser.add_argument(
        "--sort",
        action="store_true",
        help="Sort output by chromosome and position.",
    )

    parser.add_argument(
        "--collapse-sep",
        default=";",
        help="Separator used when multiple annotation values are found for one site.",
    )

    return parser.parse_args()


def read_hits_list(path: str | Path) -> list[str]:
    paths: list[str] = []

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            paths.append(line)

    return paths


def detect_separator(path: Path) -> str:
    first_data_line = ""

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith("#"):
                continue

            first_data_line = stripped
            break

    if "," in first_data_line:
        return ","

    if "\t" in first_data_line:
        return "\t"

    return r"\s+"


def normalize_colname(name: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def looks_like_header(tokens: list[str]) -> bool:
    norm = {normalize_colname(x) for x in tokens}

    return bool(
        norm & CHROM_ALIASES
        or norm & POSITION_ALIASES
        or norm & READS_ALIASES
    )


def read_table_auto(path: Path) -> pd.DataFrame:
    sep = detect_separator(path)

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        first_data_line = ""
        for line in handle:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                first_data_line = stripped
                break

    if sep == ",":
        tokens = first_data_line.split(",")
    elif sep == "\t":
        tokens = first_data_line.split("\t")
    else:
        tokens = re.split(r"\s+", first_data_line)

    header = 0 if looks_like_header(tokens) else None

    df = pd.read_csv(
        path,
        sep=sep,
        header=header,
        comment="#",
        engine="python",
    )

    if header is None:
        if df.shape[1] < 3:
            raise ValueError(
                f"{path} has no header and fewer than 3 columns. "
                "Expected Chromosome, Position, Reads."
            )

        df = df.iloc[:, :3].copy()
        df.columns = ["Chromosome", "Position", "Reads"]

    return df


def choose_column(
    df: pd.DataFrame,
    explicit: str | None,
    aliases: set[str],
    label: str,
) -> str:
    if explicit is not None:
        if explicit not in df.columns:
            raise ValueError(
                f"Explicit {label} column '{explicit}' was not found. "
                f"Available columns: {list(df.columns)}"
            )

        return explicit

    norm_to_original: dict[str, str] = {
        normalize_colname(col): col for col in df.columns
    }

    for alias in aliases:
        if alias in norm_to_original:
            return norm_to_original[alias]

    raise ValueError(
        f"Could not identify {label} column. "
        f"Available columns: {list(df.columns)}"
    )


def unique_nonempty_values(series: pd.Series) -> list[str]:
    values: list[str] = []

    for value in series:
        if pd.isna(value):
            continue

        text = str(value)

        if text == "":
            continue

        if text not in values:
            values.append(text)

    return values


def collapse_values(series: pd.Series, sep: str = ";") -> object:
    values = unique_nonempty_values(series)

    if len(values) == 0:
        return ""

    if len(values) == 1:
        return values[0]

    return sep.join(values)


def read_hit_file_full(
    path: Path,
    chrom_col: str | None = None,
    position_col: str | None = None,
    reads_col: str | None = None,
) -> tuple[pd.DataFrame, str, str, str]:
    df = read_table_auto(path)

    c_col = choose_column(df, chrom_col, CHROM_ALIASES, "chromosome")
    p_col = choose_column(df, position_col, POSITION_ALIASES, "position")
    r_col = choose_column(df, reads_col, READS_ALIASES, "reads")

    # Preserve the full original dataframe, but add internal normalized columns
    # for grouping. These internal columns are dropped before output.
    df = df.copy()

    df["_ytab_chrom"] = df[c_col].astype(str)
    df["_ytab_position"] = pd.to_numeric(df[p_col], errors="coerce")
    df["_ytab_reads"] = pd.to_numeric(df[r_col], errors="coerce")

    before = len(df)

    df = df.dropna(
        subset=["_ytab_chrom", "_ytab_position", "_ytab_reads"]
    ).copy()

    df["_ytab_position"] = df["_ytab_position"].astype(int)
    df["_ytab_reads"] = df["_ytab_reads"].round().astype(int)

    dropped = before - len(df)

    if dropped:
        print(
            f"[WARN] Dropped {dropped} malformed rows from {path}",
            file=sys.stderr,
        )

    return df, c_col, p_col, r_col


def combine_hits(
    hit_files: Iterable[Path],
    chrom_col: str | None = None,
    position_col: str | None = None,
    reads_col: str | None = None,
    keep_zero: bool = False,
    sort_output: bool = False,
    collapse_sep: str = ";",
) -> tuple[pd.DataFrame, list[dict[str, object]], dict[str, object]]:
    frames: list[pd.DataFrame] = []
    summaries: list[dict[str, object]] = []

    output_columns: list[str] | None = None
    detected_chrom_col: str | None = None
    detected_position_col: str | None = None
    detected_reads_col: str | None = None

    for path in hit_files:
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty hit file: {path}")

        df, c_col, p_col, r_col = read_hit_file_full(
            path,
            chrom_col=chrom_col,
            position_col=position_col,
            reads_col=reads_col,
        )

        if output_columns is None:
            output_columns = [
                col for col in df.columns if not col.startswith("_ytab_")
            ]
            detected_chrom_col = c_col
            detected_position_col = p_col
            detected_reads_col = r_col
        else:
            current_cols = [
                col for col in df.columns if not col.startswith("_ytab_")
            ]

            # Keep a union of all original columns across files.
            for col in current_cols:
                if col not in output_columns:
                    output_columns.append(col)

            for col in output_columns:
                if col not in df.columns and not col.startswith("_ytab_"):
                    df[col] = ""

        if not keep_zero:
            df = df[df["_ytab_reads"] > 0].copy()

        summaries.append(
            {
                "file": str(path),
                "rows_after_parse": int(len(df)),
                "unique_sites_after_parse": int(
                    df[["_ytab_chrom", "_ytab_position"]]
                    .drop_duplicates()
                    .shape[0]
                ),
                "reads_after_parse": int(df["_ytab_reads"].sum()),
                "chrom_col": c_col,
                "position_col": p_col,
                "reads_col": r_col,
            }
        )

        frames.append(df)

    if not frames:
        raise RuntimeError("No hit files were read.")

    if output_columns is None:
        raise RuntimeError("No output columns were detected.")

    if detected_chrom_col is None or detected_position_col is None or detected_reads_col is None:
        raise RuntimeError("Could not determine output column names.")

    all_hits = pd.concat(frames, ignore_index=True, sort=False)

    # Make sure all output columns exist after concat.
    for col in output_columns:
        if col not in all_hits.columns:
            all_hits[col] = ""

    group_cols = ["_ytab_chrom", "_ytab_position"]

    combined_rows: list[dict[str, object]] = []

    for (chrom, position), group in all_hits.groupby(
        group_cols,
        sort=False,
        dropna=False,
    ):
        row: dict[str, object] = {}

        for col in output_columns:
            if col == detected_chrom_col:
                row[col] = chrom
            elif col == detected_position_col:
                row[col] = int(position)
            elif col == detected_reads_col:
                row[col] = int(group["_ytab_reads"].sum())
            else:
                row[col] = collapse_values(group[col], sep=collapse_sep)

        combined_rows.append(row)

    combined = pd.DataFrame(combined_rows, columns=output_columns)

    # Ensure key numeric columns are clean.
    combined[detected_position_col] = pd.to_numeric(
        combined[detected_position_col],
        errors="coerce",
    ).astype(int)

    combined[detected_reads_col] = pd.to_numeric(
        combined[detected_reads_col],
        errors="coerce",
    ).fillna(0).round().astype(int)

    if not keep_zero:
        combined = combined[combined[detected_reads_col] > 0].copy()

    if sort_output:
        combined = combined.sort_values(
            [detected_chrom_col, detected_position_col],
            kind="mergesort",
        ).reset_index(drop=True)

    metadata = {
        "chrom_col": detected_chrom_col,
        "position_col": detected_position_col,
        "reads_col": detected_reads_col,
    }

    return combined, summaries, metadata


def main() -> None:
    args = parse_args()

    input_paths: list[str] = []

    if args.hits_txt:
        input_paths.extend(args.hits_txt)

    if args.hits_list:
        input_paths.extend(read_hits_list(args.hits_list))

    input_paths = [x for x in input_paths if str(x).strip()]

    if not input_paths:
        raise SystemExit(
            "ERROR: no input hit files provided. Use --hits-txt or --hits-list."
        )

    hit_files = [Path(x).expanduser().resolve() for x in input_paths]

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    sample_name = args.sample_name.strip()

    if not sample_name:
        raise SystemExit("ERROR: --sample-name cannot be empty.")

    if args.output_file:
        requested_output = Path(args.output_file).expanduser()
        output_hits = requested_output.resolve() if requested_output.is_absolute() else outdir / requested_output
        output_hits.parent.mkdir(parents=True, exist_ok=True)
    else:
        output_hits = outdir / f"{sample_name}_hits.txt"
    output_summary_csv = outdir / f"{sample_name}.combine_summary.csv"
    output_manifest_json = outdir / f"{sample_name}.combine_manifest.json"

    print("Combining hit files:")
    for path in hit_files:
        print(f"  {path}")

    print()
    print(f"Output directory: {outdir}")
    print(f"Combined sample:  {sample_name}")
    print(f"Output hits:      {output_hits}")
    print()

    combined, input_summaries, metadata = combine_hits(
        hit_files=hit_files,
        chrom_col=args.chrom_col,
        position_col=args.position_col,
        reads_col=args.reads_col,
        keep_zero=args.keep_zero,
        sort_output=args.sort,
        collapse_sep=args.collapse_sep,
    )

    combined.to_csv(output_hits, sep="\t", index=False)

    summary_df = pd.DataFrame(input_summaries)
    summary_df.to_csv(output_summary_csv, index=False)

    reads_col = metadata["reads_col"]

    manifest = {
        "sample_name": sample_name,
        "output_hits": str(output_hits),
        "n_input_files": len(hit_files),
        "input_files": [str(x) for x in hit_files],
        "combined_unique_sites": int(combined.shape[0]),
        "combined_reads": int(combined[reads_col].sum()),
        "output_summary_csv": str(output_summary_csv),
        "chrom_col": metadata["chrom_col"],
        "position_col": metadata["position_col"],
        "reads_col": metadata["reads_col"],
        "preserved_columns": list(combined.columns),
        "annotation_conflict_policy": (
            "For non-key/non-count columns, unique non-empty values are collapsed "
            f"with '{args.collapse_sep}'."
        ),
    }

    with open(output_manifest_json, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    print("Done.")
    print(f"Combined unique sites: {manifest['combined_unique_sites']}")
    print(f"Combined reads:        {manifest['combined_reads']}")
    print(f"Chrom column:          {manifest['chrom_col']}")
    print(f"Position column:       {manifest['position_col']}")
    print(f"Reads/count column:    {manifest['reads_col']}")
    print(f"Columns preserved:     {len(manifest['preserved_columns'])}")
    print(f"Summary CSV:           {output_summary_csv}")
    print(f"Manifest JSON:         {output_manifest_json}")


if __name__ == "__main__":
    main()
