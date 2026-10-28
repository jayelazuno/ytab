#!/usr/bin/env python3
"""Merge sequencing runs of one library into one project sample, so analyses count the library once.

A library sequenced on several lanes or runs has one row per run in the sample sheet. Comparisons and essentiality
treat each included sample as an independent library, so the runs must be summed first: this writes
<project>/create_hit_file/<library>/<library>_hits.txt (read counts summed per chromosome, strand and position), adds the
library to the sample sheet (role, condition, strain and pool taken from its runs, which must agree) and excludes the
runs, keeping their rows and QC. Then build the library's gene summary with `ytab summary --samples LIBRARY`.

  ytab merge-runs --project-config output/projects/P/config/project.yaml --library P_lib1 --runs RUN1 RUN2
"""
import argparse
import csv
from pathlib import Path
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
SAME = ("library_role", "condition", "background", "pool", "transposon", "library_design")


def resolve(path, base):
    path = Path(path)
    return path if path.is_absolute() else (base/path).resolve()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", required=True, type=Path)
    p.add_argument("--library", required=True, help="Name of the merged sample (letters, numbers, - and _)")
    p.add_argument("--runs", required=True, nargs="+", help="Sample IDs of the runs of this library")
    p.add_argument("--force", action="store_true", help="Replace an existing merged sample of this name")
    a = p.parse_args()
    cfg = yaml.safe_load(a.project_config.read_text()); base = resolve(cfg.get("repo_root", "."), a.project_config.parent)
    project = resolve(cfg["output_project_dir"], base); sheet_path = resolve(cfg.get("sample_sheet") or a.project_config.parent/"sample_sheet.csv", base)
    with open(sheet_path, newline="", encoding="utf-8") as h:rows = list(csv.DictReader(h)); fields = list(rows[0])
    by = {r["sample"]:r for r in rows}
    if len(a.runs) < 2:sys.exit("give at least two runs")
    if not a.library.replace("-", "").replace("_", "").isalnum():sys.exit("library name may contain only letters, numbers, - and _")
    missing = [r for r in a.runs if r not in by]
    if missing:sys.exit(f"not in the sample sheet: {', '.join(missing)}")
    if a.library in by and not a.force:sys.exit(f"{a.library} is already a sample; use --force to replace it")
    differ = [f for f in SAME if f in fields and len({by[r][f] for r in a.runs}) > 1]
    if differ:sys.exit(f"the runs differ in {', '.join(differ)}; they are not one library")

    merged, header = {}, None
    for run in a.runs:
        with open(project/"create_hit_file"/run/f"{run}_hits.txt", newline="") as h:
            reader = csv.reader(h, delimiter="\t"); header = next(reader)
            pos, cnt = header.index("Hit position"), header.index("Hit count")
            for row in reader:
                key = (row[0], row[1], row[pos]); kept = merged.setdefault(key, row[:])
                if kept is not row:kept[cnt] = str(int(kept[cnt])+int(row[cnt]))
    out = project/"create_hit_file"/a.library; out.mkdir(parents=True, exist_ok=True)
    with open(out/f"{a.library}_hits.txt", "w", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n"); w.writerow(header); w.writerows(merged.values())

    first = by[a.runs[0]]; runs = " ".join(a.runs)
    library = {f:"" for f in fields} | {f:first[f] for f in SAME + ("study", "culture_replicate") if f in fields}
    library.update(sample=a.library, display_name=a.library, include="True", notes=f"Library of sequencing runs {runs} (ytab merge-runs: read counts summed per site)")
    rows = [r for r in rows if r["sample"] != a.library]
    for r in rows:
        if r["sample"] in a.runs:
            r["include"] = "False"; note = f"sequencing run of library {a.library}, analysed as part of it"
            if note not in r.get("notes", ""):r["notes"] = "; ".join(x for x in (r.get("notes", ""), note) if x)
    at = max(i for i, r in enumerate(rows) if r["sample"] in a.runs)+1; rows.insert(at, library)
    tmp = sheet_path.with_suffix(".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fields, lineterminator="\n"); w.writeheader(); w.writerows(rows)
    tmp.replace(sheet_path)
    print(f"{a.library}: {len(a.runs)} runs, {len(merged)} insertion sites -> {out/(a.library+'_hits.txt')}; runs excluded in {sheet_path.name}")


if __name__ == "__main__":main()
