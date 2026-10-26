#!/usr/bin/env python3
"""Bundle reprocessed projects so others can open them in the YTAB app without rerunning anything.

A bundle holds each project folder (configuration, hit files, QC, normalization, condition comparisons, essentiality,
gene explorer and screen-planning results; no FASTQ or BAM files), the cross-species runs that use only these
projects, and MANIFEST.tsv with the size and SHA-256 of every file. Unpacked at the repository root, the projects
appear in the app's project list.

  python scripts/local/ytab_resource_bundle.py --project Gale_2020 --project Segal_2018_albicans ... --output bundle.tar.gz
  python scripts/local/ytab_resource_bundle.py ... --list     # sizes only, writes nothing
"""
import argparse
import csv
import hashlib
import io
from pathlib import Path
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[2]
SKIP_SUFFIXES = (".bam", ".bai", ".fastq", ".fastq.gz", ".fq.gz", ".sam")


def files(projects):
    out = []
    for p in projects:
        base = ROOT/"output/projects"/p
        if not (base/"config/project.yaml").exists():sys.exit(f"not a YTAB project: {base}")
        out += [f for f in base.rglob("*") if f.is_file() and not f.name.endswith(SKIP_SUFFIXES)]
    for run in sorted((ROOT/"output/cross_species").glob("*")):
        manifest = run/"manifest.json"
        if manifest.exists():
            import json
            used = {x["project"] for x in json.loads(manifest.read_text()).get("projects", [])}
            if used and used <= set(projects):out += [f for f in run.rglob("*") if f.is_file()]
    return sorted(out)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project", action="append", required=True, help="Project ID under output/projects (repeat)")
    p.add_argument("--output", type=Path); p.add_argument("--list", action="store_true")
    a = p.parse_args(); paths = files(a.project)
    total = sum(f.stat().st_size for f in paths)
    print(f"{len(paths)} files, {total/1e9:.2f} GB")
    if a.list or not a.output:return
    rows = [(str(f.relative_to(ROOT)), f.stat().st_size, sha256(f)) for f in paths]
    with tarfile.open(a.output, "w:gz") as tar:
        buf = io.StringIO(); w = csv.writer(buf, delimiter="\t", lineterminator="\n"); w.writerow(["path", "bytes", "sha256"]); w.writerows(rows)
        data = buf.getvalue().encode(); info = tarfile.TarInfo("MANIFEST.tsv"); info.size = len(data); tar.addfile(info, io.BytesIO(data))
        for f in paths:tar.add(f, arcname=str(f.relative_to(ROOT)))
    print("wrote", a.output)


if __name__ == "__main__":main()
