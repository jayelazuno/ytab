#!/usr/bin/env python3
"""Upload the dataset archives to Zenodo and record where they landed.

The archives built by `ytab_build_dataset_bundles.py` go into one Zenodo record, which gives the collection a DOI to
cite and every archive a stable URL for the app to fetch. The repository's own GitHub-Zenodo link is a different
thing: it archives releases of the code, not these files.

The deposition is left as a **draft** unless --publish is given, so the metadata can be read over on Zenodo before
anything becomes permanent. Publishing cannot be undone: a published record's files are fixed for good.

    export ZENODO_TOKEN=...                     # from https://zenodo.org/account/settings/applications/tokens/new/
    python scripts/local/ytab_upload_zenodo.py --sandbox        # rehearse against sandbox.zenodo.org
    python scripts/local/ytab_upload_zenodo.py                  # the real site, still a draft
    python scripts/local/ytab_upload_zenodo.py --deposition 123456 --publish

The token is read from the environment and never written to a file or printed. Sandbox tokens are issued separately
from real ones at sandbox.zenodo.org.

Every upload is checked: Zenodo returns each file's MD5, which is compared with the MD5 computed while sending it.
On success the record's id, DOI and per-file URLs are written to `<out-dir>/zenodo_record.json`, which is what
`resources/datasets/catalog.tsv` needs to point the Datasets page at the archives.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SITES = {"zenodo": "https://zenodo.org", "sandbox": "https://sandbox.zenodo.org"}
CHUNK = 8 << 20


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", type=Path, default=Path("output/zenodo"))
    parser.add_argument("--metadata", type=Path, default=Path("resources/datasets/zenodo_metadata.json"))
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org (a separate token)")
    parser.add_argument("--deposition", type=int, help="add to this draft instead of creating one")
    parser.add_argument("--publish", action="store_true", help="publish the draft; this is permanent")
    parser.add_argument("--dry-run", action="store_true", help="list what would be uploaded and stop")
    return parser.parse_args()


def request(method: str, url: str, token: str, data=None, headers=None, timeout=1800):
    head = {"Authorization": f"Bearer {token}"}
    head.update(headers or {})
    req = urllib.request.Request(url, data=data, method=method, headers=head)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            body = response.read()
        return json.loads(body) if body else {}
    except urllib.error.HTTPError as error:
        detail = error.read().decode(errors="replace")[:400]
        raise SystemExit(f"ERROR: {method} {url.split('?')[0]} -> {error.code} {error.reason}\n{detail}") from None


class Reader:
    """Streams a file to the request while computing its MD5, so nothing is read into memory twice."""

    def __init__(self, path: Path):
        self.handle = path.open("rb")
        self.digest = hashlib.md5()
        self.length = path.stat().st_size

    def read(self, size=-1):
        chunk = self.handle.read(CHUNK if size is None or size < 0 else size)
        self.digest.update(chunk)
        return chunk

    def close(self):
        self.handle.close()


def main() -> int:
    args = parse_args()
    root = args.repo_root.expanduser().resolve()
    site = SITES["sandbox" if args.sandbox else "zenodo"]
    out_dir = args.out_dir if args.out_dir.is_absolute() else root / args.out_dir
    index = out_dir / "bundles_index.tsv"
    if not index.is_file():
        print(f"ERROR: no index at {index}; run ytab_build_dataset_bundles.py first", file=sys.stderr)
        return 2
    with index.open(encoding="utf-8") as handle:
        archives = list(csv.DictReader(handle, delimiter="\t"))
    missing = [a["archive_name"] for a in archives if not (out_dir / a["archive_name"]).is_file()]
    if missing:
        print(f"ERROR: the index lists archives that are not here: {', '.join(missing)}", file=sys.stderr)
        return 2

    total = sum(int(a["archive_bytes"]) for a in archives)
    print(f"{len(archives)} archives, {total / 1073741824:.2f} GB -> {site}")
    if args.dry_run:
        for a in archives:
            print(f"  {a['archive_name']:44} {int(a['archive_bytes']) / 1048576:8.0f} MB")
        return 0

    token = os.environ.get("ZENODO_TOKEN", "").strip()
    if not token:
        print("ERROR: set ZENODO_TOKEN (a sandbox token when --sandbox)", file=sys.stderr)
        return 2

    metadata_path = args.metadata if args.metadata.is_absolute() else root / args.metadata
    if not metadata_path.is_file():
        print(f"ERROR: no metadata at {metadata_path}", file=sys.stderr)
        return 2
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if any("CHANGE ME" in json.dumps(v) for v in metadata.get("metadata", {}).values()):
        print(f"ERROR: {metadata_path} still has placeholder fields; fill them in first", file=sys.stderr)
        return 2

    if args.deposition:
        deposition = request("GET", f"{site}/api/deposit/depositions/{args.deposition}", token)
    else:
        deposition = request("POST", f"{site}/api/deposit/depositions", token,
                             data=json.dumps(metadata).encode(), headers={"Content-Type": "application/json"})
    bucket = deposition["links"]["bucket"]
    print(f"Draft {deposition['id']}: {deposition['links'].get('html', '')}")

    already = {f["filename"] for f in deposition.get("files", [])}
    for entry in archives:
        name = entry["archive_name"]
        if name in already:
            print(f"  {name:44} already uploaded, skipping")
            continue
        path = out_dir / name
        reader = Reader(path)
        try:
            result = request("PUT", f"{bucket}/{urllib.parse.quote(name)}", token, data=reader,
                             headers={"Content-Type": "application/octet-stream",
                                      "Content-Length": str(reader.length)})
        finally:
            reader.close()
        returned = str(result.get("checksum", "")).split(":")[-1]
        ok = returned == reader.digest.hexdigest()
        print(f"  {name:44} {reader.length / 1048576:8.0f} MB  md5 {'OK' if ok else 'MISMATCH'}")
        if not ok:
            print("ERROR: Zenodo stored a different file than was sent; not continuing", file=sys.stderr)
            return 1

    if args.publish:
        deposition = request("POST", f"{site}/api/deposit/depositions/{deposition['id']}/actions/publish", token)
        print(f"Published: {deposition.get('doi_url') or deposition.get('links', {}).get('record_html', '')}")
    else:
        print("Left as a draft. Review it on Zenodo, then rerun with --deposition "
              f"{deposition['id']} --publish when the metadata is right.")

    record = dict(site=site, deposition_id=deposition["id"], doi=deposition.get("doi") or "",
                  record_url=deposition.get("links", {}).get("record_html", ""),
                  published=bool(args.publish),
                  files={a["archive_name"]: f"{site}/record/{deposition.get('record_id') or deposition['id']}"
                                            f"/files/{urllib.parse.quote(a['archive_name'])}" for a in archives})
    (out_dir / "zenodo_record.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {out_dir / 'zenodo_record.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
