#!/usr/bin/env python3
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"));from ytab.pipeline.export_bundle import build_export_bundle
def main():
 p=argparse.ArgumentParser(description="Create a compact reproducible YTAB project export.");p.add_argument("--project-config",type=Path,required=True);p.add_argument("--output",type=Path);p.add_argument("--include-logs",action="store_true");p.add_argument("--include-browser-tracks",action="store_true");p.add_argument("--include-hit-files",action="store_true");p.add_argument("--include-bams",action="store_true");p.add_argument("--force",action="store_true");p.add_argument("--dry-run",action="store_true");a=p.parse_args()
 try:r=build_export_bundle(a.project_config,a.output,a.include_logs,a.include_browser_tracks,a.include_hit_files,a.include_bams,a.force,a.dry_run)
 except Exception as e:print(f"ERROR: {e}",file=sys.stderr);return 1
 print(f"Status: {r['status']}\nArchive: {r['archive']}")
 if a.dry_run:
  [print(f"{x['size_bytes']:>10}  {x['archive_path']}") for x in r['files']];print("Estimated bytes:",r["estimated_size_bytes"])
 else:print(f"Contents: {r['contents']}\nSHA256: {r['checksum']}\nArchive bytes: {r['size_bytes']}")
 return 0
if __name__=="__main__":raise SystemExit(main())
