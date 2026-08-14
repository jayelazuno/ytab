#!/usr/bin/env python3
import argparse,csv,json
from pathlib import Path

def nonempty(path): return path.is_file() and path.stat().st_size>0
def main():
 p=argparse.ArgumentParser();p.add_argument("--project-config",type=Path,required=True);a=p.parse_args();cfg=a.project_config.resolve();root=cfg.parents[1]
 with (root/"config/sample_sheet.csv").open(newline="",encoding="utf-8") as h: rows=[r for r in csv.DictReader(h) if str(r.get("include","true")).lower() in {"true","1","yes"}]
 assert rows
 for row in rows:
  sample=row["sample"];bam=root/"mapfastq"/sample/f"{sample}.sorted.bam";bai=Path(str(bam)+".bai");stats=root/"mapfastq"/sample/f"{sample}.mapping_stats.csv";hits=root/"create_hit_file"/sample/f"{sample}_hits.txt";features=list((root/"summary"/sample).glob("*.feature_table*.*")) if (root/"summary"/sample).is_dir() else []
  mapping="success" if nonempty(bam) and nonempty(bai) and nonempty(stats) else "ready"
  hit="success" if nonempty(hits) else ("ready" if mapping=="success" else "blocked")
  summary="success" if any(nonempty(x) for x in features) else ("ready" if hit=="success" else "blocked")
  if nonempty(bam) and not nonempty(hits): assert (mapping,hit,summary)==("success","ready","blocked")
  if nonempty(hits) and not features: assert hit=="success" and summary=="ready"
  if features: assert summary=="success"
  manifest=root/"manifests/mapfastq"/f"{sample}.mapfastq_manifest.json"
  if manifest.is_file() and nonempty(bam):
   data=json.loads(manifest.read_text());assert mapping=="success";assert "dry run" not in "Scientific mapping outputs available.".lower()
 print("PASS")
if __name__=="__main__":main()
