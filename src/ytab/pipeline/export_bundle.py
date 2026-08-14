"""Safe compact reproducibility bundle construction."""
from __future__ import annotations
import csv, hashlib, io, tarfile
from datetime import datetime, timezone
from pathlib import Path
from .project_status import load_project_status_context
def _sha(path):
 d=hashlib.sha256()
 with path.open("rb") as h:
  for chunk in iter(lambda:h.read(1024*1024),b""):d.update(chunk)
 return d.hexdigest()
def plan_export(project_config:Path,include_logs=False,include_browser_tracks=False,include_hit_files=False,include_bams=False):
 c=load_project_status_context(project_config); files=[]; project=c["project"]; export=c["export"]
 fixed=[project/"config"/n for n in ("project.yaml","sample_sheet.csv","reference_resolved.json","comparison_design.csv","final_classifier_target.txt","final_classifier_target.json")]
 files += [p for p in fixed if p.is_file()]; files += [p for p in (project/"manifests").rglob("*") if p.is_file()]
 for base in (project,export):
  for p in base.rglob("*"):
   if not p.is_file() or p.is_symlink(): continue
   rel=p.relative_to(base); low=p.name.lower(); parts=set(rel.parts)
   if "bundles" in parts or "manifests" in parts or "config" in parts: continue
   if low.endswith((".fastq",".fastq.gz",".fq",".fq.gz",".sam",".cram")): continue
   if low.endswith((".bam",".bai")) and not include_bams: continue
   if ("hit" in low and low.endswith((".txt",".tsv",".csv"))) and not include_hit_files: continue
   if "logs" in parts and not include_logs: continue
   if low.endswith((".bed",".bedgraph",".wig",".bw",".bigwig")) and not include_browser_tracks: continue
   if any(token in low for token in ("feature_table","long.raw","feature_reads_cpm","_analysis.csv")): continue
   stable=low.endswith((".csv",".tsv",".json",".yaml",".yml",".png",".svg",".html",".txt"))
   wanted=(low.startswith(("essentiality_predictions.","classifier_summary.","classifier_run_metadata.","treated_vs_parent_results","treated_vs_parent_comparison_summary","treated_vs_parent_run_metadata","normalization_","library_diagnostics","summary_stats","stats.csv","ytab_project_summary")) or low.endswith((".png",".svg")))
   if stable and wanted: files.append(p)
 root=c["root"].resolve(); unique=[]; seen=set()
 for p in files:
  resolved=p.resolve()
  if root not in resolved.parents: raise ValueError(f"Export candidate is outside repository: {p}")
  arc=str(resolved.relative_to(root))
  if ".." in Path(arc).parts or arc in seen: continue
  seen.add(arc); unique.append((resolved,arc))
 return c,sorted(unique,key=lambda x:x[1])
def build_export_bundle(project_config:Path,output=None,include_logs=False,include_browser_tracks=False,include_hit_files=False,include_bams=False,force=False,dry_run=False):
 c,files=plan_export(project_config,include_logs,include_browser_tracks,include_hit_files,include_bams); stamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"); directory=c["export"]/"bundles"; directory.mkdir(parents=True,exist_ok=True)
 archive=Path(output).resolve() if output else directory/f"ytab_{c['config']['project_id']}_{stamp}.tar.gz"; contents=archive.with_suffix("").with_suffix(".contents.csv"); checksum=archive.with_suffix(archive.suffix+".sha256")
 rows=[{"archive_path":arc,"size_bytes":p.stat().st_size,"sha256":_sha(p)} for p,arc in files]; total=sum(r["size_bytes"] for r in rows)
 if dry_run:return {"status":"dry-run","files":rows,"estimated_size_bytes":total,"archive":archive}
 if archive.exists() and not force: raise FileExistsError(f"Archive exists: {archive}")
 readme=("YTAB compact reproducibility bundle\n\nIncludes configuration, manifests, stable summaries, plots, results, and report.\n"
         "FASTQ, BAM, SAM, hit files, and normalized hit files are excluded unless explicitly requested.\n")
 with tarfile.open(archive,"w:gz") as tar:
  for p,arc in files: tar.add(p,arcname=arc,recursive=False)
  info=tarfile.TarInfo("BUNDLE_README.txt"); payload=readme.encode(); info.size=len(payload); tar.addfile(info,io.BytesIO(payload))
 with contents.open("w",newline="",encoding="utf-8") as h:
  w=csv.DictWriter(h,fieldnames=("archive_path","size_bytes","sha256")); w.writeheader(); w.writerows(rows)
 checksum.write_text(f"{_sha(archive)}  {archive.name}\n",encoding="utf-8")
 return {"status":"success","archive":archive,"contents":contents,"checksum":checksum,"file_count":len(files),"size_bytes":archive.stat().st_size}
