"""Lightweight HTML summary of existing YTAB outputs only."""
from __future__ import annotations
import csv, html, json, subprocess, webbrowser
from datetime import datetime, timezone
from pathlib import Path
from .project_status import build_project_status, load_project_status_context, write_project_status
def _rows(path,limit=None):
    if not path.is_file(): return []
    with path.open(newline="",encoding="utf-8") as h: rows=list(csv.DictReader(h))
    return rows[:limit] if limit else rows
def _table(rows):
    if not rows: return "<p>Not available.</p>"
    cols=list(rows[0]); return "<div class=scroll><table><thead><tr>"+"".join(f"<th>{html.escape(c)}</th>" for c in cols)+"</tr></thead><tbody>"+"".join("<tr>"+"".join(f"<td>{html.escape(str(r.get(c,'')))}</td>" for c in cols)+"</tr>" for r in rows)+"</tbody></table></div>"
def build_project_report(project_config:Path,output:Path|None=None,force=False,open_browser=False):
    c=load_project_status_context(project_config); status=build_project_status(project_config); write_project_status(project_config,status)
    out=Path(output).resolve() if output else c["export"]/"report"/"ytab_project_summary.html"; meta=out.with_name("ytab_project_summary_metadata.json")
    if out.exists() and not force: return {"report":out,"metadata":meta,"status":"skipped"}
    out.parent.mkdir(parents=True,exist_ok=True); cfg=c["config"]; samples=c["samples"]; parents=[x for x in samples if str(x.get("guessed_condition") or x.get("condition")).lower()=="parent"]
    treated=[x for x in samples if str(x.get("guessed_condition") or x.get("condition")).lower()=="treated"]
    ref=cfg.get("reference") or {}; project=c["project"]
    latest=lambda pattern: sorted(project.glob(pattern))[-1] if list(project.glob(pattern)) else Path("/")
    sections=[("Project overview",_table([{"project_id":cfg.get("project_id"),"species":cfg.get("species"),"samples":len(samples),"parents":len(parents),"treated":len(treated)}])),
      ("Sample design",_table([{k:v for k,v in x.items() if k in {"sample","guessed_condition","guessed_background","guessed_pool","include"}} for x in samples])),
      ("Pipeline stage status",_table(status["stages"])),("Reference resources",_table([{k:str(v) for k,v in ref.items() if not isinstance(v,(dict,list))}])),
      ("Mapping summary",_table(_rows(project/"summary"/"summary_stats.all_samples.csv"))),
      ("Library diagnostics",_table(_rows(project/"library_diagnostics"/"library_diagnostics_summary.csv"))),
      ("Normalization target evaluation",_table(_rows(project/"sample_normalization"/"normalization_target_evaluation.csv"))),
      ("Combined parent library",_table(_rows(project/"manifests"/"summary_combined"/"summary_combined_status.csv"))),
      ("Essentiality-classifier summary",_table(_rows(latest("classifier/*/classifier_summary.*.csv")))),
      ("Treated-versus-parent summary",_table(_rows(latest("treated_vs_parent/*/treated_vs_parent_comparison_summary.csv")))),
      ("Reproducibility and provenance",f"<p>Project configuration: <code>{html.escape(str(c['project_config']))}</code></p><p>Repository commit: <code>{html.escape(_git(c['root']))}</code></p>"),
      ("Warnings and limitations","".join(f"<p class=warning>{html.escape(w)}</p>" for w in status["warnings"]) or "<p>None recorded.</p>"),
      ("Stable result files",_links(c,out.parent))]
    body="".join(f"<section><h2>{html.escape(title)}</h2>{content}</section>" for title,content in sections)
    doc=f"<!doctype html><html><head><meta charset=utf-8><title>YTAB project summary</title><style>body{{font:14px system-ui;margin:2rem;max-width:1200px}}table{{border-collapse:collapse}}th,td{{border:1px solid #ccc;padding:.35rem;text-align:left}}th{{background:#eee}}.scroll{{overflow:auto}}.warning{{background:#fff4cc;padding:.7rem}}</style></head><body><h1>YTAB project summary: {html.escape(str(cfg.get('project_id')))}</h1><p>This report summarizes existing workflow outputs and does not add biological interpretation.</p>{body}</body></html>"
    out.write_text(doc,encoding="utf-8"); data={"project_id":cfg.get("project_id"),"report":str(out),"generated_at":datetime.now(timezone.utc).isoformat(),"project_config":str(c["project_config"]),"repository_commit":_git(c["root"])}; meta.write_text(json.dumps(data,indent=2)+"\n")
    if open_browser: webbrowser.open(out.as_uri())
    return {"report":out,"metadata":meta,"status":"success"}
def _git(root):
    try:return subprocess.check_output(["git","rev-parse","HEAD"],cwd=root,text=True,stderr=subprocess.DEVNULL).strip()
    except Exception:return "unavailable"
def _links(c,base):
    candidates=[]
    for pattern in ("classifier/*/essentiality_predictions.*.csv","classifier/*/classifier_summary.*.csv","treated_vs_parent/*/treated_vs_parent_results.csv","treated_vs_parent/*/treated_vs_parent_comparison_summary.csv","sample_normalization/normalization*recommendation*.csv","gene_domain_explorer/figures/*.png","gene_domain_explorer/tables/*.csv","gene_domain_explorer/manifests/*.json"):
        candidates.extend(c["project"].glob(pattern))
    return "<ul>"+"".join(f"<li><code>{html.escape(str(p))}</code></li>" for p in sorted(candidates))+"</ul>"
