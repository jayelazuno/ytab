#!/usr/bin/env python3
"""Catalogue of the public Tn-seq datasets YTAB has reprocessed, for the app's Datasets page and its search.

Each study is a row of resources/datasets/studies.tsv (citation, DOI, strain, transposon, accession, design); its runs
are the rows of pipeline/hpc/datasets/<dataset>/samples.csv whose `study` column matches. For every study and project
the catalogue adds what the runs are (roles, strains, conditions), their QC (newest `ytab qc` table) and which results
the app can open without computing anything (QC, essentiality, normalization, condition comparisons, genome browser,
library correlation, cross-species).

  ytab datasets                      # writes resources/datasets/catalog.tsv, libraries.tsv and docs/experimental_designs.md
  ytab datasets --find SRR29743984   # which study, library and analyses a run belongs to (also: a strain, drug or condition)

libraries.tsv gives each run its library (runs merged by `ytab merge-runs` share one) and every analysis that uses it, with
its role there (control, treated, reference condition, ...); experimental_designs.md lists the same per study, with the
commands that reproduce each analysis.

A new dataset appears once it has a dataset folder and a row in studies.tsv; nothing here names a study.
"""
import argparse
import csv
from collections import Counter
import json
from pathlib import Path
from statistics import median

import yaml

ROOT = Path(__file__).resolve().parents[2]
STUDIES = ROOT/"resources/datasets/studies.tsv"
DATASETS = ROOT/"pipeline/hpc/datasets"
PROJECTS = ROOT/"output/projects"

# result folder -> (app tab, file that marks a finished run)
RESULTS = [("essentiality", "Essentiality", "models/genes.tsv"), ("normalization", "Site normalization", "manifest.json"),
           ("treated_vs_parent", "Treated vs control (z-score)", "zscore_manifest.json"), ("treated_vs_control", "Treated vs control (resampling/Gumbel)", "manifest.json"),
           ("multicondition", "Multi-condition", "manifest.json"), ("genetic_interactions", "Genetic interactions", "run_manifest.json"),
           ("library_correlation", "Library correlation", "manifest.json")]


def read_tsv(path, delimiter="\t"):
    with open(path, newline="", encoding="utf-8") as h:return list(csv.DictReader(h, delimiter=delimiter))


def env(path):
    return dict(line.strip().split("=", 1) for line in open(path) if "=" in line and not line.startswith("#"))


def finished_runs(folder, marker):
    """Runs under a result folder that have their marker file (one level, or two for per-method folders)."""
    if not folder.is_dir():return []
    runs = [d for d in folder.iterdir() if d.is_dir() and ((d/marker).exists() or any((c/marker).exists() for c in d.iterdir() if c.is_dir()))]
    return sorted(r.name for r in runs)


def newest_qc(project):
    runs = sorted((PROJECTS/project/"qc").glob("qc_*/qc_table.tsv"))
    return {r["sample"]:r for r in read_tsv(runs[-1])} if runs else {}


def cross_species_runs():
    out = {}
    for m in sorted((ROOT/"output/cross_species").glob("*/manifest.json")):
        for p in json.loads(m.read_text()).get("projects", []):out.setdefault(p.get("project"), []).append(m.parent.name)
    return out


def number(value):
    try:return float(value)
    except (TypeError, ValueError):return None


def counts(values):
    return "; ".join(f"{k} {n}" for k, n in Counter(v for v in values if v).most_common())


def read_json(path):
    try:return json.loads(Path(path).read_text())
    except (OSError, ValueError):return {}


def project_sheet(project):
    """The project's samples (library rows of merged runs carry `runs`) and run -> library."""
    rows = read_tsv(PROJECTS/project/"config/sample_sheet.csv", ",") if (PROJECTS/project/"config/sample_sheet.csv").exists() else []
    library = {run:r["sample"] for r in rows for run in (r.get("runs") or "").split()}
    return rows, library


def experiments(project):
    """Every analysis of the project the app streams, with the role of each sample in it: [{kind, name, groups:
    [(label, {sample: role})], command}]. Read from the analyses' own manifests and the comparison design."""
    root, out = PROJECTS/project, []
    rel = lambda p:str(p).replace(f"{ROOT}/", "")   # paths in older manifests are absolute
    cfg = f"output/projects/{project}/config/project.yaml"
    for run in finished_runs(root/"essentiality", "models/genes.tsv"):
        m = read_json(root/"essentiality"/run/"manifest.json")
        out.append(dict(kind="Essentiality (HMM)", name=run, groups=[("libraries", {s:"library" for s in m.get("libraries", [])})],
                        command=f"ytab essentiality --project-config {cfg} --libraries {' '.join(m.get('libraries', []))} --output {rel(root/'essentiality'/run)}"))
    design = read_tsv(root/"config/comparison_design.csv", ",") if (root/"config/comparison_design.csv").exists() else []
    design = [r for r in design if r.get("include", "").lower() == "true"]
    for run in finished_runs(root/"treated_vs_control", "manifest.json"):
        m = read_json(root/"treated_vs_control"/run/"manifest.json")
        groups = [(f"{g['background']}, {g['condition']}", {**{s:"control" for s in g.get("controls", [])}, **{s:"treated" for s in g.get("treatments", [])}})
                  for g in (m.get("groups") or {}).values()]
        out.append(dict(kind=f"Treated vs control ({m.get('method_label', m.get('method', ''))}, {m.get('normalization', '')})", name=run, groups=groups,
                        command=f"ytab treated-vs-control --project-config {cfg} --method {m.get('method', 'resampling')} --normalization-run "
                                f"{rel(m.get('normalization_run', ''))} --normalization {m.get('normalization', '')} --output {rel(root/'treated_vs_control'/run)}"))
    for manifest in sorted((root/"treated_vs_parent").glob("*/zscore_manifest.json")):
        m = read_json(manifest); ids = set(m.get("comparisons", []))
        rows = [r for r in design if r["comparison_id"] in ids]
        by = {}
        for r in rows:
            g = by.setdefault(f"{r['background']}, {r['treated_condition']}", {}); g[r["parent_sample"]] = "control"; g[r["treated_sample"]] = "treated"
        for bg, b in (m.get("backgrounds") or {}).items():   # controls used only for the noise model
            extra = [s for s in b.get("controls", []) if not any(s in g for g in by.values())]
            if extra:by[f"{bg}, noise model"] = {s:"control" for s in extra}
        out.append(dict(kind=f"Treated vs control (z-score, {m.get('normalization', 'CPM')})", name=manifest.parent.name, groups=list(by.items()),
                        command=f"ytab fitness --project-config {cfg} --analysis-id {manifest.parent.name} --comparisons {','.join(sorted(ids))}"))
    for run in finished_runs(root/"multicondition", "manifest.json"):
        m = read_json(root/"multicondition"/run/"manifest.json"); ref = m.get("reference_condition", "")
        samples = {x["observation_id"].split("::", 1)[-1]:f"{x['condition']}{' (reference)' if x['condition'] == ref else ''}" for x in m.get("samples", [])}
        out.append(dict(kind=f"Multi-condition ({m.get('normalization', '')})", name=run, groups=[(", ".join(m.get("conditions", [])), samples)],
                        command=f"ytab multicondition-project --project-config {cfg} --normalized {rel(m.get('normalized', ''))} --normalization {m.get('normalization', '')} "
                                + " ".join(f"--sample {project}::{s}={c.replace(' (reference)', '')}" for s, c in samples.items())
                                + f" --reference-condition {ref} --output {rel(root/'multicondition'/run)}"))
    for run in finished_runs(root/"genetic_interactions", "run_manifest.json"):
        m = read_json(root/"genetic_interactions"/run/"run_manifest.json")
        samples = {x["observation_id"].split("::", 1)[-1]:f"{x['strain']} {x['condition']}" for x in m.get("samples", [])}
        d, norm = m.get("diagnostics", {}), m.get("normalization", {})
        ref = yaml.safe_load((root/"config/project.yaml").read_text()).get("reference", {})
        command = (f"ytab genetic-interactions --normalized {rel(Path(norm.get('input', '')).parent)} --normalization {norm.get('method', '')} "
                   f"--gff {ref.get('gff', '')} --reference {ref.get('fasta', '')} --reference-strain {d.get('reference_strain', '')} --query-strain {d.get('query_strain', '')} "
                   + " ".join(f"--sample {x['observation_id']}={x['strain']}:{x['condition']}" for x in m.get("samples", [])) + f" --output {rel(root/'genetic_interactions'/run)}")
        out.append(dict(kind="Genetic interactions", name=run, groups=[(f"{d.get('query_strain', '')} vs {d.get('reference_strain', '')}, strain x treatment", samples)], command=command))
    return out


def design_page(catalog, libraries, exps):
    """docs/experimental_designs.md: per study, its libraries (runs merged into libraries) and every analysis with the
    role of each library, and the commands that reproduce them."""
    lines = ["# Experimental designs of the public datasets", "",
             "Which runs make which library, and which libraries each analysis uses, for every reprocessed study. Written by",
             "`ytab datasets` from the sample sheets, `config/comparison_design.csv` and the analyses' manifests; search a run ID,",
             "strain, drug or condition to find where it is used. From the command line: `ytab datasets --find TEXT`.", ""]
    for c in catalog:
        libs = [l for l in libraries if l["study"] == c["study"] and l["project"] == c["project"]]
        mine = {l["library"] for l in libs if l["library"]}
        lines += [f"## {c['citation']}: {c['species']} ({c['project']})", "",
                  f"{c['strain']}; {c['transposon']}, {c['read_design']}. {c['accession']}. {c['design']}.", "",
                  "| library | runs | role | condition | strain | pool |", "|---|---|---|---|---|---|"]
        grouped = {}
        for l in libs:grouped.setdefault(l["library"] or l["sample"], []).append(l)
        for name, ls in grouped.items():
            l = ls[0]; runs = " ".join(x["sample"] for x in ls) if len(ls) > 1 or ls[0]["sample"] != name else ""
            lines.append(f"| {name}{'' if l['library'] else ' (not analysed)'} | {runs} | {l['role']} | {l['condition']} | {l['background']} | {l['pool']} |")
        used = [e for e in exps.get(c["project"], []) if any(s in mine for _, g in e["groups"] for s in g)]
        lines += ["", "**Analyses**" if used else "**Analyses:** parent libraries only (no condition analyses)", ""]
        if (PROJECTS/c["project"]/"config/analyses_workflow.yaml").exists():
            lines += [f"Rerun every condition analysis: `ytab run --workflow output/projects/{c['project']}/config/analyses_workflow.yaml`", ""]
        for e in used:
            lines.append(f"- **{e['kind']}**, `{e['name']}`")
            for label, g in e["groups"]:
                g = {s:r for s, r in g.items() if s in mine}
                if g:lines.append(f"  - {label}: " + "; ".join(s if r == "library" else f"{s} ({r})" for s, r in g.items()))
            lines.append(f"  - `{e['command']}`")
        lines.append("")
    return "\n".join(lines)


def find(libraries, catalog, text):
    words = text.lower().split(); hit = lambda row:all(w in " ".join(str(v) for v in row.values()).lower() for w in words)
    for c in (c for c in catalog if hit({k:c[k] for k in ("study", "citation", "species", "strain", "accession", "design", "conditions", "project")})):
        print(f"{c['citation']} - {c['species']}, {c['strain']} ({c['project']}): {c['runs_included']} runs; {c['in_app']}")
    for l in (l for l in libraries if hit(l)):
        print(f"  {l['sample']}: {l['role']}, {l['condition']}, strain {l['background']}, pool {l['pool']}"
              + (f"; library {l['library']}" if l["library"] and l["library"] != l["sample"] else "") + (f"; used in: {l['used_in']}" if l["used_in"] else "; not analysed"))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output", type=Path, default=ROOT/"resources/datasets", help="Folder for catalog.tsv and libraries.tsv")
    p.add_argument("--designs", type=Path, default=ROOT/"docs/experimental_designs.md", help="Experimental-design page (searchable in the app)")
    p.add_argument("--find", metavar="TEXT", nargs="+", help="Print the studies and runs matching TEXT (run ID, strain, drug, condition) and the analyses using them; writes nothing")
    a = p.parse_args()
    if a.find:
        find(read_tsv(a.output/"libraries.tsv"), read_tsv(a.output/"catalog.tsv"), " ".join(a.find)); return
    cross = cross_species_runs(); catalog, libraries = [], []
    exps = {}
    for study in read_tsv(STUDIES):
        for dataset in study["datasets"].split():
            folder = DATASETS/dataset; project = env(folder/"dataset.env")["PROJECT_ID"]
            runs = [r for r in read_tsv(folder/"samples.csv", ",") if r.get("study") == study["study"]]
            if not runs:continue
            qc = newest_qc(project); root = PROJECTS/project
            included = [r for r in runs if r.get("include", "True").strip().lower() not in ("false", "0", "no")]
            if project not in exps:exps[project] = experiments(project)
            sheet, merged = project_sheet(project); in_sheet = {r["sample"]:r for r in sheet}
            for r in runs:
                q = qc.get(r["sample"], {})
                library = merged.get(r["sample"]) or (r["sample"] if in_sheet.get(r["sample"], {}).get("include", "True") == "True" and r in included else "")
                used = [f"{e['kind']} {e['name']}" + ("" if g[library] == "library" else f" ({label}: {g[library]})") for e in exps[project] for label, g in e["groups"] if library in g]
                libraries.append(dict(study=study["study"], citation=study["citation"], project=project, sample=r["sample"], display_name=r.get("display_name", ""),
                                      role=r["library_role"], condition=r["condition"], background=r["background"], pool=r["pool"], transposon=r["transposon"],
                                      library_design=r["library_design"], included=r in included, input_reads=q.get("input_reads", ""),
                                      percent_mapped=q.get("percent_mapped", ""), insertion_sites=q.get("unique_sites", ""),
                                      sites_per_kb=q.get("insertions_per_kb", ""), top_site_share=q.get("max_share", ""), top_site=q.get("max_site", ""),
                                      library=library, used_in="; ".join(dict.fromkeys(used)), notes=r.get("notes", "")))
            density = [x for x in (number(qc.get(r["sample"], {}).get("insertions_per_kb")) for r in included) if x is not None]
            share = [(number(qc[r["sample"]].get("max_share")), qc[r["sample"]].get("max_site", "")) for r in included if number(qc.get(r["sample"], {}).get("max_share")) is not None]
            top = max(share, default=(None, ""))
            results = {label:finished_runs(root/folder_name, marker) for folder_name, label, marker in RESULTS}
            results["Genome browser"] = ["genome_browser"] if (root/"genome_browser/manifest.json").exists() else []
            results["Cross-species"] = cross.get(project, [])
            tabs = ["Library QC"] if qc else []
            tabs += [f"{label} ({len(v)})" for label, v in results.items() if v]
            catalog.append(dict(study=study["study"], citation=study["citation"], title=study["title"], journal=study["journal"], year=study["year"],
                                doi=study["doi"], pmid=study["pmid"], species=study["species"], strain=study["strain"], transposon=study["transposon"],
                                read_design=study["read_design"], accession=study["accession"], design=study["design"], project=project,
                                runs=len(runs), runs_included=len(included), roles=counts(r["library_role"] for r in included),
                                backgrounds=counts(r["background"] for r in included), conditions=counts(r["condition"] for r in included),
                                median_sites_per_kb=round(median(density), 1) if density else "", max_top_site_share=round(top[0], 3) if top[0] is not None else "",
                                max_top_site=top[1], in_app="; ".join(tabs),
                                essentiality_runs=" ".join(results["Essentiality"]), project_config=f"output/projects/{project}/config/project.yaml"))
    a.output.mkdir(parents=True, exist_ok=True)
    for name, rows in (("catalog.tsv", catalog), ("libraries.tsv", libraries)):
        with open(a.output/name, "w", newline="", encoding="utf-8") as h:
            w = csv.DictWriter(h, list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    a.designs.write_text(design_page(catalog, libraries, exps), encoding="utf-8")
    print(json.dumps(dict(studies=len({r["study"] for r in catalog}), projects=len({r["project"] for r in catalog}), runs=len(libraries),
                          output=str(a.output)), indent=1))


if __name__ == "__main__":main()
