#!/usr/bin/env python3
"""Curated C. glabrata essentiality labels from CGD (direct C. glabrata evidence) and S. cerevisiae 1:1 orthologs.

Rules are fixed here and do not look at any YTAB model output:
- C. glabrata direct (gold): CGD viable/inviable rows that do NOT cite Gale et al. 2020 (PMID 32819971).
  inviable with mutant type null or repressible -> essential; condition other than repression
  (e.g. "on media lacking riboflavin") -> conditional (not labelled). viable with mutant type null -> nonessential.
  Both essential and nonessential evidence -> conflict (not labelled).
- S. cerevisiae (silver): SGD null-mutant viable/inviable. Systematic deletion sets (haploid or homozygous
  diploid) decide when they agree; otherwise classical genetics when it agrees; mixed -> conflict.
  Transferred only through 1:1 orthogroups (exactly one C. glabrata protein and one S. cerevisiae gene).
- Label: gold where C. glabrata direct evidence exists, else silver; otherwise unlabelled with a reason.
- Gale et al. 2020 calls (CGD rows citing PMID 32819971) are kept as a comparator column, never as labels.
Gene universe and IDs: BG2 GFF (GWK60_*), CBS138 IDs (CAGL0*) from the GFF Note field.
Writes resources/species/glabrata/essentiality_labels/{labels.tsv, summary.json, sources.json}.
"""
from collections import defaultdict
from datetime import date
import csv
import gzip
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.normalization.ttr import sha256

BASE = ROOT/"resources/species/glabrata"
OUT = BASE/"essentiality_labels"; RAW = OUT/"raw"
GFF = BASE/"reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff"
ORTHO = BASE/"orthologs/C_glabrata_BG2_S_cerevisiae_orthologs.txt"
CGD = RAW/"C_glabrata_CBS138_phenotype_data.tab"; SGD = RAW/"SGD_phenotype_data.tab.gz"
GALE_PMID = "32819971"
SOURCES = dict(cgd=dict(url="http://www.candidagenome.org/download/phenotype/C_glabrata_CBS138_phenotype_data.tab", retrieved="2026-10-19"),
               sgd=dict(url="https://downloads.yeastgenome.org/curation/literature/phenotype_data.tab", retrieved="2026-10-19",
                        stored="gzip of the downloaded file"))


def bg2_genes():
    genes = {}
    for line in open(GFF):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if a.get("gene_biotype") != "protein_coding":continue
        note = a.get("Note", ""); cagl = re.match(r"(CAGL0[A-M]\d+g)~", note)
        genes[a["locus_tag"]] = dict(cglab_id=a["locus_tag"], name=a.get("Name", ""), cagl_id=cagl.group(1) if cagl else "",
                                     chromosome=f[0], start=int(f[3]), end=int(f[4]))
    return genes


def one_to_one():
    """cglab_id -> S. cerevisiae systematic name for 1:1 orthogroups; also orthogroup type per cglab_id."""
    groups = defaultdict(lambda:dict(cg=set(), sc=set(), rows=[]))
    for r in csv.DictReader(open(ORTHO), delimiter="\t"):
        g = groups[r["Orthogroup"]]
        g["cg"] |= set(re.findall(r"QNG\d+\.\d", r["C_glabrata_raw"] or ""))
        g["sc"] |= set(re.findall(r"SGDID_(S\d+)", r["S_cerevisiae_raw"] or ""))
        g["rows"].append(r)
    mapping, kind = {}, {}
    for og, g in groups.items():
        k = "1:1" if len(g["cg"]) == 1 and len(g["sc"]) == 1 else f"{len(g['cg'])}:{len(g['sc'])}"
        for r in g["rows"]:
            kind[r["cglab_id"]] = (k, og, r["scer_id"], r["scer_gene_name"])
            if k == "1:1":mapping[r["cglab_id"]] = r["scer_id"]
    return mapping, kind


def sgd_essentiality():
    calls = defaultdict(lambda:defaultdict(set))
    with gzip.open(SGD, "rt") as f:
        for r in csv.reader(f, delimiter="\t"):
            if len(r) < 10 or r[6] != "null" or r[9] not in ("viable", "inviable"):continue
            source = "systematic" if "systematic mutation set" in r[5] else "classical" if "classical genetics" in r[5] else None
            if source:calls[r[0]][source].add(r[9])
    out = {}
    for gene, c in calls.items():
        for source in ("systematic", "classical"):
            if c[source]:
                out[gene] = ("conflict" if len(c[source]) > 1 else "essential" if c[source] == {"inviable"} else "nonessential", source); break
    return out


def cgd_evidence():
    direct, gale = defaultdict(list), {}
    for r in csv.DictReader(open(CGD), delimiter="\t"):
        if r["Phenotype"] not in ("viable", "inviable"):continue
        pmid = (re.search(r"PMID: (\d+)", r["Reference"]) or [None, ""])[1]
        if pmid == GALE_PMID:
            gale.setdefault(r["Feature Name"], set()).add(r["Phenotype"]); continue
        cond = r["Condition"].strip()
        if r["Phenotype"] == "inviable" and r["Mutant Type"] in ("null", "repressible"):
            e = "conditional" if cond and "repress" not in cond.lower() else "essential"
        elif r["Phenotype"] == "viable" and r["Mutant Type"] == "null":e = "nonessential"
        else:continue
        direct[r["Feature Name"]].append(dict(evidence=e, pmid=pmid, strain=r["Strain Background"], condition=cond, mutant=r["Mutant Type"]))
    resolved = {}
    for gene, ev in direct.items():
        kinds = {x["evidence"] for x in ev}
        call = ("conflict" if {"essential", "nonessential"} <= kinds else "essential" if "essential" in kinds else
                "nonessential" if kinds == {"nonessential"} else "conditional" if "nonessential" not in kinds else "nonessential (conditional inviable)")
        resolved[gene] = (call, ";".join(sorted({x["pmid"] for x in ev})), ";".join(sorted({x["strain"] for x in ev if x["strain"]})),
                          ";".join(sorted({x["condition"] for x in ev if x["condition"]})))
    gale = {g:("conflict" if len(p) > 1 else "essential" if p == {"inviable"} else "nonessential") for g, p in gale.items()}
    return resolved, gale


def main():
    genes = bg2_genes(); mapping, kind = one_to_one(); sc = sgd_essentiality(); direct, gale = cgd_evidence()
    rows = []
    for gid, g in genes.items():
        k, og, scer, scer_name = kind.get(gid, ("no orthogroup", "", "", ""))
        sc_call, sc_source = sc.get(mapping.get(gid, ""), ("", "")) if gid in mapping else ("", "")
        cg_call, pmids, strains, conditions = direct.get(g["cagl_id"], ("", "", "", ""))
        if cg_call in ("essential", "nonessential"):label, tier, reason = cg_call, "gold", ""
        elif cg_call == "nonessential (conditional inviable)":label, tier, reason = "nonessential", "gold", "inviable only under a nutrient condition"
        elif cg_call in ("conflict", "conditional"):label, tier, reason = "", "", f"C. glabrata direct evidence {cg_call}"
        elif gid not in mapping:label, tier, reason = "", "", f"orthology {k}"
        elif sc_call in ("essential", "nonessential"):label, tier, reason = sc_call, "silver", ""
        else:label, tier, reason = "", "", "S. cerevisiae conflict" if sc_call == "conflict" else "S. cerevisiae ortholog without viability data"
        rows.append(g | dict(orthology=k, orthogroup=og, scer_id=scer if k == "1:1" else "", scer_name=scer_name if k == "1:1" else "",
                             sc_call=sc_call, sc_source=sc_source, cg_direct_call=cg_call, cg_direct_pmids=pmids, cg_direct_strains=strains,
                             cg_direct_conditions=conditions, label=label, tier=tier, unlabelled_reason=reason,
                             gale2020_comparator=gale.get(g["cagl_id"], "")))
    keys = list(rows[0])
    with open(OUT/"labels.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, keys, delimiter="\t"); w.writeheader(); w.writerows(rows)
    count = lambda pred:sum(1 for r in rows if pred(r))
    both = [r for r in rows if r["cg_direct_call"] in ("essential", "nonessential") and r["sc_call"] in ("essential", "nonessential")]
    labelled = [r for r in rows if r["label"] and r["gale2020_comparator"] in ("essential", "nonessential")]
    summary = dict(
        genes=len(rows), with_cagl_id=count(lambda r:r["cagl_id"]), one_to_one=count(lambda r:r["orthology"] == "1:1"),
        labels={t:{l:count(lambda r, t=t, l=l:r["tier"] == t and r["label"] == l) for l in ("essential", "nonessential")} for t in ("gold", "silver")},
        unlabelled=dict(sorted({x:count(lambda r, x=x:not r["label"] and r["unlabelled_reason"] == x)
                                for x in {r["unlabelled_reason"] for r in rows if not r["label"]}}.items(), key=lambda kv:-kv[1])),
        gold_vs_sc_ortholog=dict(genes=len(both), agree=sum(r["cg_direct_call"] == r["sc_call"] for r in both),
                                 discordant=[f"{r['name']} (Cg {r['cg_direct_call']}, Sc {r['sc_call']})" for r in both if r["cg_direct_call"] != r["sc_call"]]),
        comparator_only_gale2020=dict(labelled_genes_with_gale_call=len(labelled),
            agree=sum(r["label"] == r["gale2020_comparator"] for r in labelled),
            by_tier={t:dict(n=sum(r["tier"] == t for r in labelled), agree=sum(r["tier"] == t and r["label"] == r["gale2020_comparator"] for r in labelled))
                     for t in ("gold", "silver")}))
    (OUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    sources = SOURCES | dict(files={p.name:sha256(p) for p in (CGD, SGD, ORTHO, GFF)}, built=str(date.today()), script_sha256=sha256(Path(__file__)))
    (OUT/"sources.json").write_text(json.dumps(sources, indent=2)+"\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":main()
