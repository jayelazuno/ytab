#!/usr/bin/env python3
"""Claims check, biology: do the genes each method calls make sense? Per study and direction, the genes called by both
the published analysis and YTAB, by YTAB only and by the paper only, each tested for GO biological-process enrichment
against the genes both tested (hypergeometric, annotations propagated, BH over terms of 5-500 genes;
scripts/paper/go_enrichment.py). S. cerevisiae GO (org.Sc.sgd.db, exported to output/claims_check/go/) through 1:1
orthologs for C. glabrata and C. neoformans; PomBase GO for S. pombe. C. auris has no ortholog table here.

A set that is mostly noise has few or no enriched terms; a set of real responders is enriched for the pathways the
drug or condition acts on. Writes output/claims_check/biology_sets.tsv (set sizes, top terms) and biology_terms.tsv.
"""
from __future__ import annotations
import importlib.util, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]; OUT = ROOT/"output/claims_check"; REPORT = ROOT/"output/comparison/report"
sys.path.insert(0, str(ROOT/"src"))
from ytab.fitness import agreement as ag  # noqa: E402
_spec = importlib.util.spec_from_file_location("go", ROOT/"scripts/paper/go_enrichment.py")
go = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(go)
SP = ROOT/"resources/species"


def scer_annotations():
    a = pd.read_csv(OUT/"go/scer_go_bp_propagated.tsv", sep="\t")
    return a.rename(columns={"ORF":"gene", "GOALL":"term", "TERM":"name"}).assign(aspect="P")[["gene", "term", "aspect", "name"]].drop_duplicates(["gene", "term"])


def pombe_annotations():
    a = go.annotations(SP/"pombe/annotation/go/gene_association.pombase.gz", SP/"pombe/annotation/go/go-basic.obo")
    return a[a.aspect == "P"]


def to_scer(species):
    """Gene ID -> S. cerevisiae ORF (1:1 orthologs only)."""
    if species == "glabrata":
        l = pd.read_csv(SP/"glabrata/essentiality_labels/labels.tsv", sep="\t", dtype=str)
        l = l[l.orthology.eq("1:1") & l.scer_id.notna()] if "orthology" in l else l[l.scer_id.notna()]
        return dict(zip(l.cglab_id, l.scer_id))
    if species == "glabrata_CBS138":
        a = pd.read_csv(SP/"glabrata/strains/CBS138/annotation/gene_annotation.tsv", sep="\t", dtype=str)
        a = a[a.scer_relationship.eq("1:1") & a.scer_gene_id.notna()]
        return dict(zip(a.gene_id, a.scer_gene_id))
    if species == "neoformans":
        l = pd.read_csv(SP/"neoformans/essentiality_labels/labels.tsv", sep="\t", dtype=str)
        l = l[l.scer_ortholog.notna()]
        return dict(zip(l.gene_id, l.scer_ortholog))
    raise ValueError(species)


def sets_from_calls(ours: pd.DataFrame, theirs: pd.DataFrame, direction):
    """ours / theirs: DataFrames with ORF, call (tested genes only). Gene-level sets: called in any condition."""
    tested = set(ours.ORF) & set(theirs.ORF)
    o = set(ours.ORF[ours.call.eq(direction)]) & tested; t = set(theirs.ORF[theirs.call.eq(direction)]) & tested
    return tested, {"both":o & t, "YTAB only":o - t, "published only":t - o}


def study_inputs():
    """(study, species, direction, tested, sets, note) for every study with gene-level published calls."""
    out = []
    def report_pair(dataset, keep_collapsed=True):
        pub = ag.read_published(ROOT/"output/comparison"/dataset/"published.tsv")
        res = ag.read_method(REPORT/dataset/"none/resampling_betageom", "resampling")
        if not keep_collapsed:
            c = pd.read_csv(REPORT/"tables"/f"{dataset}_collapsed_conditions.tsv", sep="\t"); gone = set(zip(c.background[c.collapsed], c.condition[c.collapsed]))
            res = res[[(b, k) not in gone for b, k in zip(res.background, res.condition)]]
            pub = pub[[(b, k) not in gone for b, k in zip(pub.background, pub.condition)]] if "background" in pub else pub
        m = res[res.tested].merge(pub, on=["ORF", "background", "condition"], suffixes=("", "_pub"))
        return m[["ORF", "call"]], m[["ORF", "call_pub"]].rename(columns={"call_pub":"call"})
    o, t = report_pair("Gale_2020")
    for d in ("depleted", "enriched"):out.append(("Gale 2020 (FLU)", "glabrata", d, *sets_from_calls(o, t, d), "resampling BGC vs published |z| >= 3"))
    o, t = report_pair("Nickels_2025_glabrata", keep_collapsed=False)
    o = o[o.ORF.isin(o.ORF)];
    for d in ("depleted", "enriched"):out.append(("Nickels 2025 (echinocandins)", "glabrata", d, *sets_from_calls(o, t, d), "pooled over 22 non-collapsed conditions"))
    # Billmyre 2025: their Bonferroni hits (S1 Table) vs YTAB resampling
    s = pd.read_excel(OUT/"papers/billmyre2025/supp/pbio.3003184.s009.xlsx", sheet_name="Supplemental Table 1")
    fc = "Mean(Log10(Fluconazole/DMSO frequency)) in Gene Body"; pv = "P value of Mean Fold Change in Gene Body Compared to Intergenic Inserts"
    s[pv] = pd.to_numeric(s[pv], errors="coerce"); s = s[pd.to_numeric(s["Insert sites in Gene Body"], errors="coerce") >= 5]
    s["call"] = np.where(s[pv]*len(s) < .05, np.where(s[fc] > 0, "enriched", "depleted"), "none")
    g = pd.read_csv(OUT/"Billmyre_2025_neoformans_cc/resampling_aBGC/genes.tsv", sep="\t", low_memory=False); g = g[g.status == "tested"]
    for d in ("depleted", "enriched"):out.append(("Billmyre 2025 (FLU)", "neoformans", d, *sets_from_calls(g[["ORF", "call"]], s.rename(columns={"Gene":"ORF"})[["ORF", "call"]], d),
                                                  "one of two DMSO libraries (the other lacks a hit file)"))
    # Billmyre 2022: 532 required genes vs YTAB depleted after sex (long outgrowth)
    r = pd.read_excel(OUT/"papers/billmyre2022/supp/pgen.1010462.s018.xlsx", sheet_name="required_genes")
    g = pd.read_csv(OUT/"Billmyre_2022_pombe/resampling_aBGC/genes.tsv", sep="\t", low_memory=False); g["ORF"] = g.ORF.str.removeprefix("gene:")
    g = g[(g.status == "tested") & (g.condition == "meiotic_longOutgrowth")]
    t = pd.DataFrame(dict(ORF=g.ORF, call=np.where(g.ORF.isin(set(r.GeneID.astype(str))), "depleted", "none")))
    out.append(("Billmyre 2022 (sex)", "pombe", "depleted", *sets_from_calls(g[["ORF", "call"]], t, "depleted"), "long outgrowth"))
    # Chow 2024: their full list is not deposited; their rule applied to YTAB's p (uncorrected) vs FDR calls
    g = pd.read_csv(OUT/"Chow_2024_glabrata/resampling_aBGC/genes.tsv", sep="\t", low_memory=False); g = g[g.status == "tested"]
    rule = g.assign(call=np.where((g.p_value < .05) & (g.log2FC > 1), "enriched", "none"))
    out.append(("Chow 2024 (FLU)", "glabrata_CBS138", "enriched", *sets_from_calls(g[["ORF", "call"]], rule[["ORF", "call"]], "enriched"),
                "published list not deposited: their rule (uncorrected p < 0.05, log2FC > 1) on YTAB's p"))
    return out


def main():
    ann = {"scer":scer_annotations(), "pombe":pombe_annotations()}
    rows, terms = [], []
    for study, species, direction, tested, sets, note in study_inputs():
        if species == "pombe":a, mapper = ann["pombe"], (lambda x:x)
        else:
            m = to_scer(species); a = ann["scer"]; mapper = m.get
        universe = {mapper(g) for g in tested} - {None}
        for name, genes in sets.items():
            mapped = {mapper(g) for g in genes} - {None}
            e = go.enrich(mapped, universe, a) if len(mapped) >= 5 else pd.DataFrame()
            sig = e[e.q_bh < .05] if not e.empty else e
            rows.append(dict(study=study, direction=direction, set=name, genes=len(genes), mapped=len(mapped), universe=len(universe),
                             enriched_terms=len(sig), top_terms="; ".join(f"{r.name} ({r.k}/{r.K}, q={r.q_bh:.1e})" for r in sig.head(6).itertuples()) if len(sig) else "",
                             note=note))
            if len(sig):terms.append(sig.head(25).assign(study=study, direction=direction, set=name))
    pd.DataFrame(rows).to_csv(OUT/"biology_sets.tsv", sep="\t", index=False)
    if terms:pd.concat(terms).to_csv(OUT/"biology_terms.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 160):
        print(pd.DataFrame(rows)[["study", "direction", "set", "genes", "mapped", "enriched_terms", "top_terms"]].to_string(index=False))


if __name__ == "__main__":main()
