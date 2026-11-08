#!/usr/bin/env python3
"""Claims check, essentiality: each study's published essential / nonessential calls, and YTAB's label-free HMM on the
same libraries, both scored against independent evidence (the curated labels: deletion or repressible-allele evidence
in the species itself = gold; S. cerevisiae null-mutant viability through 1:1 orthologs = silver; transposon calls
never label). Untracked study; writes output/claims_check/essentiality_*.tsv.

Per study and caller, on the genes both callers scored:
  ES called, of which contradicted (independent label nonessential) per tier;
  NE called, of which contradicted (independent label essential) per tier;
  agreement of the published calls with the HMM (Cohen's kappa) and library density (sites per coding kb).
Published calls come from the label tables' comparator columns (Gale 2020: CGD; Segal 2018: CGD; Billmyre 2025: S1
Table), which the label curation keeps out of the labels themselves.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]; OUT = ROOT/"output/claims_check"
SP = ROOT/"resources/species"
STUDIES = [   # study, labels, gene-id column, published column, HMM genes.tsv, HMM id column in labels, density (sites/coding kb), density source
    ("Gale 2020 (C. glabrata, Hermes)", SP/"glabrata/essentiality_labels/labels.tsv", "cglab_id", "gale2020_comparator",
     ROOT/"output/projects/Gale_2020/essentiality/essentiality_20261025_131741/models/genes.tsv", 20.8, "S7 Table (R16)"),
    ("Segal 2018 (C. albicans, Ac/Ds)", SP/"albicans/essentiality_labels/labels.tsv", "gene_id", "transposon_comparator",
     ROOT/"output/projects/Segal_2018_albicans/essentiality/essentiality_20261025_131741/models/genes.tsv", 14.7, "S7 Table (R16)"),
    ("Billmyre 2025 (C. neoformans, Ac/Ds)", SP/"neoformans/essentiality_labels/labels.tsv", "gene_id", "billmyre2025_comparator",
     OUT/"Billmyre_2025_neoformans_cc/essentiality/models/genes.tsv", 19.8, "R26 (one of two parents; the other lacks a hit file)"),
]


def norm_call(x):
    x = str(x).strip().lower()
    return "ES" if x in ("essential", "es", "ess") else "NE" if x in ("nonessential", "ne", "ness", "non-essential") else None


def kappa(a, b):
    t = pd.crosstab(a, b); n = t.values.sum()
    if n == 0 or t.shape != (2, 2):return float("nan")
    po = t.values.trace()/n; pe = sum(t.sum(axis=1)[k]*t.sum(axis=0)[k] for k in t.index if k in t.columns)/n**2
    return (po-pe)/(1-pe)


def score(calls, labels):
    rows = {}
    for c in ("ES", "NE"):
        sub = labels[calls == c]; opposite = "nonessential" if c == "ES" else "essential"
        rows[f"{c}_called"] = int((calls == c).sum())
        for tier in ("gold", "silver"):
            lab = sub[sub.tier == tier]
            rows[f"{c}_{tier}_labelled"] = len(lab); rows[f"{c}_{tier}_contradicted"] = int((lab.label == opposite).sum())
    return rows


def main():
    out = []
    for study, lab_path, idcol, pubcol, hmm_path, density, dsrc in STUDIES:
        lab = pd.read_csv(lab_path, sep="\t", dtype=str)
        lab["published"] = lab[pubcol].map(norm_call)
        if hmm_path.exists():
            h = pd.read_csv(hmm_path, sep="\t", dtype=str)
            h["hmm"] = h.hmm_call.map(lambda x:"ES" if x == "ES" else "NE" if x in ("NE", "GA") else None)   # GD (growth defect) left out, as in R16
            lab = lab.merge(h[["ORF", "hmm"]], left_on=idcol, right_on="ORF", how="left")
        else:
            lab["hmm"] = None
        both = lab[lab.published.notna() & lab.hmm.notna()]
        for caller in ("published", "hmm"):
            out.append(dict(study=study, caller=caller, density=density, density_source=dsrc, genes_scored_by_both=len(both),
                            **score(both[caller], both), kappa_published_vs_hmm=round(kappa(both.published, both.hmm), 3) if len(both) else None))
    t = pd.DataFrame(out)
    for c in ("ES", "NE"):
        for tier in ("gold", "silver"):
            t[f"{c}_{tier}_rate"] = (t[f"{c}_{tier}_contradicted"]/t[f"{c}_{tier}_labelled"]).round(3)
    t.to_csv(OUT/"essentiality_summary.tsv", sep="\t", index=False)
    cols = ["study", "caller", "density", "genes_scored_by_both", "ES_called", "ES_gold_labelled", "ES_gold_contradicted", "ES_silver_labelled",
            "ES_silver_contradicted", "NE_called", "NE_gold_labelled", "NE_gold_contradicted", "NE_silver_contradicted", "kappa_published_vs_hmm"]
    print(t[cols].to_string(index=False))


if __name__ == "__main__":main()
