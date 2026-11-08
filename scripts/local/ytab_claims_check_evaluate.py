#!/usr/bin/env python3
"""Claims check, step 3: each published claim (output/claims_check/claims.tsv) against YTAB's calls on the same data.

claims.tsv columns: study, project, conditions (regex on the condition names; ".*" for all), gene (published name),
direction (enriched | depleted | any), published_rule (how the paper called it), validated (yes if confirmed outside Tn-seq,
e.g. a deletion mutant), source (where in the paper). Gene names are resolved through the run's own names and
descriptions, the species annotation (gene_name, scer_gene_name) and, for C. glabrata, the S. cerevisiae ortholog
table. YTAB runs: output/claims_check/<project>/ (step 1) or, for Gale 2020 and Nickels 2025, the method-comparison
report runs (aBGC, all sites).

Per claim and condition: resampling (aBGC) call, q, log2FC and the gene's sites and control/treated reads; z-score
(aBGC) call and z where it ran. Verdict:
  supported            resampling calls it in the published direction (q < 0.05)
  direction only       same direction, not significant: the data point the same way without enough evidence
  data-limited         fewer than MIN_SITES sites, or no control reads: no test can support it on these libraries
  opposite             significant in the other direction
  not found            the gene is not in the run (name not resolved, or untested)
Output: output/claims_check/claims_evaluated.tsv and a per-study summary.
"""
from __future__ import annotations
import re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]; OUT = ROOT/"output/claims_check"


def config(project):
    """Project config: the project's own, or a study-local copy (output/claims_check/projects/; e.g. libraries left out)."""
    own = ROOT/"output/projects"/project/"config/project.yaml"
    return own if own.exists() else OUT/"projects"/project/"config/project.yaml"
MIN_SITES = 3
REPORT = ROOT/"output/comparison/report"
RUNS = {"Gale_2020": REPORT/"Gale_2020/none", "Nickels_2025_glabrata": REPORT/"Nickels_2025_glabrata/none"}
ANNOTATION = {"glabrata_CBS138": ROOT/"resources/species/glabrata/strains/CBS138/annotation/gene_annotation.tsv",
              "cerevisiae": ROOT/"resources/species/cerevisiae/annotation/gene_annotation.tsv",
              "pombe": ROOT/"resources/species/pombe/annotation/gene_annotation.tsv",
              "albicans": ROOT/"resources/species/albicans/annotation/gene_annotation.tsv"}
ORTHOLOGS = {"glabrata_CBS138": ROOT/"resources/species/glabrata/strains/CBS138/orthologs/C_glabrata_CBS138_S_cerevisiae_orthologs.txt",
             "glabrata": ROOT/"resources/species/glabrata/orthologs/C_glabrata_BG2_S_cerevisiae_orthologs.txt"}


def species(project):
    import yaml
    cfg = yaml.safe_load(open(config(project)))
    return str(cfg.get("species") or cfg.get("reference", {}).get("species") or "")


def aliases(project, genes):
    """ORF -> set of lower-case names it may be published under."""
    out = {o:{str(o).lower(), str(n).lower()} for o, n in zip(genes.ORF, genes.Name)}
    for o, d in zip(genes.ORF, genes.get("Description", pd.Series([""]*len(genes)))):
        for w in re.findall(r"\b[A-Z][A-Z0-9]{2,7}\b", str(d)):out[o].add(w.lower())
    sp = species(project)
    for o, extra in gff_names(project).items():   # GFF names and old locus tags (papers cite either)
        if o in out:out[o] |= extra
    for key, path in ANNOTATION.items():
        if key in sp and path.exists():
            a = pd.read_csv(path, sep="\t", dtype=str).fillna("")
            for _, r in a.iterrows():
                for o in (r.get("gene_id"), r.get("systematic_id")):
                    if o in out:out[o] |= {x.lower() for x in (r.get("gene_name", ""), r.get("scer_gene_name", "")) if x}
    for key, path in ORTHOLOGS.items():
        if key in sp and path.exists():
            o = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str).fillna("")
            for _, r in o.iterrows():
                for col in (0, 1):
                    if r[col] in out and r[4]:out[r[col]].add(r[4].lower())
    return out


def gff_names(project):
    """locus tag / ID -> {gene name, old locus tag and its numeric suffix} from the project's GFF."""
    import yaml
    cfg = yaml.safe_load(open(config(project))); gff = Path(cfg["reference"]["gff"])
    gff = gff if gff.is_absolute() else ROOT/gff; out = {}
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        names = {a.get(k, "").lower() for k in ("Name", "gene", "old_locus_tag", "locus_tag") if a.get(k)}
        names |= {n.split("_")[-1] for n in names if "_" in n}
        for k in ("locus_tag", "ID", "Name"):
            if a.get(k):out.setdefault(a[k].removeprefix("gene-"), set()).update(names)
    return out


def run_dir(project):
    return RUNS.get(project, OUT/project)


def load(project):
    base = run_dir(project); res = pd.read_csv(base/"resampling_aBGC/genes.tsv", sep="\t", low_memory=False)
    res["ORF"] = res.ORF.astype(str).str.removeprefix("gene:")
    z = None; zt = base/"zscore_aBGC/tables/treated_vs_parent.combined.log2fc_z.csv"
    zt = zt if zt.exists() else zt.with_suffix(".csv.gz")
    if zt.exists():
        z = pd.read_csv(zt).rename(columns={"feature_id":"ORF", "treated_condition":"condition"})[["ORF", "background", "condition", "z", "call"]]
    return res, z


def design(project):
    """(background, condition) -> (control libraries, treated libraries) from the comparison design."""
    import yaml
    cfg = yaml.safe_load(open(config(project))); base = Path(cfg.get("output_project_dir") or f"output/projects/{project}")
    base = base if base.is_absolute() else ROOT/base
    d = pd.read_csv(base/"config/comparison_design.csv", dtype=str)
    d = d[d.include.str.strip().str.lower().isin(["1", "true", "yes", "y"])]
    return {k:(g.parent_sample.nunique(), g.treated_sample.nunique()) for k, g in d.groupby(["background", "treated_condition"])}


def min_p(sites, controls, treated, permutations=10000):
    """Smallest two-sided p-value site-restricted resampling can reach: every site's labels can be assigned in
    C(J+K, K) ways, so no outcome is rarer than 1 / prod C(J+K, K); Monte Carlo cannot go below 1/(permutations+1)."""
    from math import comb, log10
    if not sites or not controls or not treated:return 1.0
    exact = -sites*log10(comb(controls+treated, treated))
    return min(1.0, max(2*10**exact, 1/(permutations+1)))


def verdict(r, direction):
    if r is None:return "not found"
    if r["call"] == direction or (direction == "any" and r["call"] in ("enriched", "depleted")):return "supported"
    if r["call"] in ("enriched", "depleted"):return "opposite"
    if r["n_sites"] < MIN_SITES or r["sum_control"] <= 0:return "data-limited"
    if direction == "any":return "not supported"
    same = (r["log2FC"] > 0) == (direction == "enriched")
    return "direction only" if same else "not supported"


def main():
    claims = pd.read_csv(OUT/"claims.tsv", sep="\t", dtype=str).fillna("")
    rows, cache = [], {}
    for c in claims.itertuples():
        if c.project not in cache:
            res, z = load(c.project); cache[c.project] = (res, z, aliases(c.project, res.drop_duplicates("ORF")))
        res, z, names = cache[c.project]
        orfs = {o for o, a in names.items() if c.gene.lower() in a}
        sub = res[res.condition.str.fullmatch(c.conditions) & res.ORF.isin(orfs)]
        conds = sorted(res.condition[res.condition.str.fullmatch(c.conditions)].unique())
        libs = design(c.project)
        for cond in conds:
            hit = sub[sub.condition == cond]
            r = hit.sort_values("q_bh").iloc[0] if len(hit) else None
            called = res[(res.condition == cond) & res.call.isin(["enriched", "depleted"])]
            cutoff = float(called.p_value.max()) if len(called) else 0.   # the p-value the realized FDR cut needed
            j, k = libs.get((r["background"], cond), (0, 0)) if r is not None else (0, 0)
            floor = min_p(int(r["n_sites"] or 0), j, k) if r is not None else None
            powered = None if r is None else (floor <= cutoff)
            zr = None if z is None or r is None else z[(z.condition == cond) & (z.ORF == r["ORF"])]
            rows.append(dict(study=c.study, project=c.project, condition=cond, gene=c.gene, direction=c.direction, validated=c.validated,
                             published_rule=c.published_rule, ORF=None if r is None else r["ORF"],
                             sites=None if r is None else r["n_sites"], control_reads=None if r is None else round(r["sum_control"], 1),
                             treated_reads=None if r is None else round(r["sum_treatment"], 1), log2FC=None if r is None else round(r["log2FC"], 2),
                             q=None if r is None else r["q_bh"], call=None if r is None else r["call"],
                             zscore_call=None if zr is None or zr.empty else zr.call.iloc[0], z=None if zr is None or zr.empty else round(zr.z.iloc[0], 2),
                             libraries=f"{j} vs {k}", min_p=None if floor is None else float(f"{floor:.2g}"), fdr_cutoff_p=float(f"{cutoff:.2g}"),
                             powered=powered, verdict=verdict(r, c.direction) if powered is not False or r["call"] == c.direction else "not powered",
                             source=c.source))
    e = pd.DataFrame(rows); e.to_csv(OUT/"claims_evaluated.tsv", sep="\t", index=False)
    s = e.groupby(["study", "verdict"]).size().unstack(fill_value=0); s.to_csv(OUT/"claims_summary.tsv", sep="\t")
    print(s.to_string())


if __name__ == "__main__":main()
