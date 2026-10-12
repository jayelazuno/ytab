"""Screen Coverage inputs from the accepted, immutable parent-union audit.

No essentiality classifier is consulted. Representation and the two sequencing
abundance proxies have deliberately separate denominators.
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from .library import INFORMATIVE_DEFINITION, is_parent_sample

BASES = ("sites", "h2o2", "zn")
PROJECTS = {"h2o2": "H2O2_screen_v1", "zn": "Zn_toxicity_screen"}
MODEL_LABELS = {"sites": "Equal-insertion union", "h2o2": "H2O2 abundance sensitivity", "zn": "Zn abundance sensitivity"}


def load_audited_library(audit_dir, selected_sources=None, universe="informative", basis="sites", conditioning="physical"):
    audit_dir = Path(audit_dir)
    paths = [audit_dir / name for name in ("parent_union_audit.json", "parent_union_gene_inventory.csv", "parent_union_site_provenance.csv")]
    audit = json.loads(paths[0].read_text())
    allowed = {f"{r['project']}::{r['sample']}" for r in audit["samples"]}
    if any(not is_parent_sample(r) for r in audit["samples"]):
        raise ValueError("Audit contains a non-parent sample")
    chosen = allowed if selected_sources is None else set(selected_sources)
    if not chosen or not chosen <= allowed:
        raise ValueError("Select parent sources from the accepted audit only")
    if universe not in {"all", "observed", "informative", "alleles2", "alleles3"} or basis not in BASES or conditioning not in {"physical", "conditional"}:
        raise ValueError("Invalid Screen Coverage settings")
    with paths[1].open() as f:
        all_genes = list(csv.DictReader(f))
    genes = {r["gene_id"]: r for r in all_genes if r["biotype"] == "protein_coding"}
    # Lock to the accepted provenance, rather than rediscovering a denominator in the UI.
    if len(genes) != audit["gene_annotation"]["protein_coding_genes"]:
        raise ValueError("Gene inventory disagrees with accepted annotation audit")
    coordinates, gids, reads = [], [], []
    seen = set()
    by_project_sites = {b:set() for b in PROJECTS}
    by_project_genes = {b:set() for b in PROJECTS}
    with paths[2].open() as f:
        for r in csv.DictReader(f):
            coordinate = (r["chromosome"], int(r["position"]))
            if coordinate in seen:
                raise ValueError("Duplicate coordinate in accepted site audit")
            seen.add(coordinate)
            support = json.loads(r["source_read_support_json"])
            if not set(support) <= allowed:
                raise ValueError("Unapproved source in site audit")
            support = {s:float(v) for s,v in support.items() if s in chosen}
            if not support:
                continue
            if any(not np.isfinite(v) or v <= 0 for v in support.values()):
                raise ValueError("Invalid read support in site audit")
            gid = r["unambiguous_CDS_locus"] if r["unambiguous_CDS_locus"] in genes else ""
            values = []
            for b,project in PROJECTS.items():
                value = sum(v for s,v in support.items() if s.split("::")[0] == project)
                values.append(value)
                if value:
                    by_project_sites[b].add(coordinate)
                    if gid:
                        by_project_genes[b].add(gid)
            coordinates.append(coordinate);gids.append(gid);reads.append(values)
    if len(seen) != audit["union"]["unique_sites"]:
        raise ValueError("Site inventory disagrees with accepted union audit")
    reads = np.asarray(reads)
    m = len(coordinates)
    counts = {gid:0 for gid in genes}
    for gid in gids:
        if gid:
            counts[gid] += 1
    # ORF-span representation is used only for the explicitly named observed universe.
    pos = {}
    for chrom in {c for c,_ in coordinates}:
        pos[chrom] = np.array(sorted(p for c,p in coordinates if c==chrom))
    observed = set()
    for gid,g in genes.items():
        pp = pos.get(g["chromosome"],np.array([]))
        if np.searchsorted(pp,int(g["end"]),side="right") > np.searchsorted(pp,int(g["start"])):
            observed.add(gid)
    universes = dict(all=set(genes),observed=observed,informative={g for g,n in counts.items() if n>=1},
                     alleles2={g for g,n in counts.items() if n>=2},alleles3={g for g,n in counts.items() if n>=3})
    if chosen == allowed:
        expected = audit["union"]
        checks = [(m,expected["unique_sites"]),(sum(counts.values()),expected["unambiguous_protein_coding_sites"]),
                  (len(observed),expected["protein_genes_with_observed_ORF_insertion"])]
        checks += [(len(universes[key]),expected[name]) for key,name in
                   [("informative","protein_genes_with_coding_insertion"),("alleles2","protein_genes_with_ge2_sites"),("alleles3","protein_genes_with_ge3_sites")]]
        if any(a!=b for a,b in checks):
            raise ValueError("Derived representation does not match the accepted parent-union audit")
    selected = sorted(universes[universe])
    if not selected:
        raise ValueError("Selected represented gene universe is empty")
    index = {g:i for i,g in enumerate(selected)}
    mask = np.array([gid in index for gid in gids])
    allele_genes = np.array([index[gid] for gid in gids if gid in index],dtype=int)
    # Reads are summed within each experiment only. This is sequencing-depth-
    # weighted evidence, NOT measured mixing proportions or viable frequencies.
    totals = reads.sum(axis=0)
    site_probabilities = {"sites":np.full(m,1/m)}
    for j,b in enumerate(PROJECTS):
        site_probabilities[b] = reads[:,j]/totals[j] if totals[j] else np.zeros(m)
    allele_p, gene_p = {}, {}
    for b in BASES:
        p = site_probabilities[b][mask]
        if conditioning == "conditional" and p.sum():
            p = p/p.sum()
        allele_p[b] = p
        gene_p[b] = np.bincount(allele_genes,weights=p,minlength=len(selected))
    rows = []
    for gid in selected:
        g=genes[gid];i=index[gid]
        rows.append(dict(gene_id=gid,common_name="",chromosome=g["chromosome"],start=int(g["start"]),end=int(g["end"]),
                         length=int(g["coding_length"]),essentiality="Not used",informative_sites=counts[gid],
                         **{f"p_{b}":float(gene_p[b][i]) for b in BASES},p=float(gene_p[basis][i])))
    for j,b in enumerate(PROJECTS):
        values = np.bincount(allele_genes,weights=reads[:,j][mask],minlength=len(selected))
        for row,v in zip(rows,values):
            row[f"read_abundance_{b}"] = float(v)
    h,z = by_project_sites["h2o2"],by_project_sites["zn"]
    hg,zg = by_project_genes["h2o2"],by_project_genes["zn"]
    representation = dict(annotated=len(genes),represented=len(universes["informative"]),ge2=len(universes["alleles2"]),
        ge3=len(universes["alleles3"]),missing=len(genes)-len(universes["informative"]),reference_fraction=len(universes["informative"])/len(genes),
        h2o2_genes=len(hg),zn_genes=len(zg),shared_genes=len(hg&zg),h2o2_only_genes=len(hg-zg),zn_only_genes=len(zg-hg),
        shared_sites=len(h&z),h2o2_only_sites=len(h-z),zn_only_sites=len(z-h))
    return dict(genes=rows,universes={k:len(v) for k,v in universes.items()},representation=representation,
        annotated_genes=len(genes),informative_genes=len(universes["informative"]),observed_genes=len(observed),
        informative_sites=sum(counts.values()),unique_sites=m,selected_informative_sites=len(allele_genes),
        rarefaction_other_sites=m-len(allele_genes),pools=[r for r in audit["samples"] if f"{r['project']}::{r['sample']}" in chosen],
        basis=basis,conditioning=conditioning,parents_only=True,classifier="Not used",informative_definition=INFORMATIVE_DEFINITION,
        sources=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],
        audit_matches=chosen==allowed,abundance_total_reads={b:float(totals[j]) for j,b in enumerate(PROJECTS)},
        _allele_genes=allele_genes,_allele_p=allele_p,
        limitations=[
            "The accepted mapping-reference annotation defines protein-coding genes; its 14 pseudogenes are excluded. Essentiality predictions never define the default denominator.",
            "Parent-union representation is evidence of observed sites, not a physical mixture or viable-colony frequency distribution.",
            "The equal-insertion model assumes every union coordinate is an equally recoverable insertion genotype. Neither equal abundance nor one insertion per clone has been measured.",
            "H2O2 and Zn sequencing counts are summed and normalized separately within each experiment. No cross-experiment averaging is performed. Within-experiment pooling is sequencing-depth-weighted, not measured physical mixing.",
            "H2O2 is the primary abundance sensitivity because its parents were thawed and sequenced; Zn outgrowth is secondary sensitivity. These experimental histories were provided by the user.",
            "Read support is affected by PCR, extraction, mapping, growth and sampling; it is not a viable-colony count. Zero observed reads receive zero model probability without imputation.",
            "CDS insertions are a disruption proxy, not verified knockouts. Ambiguous coding assignments and pseudogene insertions remain non-target sampling mass.",
            "Colony redundancy permits recovery of the same allele; distinct-allele simulations count each genomic coordinate once. Coordinates do not establish clone identity, independence or null function.",
            "Site rarefaction samples a finite catalogue without replacement, whereas colony simulations sample with replacement. They answer different questions."])
