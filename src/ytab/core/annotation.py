"""Genome annotation and gene-to-site mapping, shared by every YTAB analysis that works per gene.

`annotations` reads a GFF3 into each gene's own CDS blocks; `gene_sites` maps a coordinate catalogue onto those
blocks, with the N/C-terminal trimming and stop-codon exclusion the callers ask for. Both are method-agnostic:
the normalization, QC, screen-coverage, fitness and essentiality code all read genes through them."""

import re
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import unquote

import numpy as np


def annotations(path):
    """Resolve CDS parents through transcript ancestry, using each gene's own CDS."""
    records = {}
    ambiguous = set()
    cds = []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("##FASTA"):break
            if line.startswith("#") or not line.strip():continue
            fields=line.rstrip().split("\t")
            if len(fields)!=9:raise ValueError("Expected nine-column GFF3")
            attrs={unquote(k):unquote(v) for k,v in
                   (entry.split("=",1) for entry in fields[8].split(";") if "=" in entry)}
            rec=dict(chrom=fields[0],kind=fields[2],start=int(fields[3]),end=int(fields[4]),
                     strand=fields[6],attrs=attrs)
            if rec["start"]<1 or rec["end"]<rec["start"]:raise ValueError("Invalid inclusive GFF coordinates")
            if fields[2]=="CDS":cds.append(rec)
            elif attrs.get("ID"):
                # A feature split across several rows carries one ID on each (RefSeq does this for a misc_feature
                # interrupted by another annotation). That is ambiguous only if a CDS resolves through it, so the
                # error is deferred to ancestor() rather than rejecting an otherwise valid annotation.
                if attrs["ID"] in records:ambiguous.add(attrs["ID"])
                else:records[attrs["ID"]]=rec
    genes={}
    def ancestor(identifier,seen):
        if identifier in seen:raise ValueError("Cyclic GFF ancestry")
        if identifier in ambiguous:raise ValueError(f"Duplicate non-CDS GFF ID in CDS ancestry: {identifier}")
        if identifier not in records:raise ValueError(f"Unresolved GFF CDS parent: {identifier}")
        record=records[identifier]
        if record["kind"] in {"gene","pseudogene"} or not record["attrs"].get("Parent"):
            return [identifier]
        return [g for p in record["attrs"]["Parent"].split(",") for g in ancestor(p,seen|{identifier})]
    for record in cds:
        parents=record["attrs"].get("Parent","").split(",")
        if not all(parents):raise ValueError("CDS must have explicit GFF Parent")
        for gid in {g for p in parents for g in ancestor(p,set())}:
            parent=records[gid]
            if parent["chrom"]!=record["chrom"] or parent["strand"]!=record["strand"]:
                raise ValueError("CDS and gene disagree on chromosome/strand")
            identifier=parent["attrs"].get("locus_tag",gid)
            gene=genes.setdefault(gid,dict(id=identifier,gff_id=gid,name=parent["attrs"].get("gene") or parent["attrs"].get("Gene") or parent["attrs"].get("Name",identifier),
                description=parent["attrs"].get("description",parent["attrs"].get("Note","")),
                chrom=parent["chrom"],strand=parent["strand"],blocks=[]))
            gene["blocks"].append((record["start"],record["end"]))
    if not genes:raise ValueError("No protein-coding CDS groups found in GFF3")
    if len({g["id"] for g in genes.values()})!=len(genes):raise ValueError("Duplicate gene identifiers")
    for gene in genes.values():
        merged=[]
        for a,b in sorted(set(gene["blocks"])):
            if merged and a<=merged[-1][1]+1:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
            else:merged.append((a,b))
        gene["blocks"]=merged
    return sorted(genes.values(),key=lambda g:g["id"])


def gene_sites(genes, sites, settings):
    """Catalogue indices in each gene's own CDS union (N/C-terminal trimming and last 3 bases excluded, as in TRANSIT's Gumbel method)."""
    by_chrom = {}
    for i, (chrom, p) in enumerate(sites):by_chrom.setdefault(chrom, []).append((p, i))
    chroms = {}
    for chrom, entries in by_chrom.items():
        entries.sort(); positions, indices = np.asarray(entries, dtype=np.int64).T; chroms[chrom] = (positions, indices)
    result = []
    for gene in genes:
        positions, indices = chroms.get(gene["chrom"], (np.array([], dtype=np.int64), np.array([], dtype=np.int64)))
        selected = np.asarray([k for a, b in gene["blocks"]
                               for k in range(np.searchsorted(positions, a), np.searchsorted(positions, b, side="right"))], dtype=int)
        p = positions[selected]; idx = indices[selected]
        lengths = np.array([b-a+1 for a, b in gene["blocks"]]); length = int(lengths.sum())
        offset = np.zeros(len(p), dtype=int); cumulative = 0
        for (a, b), width in zip(gene["blocks"], lengths):
            mask = (p >= a) & (p <= b); offset[mask] = cumulative+p[mask]-a; cumulative += int(width)
        oriented = offset if gene["strand"] != "-" else length-1-offset
        keep = (oriented >= length*settings.nterm/100) & (oriented < length*(1-settings.cterm/100))
        if settings.ignore_stop:keep &= oriented < length-3
        result.append(idx[keep])
    return result


# Default per-gene site rules: the whole CDS, minus the stop codon. `gene_sites` takes any object with these
# three attributes, so a method with its own trimming passes its own Settings instead.
SITE_RULES = SimpleNamespace(nterm=0., cterm=0., ignore_stop=True)


def gene_key(gene_id):
    """Comparable gene ID: annotation prefixes such as 'gene:' (Ensembl/PomBase GFF3) removed."""
    return re.sub(r"^(gene[:-])", "", gene_id)
