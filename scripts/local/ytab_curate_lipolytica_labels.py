#!/usr/bin/env python3
"""Curate Y. lipolytica essentiality labels that are independent of transposon data.

Why not use the published calls. Patterson et al. 2018 (Metab Eng 48:184-196) classify about 22% of genes essential
from a Poisson analysis of Hermes insertion density - computed on PRJNA356612, **the same 33 sequencing runs R26
maps**. Scoring a YTAB model against them would compare two analyses of one dataset and could say nothing about
accuracy, so they are a comparator at most and never labels.

Evidence used instead:
  silver both labels: the viability of the 1:1 S. cerevisiae ortholog in SGD, from this repository's own curated
         cerevisiae labels (gold tier only).

There is no gold tier for this species: no Y. lipolytica deletion collection or knockout set was obtained.

Ortholog evidence. Unlike C. neoformans and K. phaffii, no published ortholog table was available - Patterson et al.
2018 has no PubMed Central copy. The CLIB89 RefSeq annotation (GCF_001761485.1) instead carries a per-CDS `Note`
of the form "similar to uniprot|<acc> Saccharomyces cerevisiae <systematic name> <gene>", inherited from the
Genolevures annotation. Only the two stronger qualifiers are used, "similar to" and "highly similar to"; "weakly
similar to" and "some similarities with" are excluded. A pair is kept when the Y. lipolytica gene carries exactly
one such name and that name is carried by no other gene.

This is a similarity annotation, not curated orthology, and it is weaker evidence than the ortholog tables used for
the other two species. It is recorded as such in the README and in `sources.json`.

  python scripts/local/ytab_curate_lipolytica_labels.py
"""
import csv
import collections
import hashlib
import json
import re
import urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LABELS_DIR = ROOT/"resources/species/lipolytica/essentiality_labels"
SCER_LABELS = ROOT/"resources/species/cerevisiae/essentiality_labels/labels.tsv"
GFF = ROOT/"resources/species/lipolytica/reference_genome/GCF_001761485.1_ASM176148v1_genomic.gff"
STRONG = ("similar to", "highly similar to")
NOTE = re.compile(r"(highly similar to|similar to|weakly similar to|some similarities with)\s+uniprot\|(\S+)\s+"
                  r"Saccharomyces cerevisiae\s+(Y[A-Z]{2}\d{3}[WC](?:-[A-Z])?)", re.I)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cerevisiae_labels():
    out = {}
    with SCER_LABELS.open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["tier"] == "gold" and r["label"] in ("essential", "nonessential"):out[r["gene_id"]] = r["label"]
    return out


def annotation():
    """({locus tag: gene name} for genes, [(locus tag, systematic name, qualifier)] from the CDS similarity notes)."""
    genes, notes = {}, []
    for line in GFF.open():
        if line.startswith("#"):continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "gene":
            tag = attrs.get("locus_tag") or attrs.get("ID", "").removeprefix("gene-")
            if tag:genes[tag] = attrs.get("gene", "")
        elif f[2] == "CDS" and attrs.get("locus_tag"):
            for m in NOTE.finditer(urllib.parse.unquote(attrs.get("Note", ""))):
                notes.append((attrs["locus_tag"], m.group(3).upper(), m.group(1).lower()))
    return genes, notes


def main():
    scer = cerevisiae_labels()
    genes, notes = annotation()
    strong = [(tag, name) for tag, name, qualifier in notes if qualifier in STRONG]
    per_gene, per_ortholog = collections.Counter(t for t, _ in strong), collections.Counter(n for _, n in strong)
    one_to_one = {t:n for t, n in strong if per_gene[t] == 1 and per_ortholog[n] == 1}
    silver = {t:n for t, n in one_to_one.items() if n in scer}
    qualifiers = {tag:qualifier for tag, _, qualifier in notes}

    rows, counts = [], collections.Counter()
    for gene in sorted(genes):
        ortholog = silver.get(gene, "")
        label = scer[ortholog] if ortholog else ""
        reason = "" if label else (
            "no S. cerevisiae similarity note" if gene not in {t for t, _, _ in notes} else
            "only a weak similarity note" if gene not in {t for t, _ in strong} else
            "similarity is not 1:1" if gene not in one_to_one else
            "the 1:1 match has no curated SGD viability")
        rows.append(dict(gene_id=gene, name=genes[gene], label=label, tier="silver" if label else "",
                         unlabelled_reason=reason, scer_ortholog=ortholog,
                         similarity_qualifier=qualifiers.get(gene, "")))
        counts[("silver" if label else "unlabelled", label or "none")] += 1

    LABELS_DIR.mkdir(parents=True, exist_ok=True)
    with (LABELS_DIR/"labels.tsv").open("w") as handle:
        w = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    summary = dict(
        genes=len(rows), by_tier={f"{t}:{l}":n for (t, l), n in sorted(counts.items())},
        gold_source=None, silver_source="1:1 S. cerevisiae ortholog viability (SGD) through the CLIB89 annotation's "
                                        "'similar to' / 'highly similar to' notes",
        similarity_notes=len(notes), strong_notes=len(strong), one_to_one=len(one_to_one), silver_labels=len(silver),
        excluded_qualifiers=["weakly similar to", "some similarities with"],
        caveat="Similarity annotation, not curated orthology; weaker evidence than the ortholog tables used for "
               "C. neoformans and K. phaffii. No gold tier: no Y. lipolytica knockout set was obtained.")
    (LABELS_DIR/"summary.json").write_text(json.dumps(summary, indent=1)+"\n")
    (LABELS_DIR/"sources.json").write_text(json.dumps(dict(
        reference_gff=dict(path=str(GFF.relative_to(ROOT)), sha256=sha(GFF),
                           note="CLIB89 RefSeq annotation; the S. cerevisiae similarity notes are inherited from the "
                                "Genolevures annotation of Y. lipolytica"),
        cerevisiae_labels=dict(path=str(SCER_LABELS.relative_to(ROOT)), sha256=sha(SCER_LABELS)),
        patterson_2018=dict(citation="Patterson K et al. (2018) Functional genomics for the oleaginous yeast Yarrowia "
                                     "lipolytica. Metab Eng 48:184-196", pmid="29792930",
                            use="comparator only: their essential calls come from Hermes insertion density on "
                                "PRJNA356612, the same runs R26 maps"),
        script=dict(path=str(Path(__file__).relative_to(ROOT)), sha256=sha(__file__))), indent=1)+"\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
