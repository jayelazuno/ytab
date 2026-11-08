#!/usr/bin/env python3
"""Curate C. neoformans essentiality labels that are independent of transposon data.

Why not use the published calls. Billmyre et al. 2025 (PLoS Biol 23:e3003184) classify every gene ESS/NESS/UNK, but
their own S1 Table description states the classification is "the mean from 0-1 of essentiality prediction score from
100 random forest model runs" fitted to insertion features of the same libraries YTAB analyses. Scoring a YTAB model
against them would measure agreement between two pipelines on shared data. They are therefore kept as
`billmyre2025_comparator` and are never labels - the rule the C. glabrata labels already follow for Gale et al. 2020
(resources/species/glabrata/essentiality_labels/README.md).

Evidence used instead, in order of precedence:
  gold   nonessential: a viable deletion strain exists in the Madhani 2024 deletion collection (S1 Table column
         "In Madhani Deletion Collection" == True). Direct C. neoformans evidence, not derived from transposons.
         Absence from the collection is not evidence of essentiality - a gene may simply not have been attempted -
         so False yields no label.
  silver both labels: the viability of the 1:1 S. cerevisiae ortholog in SGD, taken from this repository's own
         curated cerevisiae labels (gold tier only). Orthologs come from Billmyre's S4 Table, which lists one
         S. cerevisiae accession per gene; a pair is 1:1 when the C. neoformans gene has exactly one accession and
         that accession is listed for no other C. neoformans gene. SGD accessions are resolved to systematic names
         through SGD_features.tab.

Where the two disagree the gold label wins and the silver one is recorded in `conflict`, as for C. glabrata.

  python scripts/local/ytab_curate_neoformans_labels.py                 # writes labels.tsv, sources.json, summary.json

Outputs go to resources/species/neoformans/essentiality_labels/. The raw supplements must already be in raw/ (they
are downloaded once by hand and kept with their hashes, as the other species' raw inputs are).
"""
import csv
import collections
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LABELS_DIR = ROOT/"resources/species/neoformans/essentiality_labels"
RAW = LABELS_DIR/"raw/billmyre_2025_plosbiol"
S1 = RAW/"S1 Table essentiality predictions and fluconazole TN-seq.xlsx"
S4 = RAW/"S4 Table cross-species essentiality and orthologs.xlsx"
SCER_LABELS = ROOT/"resources/species/cerevisiae/essentiality_labels/labels.tsv"
SGD_FEATURES = ROOT/"resources/species/cerevisiae/reference_genome/SGD_features.tab"
GFF = ROOT/"resources/species/neoformans/reference_genome/GCF_000149245.1_CNA3_genomic.gff"


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sheet(path, name):
    import openpyxl
    ws = openpyxl.load_workbook(path, read_only=True)[name]
    rows = ws.iter_rows(values_only=True); header = list(next(rows))
    return [dict(zip(header, r)) for r in rows]


def cerevisiae_labels():
    """{systematic name: essential|nonessential} from this repository's curated gold cerevisiae labels."""
    out = {}
    with SCER_LABELS.open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["tier"] == "gold" and r["label"] in ("essential", "nonessential"):out[r["gene_id"]] = r["label"]
    return out


def sgd_accessions():
    """{SGD accession: systematic name}, from the SGD feature table already in the repository."""
    out = {}
    for row in csv.reader(SGD_FEATURES.open(), delimiter="\t"):
        if len(row) > 3 and row[0].startswith("S0") and row[3]:out[row[0]] = row[3]
    return out


def annotated_genes():
    """{locus tag: gene name} for the protein-coding genes of the reference the libraries are mapped to."""
    names = {}
    for line in GFF.open():
        if line.startswith("#"):continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        tag = attrs.get("locus_tag") or attrs.get("ID", "").removeprefix("gene-")
        if tag:names[tag] = attrs.get("gene", "")
    return names


def main():
    s1, s4 = sheet(S1, "Supplemental Table 1"), sheet(S4, "Supplemental Table 4")
    scer, accession, annotated = cerevisiae_labels(), sgd_accessions(), annotated_genes()

    # silver: viability of the 1:1 S. cerevisiae ortholog
    pairs = [(r["Gene"], r["S. cerevisiae ortholog"]) for r in s4 if r["Gene"] and r["S. cerevisiae ortholog"]]
    per_gene, per_ortholog = collections.Counter(g for g, _ in pairs), collections.Counter(o for _, o in pairs)
    one_to_one = {g:o for g, o in pairs if per_gene[g] == 1 and per_ortholog[o] == 1}
    silver = {g:(scer[accession[o]], accession[o]) for g, o in one_to_one.items()
              if o in accession and accession[o] in scer}

    # gold: a viable deletion strain exists
    madhani = {r["Gene"]:str(r["In Madhani Deletion Collection"]) == "True" for r in s1}
    comparator = {r["Gene"]:{"ESS":"essential", "NESS":"nonessential"}.get(r["Essentiality Classification"], "")
                  for r in s1}

    rows, counts = [], collections.Counter()
    for gene in sorted(set(s1_gene for s1_gene in madhani) | set(annotated)):
        gold = "nonessential" if madhani.get(gene) else ""
        silver_label, ortholog = silver.get(gene, ("", ""))
        label, tier = (gold, "gold") if gold else (silver_label, "silver" if silver_label else "")
        conflict = silver_label if gold and silver_label and silver_label != gold else ""
        reason = "" if label else (
            "no 1:1 S. cerevisiae ortholog in Billmyre S4 Table" if gene not in one_to_one else
            "1:1 ortholog has no curated SGD viability")
        rows.append(dict(gene_id=gene, name=annotated.get(gene, ""), label=label, tier=tier,
                         unlabelled_reason=reason, scer_ortholog=ortholog, conflict=conflict,
                         in_madhani_collection=str(madhani.get(gene, "")).lower(),
                         billmyre2025_comparator=comparator.get(gene, ""),
                         in_reference_annotation=str(gene in annotated).lower()))
        counts[(tier or "unlabelled", label or "none")] += 1

    LABELS_DIR.mkdir(parents=True, exist_ok=True)
    with (LABELS_DIR/"labels.tsv").open("w") as handle:
        w = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    conflicts = [r for r in rows if r["conflict"]]
    agreement = collections.Counter((r["label"], r["billmyre2025_comparator"]) for r in rows
                                    if r["label"] and r["billmyre2025_comparator"])
    agreed = sum(n for (a, b), n in agreement.items() if a == b)
    summary = dict(
        genes=len(rows), reference_annotation=len(annotated),
        by_tier={f"{t}:{l}":n for (t, l), n in sorted(counts.items())},
        silver_source="1:1 S. cerevisiae ortholog viability (SGD, via this repository's curated cerevisiae labels)",
        gold_source="viable deletion strain in the Madhani 2024 collection (Billmyre S1 Table)",
        one_to_one_orthologs=len(one_to_one), silver_labels=len(silver),
        conflicts_gold_over_silver=len(conflicts),
        comparator_agreement=dict(labelled_and_called=sum(agreement.values()), agree=agreed,
                                  fraction=agreed/sum(agreement.values()) if agreement else None,
                                  table={f"{a} vs {b}":n for (a, b), n in sorted(agreement.items())}))
    (LABELS_DIR/"summary.json").write_text(json.dumps(summary, indent=1)+"\n")
    (LABELS_DIR/"sources.json").write_text(json.dumps(dict(
        billmyre_2025=dict(citation="Billmyre RB et al. (2025) Landscape of essential growth and fluconazole-resistance "
                                    "genes in the human fungal pathogen Cryptococcus neoformans. PLoS Biol 23:e3003184",
                           doi="10.1371/journal.pbio.3003184",
                           files={p.name:sha(p) for p in sorted(RAW.glob("*.xlsx"))}),
        cerevisiae_labels=dict(path=str(SCER_LABELS.relative_to(ROOT)), sha256=sha(SCER_LABELS)),
        sgd_features=dict(path=str(SGD_FEATURES.relative_to(ROOT)), sha256=sha(SGD_FEATURES)),
        reference_gff=dict(path=str(GFF.relative_to(ROOT)), sha256=sha(GFF)),
        script=dict(path=str(Path(__file__).relative_to(ROOT)), sha256=sha(__file__))), indent=1)+"\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
