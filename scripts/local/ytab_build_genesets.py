#!/usr/bin/env python3
"""Build a GMT of gene sets from GO annotations, for the comet plot's overlays.

Sets are written on the gene names of the species the annotations describe. The app matches a gene to a set by
its identifier, its own name **or its ortholog's name**, so annotations for a well-studied reference species
cover a sparsely annotated one through its orthologs. That indirection is what makes this usable: C. glabrata's
own GO (CGD) is 45% root "unknown" terms and holds one mitochondrial-ribosome gene, while S. cerevisiae's is
2.5% root and holds 47.

**Annotations are propagated up the ontology.** A GAF records only the term a curator chose, not its ancestors,
so a gene annotated to "ergosterol biosynthetic process" never appears under "steroid biosynthetic process"
unless the is_a and part_of edges are followed. Without that step a set is silently undercounted -- the same
containment problem that makes a GO enrichment report nested restatements of one finding.

  python ytab_build_genesets.py GAF OBO OUT.gmt [--sets sets.tsv]

The set definitions are data, not code: --sets takes a TSV of `set_name<TAB>GO term<TAB>description`. With no
--sets the built-in list below is used, which covers the families a Tn-seq drug screen usually wants.
"""
import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

# name, GO term, description. Chosen for drug and fitness screens; override with --sets.
DEFAULT_SETS = [
    ("ERG",            "GO:0006696", "Ergosterol biosynthetic process"),
    ("PDR_CDR_efflux", "GO:0140359", "ABC-type transmembrane transporter activity"),
    ("ribosome_40S",   "GO:0022627", "Cytosolic small ribosomal subunit"),
    ("ribosome_60S",   "GO:0022625", "Cytosolic large ribosomal subunit"),
    ("mito_ribosome",  "GO:0005761", "Mitochondrial ribosome"),
    ("OXPHOS",         "GO:0006119", "Oxidative phosphorylation"),
    ("autophagy",      "GO:0006914", "Autophagy"),
    ("cell_wall",      "GO:0071555", "Cell wall organization"),
    ("DNA_repair",     "GO:0006281", "DNA repair"),
]

# Deliberately not a default: GO:0005739 (mitochondrion) is a compartment, not a response, and it covers 1,070
# of C. glabrata's 5,234 genes -- 20% of the comet, enough to obscure the body the plot is read against. Every
# family above stays under 7%, and the ones a drug screen actually overlays stay under 2%. Add it with --sets
# if you want it.

# Only these edges carry annotations upward. is_a is subtype; part_of is containment. regulates does NOT:
# a gene regulating a process is not in it, and following those edges sweeps regulators into every set.
ANCESTOR_EDGES = ("is_a", "part_of")


def read_obo(path):
    """term -> set of its parents, over is_a and part_of. Obsolete terms are dropped."""
    parents, alias = defaultdict(set), {}
    term = None
    obsolete = False
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "[Term]":
                term, obsolete = None, False
            elif line.startswith("id: GO:"):
                term = line[4:].strip()
            elif line.startswith("is_obsolete: true"):
                obsolete = True
                if term:
                    parents.pop(term, None)
            elif term and not obsolete:
                if line.startswith("is_a: GO:"):
                    parents[term].add(line[6:].split("!")[0].strip())
                elif line.startswith("relationship: part_of GO:"):
                    parents[term].add(line[len("relationship: part_of "):].split("!")[0].strip())
                elif line.startswith("alt_id: GO:"):
                    alias[line[len("alt_id: "):].strip()] = term
    return parents, alias


def descendants(parents, roots):
    """Every term whose ancestry reaches one of `roots` -- the set a gene must be annotated to, directly or not."""
    children = defaultdict(set)
    for term, ups in parents.items():
        for up in ups:
            children[up].add(term)
    out, stack = set(roots), list(roots)
    while stack:
        for child in children.get(stack.pop(), ()):
            if child not in out:
                out.add(child)
                stack.append(child)
    return out


def read_gaf(path):
    """GO term -> set of gene symbols. NOT qualifiers are skipped: they assert the gene is *not* in the term."""
    by_term = defaultdict(set)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                continue
            symbol, qualifier, term = f[2].strip(), f[3].strip(), f[4].strip()
            if not symbol or not term or "NOT" in qualifier.upper():
                continue
            by_term[term].add(symbol)
    return by_term


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("gaf", type=Path, help="GO annotation file (.gaf or .gaf.gz)")
    p.add_argument("obo", type=Path, help="ontology (go-basic.obo)")
    p.add_argument("output", type=Path)
    p.add_argument("--sets", type=Path, help="TSV: set name, GO term, description")
    p.add_argument("--min-genes", type=int, default=3, help="skip a set with fewer members")
    a = p.parse_args(argv)

    if a.sets:
        rows = []
        for line in a.sets.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = [x.strip() for x in line.split("\t")]
            if len(parts) < 2:
                raise SystemExit(f"--sets needs name and GO term per line: {line[:60]!r}")
            rows.append((parts[0], parts[1], parts[2] if len(parts) > 2 else parts[0]))
    else:
        rows = DEFAULT_SETS

    parents, alias = read_obo(a.obo)
    by_term = read_gaf(a.gaf)
    print(f"ontology: {len(parents)} terms   annotations: {sum(len(v) for v in by_term.values())}")

    lines = [
        f"# Gene sets from {a.gaf.name} over {a.obo.name}, one line per set: name, description, members.",
        "# Members are every gene annotated to the term OR to any of its descendants through is_a and part_of,",
        "# because a GAF records only the term a curator chose. Regenerate with scripts/local/ytab_build_genesets.py.",
    ]
    report = []
    for name, term, description in rows:
        term = alias.get(term, term)
        wanted = descendants(parents, {term})
        members = sorted(set().union(*(by_term.get(t, set()) for t in wanted)) if wanted else set())
        direct = len(by_term.get(term, set()))
        if len(members) < a.min_genes:
            report.append(f"  {name:16s} {len(members):4d}  (skipped, fewer than {a.min_genes})")
            continue
        lines.append("\t".join([name, f"{description} ({term})", *members]))
        report.append(f"  {name:16s} {len(members):4d} genes  ({direct} annotated to {term} itself)")

    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text("\n".join(lines) + "\n")
    print(f"{a.output}  ({len(lines) - 3} sets)")
    print("\n".join(report))


if __name__ == "__main__":
    raise SystemExit(main())
