#!/usr/bin/env python3
"""One gene across species: each species' insertion profile over the gene, from a representative library set.

This answers the species-level question ("is this gene essential in these yeasts?") rather than the data-level one
("what do my libraries say about it?"), so each species contributes one library set, not all of them:

  representative   among a species' essentiality runs, the one whose libraries reach the highest pooled insertion
                   density (sites per coding kb, median gene), provided it reaches the density an essentiality call
                   needs at all. Density is the property that decides whether a call means anything, and the
                   deepest set is usually the unselected wild-type pool, so the rule needs no list of strain names.
                   A species with nothing dense enough is reported as such instead of being given a call.
  search           a union over every represented species' own identifiers, names and aliases, not a route
                   through one of them: a gene present in a species and named exactly as asked for is found
                   there even when no ortholog table connects it to S. cerevisiae.
  orthologs        through the S. cerevisiae hub, whole families rather than one-to-one pairs: where a family has
                   expanded in a species, every member is profiled, because the expansion is usually the point of
                   the comparison. A species with no counterpart is reported as having none, not dropped. Each
                   gene records whether it was reached as an ortholog of the query or only by carrying its name.
  profile          reads at each insertion site inside the gene's coding blocks plus flanks, summed over the
                   representative's libraries, with each site placed by its fraction along the coding sequence so
                   that genes of different length and orientation can be drawn on one axis.

  ytab gene-profiles --gene ATH1 --output output/gene_profiles/ATH1
  ytab gene-profiles --gene CAGL0K05137g --output DIR --min-density 18
"""
import argparse
import csv
import json
from pathlib import Path
from urllib.parse import unquote
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
sys.path.append(str(ROOT/"scripts/local"))
from ytab.essentiality.control_conditioned_gumbel import annotations          # noqa: E402
from ytab.qc.library_qc import HMM_SITES_PER_CODING_KB, HYBRID_SITES_PER_CODING_KB  # noqa: E402
import ytab_cross_species as xs                                               # noqa: E402

FLANK = .25


def projects():
    """Every project configuration, with its species."""
    for config in sorted((ROOT/"output/projects").glob("*/config/project.yaml")):
        cfg = yaml.safe_load(config.read_text())
        if cfg.get("project_id"):
            yield config, cfg


def readiness(project_root):
    """Pooled density per strain, from the newest QC run (`ytab qc` writes it)."""
    files = sorted((project_root/"qc").glob("*/essentiality_readiness.tsv"),
                   key=lambda p: p.stat().st_mtime, reverse=True)
    if not files:
        return []
    with files[0].open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def candidates():
    """Every essentiality run that could represent its species, with the density of the libraries it was fitted on."""
    out = []
    for config, cfg in projects():
        project_root = ROOT/"output/projects"/cfg["project_id"]
        rows = readiness(project_root)
        for manifest_path in (project_root/"essentiality").glob("*/manifest.json"):
            run = manifest_path.parent
            if not (run/"models/genes.tsv").is_file():
                continue
            manifest = json.loads(manifest_path.read_text())
            libraries = [str(x) for x in (manifest.get("libraries") or [])]
            match = next((r for r in rows if set((r.get("libraries") or "").split()) & set(libraries)), None)
            density = float(match["sites_per_coding_kb"]) if match and match.get("sites_per_coding_kb") else None
            out.append(dict(species=cfg.get("species", ""), project=cfg["project_id"], project_config=config,
                            run=run, libraries=libraries, background=(match or {}).get("background", ""),
                            sites_per_coding_kb=density))
    return out


def densest_libraries():
    """The densest unselected library set each species has, whether or not an essentiality run was ever fitted.

    A species whose libraries are all too sparse never gets a run, so looking only at runs would leave it out of
    the comparison altogether — and silence reads as "no data" when the honest answer is "these libraries cannot
    answer this". The QC readiness table covers every project, so it is the one that decides who appears.
    """
    best = {}
    for _, cfg in projects():
        species = cfg.get("species", "")
        if not species:
            continue
        for row in readiness(ROOT/"output/projects"/cfg["project_id"]):
            value = row.get("sites_per_coding_kb")
            if not value:
                continue
            density = float(value)
            if species not in best or density > best[species][1]:
                best[species] = (cfg["project_id"], density, row.get("background", ""))
    return best


def representatives(min_density):
    """One run per species: the densest, and only if it reaches the density a call needs."""
    best, rejected = {}, {}
    for entry in candidates():
        species = entry["species"]
        if not species:
            continue
        density = entry["sites_per_coding_kb"]
        if density is None:
            continue
        if density >= min_density:
            if species not in best or density > best[species]["sites_per_coding_kb"]:
                best[species] = entry
        elif species not in rejected or density > rejected[species]["sites_per_coding_kb"]:
            rejected[species] = entry
    for species, (project, density, background) in densest_libraries().items():
        if species in best or species in rejected:
            continue
        rejected[species] = dict(species=species, project=project, run=None, libraries=[],
                                 background=background, sites_per_coding_kb=density)
    for species in list(rejected):
        if species in best:
            rejected.pop(species)
    return best, rejected


def ortholog_families(project_config):
    """({hub gene: [this species' genes]}, {this species' gene: [hub genes]}) from the project's ortholog table.

    Whole families, not one-to-one pairs. A gene family that expanded in one species has several counterparts
    for one S. cerevisiae gene, and the expansion is usually the interesting part of the comparison: dropping
    those genes would report the species as having no ortholog when what it has is more than one. Both empty
    when the project has no ortholog table, which is how the hub species itself is recognised.
    """
    cfg = yaml.safe_load(Path(project_config).read_text())
    path = (cfg.get("reference") or {}).get("orthology_file")
    if not path:
        return {}, {}
    path = Path(path)
    path = path if path.is_absolute() else ROOT/path
    if not path.is_file():
        return {}, {}
    families, forward = {}, {}
    for gene, hub in xs.read_ortholog_pairs(path):
        families.setdefault(hub, []).append(gene)
        forward.setdefault(gene, []).append(hub)
    return ({hub: sorted(set(genes)) for hub, genes in families.items()},
            {gene: sorted(set(hubs)) for gene, hubs in forward.items()})


ORTHOGROUP_DIR = ROOT/"resources/comparative/orthology"
ORTHOGROUP_TABLE = ORTHOGROUP_DIR/"orthogroups_long.csv"

def canonical_gene_id(value):
    """A gene id as the rest of the toolkit spells it: without the GFF's `gene:` or `gene-` decoration.

    PomBase writes a gene's ID as `gene:SPBC32H8.12c` while its hit files, calls and annotation blocks all use
    the bare id. Indexing the decorated form made a second, unresolvable copy of the gene: it matched a query,
    was offered as a counterpart, then failed to find its own coding blocks and drew an empty panel beside the
    real one.
    """
    return str(value).split(":", 1)[-1].replace("gene-", "").strip()


_SEARCH_INDEX = {}
_SYNONYMS = None
_ORTHOGROUPS = None


def id_synonyms():
    """{lowercased identifier: other identifiers for the same gene}, from the id maps beside the orthogroup table.

    Resources name the same gene in different ways - a C. glabrata gene is GWK60_L07469 in the BG2 annotation and
    CAGL0L07766g in the orthogroup table - so one has to be translated into the other before the two can be
    joined. Rather than a translation per species, every mapping file in the directory is read as rows of
    interchangeable identifiers, whatever its columns are called: a species whose annotation already uses the
    orthogroup table's names needs no file, and one that does not is handled by dropping a map in beside these.
    """
    global _SYNONYMS
    if _SYNONYMS is not None:
        return _SYNONYMS
    _SYNONYMS = {}
    for path in sorted(ORTHOGROUP_DIR.glob("*.csv")):
        if path == ORTHOGROUP_TABLE:
            continue
        with path.open() as handle:
            for row in csv.DictReader(handle):
                # Identifiers only: a descriptive column holds prose, which would make unrelated genes synonymous.
                ids = [str(v).strip() for v in row.values()
                       if v and len(str(v).strip()) > 3 and " " not in str(v).strip()]
                ids = list(dict.fromkeys(ids))
                for one in ids:
                    _SYNONYMS.setdefault(one.lower(), set()).update(x for x in ids if x != one)
    return _SYNONYMS


def orthogroups():
    """({hub gene: orthogroup ids}, {orthogroup id: [(species, gene id, gene name)]}) from the OrthoFinder table.

    An orthogroup is a gene family, not an ortholog assignment: two genes sharing one means they are related, not
    that either is the other's counterpart. It is read anyway because a curated pairwise table often cannot say a
    family exists at all - CGD's C. albicans table is one-to-one by construction, so the 824 S. cerevisiae genes
    with more than one albicans relative are invisible in it - and because this is additive, a gap in the
    orthogroup table costs an extra edge rather than a wrong answer.
    """
    global _ORTHOGROUPS
    if _ORTHOGROUPS is not None:
        return _ORTHOGROUPS
    by_hub, members = {}, {}
    if not ORTHOGROUP_TABLE.is_file():
        _ORTHOGROUPS = (by_hub, members)
        return _ORTHOGROUPS
    with ORTHOGROUP_TABLE.open() as handle:
        for row in csv.DictReader(handle):
            group = str(row.get("orthogroup_id") or "").strip()
            gene = str(row.get("gene_id") or "").strip()
            if not group or not gene:
                continue
            name = str(row.get("gene_name") or "").strip()
            members.setdefault(group, []).append((str(row.get("species") or "").strip(), gene, name))
            # The hub is keyed on whatever the table calls the S. cerevisiae gene, matched later against the hub
            # id the ortholog tables use; a row whose species is blank is still indexed, because it is a member
            # of the family even when the table forgot to label it.
            for token in (gene, name):
                if token:
                    by_hub.setdefault(token.lower(), set()).add(group)
    _ORTHOGROUPS = (by_hub, members)
    return _ORTHOGROUPS


def orthogroup_counterparts(project_config, species, hubs):
    """This species' own gene ids that share an orthogroup with any of the hub genes.

    Only the rows the table attributes to this species, plus the rows it leaves unattributed: gene names repeat
    across species, so resolving every member against one species' annotation would let another species' HXT1
    match this one's. An unattributed row is still tried, because it is a member of the family whose species
    label is merely missing, and it is kept only if this annotation recognises it.

    Each member is translated into the identifier this species' annotation actually uses, through that
    annotation's own index of ids, names and aliases, widened by the shared id maps. A member that cannot be
    translated is skipped rather than guessed at. A species the table does not cover simply yields nothing,
    which is how the species with no orthogroup data fall back to the other routes.
    """
    by_hub, members = orthogroups()
    if not by_hub:
        return set()
    groups = set()
    for hub in hubs:
        groups.update(by_hub.get(str(hub).lower(), ()))
    if not groups:
        return set()
    wanted = str(species).strip().lower()
    index, synonyms, found = search_index(project_config), id_synonyms(), set()
    for group in groups:
        for row_species, gene, name in members.get(group, ()):
            if row_species and row_species.strip().lower() != wanted:
                continue
            tokens = [gene, name]
            tokens.extend(synonyms.get(gene.lower(), ()))
            if name:
                tokens.extend(synonyms.get(name.lower(), ()))
            for token in tokens:
                resolved = index.get(str(token).lower())
                if resolved:
                    found.add(resolved)
                    break
    return found


def search_index(project_config):
    """{lowercased token: gene id} for one species: its gene ids, current names and every recorded alias.

    A curated annotation records a gene's other names beside its current one - SGD gives ADY2 `Alias=ADY2,ATO1`,
    and CGD does the same - and those are often the names people search with, because the literature uses them.
    Read straight from the GFF's gene rows rather than through the annotation parser, which keeps only the
    current name. Cached per process: a union search reads every represented species' annotation, and the same
    query would otherwise reread each one.
    """
    key = str(project_config)
    if key in _SEARCH_INDEX:
        return _SEARCH_INDEX[key]
    cfg = yaml.safe_load(Path(project_config).read_text())
    gff = Path((cfg.get("reference") or {})["gff"])
    gff = gff if gff.is_absolute() else ROOT/gff
    index = {}
    with gff.open() as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "pseudogene"}:
                continue
            attrs = {}
            for entry in fields[8].split(";"):
                if "=" in entry:
                    key, value = entry.split("=", 1)
                    attrs[unquote(key)] = unquote(value)
            gene_id = canonical_gene_id(attrs.get("locus_tag") or attrs.get("ID") or "")
            if not gene_id:
                continue
            names = [gene_id, attrs.get("ID", "")]
            for field in ("Alias", "gene_synonym", "Gene_synonym", "synonym", "gene", "Gene", "Name"):
                names.extend(part.strip() for part in str(attrs.get(field, "")).split(",") if part.strip())
            for name in names:
                if name:
                    index.setdefault(name.lower(), gene_id)
    _SEARCH_INDEX[key] = index
    return index


def resolve_query(query, chosen):
    """({species: {gene ids named by the query}}, [hub genes they imply]) over every represented species.

    A union over the species rather than a route through one of them. Routing a query to S. cerevisiae first and
    asking every other species for counterparts can only find what the ortholog tables connect, and that silently
    loses genes which are present, named exactly as asked for, and simply not mapped: C. albicans carries a gene
    called ATO1, and a query for ATO1 reported the species as having no counterpart because the table has no edge
    from it to S. cerevisiae ADY2. Searching every species' own identifiers, names and aliases finds the gene
    wherever it is; the hub genes it maps to then still bring in the orthologs, so the two are additive.
    """
    wanted = query.strip().lower()
    seeds, hubs = {}, set()
    for species, entry in chosen.items():
        gene = search_index(entry["project_config"]).get(wanted)
        if not gene:
            continue
        seeds.setdefault(species, set()).add(gene)
        _, forward = ortholog_families(entry["project_config"])
        if not forward:
            # No ortholog table means this project is the S. cerevisiae hub itself, so the gene is its own hub.
            hubs.add(gene)
        else:
            hubs.update(forward.get(gene, []))
    return seeds, sorted(hubs)


def comparison_genes(project_config, species, hubs, seeded):
    """({gene id: how it was found}) for one species, from three routes unioned in order of what each asserts.

    The routes are kept apart because they are not equally strong evidence, and the strongest one that reached a
    gene is the one reported:

      ortholog    the species' own curated ortholog table says this gene corresponds to the query. Strongest,
                  and authoritative where it exists: for the distant species it also has the best coverage.
      orthogroup  OrthoFinder puts the two in one gene family. That is relatedness, not correspondence, and it
                  is the only route that can report a family a one-to-one table cannot express.
      name        the gene carries the name asked for and nothing connects it to the rest. Weakest, and often
                  the only route to a lineage-specific expansion, which keeps the name without a counterpart.

    Unioned rather than ranked against each other, so adding a route can only add genes: a gap in any one source
    costs an edge, never a wrong or a lost answer.
    """
    families, _ = ortholog_families(project_config)
    found = {}
    for hub in hubs:
        for gene in ([hub] if not families else families.get(hub, [])):
            found[gene] = "ortholog"
    for gene in sorted(orthogroup_counterparts(project_config, species, hubs)):
        found.setdefault(gene, "orthogroup")
    for gene in sorted(seeded):
        found.setdefault(gene, "name")
    return found


def family_note(row):
    """'2 of 3' where the family has expanded in this species, blank where the gene has a single counterpart."""
    n = int(row.get("n_orthologs") or 0)
    return f"{row.get('ortholog_index')} of {n}" if n > 1 else ""


def gene_blocks(project_config, gene_id):
    """The gene's CDS blocks in its own annotation, by systematic id or GFF id."""
    cfg = yaml.safe_load(Path(project_config).read_text())
    gff = Path((cfg.get("reference") or {})["gff"])
    gff = gff if gff.is_absolute() else ROOT/gff
    for gene in annotations(gff):
        ids = {str(gene.get("id", "")).split(":", 1)[-1].replace("gene-", ""),
               str(gene.get("gff_id", "")).split(":", 1)[-1].replace("gene-", "")}
        if gene_id in ids and gene["blocks"]:
            return gene
    return None


def profile_many(entry, genes):
    """{gene id: (rows, coding length)} for several genes at once, in one pass over each library's hit file.

    Taken together rather than one gene at a time because a species where the family has expanded would
    otherwise re-read the same multi-million-row files once per family member, which is the whole cost of the
    query: the windows are few and checking a position against all of them is far cheaper than another pass.
    """
    specs = []
    for gene_id, gene in genes:
        length = sum(b-a+1 for a, b in gene["blocks"])
        lo, hi = gene["blocks"][0][0], gene["blocks"][-1][1]
        # The flank is a fraction of the CODING length, not of the genomic span: for a gene with introns the span is
        # the wider of the two, and padding on it collects sites past the window the axis claims to show.
        pad = int(FLANK*length)
        specs.append(dict(gene_id=gene_id, gene=gene, length=length, lo=lo, hi=hi,
                          window=(lo-pad, hi+pad), sites={}))
    by_chrom = {}
    for spec in specs:
        by_chrom.setdefault(spec["gene"]["chrom"], []).append(spec)
    project_root = ROOT/"output/projects"/entry["project"]
    for library in entry["libraries"]:
        hits = project_root/"create_hit_file"/library/f"{library}_hits.txt"
        if not hits.is_file():
            continue
        with hits.open() as handle:
            # csv.reader with column indices rather than DictReader: these files run to millions of rows and
            # building a dict per row is most of the cost of answering a query about one gene.
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if not header:
                continue
            chrom_at, pos_at, count_at = (header.index(k) for k in ("Chromosome", "Hit position", "Hit count"))
            for row in reader:
                group = by_chrom.get(row[chrom_at])
                if not group:
                    continue
                position = int(row[pos_at])
                for spec in group:
                    if spec["window"][0] <= position <= spec["window"][1]:
                        spec["sites"][position] = spec["sites"].get(position, 0)+int(row[count_at])
    return {spec["gene_id"]: (_place(spec), spec["length"]) for spec in specs}


def _place(spec):
    """The collected sites as rows, each placed by its fraction along the coding sequence."""
    gene, length, lo, hi = spec["gene"], spec["length"], spec["lo"], spec["hi"]
    rows = []
    for position, reads in sorted(spec["sites"].items()):
        if position < lo:
            offset = position-lo
        elif position > hi:
            offset = length-1+(position-hi)
        else:
            offset = sum(min(position, end)-start+1 for start, end in gene["blocks"] if start <= position)-1
        fraction = offset/length if gene["strand"] != "-" else (length-1-offset)/length
        rows.append(dict(position=position, cds_fraction=round(fraction, 6),
                         in_cds=0 <= fraction <= 1, reads=reads))
    return rows


def gene_calls(query, chosen, too_sparse):
    """One row per (query, species, counterpart): what that species' representative calls it, and how it matched.

    The calls alone, without the insertion profiles. A list of genes is asked the same question as a single gene
    - which species call it essential - but answering it does not need the hit files, only the run's gene table,
    so a group of a hundred genes costs about what one profiled gene costs to resolve.
    """
    seeds, hubs = resolve_query(query, chosen)
    rows = []
    if not seeds:
        for species, entry in sorted(chosen.items()):
            rows.append(dict(query=query, species=species, label=xs.species_label(species),
                             project=entry["project"], gene="", gene_name="", matched_by="", call="", P_ES="",
                             sites_per_coding_kb=round(entry["sites_per_coding_kb"], 2),
                             represented=1, status="not found in any species"))
        return rows
    for species, entry in sorted(chosen.items()):
        found = comparison_genes(entry["project_config"], species, hubs, seeds.get(species, set()))
        calls = xs.hmm_calls(entry["run"])
        base = dict(query=query, species=species, label=xs.species_label(species), project=entry["project"],
                    sites_per_coding_kb=round(entry["sites_per_coding_kb"], 2), represented=1)
        if not found:
            rows.append(dict(base, gene="", gene_name="", matched_by="", call="", P_ES="",
                             status="no counterpart"))
            continue
        for gene_id in sorted(found):
            row = dict(base, gene=gene_id, gene_name="", matched_by=found[gene_id],
                       call="", P_ES="", status="no call")
            if gene_id in calls:
                p_es, essential, name = calls[gene_id]
                # The name the species itself gives the gene. Without it a panel reads as a bare systematic id,
                # and a genuine paralog - YGR088W beside YDR256C - is indistinguishable from the same gene twice.
                row["gene_name"] = "" if str(name) == gene_id else str(name)
                row["call"] = "ES" if essential else "not ES"
                row["P_ES"] = round(p_es, 4)
                row["status"] = "called"
            rows.append(row)
    for species, entry in sorted(too_sparse.items()):
        rows.append(dict(query=query, species=species, label=xs.species_label(species), project=entry["project"],
                         gene="", gene_name="", matched_by="", call="", P_ES="",
                         sites_per_coding_kb=round(entry["sites_per_coding_kb"], 2), represented=0,
                         status="libraries cannot support an essentiality call"))
    return rows


def run_group(args, parser):
    """A list of genes, answered exactly as one gene is: the same representatives, the same three matching routes
    and the same insertion profiles, so a group reads like a stack of Comparative View figures.

    Profiling the whole list costs little more than profiling one gene of it. The expensive part is reading each
    representative library's hit file, and that is one pass whatever is being looked for, so every counterpart of
    every gene in the group is collected in the same pass rather than a pass per gene.
    """
    queries = [line.strip() for line in Path(args.genes).read_text().splitlines() if line.strip()]
    queries = list(dict.fromkeys(queries))
    if not queries:
        parser.exit(1, f"ERROR: {args.genes} holds no gene names.\n")
    chosen, too_sparse = representatives(args.min_density)
    if not chosen:
        parser.exit(1, "ERROR: no project reaches the density an essentiality call needs; run `ytab qc` first.\n")

    rows, wanted = [], {}
    for query in queries:
        for row in gene_calls(query, chosen, too_sparse):
            rows.append(row)
            if row.get("gene") and row.get("represented") == 1:
                wanted.setdefault(row["species"], {}).setdefault(row["gene"], []).append(query)

    profiles = []
    if not args.no_profiles:
        for species, entry in sorted(chosen.items()):
            genes = wanted.get(species) or {}
            located = {g: gene_blocks(entry["project_config"], g) for g in sorted(genes)}
            pairs = [(g, b) for g, b in located.items() if b is not None]
            if not pairs:
                continue
            drawn = profile_many(entry, pairs)
            for gene_id, (sites, length) in drawn.items():
                for query in genes[gene_id]:
                    for item in sites:
                        profiles.append(dict(query=query, species=species, label=xs.species_label(species),
                                             project=entry["project"], gene=gene_id, cds_length=length, **item))
        # Whether a gene was profiled is part of the answer: an empty panel because the libraries hold no
        # insertion there is not the same as one because the gene was never looked for.
        drawn_keys = {(r["query"], r["species"], r["gene"]) for r in profiles}
        for row in rows:
            if row.get("gene") and row.get("represented") == 1:
                row["profiled"] = 1 if (row["query"], row["species"], row["gene"]) in drawn_keys else 0

    args.output.mkdir(parents=True, exist_ok=True)
    fields = ["query", "species", "label", "project", "gene", "gene_name", "matched_by", "call", "P_ES",
              "sites_per_coding_kb", "represented", "profiled", "status"]
    xs.write(args.output/"group_calls.tsv", [{k: r.get(k, "") for k in fields} for r in rows], fields)
    if profiles:
        xs.write(args.output/"profiles.tsv", profiles)
    resolved = {r["query"] for r in rows if r.get("gene")}
    manifest = dict(status="success", genes=len(queries), resolved=len(resolved),
                    unresolved=sorted(set(queries) - resolved), min_density=args.min_density,
                    flank_fraction=FLANK, sites=len(profiles),
                    species=sorted({r["species"] for r in rows if r.get("represented") == 1}),
                    representatives={s: dict(project=e["project"], run=e["run"].name,
                                             libraries=e["libraries"],
                                             sites_per_coding_kb=round(e["sites_per_coding_kb"], 2))
                                     for s, e in sorted(chosen.items())},
                    files=dict(group_calls="group_calls.tsv",
                               profiles="profiles.tsv" if profiles else ""))
    (args.output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps(manifest, indent=2) if args.json else
          f"{len(resolved)} of {len(queries)} genes resolved, {len(profiles)} sites, into {args.output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gene", help="Gene name or systematic id, in any species YTAB has reprocessed")
    parser.add_argument("--genes", type=Path,
                        help="A file of gene names, one per line: report each one's call in every species "
                             "without profiling it, which is what a gene group needs and is far cheaper")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-density", type=float, default=HYBRID_SITES_PER_CODING_KB,
                        help=f"Pooled sites per coding kb a species must reach to be given a call "
                             f"(default {HYBRID_SITES_PER_CODING_KB:g}; the label-free HMM asks for {HMM_SITES_PER_CODING_KB:g})")
    parser.add_argument("--no-profiles", action="store_true",
                        help="With --genes, report the calls only and skip reading the hit files")
    parser.add_argument("--force", action="store_true", help="Rebuild even when the cached answer still applies")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    if bool(args.gene) == bool(args.genes):
        parser.exit(2, "ERROR: give either --gene (one gene, profiled) or --genes (a list, calls only).\n")
    if args.genes:
        return run_group(args, parser)

    chosen, too_sparse = representatives(args.min_density)
    if not chosen:
        parser.exit(1, "ERROR: no project reaches the density an essentiality call needs; run `ytab qc` first.\n")
    seeds, hubs = resolve_query(args.gene, chosen)
    if not seeds:
        parser.exit(1, f"ERROR: {args.gene} did not match a gene, name or alias in any representative species.\n")

    # Answering one gene means scanning each representative library's hit file, so a finished answer is reused
    # until the inputs change: what the query resolved to, the density rule, and which run represents each species.
    resources = sorted((p.name, p.stat().st_mtime_ns) for p in ORTHOGROUP_DIR.glob("*.csv")) \
        if ORTHOGROUP_DIR.is_dir() else []
    signature = json.dumps(dict(hubs=hubs, seeds={k: sorted(v) for k, v in sorted(seeds.items())},
                                min_density=args.min_density, flank=FLANK, resources=resources,
                                runs={s: [e["run"].name, e["libraries"]] for s, e in sorted(chosen.items())}),
                           sort_keys=True)
    done = args.output/"manifest.json"
    if not args.force and done.is_file():
        previous = json.loads(done.read_text())
        if previous.get("signature") == signature:
            print(json.dumps(previous, indent=2) if args.json else
                  f"reused {args.output} ({previous.get('sites', 0)} sites; --force to rebuild)")
            return

    args.output.mkdir(parents=True, exist_ok=True)
    profiles, summary = [], []
    for species, entry in sorted(chosen.items()):
        found = comparison_genes(entry["project_config"], species, hubs, seeds.get(species, set()))
        gene_ids = sorted(found)
        calls = xs.hmm_calls(entry["run"])

        def base(gene_id="", n=0, index=0):
            return dict(species=species, label=xs.species_label(species), project=entry["project"],
                        run=entry["run"].name, libraries=" ".join(entry["libraries"]),
                        n_libraries=len(entry["libraries"]), background=entry["background"],
                        sites_per_coding_kb=round(entry["sites_per_coding_kb"], 2),
                        gene=gene_id, gene_name="", n_orthologs=n, ortholog_index=index, represented=1,
                        matched_by=found.get(gene_id, ""), call="", P_ES="", status="")

        # A species with no counterpart is still reported as a row of its own, so that the display can show it as
        # an empty panel: "this species has no ortholog" and "this species was not looked at" are different
        # answers, and a species dropped from the table is indistinguishable from one never considered.
        if not gene_ids:
            row = base(); row["status"] = "no ortholog"
            summary.append(row); continue
        # The whole family's blocks are looked up first, then profiled in one pass over the libraries.
        located = {gene_id: gene_blocks(entry["project_config"], gene_id) for gene_id in gene_ids}
        drawn = profile_many(entry, [(gene_id, gene) for gene_id, gene in located.items() if gene is not None])
        for index, gene_id in enumerate(gene_ids, start=1):
            row = base(gene_id, len(gene_ids), index)
            if gene_id in calls:
                p_es, essential, name = calls[gene_id]
                row["call"], row["P_ES"] = ("ES" if essential else "not ES"), round(p_es, 4)
                row["gene_name"] = "" if str(name) == gene_id else str(name)
            if located[gene_id] is None:
                row["status"] = "not in this species' annotation"
                summary.append(row); continue
            rows, length = drawn[gene_id]
            row["status"] = "profiled" if rows else "no insertion in the gene or its flanks"
            row["cds_length"] = length
            row["sites"] = len(rows)
            summary.append(row)
            for item in rows:
                profiles.append(dict(species=species, label=xs.species_label(species), project=entry["project"],
                                     gene=gene_id, n_orthologs=len(gene_ids), ortholog_index=index,
                                     matched_by=found.get(gene_id, ""), cds_length=length, **item))
    for species, entry in sorted(too_sparse.items()):
        # Two different answers that an empty column would blur: the libraries cannot support a call, or they
        # could but nobody has fitted one yet.
        density = entry["sites_per_coding_kb"]
        status = (f"no essentiality run yet ({density:.1f} sites/coding kb would support one)"
                  if entry["run"] is None and density >= args.min_density
                  else f"no library dense enough to call essentiality (needs {args.min_density:g})")
        summary.append(dict(species=species, label=xs.species_label(species), project=entry["project"],
                            run=entry["run"].name if entry["run"] is not None else "",
                            libraries=" ".join(entry["libraries"]), n_libraries=len(entry["libraries"]),
                            background=entry["background"], sites_per_coding_kb=round(density, 2),
                            gene="", gene_name="", n_orthologs=0, ortholog_index=0, represented=0,
                            matched_by="", call="", P_ES="", status=status))

    fields = ["species", "label", "project", "run", "libraries", "n_libraries", "background", "sites_per_coding_kb",
              "gene", "gene_name", "n_orthologs", "ortholog_index", "represented", "matched_by", "call",
              "P_ES", "cds_length", "sites", "status"]
    xs.write(args.output/"representatives.tsv", [{k: r.get(k, "") for k in fields} for r in summary], fields)
    if profiles:
        xs.write(args.output/"profiles.tsv", profiles)
    manifest = dict(status="success", signature=signature, query=args.gene, hub_genes=hubs,
                    found_in=sorted(seeds), seeds={k: sorted(v) for k, v in sorted(seeds.items())},
                    min_density=args.min_density, flank_fraction=FLANK,
                    species_profiled=len({r["species"] for r in summary if r.get("status") == "profiled"}),
                    genes_profiled=sum(1 for r in summary if r.get("status") == "profiled"),
                    species_without_ortholog=len({r["species"] for r in summary if r.get("status") == "no ortholog"}),
                    expanded_families={r["species"]: r["n_orthologs"] for r in summary
                                       if int(r.get("n_orthologs") or 0) > 1},
                    species_too_sparse=len(too_sparse), sites=len(profiles),
                    files=dict(profiles="profiles.tsv", representatives="representatives.tsv"))
    (args.output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps(manifest, indent=2) if args.json else
          "\n".join(f"{r['species']:18s} {r.get('gene',''):16s} "
                    f"{family_note(r):8s} {r.get('matched_by',''):9s} {r.get('call',''):7s} "
                    f"{r.get('sites_per_coding_kb','')!s:>7s} sites/kb  {r['status']}"
                    for r in summary))


if __name__ == "__main__":
    main()
