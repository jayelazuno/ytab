"""Read-only adapter for YTAB's raw hits, CDS annotation and classifier calls."""
from collections import defaultdict
import csv
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import yaml

from .statistics import empirical_gene_probabilities


INFORMATIVE_DEFINITION = (
    "Insertion position within annotated CDS blocks (1-based inclusive), reusing "
    "the coding-exon definition in YTAB GlabrataFeatureDB/SummaryTable. This is a "
    "coding-insertion proxy for informative disruption, not a validated knockout "
    "call. No first-100-codon or central-ORF cutoff is defined by this pipeline."
)


def _path(value, root):
    p = Path(value)
    return p if p.is_absolute() else root / p


def _coding_genes_by_parent(gff):
    """CDS union per gene for GFF3 files whose CDS name their gene only through Parent (CDS -> transcript -> gene)."""
    from ytab.essentiality.control_conditioned_gumbel import annotations
    genes = {}
    for g in annotations(gff):
        blocks = sorted(g["blocks"])
        genes[g["id"]] = dict(gene_id=g["id"], common_name=g["name"] if g["name"] != g["id"] else "", chromosome=g["chrom"], strand=g["strand"],
                              blocks=blocks, essentiality="Unknown", length=sum(b-a+1 for a, b in blocks), start=blocks[0][0], end=blocks[-1][1])
    if not genes:
        raise ValueError("No protein-coding CDS annotations available")
    return genes


def read_coding_genes(gff):
    """CDS union per locus, excluding RNA-only features and pseudogenes without CDS."""
    genes = {}
    with Path(gff).open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            gid = attrs.get("locus_tag") or attrs.get("gene_id")
            if not gid:  # CDS linked to its gene only through Parent (Ensembl-style GFF3): resolve the ancestry instead
                return _coding_genes_by_parent(gff)
            gene = genes.setdefault(gid, dict(gene_id=gid, common_name=unquote(attrs.get("gene", "")),
                chromosome=fields[0], strand=fields[6], blocks=[], essentiality="Unknown"))
            if gene["chromosome"] != fields[0] or gene["strand"] != fields[6]:
                raise ValueError(f"Inconsistent coordinates for {gid}")
            gene["blocks"].append((int(fields[3]), int(fields[4])))
    for gene in genes.values():
        merged = []
        for a, b in sorted(set(gene["blocks"])):
            if merged and a <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
            else:
                merged.append((a, b))
        gene["blocks"] = merged
        gene["length"] = sum(b - a + 1 for a, b in merged)
        gene["start"], gene["end"] = merged[0][0], merged[-1][1]
    if not genes:
        raise ValueError("No protein-coding CDS annotations available")
    return genes


def _classify(genes, path):
    if path is None:
        return "No classifier available; unknown genes are excluded from non-essential universes"
    with path.open() as handle:
        reader = csv.DictReader(handle)
        labels = [c for c in reader.fieldnames if "ess. for fpr" in c.lower()]
        if len(labels) != 1 or "Standard name" not in reader.fieldnames:
            raise ValueError("Classifier must contain Standard name and one essentiality verdict column")
        seen = set()
        for row in reader:
            gid = row["Standard name"]
            if gid in seen:
                raise ValueError(f"Duplicate classifier gene {gid}")
            seen.add(gid)
            if gid in genes:
                genes[gid]["essentiality"] = {"yes": "Essential", "no": "Non-essential"}.get(row[labels[0]].strip().lower(), "Unknown")
    return labels[0]


def is_parent_sample(row):
    """Explicit treatment metadata overrides misleading 'parent' tokens in names."""
    values = [str(row.get(k, "")).strip().lower() for k in
              ("control_or_treated", "library_role", "condition", "guessed_condition", "fitness_role", "treatment")]
    if any(re.search(r"\btreated\b|\btreatment\b", v) for v in values):
        return False
    if any(v in {"parent", "control", "mock", "untreated", "parent/control", "parental"} for v in values):
        return True
    name = row.get("sample", "").lower()
    return not re.search(r"(?:^|[-_ ])treated(?:$|[-_ ])", name) and bool(
        re.search(r"(?:^|[-_ ])(?:parent|control|mock)(?:$|[-_ ])", name))


def load_library(project_config, samples, classifier=None, weights=None, universe="informative", basis="sites", conditioning="physical",
                 parent_sources=None, catalogue_mode="union"):
    config = Path(project_config).resolve()
    root = next((p for p in config.parents if (p / "src/ytab").is_dir()), None)
    if root is None:
        raise ValueError("Project must belong to a YTAB repository")
    cfg = yaml.safe_load(config.read_text())
    project = _path(cfg["output_project_dir"], root)
    ref_file = project / "config/reference_resolved.json"
    ref = json.loads(ref_file.read_text()) if ref_file.exists() else cfg["reference"]
    gff = _path(ref["gff"], root)
    genes = read_coding_genes(gff)
    sources = [config, gff]
    if ref_file.exists():
        sources.append(ref_file)
    if classifier:
        cp = Path(classifier).resolve()
        if not cp.is_relative_to((project / "classifier").resolve()) or not cp.is_file():
            raise ValueError("Select an existing classifier result from this project")
    else:
        cp = None
    label = _classify(genes, cp)
    if cp:
        sources.append(cp)
    available, selected_samples = {}, []
    for source in parent_sources or [dict(project=str(config), samples=samples)]:
        source_config = Path(source["project"]).resolve()
        source_cfg = yaml.safe_load(source_config.read_text())
        source_project = _path(source_cfg["output_project_dir"], root)
        source_ref_file = source_project / "config/reference_resolved.json"
        source_ref = json.loads(source_ref_file.read_text()) if source_ref_file.exists() else source_cfg["reference"]
        source_gff = _path(source_ref["gff"], root)
        # Coordinate union is valid only in the same annotated reference space.
        if hashlib.sha256(source_gff.read_bytes()).digest() != hashlib.sha256(gff.read_bytes()).digest():
            raise ValueError("Parent union requires matching reference GFF coordinates; incompatible projects selected")
        sheet = _path(source_cfg["sample_sheet"], root)
        sources.extend([source_config, source_gff, sheet])
        if source_ref_file.exists():
            sources.append(source_ref_file)
        with sheet.open() as handle:
            rows = {r["sample"]: r for r in csv.DictReader(handle)
                    if r.get("include", "true").lower() in {"true", "yes", "1"} and is_parent_sample(r)}
        for sample in source["samples"]:
            if sample not in rows:
                raise ValueError(f"Screen Coverage accepts included parent/control samples only: {sample}")
            key = f"{source_project.name}::{sample}"
            if key in selected_samples:
                raise ValueError("Select each parent sample only once")
            selected_samples.append(key)
            available[key] = (rows[sample], source_project)
    samples = selected_samples
    if not samples:
        raise ValueError("Choose one or more included parent/control samples")
    if catalogue_mode not in {"union", "mixture"}:
        raise ValueError("Unknown parent catalogue mode")
    if weights is None:
        weights = [1.] * len(samples)
    mix = np.asarray(weights, dtype=float)
    if len(mix) != len(samples) or np.any(~np.isfinite(mix)) or np.any(mix <= 0):
        raise ValueError("Provide one positive physical mixture weight per selected sample")
    mix /= mix.sum()
    if universe not in {"all", "nonessential", "observed", "informative"} or basis not in {"reads", "sites"} or conditioning not in {"physical", "conditional"}:
        raise ValueError("Invalid coverage settings")
    chrom_genes = defaultdict(list)
    for g in genes.values():
        chrom_genes[g["chromosome"]].append(g)
    union_sites, informative_sites, observed = set(), defaultdict(set), defaultdict(set)
    abundances = defaultdict(float)
    gene_p = {b: defaultdict(float) for b in ("reads", "sites")}
    pool_summary, assignments = [], {}
    for sample, fraction in zip(samples, mix):
        row, source_project = available[sample]
        sample_name = row["sample"]
        hit_file = _path(row.get("hit_file", ""), root)
        if not hit_file.is_file():
            hit_file = source_project / "create_hit_file" / sample_name / f"{sample_name}_hits.txt"
        sources.append(hit_file)
        sites = defaultdict(float)
        with hit_file.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"Chromosome", "Hit position", "Hit count"}.issubset(reader.fieldnames):
                raise ValueError("Raw insertion file lacks coordinates or read abundance")
            for hit in reader:
                pos, reads = int(hit["Hit position"]), float(hit["Hit count"])
                if pos < 1 or not np.isfinite(reads) or reads < 0:
                    raise ValueError("Invalid raw insertion coordinates/read abundance")
                if reads > 0:
                    sites[(hit["Chromosome"], pos)] += reads
        if not sites:
            raise ValueError(f"No positive-abundance insertions for {sample}")
        # Sorted coordinates plus searchsorted avoids scanning all genes per hit.
        per_gene = defaultdict(lambda: [0, 0.])
        ambiguous = set()
        for chrom in {x[0] for x in sites}:
            positions = np.array(sorted(pos for cc, pos in sites if cc == chrom))
            hits_by_pos = defaultdict(set)
            for g in chrom_genes[chrom]:
                for pos in positions[np.searchsorted(positions, g["start"]):np.searchsorted(positions, g["end"], side="right")]:
                    observed[g["gene_id"]].add((chrom, int(pos)))
                for a, b in g["blocks"]:
                    for pos in positions[np.searchsorted(positions, a):np.searchsorted(positions, b, side="right")]:
                        hits_by_pos[int(pos)].add(g["gene_id"])
            for pos, gids in hits_by_pos.items():
                key = (chrom, pos)
                if len(gids) != 1:
                    ambiguous.add(key)
                    continue
                gid = next(iter(gids))
                assignments[key] = gid
                informative_sites[gid].add(key)
                abundances[gid] += sites[key]
                per_gene[gid][0] += 1
                per_gene[gid][1] += sites[key]
        total_reads = sum(sites.values())
        for gid, (count, reads) in per_gene.items():
            gene_p["sites"][gid] += fraction * count / len(sites)
            gene_p["reads"][gid] += fraction * reads / total_reads
        union_sites.update(sites)
        pool_summary.append(dict(sample=sample_name, project=source_project.name, pool=row.get("pool", ""), role="parent/control",
            background=row.get("background", ""), mixture_fraction=float(fraction), unique_sites=len(sites),
            read_abundance=total_reads, coding_sites=sum(v[0] for v in per_gene.values()), ambiguous_sites=len(ambiguous)))
    if catalogue_mode == "union":
        gene_p["sites"] = defaultdict(float, {gid: len(sites) / len(union_sites)
                                             for gid, sites in informative_sites.items()})
    nonessential = {gid for gid, g in genes.items() if g["essentiality"] == "Non-essential"}
    universes = dict(all=set(genes), nonessential=nonessential,
        observed=nonessential & set(observed), informative=nonessential & set(informative_sites))
    selected = sorted(universes[universe])
    if not selected:
        raise ValueError("Chosen universe is empty; select a classifier result or All annotated protein-coding genes")
    rows = []
    for gid in selected:
        g = genes[gid]
        rows.append({k: g[k] for k in ("gene_id", "common_name", "chromosome", "start", "end", "strand", "length", "essentiality")})
        rows[-1].update(observed_sites=len(observed[gid]), informative_sites=len(informative_sites[gid]),
            read_abundance=abundances[gid], p_sites=gene_p["sites"][gid], p_reads=gene_p["reads"][gid])
    for b in ("sites", "reads"):
        p = np.array([g[f"p_{b}"] for g in rows])
        if conditioning == "conditional":
            p = empirical_gene_probabilities(p)
        for g, prob in zip(rows, p):
            g[f"p_{b}"] = float(prob)
            if b == basis:
                g["p"] = float(prob)
    fingerprints = []
    for path in dict.fromkeys(sources):
        fingerprints.append(dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    return dict(genes=rows, pools=pool_summary, universes={k: len(v) for k, v in universes.items()},
        annotated_genes=len(genes), essential_genes=sum(g["essentiality"] == "Essential" for g in genes.values()),
        unknown_genes=sum(g["essentiality"] == "Unknown" for g in genes.values()),
        nonessential_without_informative=sum(not informative_sites[gid] for gid in nonessential),
        observed_genes=sum(bool(v) for v in observed.values()),
        informative_genes=sum(bool(v) for v in informative_sites.values()),
        unique_sites=len(union_sites), informative_sites=sum(len(v) for v in informative_sites.values()),
        selected_informative_sites=sum(g["informative_sites"] for g in rows),
        catalogue_mode=catalogue_mode, parents_only=True,
        repeated_site_observations=sum(x["unique_sites"] for x in pool_summary)-len(union_sites),
        rarefaction_other_sites=len(union_sites) - sum(g["informative_sites"] for g in rows),
        classifier=str(cp) if cp else "Unavailable", classifier_verdict=label,
        informative_definition=INFORMATIVE_DEFINITION, basis=basis, conditioning=conditioning,
        sources=fingerprints, read_abundance_usable=True,
        limitations=[
            "Coding insertions are candidate disruptive alleles; functional knockout is unverified.",
            "Essentiality calls are project classifier predictions, conditional on growth conditions; Unknown is not Non-essential.",
            "Read fractions are estimates affected by PCR, extraction, mapping, growth and sequencing bottlenecks; viable cell frequencies are unavailable.",
            "Raw mapped positive-read sites are used; no additional mappability exclusion mask is inferred.",
            "Overlapping CDS assignments are excluded from gene recovery and retained as non-target sampling mass.",
            "The default site model uses a distinct-coordinate UNION of selected parents, not an intersection or a sum of duplicate coordinates. Treated samples are excluded.",
            "A catalogue union maximizes observed insertion representation but is not evidence of a physically mixed population. Colony estimates assume each union site is an equally recoverable mutant.",
            "Read-weighted results are a sensitivity analysis: normalize each parent's read distribution then average using the stated weights (equal by default). This does not measure physical union abundance. Summed read abundance is reported only as sequencing support.",
            "Redundancy counts sampled colonies and can include repeated recovery of the same allele; it does not guarantee independent disruptive alleles.",
            "Rarefaction counts distinct genomic sites without replacement; it is not sequencing-read rarefaction or viable-colony sampling."])
