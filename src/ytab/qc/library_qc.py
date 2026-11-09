"""Per-library QC table for a project: every QC metric in one table, flagged and commented, plus plot data.

Collects, for each included library, the read and mapping metrics (read cleanup, MapFastq), the Tn-seq count
statistics computed from the hit file (tnseq_stats.statistics, over the project's empirical site catalogue), the
strand correlation, the SummaryTable gene-level metrics and the saved library diagnostics (MidLC, jackpots).

Flags (green / yellow / red; grey = not flagged) never use invented cut-offs:
- "rule of thumb": only where published Tn-seq QC guidance is scale-free and carries over to these libraries (skewness of the
  positive counts: above 30 often signals a problem, below 50 is desired; see docs/library_qc.md for the source).
- "project-relative": a library is compared with the project's libraries of the same role (parents with parents,
  treated with treated) by robust z = (x - median) / (1.4826 MAD), in the metric's bad direction only; yellow beyond
  2.5 and red beyond 3 (Leys et al. 2013, J Exp Soc Psychol 49:764, "moderately" and "very" conservative). Needs at
  least 4 libraries of the role; otherwise the flag is grey with that reason. The spread is at least 10% of the role
  median (MIN_RELATIVE_SPREAD, a YTAB choice), so near-identical libraries are not flagged for trivial differences.
- descriptive metrics (e.g. kurtosis, not very meaningful on its own) are not flagged.
Each flag carries a comment and a normalization hint where one applies (skewed or jackpot-dominated libraries: BGC;
isolated outliers: TTR is robust to them).
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import gzip
import json
import math
from pathlib import Path
import re

import numpy as np

from ytab.qc.tnseq_stats import statistics

ROBUST_YELLOW, ROBUST_RED, MIN_ROLE_LIBRARIES = 2.5, 3.0, 4
# Spread floor: with few, near-identical libraries the MAD is tiny and trivial differences would be flagged; a
# difference smaller than this fraction of the role median is never flagged. A YTAB choice, not a published value.
MIN_RELATIVE_SPREAD = 0.10
LEYS = "Compared with the project's libraries of the same role (robust z; Leys et al. 2013)"
RULE_OF_THUMB = "Tn-seq rule of thumb (skewness)"
BGC_HINT = "Consider BGC (betageom): it suits skewed libraries dominated by a few high counts."
TTR_HINT = "TTR is robust to isolated outlier sites; BGC if many sites are affected."

# key, theme, label, bad direction ("low", "high", "both" or None), rule ("rule_of_thumb", "relative", None), comment
METRICS = [
    ("input_reads", "reads", "Reads sequenced", "low", "relative", "Reads entering read cleanup (the junction read)."),
    ("retained_fraction", "reads", "Reads retained after cleanup", "low", "relative", "Fraction of reads kept after transposon-end and adaptor trimming and the length filter."),
    ("adapter_fraction", "reads", "Reads with 3' adaptor", "high", "relative", "Fraction of reads with adaptor at the 3' end (short fragments or primer dimers)."),
    ("mapped_reads", "reads", "Mapped reads", "low", "relative", "Primary alignments that map to the reference."),
    ("percent_mapped", "reads", "Mapped (%)", "low", "relative", "Primary alignments that map to the reference."),
    ("percent_mapq", "reads", "Mapped with MAPQ above threshold (%)", "low", "relative", "Uniquely placed alignments; low values mean repeats or the wrong reference."),
    ("mean_mapq", "reads", "Mean MAPQ", "low", "relative", "Mean mapping quality of mapped primary alignments."),
    ("percent_donor", "reads", "Reads from the transposon donor (%)", "high", "relative", "Reads that start at the transposon end and run into the donor plasmid or donor locus rather than the genome: the donor was not lost after induction. Counted against the donor sequence when the project names one (mapping.donor_fasta), otherwise estimated from recurring unmapped flanks (a donor junction is one fixed sequence). These reads carry no insertion; high values waste depth and mean the donor persists in the pool."),
    ("percent_donor_shared", "reads", "Donor reads placed on the genome (%)", "high", "relative", "Donor reads whose flank the donor shares with the host genome (a promoter, marker or CEN/ARS fragment on the plasmid): they map to that locus and look like insertions there. Only counted when the donor sequence is known."),
    ("total_cts", "insertions", "Reads at insertion sites", "low", "relative", "Reads supporting called insertions."),
    ("unique_sites", "insertions", "Insertion sites", "low", "relative", "Distinct insertion positions (strands combined)."),
    ("insertions_per_kb", "insertions", "Insertions per kb of genome", "low", "relative", "Insertion sites per kilobase of reference: how densely the genome is covered (Hermes has no fixed set of possible sites)."),
    ("density", "insertions", "Catalogue density", "low", "relative", "Fraction of the project's catalogue (positions with an insertion in any library) that this library hits. Hermes has no fixed set of possible sites, so there is no absolute target."),
    ("mean_ct", "insertions", "Mean reads per catalogue site", None, None, "Mean reads per catalogue position, empty positions included."),
    ("NZmean", "insertions", "Mean reads per occupied site", "both", "relative", "Mean reads per position with an insertion; too high or too low can indicate problems."),
    ("NZmedian", "insertions", "Median reads per occupied site", "both", "relative", "Median reads per position with an insertion; robust to spikes."),
    ("max_ct", "insertions", "Max reads at a site", None, None, "Largest read count at one position; depends on depth, see the top site's share."),
    ("max_share", "insertions", "Top site's share of reads", "high", "relative", "Fraction of all reads at the single largest site (a PCR jackpot or strong selection)."),
    ("top1pct_share", "insertions", "Top 1% of sites' share of reads", "high", "relative", "Fraction of reads at the top 1% of sites; high values mean a few sites dominate."),
    ("skewness", "insertions", "Skewness", "high", "rule_of_thumb", "Skewness of the read counts at occupied sites: above 30 often signals a problem, below 50 is desired (can be higher under strong selection)."),
    ("kurtosis", "insertions", "Kurtosis", None, None, "Excess kurtosis of the read counts at occupied sites; not very meaningful on its own."),
    ("pickands_tail_index", "insertions", "Tail index (Pickands)", "high", "relative", "Heaviness of the upper tail of the read counts; higher means a few sites carry more of the reads."),
    ("fr_corr", "insertions", "Strand correlation (W vs C)", None, None, "Correlation of read counts on the two strands at a position. For Hermes and Ac/Ds the two orientations are separate insertions; for piggyBac they are the two ends of one occupied TTAA. Shown, not flagged."),
    ("sites_per_coding_kb", "genes", "Insertion sites per coding kb (median gene)", "low", "relative", "Distinct insertion positions inside a gene's coding blocks per kb of coding sequence, median over genes. This is the density the essentiality recommendation is stated in (docs/methods/00_recommendations.md): about 36 for the label-free HMM, about 18 when curated labels are available, reached by pooling independent inductions. Insertions per kb of genome is not the same measure, because it counts the whole reference rather than coding sequence."),
    ("percent_features_hit", "genes", "Genes hit (%)", "low", "relative", "SummaryTable: features with at least one insertion."),
    ("mean_hits_per_feature", "genes", "Insertions per gene", "low", "relative", "SummaryTable: mean insertion sites per feature."),
    ("percent_hits_in_features", "genes", "Insertions in genes (%)", "both", "relative", "SummaryTable: share of insertion sites inside features."),
    ("neighborhood_spearman", "genes", "Insertion-neighbourhood correlation", "low", "relative", "SummaryTable: Spearman correlation between a gene's insertions and its neighbourhood's."),
    ("one_site_genes", "genes", "Genes led by one site (%)", "high", "relative", "Of the genes with at least 3 insertion sites, the share in which the largest site holds at least half of the gene's reads: a jackpot, an insertion that does not disrupt the gene or a selected site, rather than many independent mutants. High values make gene totals depend on single sites."),
    ("midlc_est", "saturation", "MidLC estimate", None, None, "Library diagnostics: reads at which half the discoverable sites are found."),
    ("depth_ratio", "saturation", "Depth / MidLC", None, None, "Reads at insertion sites divided by the MidLC estimate (above 1: sequenced past half-saturation)."),
    ("jackpot_frac_reads", "saturation", "Reads in jackpot sites", "high", "relative", "Library diagnostics: fraction of reads in the top jackpot sites."),
    ("n_jackpot_sites", "saturation", "Jackpot sites", "high", "relative", "Library diagnostics: number of sites flagged as jackpots."),
]
TOP_SITE_HINT = "Consider a top-site treatment before normalizing (Normalization → Top sites: winsorize or drop the top site); see Genes with the most reads."
HINTS = {"skewness":BGC_HINT, "max_share":TTR_HINT, "top1pct_share":BGC_HINT, "jackpot_frac_reads":TTR_HINT,
         "pickands_tail_index":BGC_HINT, "one_site_genes":TOP_SITE_HINT}
ONE_SITE_MIN_SITES, ONE_SITE_SHARE = 3, .5
THEMES = {"reads":"Reads & Mapping", "insertions":"Insertions", "genes":"Genes", "saturation":"Saturation"}


def _float(value):
    try:
        value = float(str(value).strip().rstrip("%"))
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def _resolve(value, root: Path) -> Path | None:
    if not value:return None
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else root/path


def genome_length(fasta: Path) -> int:
    fai = Path(str(fasta)+".fai")
    if fai.is_file():return sum(int(line.split("\t")[1]) for line in open(fai) if line.strip())
    total = 0
    with open(fasta) as handle:
        for line in handle:
            if not line.startswith(">"):total += len(line.strip())
    return total


def coding_genes(fasta: Path, gff: Path) -> list[dict]:
    """Annotated genes with CDS blocks, on the chromosomes the reads were mapped to."""
    from ytab.essentiality.control_conditioned_gumbel import annotations
    chroms = {line[1:].split()[0] for line in open(fasta) if line.startswith(">")} if fasta and fasta.is_file() else set()
    return [g for g in annotations(gff) if g["blocks"] and (not chroms or g["chrom"] in chroms)]


TOP_GENES = 30


def gene_reads(sample: str, sites: dict, genes: list[dict]) -> list[dict]:
    """Every gene with reads in one library (reads at sites in the gene's coding sequence), with its share of the
    library's reads, its insertion sites, its largest site's reads and that site's share of the gene's reads."""
    from ytab.essentiality.control_conditioned_gumbel import gene_sites
    from ytab.essentiality.library_correlation import SITE_RULES, gene_key
    if not genes or not sites:return []
    keys = list(sites); reads = np.fromiter((sites[k] for k in keys), dtype=float, count=len(keys)); total = float(reads.sum())
    rows = []
    for g, idx in zip(genes, gene_sites(genes, keys, SITE_RULES)):
        if not len(idx):continue
        r = reads[idx]; gene_reads = float(r.sum()); top = int(np.argmax(r))
        rows.append(dict(sample=sample, feature_id=gene_key(g["id"]), name=g.get("name") or "", chrom=g["chrom"], reads=gene_reads,
                         library_share=gene_reads/total if total else None, sites=int(len(idx)), top_site=f"{keys[idx[top]][0]}:{keys[idx[top]][1]}",
                         top_site_reads=float(r[top]), top_site_share=float(r[top])/gene_reads if gene_reads else None,
                         reads_without_top_site=gene_reads-float(r[top])))
    return rows


def top_genes(sample: str, sites: dict, genes: list[dict], n: int = TOP_GENES, rows: list[dict] | None = None) -> list[dict]:
    """The n genes holding the most reads of one library (gene_reads rows, largest first, ranked). A gene whose reads
    come mostly from one site is a jackpot, a selection marker or the donor locus rather than a gene many independent
    mutants support."""
    rows = sorted(gene_reads(sample, sites, genes) if rows is None else rows, key=lambda x:-x["reads"])
    return [dict(x, rank=i+1) for i, x in enumerate(rows[:n])]


def one_site_genes(rows: list[dict]) -> float | None:
    """Percent of the genes with >= 3 insertion sites whose largest site holds >= half of the gene's reads."""
    judged = [r for r in rows if r["sites"] >= ONE_SITE_MIN_SITES and r["reads"] > 0]
    return 100*sum(r["top_site_share"] >= ONE_SITE_SHARE for r in judged)/len(judged) if judged else None


TOP_GENE_COLUMNS = ["sample", "rank", "feature_id", "name", "chrom", "reads", "library_share", "sites", "top_site", "top_site_reads", "top_site_share",
                    "reads_without_top_site"]


def refresh_top_genes(run: Path) -> Path:
    """Add plots/top_genes.tsv.gz and the "Genes led by one site" metric (table column, flags, metrics list, Excel) to
    an earlier QC run, from the hit files and GFF of its project config."""
    from ytab.pipeline.mapfastq_runner import load_project_for_mapping
    run = Path(run).resolve(); manifest = json.loads((run/"manifest.json").read_text())
    config = load_project_for_mapping(Path(manifest["project_config"])); root = Path(config["_repo_root"])
    if not (root/"output").is_dir():root = Path(__file__).resolve().parents[3]
    project = _resolve(config.get("output_project_dir"), root)
    fasta = _resolve((config.get("reference") or {}).get("fasta"), root); gff = _resolve((config.get("reference") or {}).get("gff"), root)
    genes = coding_genes(fasta, gff) if gff and gff.is_file() else []
    rows, share = [], {}
    for s in config.get("samples", []):
        if str(s.get("include", True)).lower() not in {"true", "1", "yes"}:continue
        name = str(s["sample"]); path = _resolve(s.get("hit_file"), root) if s.get("hit_file") else project/"create_hit_file"/name/f"{name}_hits.txt"
        if path and path.is_file():
            per_gene = gene_reads(name, read_hits(path)[0], genes); rows += top_genes(name, None, genes, rows=per_gene); share[name] = one_site_genes(per_gene)
    out = run/"plots/top_genes.tsv.gz"; _write(out, TOP_GENE_COLUMNS, rows)
    table = list(csv.DictReader(open(run/"qc_table.tsv"), delimiter="\t"))
    if table:
        key, theme, label, direction, rule, comment = next(m for m in METRICS if m[0] == "one_site_genes")
        for r in table:r[key] = share.get(r["sample"])
        columns = list(table[0]) if key in table[0] else list(table[0])+[key]
        _write(run/"qc_table.tsv", columns, table)
        flags = [f for f in csv.DictReader(open(run/"qc_flags.tsv"), delimiter="\t") if f["metric"] != key]
        values = {r["sample"]:_float(r.get(key)) for r in table}; roles = {r["sample"]:r["role"] for r in table}
        for s, (colour, reason) in flag(values, roles, key, direction, rule).items():
            flags.append(dict(sample=s, metric=key, label=label, theme=THEMES[theme], value=values[s], flag=colour, reason=reason, source=LEYS,
                              comment=comment, normalization_hint=HINTS.get(key, "") if colour in {"yellow", "red"} else ""))
        _write(run/"qc_flags.tsv", list(flags[0]), flags)
        _write(run/"metrics.tsv", ["metric", "theme", "label", "bad_direction", "rule", "comment"],
               [dict(metric=k, theme=THEMES[t], label=l, bad_direction=d or "", rule=r or "descriptive", comment=c) for k, t, l, d, r, c in METRICS])
        write_xlsx(run/"qc_table.xlsx", table, flags)
        manifest["flag_counts"] = {c:sum(f["flag"] == c for f in flags) for c in ("red", "yellow", "green", "grey")}
    manifest.setdefault("files", {})["top_genes"] = "plots/top_genes.tsv.gz"; (run/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return out


def coding_density(sites, genes: list[dict]) -> float | None:
    """Median over genes of distinct insertion positions per kb of coding sequence.

    `sites` is any collection of (chromosome, position) keys, so it takes both one library's site counts and the
    union of several libraries' sites.

    The measure the essentiality recommendation is written in, so that the QC table can be read against it
    directly. It is the median gene rather than the mean because a handful of very densely hit genes would
    otherwise make a library look usable when most of its genes carry nothing.
    """
    if not genes:return None
    positions = {}
    for chrom, position in sites:positions.setdefault(chrom, []).append(position)
    arrays = {c:np.array(sorted(p)) for c, p in positions.items()}
    per_kb = []
    for gene in genes:
        found = arrays.get(gene["chrom"]); length = sum(end-start+1 for start, end in gene["blocks"])
        if not length:continue
        n = sum(int(np.searchsorted(found, end, "right")-np.searchsorted(found, start, "left"))
                for start, end in gene["blocks"]) if found is not None else 0
        per_kb.append(n/length*1000)
    return float(np.median(per_kb)) if per_kb else None


# How dense a library set has to be before an essentiality call is worth making (docs/methods/00_recommendations.md).
HMM_SITES_PER_CODING_KB = 36.
HYBRID_SITES_PER_CODING_KB = 18.


def essentiality_readiness(rows: list[dict], hits: dict, genes: list[dict]) -> list[dict]:
    """Per strain, the density an essentiality run would actually have, and whether that clears the recommendation.

    An HMM is fitted on the unselected (parent, else control) libraries of one strain together, so the number that
    decides whether the call is worth making is their pooled density, not any single library's: pooling independent
    inductions is how a sparse set reaches a usable density, and summing the per-library figures would overstate it
    because the sets of occupied sites overlap. The app reads this rather than recomputing it.
    """
    out = []
    for background in dict.fromkeys(r["background"] for r in rows):
        same = [r for r in rows if r["background"] == background]
        parents = [r["sample"] for r in same if r["role"] == "parent"]
        unselected = parents or [r["sample"] for r in same if r["role"] == "control"]
        usable = [s for s in unselected if s in hits]
        pooled = set().union(*(set(hits[s][0]) for s in usable)) if usable else set()
        density = coding_density(pooled, genes) if usable else None
        out.append(dict(background=background, role="parent" if parents else "control", libraries=" ".join(usable),
                        n_libraries=len(usable), sites_per_coding_kb=density,
                        supports_hmm=density is not None and density >= HMM_SITES_PER_CODING_KB,
                        supports_hybrid=density is not None and density >= HYBRID_SITES_PER_CODING_KB,
                        hmm_target=HMM_SITES_PER_CODING_KB, hybrid_target=HYBRID_SITES_PER_CODING_KB))
    return out


def read_hits(path: Path) -> tuple[dict, dict]:
    """Site -> reads (strands combined) and site -> (W, C) reads, from a CreateHitFile hit file."""
    sites, strands = {}, {}
    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t"); header = next(reader)
        chrom, source, pos, count = (header.index(k) for k in ("Chromosome", "Source", "Hit position", "Hit count"))
        for row in reader:
            key = (row[chrom], int(row[pos])); n = int(row[count])
            sites[key] = sites.get(key, 0)+n
            w, c = strands.get(key, (0, 0)); strands[key] = (w+n, c) if row[source] == "W" else (w, c+n)
    return sites, strands


def insertion_metrics(sites: dict, strands: dict, universe: int, genome_bp: int) -> dict:
    counts = np.fromiter(sites.values(), dtype=float)
    out = statistics(counts, universe) if len(counts) else {}
    total = counts.sum() if len(counts) else 0.
    ordered = np.sort(counts)[::-1]
    out.update(unique_sites=len(counts), insertions_per_kb=len(counts)/(genome_bp/1000) if genome_bp else None,
               max_share=float(ordered[0]/total) if total else None,
               top1pct_share=float(ordered[:max(1, len(ordered)//100)].sum()/total) if total else None)
    if out.get("max_ct") is not None:
        out["max_site"] = "%s:%d" % max(sites, key=sites.get)
    w, c = (np.array(x, dtype=float) for x in zip(*strands.values())) if strands else (np.array([]), np.array([]))
    out["fr_corr"] = float(np.corrcoef(w, c)[0, 1]) if len(w) > 2 and w.std() > 0 and c.std() > 0 else None
    return out


def plot_data(sample: str, counts: np.ndarray, n_points: int = 400) -> dict[str, list]:
    """Count distribution plots: positive-count histogram (log bins), QQ against a geometric fit, ranked counts, and the
    Lorenz curve (share of reads held by the top fraction of occupied sites; log-spaced fractions plus the top site)."""
    x = np.sort(counts[counts > 0])
    if not len(x):return {"histogram":[], "qq":[], "ranked":[], "lorenz":[]}
    edges = np.unique(np.floor(np.logspace(0, math.log10(x.max()+1), 40)).astype(int))
    hist, _ = np.histogram(x, bins=np.append(edges, x.max()+1))
    histogram = [dict(sample=sample, bin_start=int(a), bin_end=int(b)-1, sites=int(h)) for a, b, h in zip(edges, np.append(edges[1:], x.max()+1), hist)]
    # Geometric on {1, 2, ...} with the sample's NZmean (method of moments): p = 1/NZmean.
    p = 1/x.mean(); probs = (np.arange(n_points)+.5)/n_points
    observed = np.quantile(x, probs); expected = np.maximum(1, np.ceil(np.log1p(-probs)/np.log1p(-p))) if p < 1 else np.ones(n_points)
    qq = [dict(sample=sample, probability=round(float(q), 5), observed=float(o), geometric=float(e)) for q, o, e in zip(probs, observed, expected)]
    ranks = np.unique(np.round(np.logspace(0, math.log10(len(x)), n_points)).astype(int))
    ranked = [dict(sample=sample, rank=int(r), reads=float(x[::-1][r-1])) for r in ranks]
    cum = np.cumsum(x[::-1])/x.sum(); tops = np.unique(np.maximum(1, np.round(np.geomspace(1, len(x), n_points)).astype(int)))
    lorenz = [dict(sample=sample, fraction_sites=round(int(k)/len(x), 8), share_reads=round(float(cum[k-1]), 6)) for k in tops]
    return {"histogram":histogram, "qq":qq, "ranked":ranked, "lorenz":lorenz}


def _cleanup(path: Path) -> dict:
    if not path.is_file():return {}
    data = json.loads(path.read_text()); mate = (data.get("mates") or [{}])[0]
    reads = _float(mate.get("input_reads"))
    return dict(input_reads=reads, retained_fraction=(_float(data.get("retained_reads_or_pairs")) or 0)/reads if reads else None,
                adapter_fraction=(_float(mate.get("adapter_3p_reads")) or 0)/reads if reads else None)


def _mapping_log(path: Path) -> dict:
    """Read counts from the mapping log when there is no preprocessing JSON (runs made on the cluster): reads given to
    the first step (transposon or adapter removal, or alignment), reads that reached the alignment, reads with the
    3' adaptor."""
    if not path.is_file():return {}
    text = path.read_text(errors="replace")
    counts = [int(n) for n in re.findall(r"^(\d+) reads; of these:", text, flags=re.M)]
    if not counts:return {}
    align = re.search(r"=== Sequence alignment ===\s*\n(\d+) reads; of these:", text)
    adapter = re.search(r"=== Adapter removal ===\s*\n\d+ reads; of these:\s*\n\s*(\d+) \([\d.]+%\) contained the adapter", text)
    reads = float(counts[0])
    return dict(input_reads=reads, retained_fraction=int(align.group(1))/reads if align else None,
                adapter_fraction=int(adapter.group(1))/reads if adapter else None)


def _mapping(path: Path) -> dict:
    if not path.is_file():return {}
    row = next(csv.DictReader(open(path)), {})
    return dict(mapped_reads=_float(row.get("primary_mapped")), total_records=_float(row.get("total_records")), percent_mapped=_float(row.get("percent_mapped")), percent_mapq=_float(row.get("percent_mapq_ge_threshold")),
                mean_mapq=_float(row.get("avg_mapq_mapped_primary")))


def _donor(path: Path) -> dict:
    """Donor reads (ytab.mapping.donor_reads): counted against the donor sequence when one was given, otherwise the
    recurring-unmapped-flank estimate; `donor_source` says which."""
    if not path.is_file():return {}
    row = next(csv.DictReader(open(path)), {})
    counted = _float(row.get("percent_donor"))
    return dict(percent_donor=counted if counted is not None else _float(row.get("percent_recurrent_unmapped_flank")),
                percent_donor_shared=_float(row.get("percent_donor_shared")),
                donor_source=row.get("donor") or ("estimate (recurring unmapped flanks)" if row.get("percent_recurrent_unmapped_flank") else ""))


def _summary(project: Path, sample: str, combined: dict) -> dict:
    row = combined.get(sample)
    if row is None:
        path = project/"summary"/sample/"stats.csv"
        if path.is_file():
            raw = next(csv.DictReader(open(path)), {})
            row = {"percent_features_hit":raw.get("% of features hit"), "mean_hits_per_feature":raw.get("Mean hits per feature"),
                   "percent_hits_in_features":raw.get("% of hits in features")}
    out = {k:_float((row or {}).get(k)) for k in ("percent_features_hit", "mean_hits_per_feature", "percent_hits_in_features")}
    corr = project/"summary"/sample/"insertion_vs_neighborhood_correlations.txt"
    match = re.search(r"Spearman: \S+\(statistic=(?:np\.float64\()?([-0-9.eE]+)", corr.read_text()) if corr.is_file() else None
    out["neighborhood_spearman"] = float(match.group(1)) if match else None
    return out


def _diagnostics(project: Path) -> dict:
    """Saved library diagnostics per sample: the project-level summary, else the newest run's."""
    root = project/"library_diagnostics"
    candidates = [root/"library_diagnostics.summary.csv", *sorted(root.glob("runs/*/library_diagnostics.summary.csv"), key=lambda p:p.stat().st_mtime, reverse=True)]
    rows = {}
    for path in candidates:
        if path.is_file():
            for row in csv.DictReader(open(path)):rows.setdefault(row["sample"], row)
    return {s:dict(midlc_est=_float(r.get("midlc_est")), depth_ratio=_float(r.get("depth_ratio_R_over_midlc")),
                   jackpot_frac_reads=_float(r.get("jackpot_frac_reads")), n_jackpot_sites=_float(r.get("n_jackpot_sites")))
            for s, r in rows.items()}


def flag(values: dict[str, float | None], roles: dict[str, str], key: str, direction, rule) -> dict[str, tuple[str, str]]:
    """sample -> (colour, reason) for one metric."""
    out = {}
    if rule == "rule_of_thumb":  # skewness
        for s, v in values.items():
            out[s] = ("grey", "not available") if v is None else ("red", "above 50 (below 50 is desired)") if v > 50 else \
                     ("yellow", "above 30 (often signals a problem)") if v > 30 else ("green", "within the rule of thumb (30)")
        return out
    if rule != "relative":
        return {s:("grey", "descriptive; not flagged") for s in values}
    for role in set(roles.values()):
        members = [s for s in values if roles[s] == role]
        known = np.array([values[s] for s in members if values[s] is not None], dtype=float)
        median = float(np.median(known)) if len(known) else None
        scale = max(1.4826*float(np.median(np.abs(known-median))), MIN_RELATIVE_SPREAD*abs(median)) if len(known) else 0.
        for s in members:
            v = values[s]
            if v is None:out[s] = ("grey", "not available"); continue
            if len(known) < MIN_ROLE_LIBRARIES:out[s] = ("grey", f"fewer than {MIN_ROLE_LIBRARIES} {role} libraries to compare"); continue
            if scale == 0:out[s] = ("green", f"{role} libraries all zero"); continue
            z = (v-median)/scale; bad = -z if direction == "low" else z if direction == "high" else abs(z)
            where = "below" if z < 0 else "above"
            out[s] = ("red" if bad > ROBUST_RED else "yellow" if bad > ROBUST_YELLOW else "green",
                      f"robust z {z:+.1f} vs {len(known)} {role} libraries (median {median:.4g})" + ("" if bad <= ROBUST_YELLOW else f"; {where} the others"))
    return out


def build(project_config: Path, output: Path) -> dict:
    from ytab.pipeline.mapfastq_runner import load_project_for_mapping
    config_path = Path(project_config).resolve(); config = load_project_for_mapping(config_path)
    root = Path(config["_repo_root"])
    if not (root/"output").is_dir():root = Path(__file__).resolve().parents[3]  # a config written on another machine
    project = _resolve(config.get("output_project_dir"), root)
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):raise ValueError("Choose a new, empty output directory")
    samples = [s for s in config.get("samples", []) if str(s.get("include", True)).lower() in {"true", "1", "yes"}]
    if not samples:raise ValueError("No included samples")
    fasta = _resolve((config.get("reference") or {}).get("fasta"), root)
    genome_bp = genome_length(fasta) if fasta and fasta.is_file() else 0
    gff = _resolve((config.get("reference") or {}).get("gff"), root)
    genes = coding_genes(fasta, gff) if gff and gff.is_file() else []
    combined = {}
    combined_path = project/"summary/summary_stats.all_samples.csv"
    if combined_path.is_file():combined = {r["sample"]:r for r in csv.DictReader(open(combined_path))}
    diagnostics = _diagnostics(project)
    hits, missing = {}, []
    for s in samples:
        name = str(s["sample"]); path = _resolve(s.get("hit_file"), root) if s.get("hit_file") else project/"create_hit_file"/name/f"{name}_hits.txt"
        if path and path.is_file():hits[name] = read_hits(path)
        else:missing.append(name)
    universe = len(set().union(*(set(h[0]) for h in hits.values()))) if hits else 0
    rows, plots, top = [], {"histogram":[], "qq":[], "ranked":[], "lorenz":[]}, []
    for s in samples:
        name = str(s["sample"]); role = str(s.get("library_role") or s.get("condition") or s.get("guessed_condition") or "library").strip().lower()
        row = dict(sample=name, display_name=s.get("display_name") or name, role=role, pool=s.get("pool") or s.get("guessed_pool") or "",
                   background=s.get("background") or s.get("guessed_background") or "")
        # The cleanup report was named after Hermes before the step was understood to apply to any transposon;
        # runs made then still carry that name.
        cleanup = _cleanup(project/"mapfastq"/name/f"{name}.read_cleanup.json") \
            or _cleanup(project/"mapfastq"/name/f"{name}.hermes_preprocessing.json")
        row |= (cleanup or _mapping_log(project/"mapfastq"/name/f"{name}_log.txt")) \
            | _mapping(project/"mapfastq"/name/f"{name}.mapping_stats.csv") | _donor(project/"mapfastq"/name/f"{name}.donor_reads.csv")
        if name in hits:
            row |= insertion_metrics(*hits[name], universe, genome_bp)
            row["sites_per_coding_kb"] = coding_density(hits[name][0], genes)
            data = plot_data(name, np.fromiter(hits[name][0].values(), dtype=float))
            for k in plots:plots[k] += data[k]
            per_gene = gene_reads(name, hits[name][0], genes)
            top += top_genes(name, hits[name][0], genes, rows=per_gene); row["one_site_genes"] = one_site_genes(per_gene)
        row |= _summary(project, name, combined) | diagnostics.get(name, {})
        rows.append(row)
    flags = []
    roles = {r["sample"]:r["role"] for r in rows}
    for key, theme, label, direction, rule, comment in METRICS:
        values = {r["sample"]:_float(r.get(key)) for r in rows}
        for s, (colour, reason) in flag(values, roles, key, direction, rule).items():
            flags.append(dict(sample=s, metric=key, label=label, theme=THEMES[theme], value=values[s], flag=colour, reason=reason,
                              source=RULE_OF_THUMB if rule == "rule_of_thumb" else LEYS if rule == "relative" else "descriptive",
                              comment=comment, normalization_hint=HINTS.get(key, "") if colour in {"yellow", "red"} else ""))
    readiness = essentiality_readiness(rows, hits, genes)
    output.mkdir(parents=True)
    keys = ["sample", "display_name", "role", "pool", "background"]+[m[0] for m in METRICS]+["max_site", "pti_valid_terms", "catalogue_sites", "occupied_sites", "stats_notes"]
    _write(output/"qc_table.tsv", keys, rows)
    _write(output/"essentiality_readiness.tsv", list(readiness[0]) if readiness else ["background"], readiness)
    _write(output/"qc_flags.tsv", list(flags[0]), flags)
    _write(output/"metrics.tsv", ["metric", "theme", "label", "bad_direction", "rule", "comment"],
           [dict(metric=k, theme=THEMES[t], label=l, bad_direction=d or "", rule=r or "descriptive", comment=c) for k, t, l, d, r, c in METRICS])
    (output/"plots").mkdir()
    for name, data in plots.items():
        _write(output/"plots"/f"{name}.tsv.gz", list(data[0]) if data else ["sample"], data)
    _write(output/"plots/top_genes.tsv.gz", TOP_GENE_COLUMNS, top)
    write_xlsx(output/"qc_table.xlsx", rows, flags)
    manifest = dict(status="success", project=config.get("project_id"), project_config=str(config_path), libraries=len(rows),
                    libraries_without_hits=missing, catalogue_sites=universe, genome_bp=genome_bp,
                    flag_counts={c:sum(f["flag"] == c for f in flags) for c in ("red", "yellow", "green", "grey")},
                    rules=dict(rule_of_thumb=RULE_OF_THUMB, relative=LEYS, yellow=ROBUST_YELLOW, red=ROBUST_RED, min_role_libraries=MIN_ROLE_LIBRARIES, min_relative_spread=MIN_RELATIVE_SPREAD,
                               comparison="within role (parents with parents, treated with treated)"),
                    created=datetime.now(timezone.utc).isoformat(),
                    files=dict(table="qc_table.tsv", flags="qc_flags.tsv", metrics="metrics.tsv", excel="qc_table.xlsx",
                               histogram="plots/histogram.tsv.gz", qq="plots/qq.tsv.gz", ranked="plots/ranked.tsv.gz",
                               lorenz="plots/lorenz.tsv.gz", top_genes="plots/top_genes.tsv.gz", essentiality_readiness="essentiality_readiness.tsv"))
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def _write(path: Path, columns: list[str], rows: list[dict]) -> None:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="ignore", lineterminator="\n"); writer.writeheader()
        for row in rows:writer.writerow({k:("" if row.get(k) is None else row.get(k)) for k in columns})


def write_xlsx(path: Path, rows: list[dict], flags: list[dict]) -> None:
    """One row per library, one column per metric; cells coloured by flag, with the flag's reason and comment as a note."""
    from openpyxl import Workbook
    from openpyxl.comments import Comment
    from openpyxl.styles import Font, PatternFill
    fills = {"red":"F4C7C3", "yellow":"FCE8B2", "green":"B7E1CD"}
    lookup = {(f["sample"], f["metric"]):f for f in flags}
    book = Workbook(); sheet = book.active; sheet.title = "QC"
    header = ["Library", "Role"]+[m[2] for m in METRICS]
    sheet.append(header)
    for cell in sheet[1]:cell.font = Font(bold=True)
    for j, (key, *_rest) in enumerate(METRICS, start=3):
        sheet.cell(1, j).comment = Comment(f"{_rest[4]}", "YTAB")
    for i, row in enumerate(rows, start=2):
        sheet.cell(i, 1, row["display_name"]); sheet.cell(i, 2, row["role"])
        for j, (key, *_rest) in enumerate(METRICS, start=3):
            f = lookup[(row["sample"], key)]; cell = sheet.cell(i, j, f["value"])
            if f["flag"] in fills:cell.fill = PatternFill("solid", fgColor=fills[f["flag"]])
            note = f"{f['flag']}: {f['reason']}" + (f"\n{f['normalization_hint']}" if f["normalization_hint"] else "")
            cell.comment = Comment(note, "YTAB")
    legend = book.create_sheet("Read me")
    for line in (["Colours", "green: within the cutoff; yellow: worth a look; red: outside the cutoff; no colour: not flagged"],
                 ["Rule of thumb", RULE_OF_THUMB+": above 30 yellow, above 50 red"],
                 ["Project-relative", LEYS+f"; compared within role; needs {MIN_ROLE_LIBRARIES}+ libraries of the role; spread at least {MIN_RELATIVE_SPREAD:.0%} of the median"],
                 ["Per cell", "hover a cell for its flag reason and, where one applies, a normalization hint"]):
        legend.append(line)
    flat = book.create_sheet("Flags")
    flat.append(list(flags[0]))
    for f in flags:flat.append([f[k] for k in flags[0]])
    book.save(path)
