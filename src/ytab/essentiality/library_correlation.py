"""Correlation of Hermes libraries across pools, strains, conditions and species, generalized from TRANSIT's corrplot (DeJesus et al. 2015).

Each library is summarized per gene by insertion density: distinct insertion sites (>= min_reads reads) in
the CDS (last 3 bp excluded) per kb. Optionally "reads": reads per kb. Libraries are
compared by Spearman correlation across genes:
- same species: genes matched by ID;
- different species: genes matched through 1:1 ortholog groups (a table of group, species, gene; a group
  is used only if it has exactly one gene in each of the two species).
Only genes on chromosomes carrying insertions in both libraries are compared. Genes on chromosomes absent
from the alignment reference (e.g. the B haplotype of a diploid C. albicans annotation when reads are
mapped to HapA) would be zero in every library and inflate the correlation.
Spearman is invariant to any per-library monotone scaling, so no normalization is needed and libraries of
different depth compare directly. Insertion density is the default because it reflects essentiality and
target accessibility, which are shared between independent libraries. Reads mostly reflect clone
abundance, which is specific to one pool.

Every pair is labelled from its metadata:
- same pool (shared clones);
- same strain, different pool (independent libraries);
- different strain;
- different species.
Samples split from one pool share their clones, so their correlation measures library structure
rather than biology; such pairs are reported but are not evidence of biological similarity.

Clone-sharing evidence (same-species pairs): site overlap (Jaccard) depends on sequencing depth and
complexity, so the primary check is the Spearman correlation of log reads at sites detected in both
libraries. Split cultures share founder abundances; independent pools share only hotspot positions.
Calibration on Gale et al. (2020):
- same pool: 0.89-0.99;
- independent pools: 0.42-0.43;
- different strains: 0.33.
Labels: >= 0.7 "shared clones", <= 0.5 "no shared clone abundance", otherwise "ambiguous". A pair whose
evidence contradicts its declared relation is flagged. The thresholds come from one dataset and are
indicative.
"""
from dataclasses import asdict, dataclass
import csv
import gzip
import json
from pathlib import Path
import re
from types import SimpleNamespace

import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.stats import spearmanr

from ytab.essentiality.control_conditioned_gumbel import annotations, gene_sites

REVISION = "library_correlation_v1"
SHARED, INDEPENDENT = .7, .5   # shared-site read correlation thresholds (Gale calibration, indicative)
SITE_RULES = SimpleNamespace(nterm=0., cterm=0., ignore_stop=True)
PAIR_TYPES = ("same pool (shared clones)", "same strain, different pool", "different strain", "different species", "unknown")


@dataclass
class Sample:
    name: str
    hits: str
    species: str
    gff: str
    strain: str = ""
    pool: str = ""
    condition: str = ""
    source: str = ""


def gene_key(gene_id):
    """Comparable gene ID: annotation prefixes such as 'gene:' (Ensembl/PomBase GFF3) removed."""
    return re.sub(r"^(gene[:-])", "", gene_id)


def read_hits(path, min_reads=1):
    """{chromosome: (sorted positions, reads)} from a YTAB *_hits.txt or a table with chromosome/position/reads columns.
    Counts at one position (both strands) are summed."""
    path = Path(path); text = path.read_text().splitlines()
    delimiter = "\t" if "\t" in text[0] else ","
    rows = csv.DictReader(text, delimiter=delimiter); fields = {f.lower().strip():f for f in rows.fieldnames}
    pick = lambda *names:next((fields[n] for n in names if n in fields), None)
    chrom, pos, reads = pick("chromosome", "chrom", "chr"), pick("hit position", "position", "pos"), pick("hit count", "reads", "count", "counts")
    if not (chrom and pos and reads):raise ValueError(f"{path}: need chromosome, position and read-count columns (YTAB hit file or chromosome/position/reads)")
    totals = {}
    for r in rows:
        key = (r[chrom], int(float(r[pos]))); totals[key] = totals.get(key, 0)+int(float(r[reads]))
    out = {}
    for (c, p), n in totals.items():
        if n >= min_reads:out.setdefault(c, []).append((p, n))
    return {c:(np.array([p for p, _ in v]), np.array([n for _, n in v], dtype=float)) for c, v in ((c, sorted(v)) for c, v in out.items())}


def gene_metrics(genes, hits):
    """Per gene: CDS length, insertion sites, reads, density (sites per kb) and reads per kb."""
    sites = [(c, int(p)) for c, (ps, _) in sorted(hits.items()) for p in ps]
    reads = np.concatenate([r for _, (_, r) in sorted(hits.items())]) if hits else np.zeros(0)
    idx = gene_sites(genes, sites, SITE_RULES)
    rows = {}
    for g, i in zip(genes, idx):
        length = sum(b-a+1 for a, b in g["blocks"])
        kb = max(length-3, 1)/1000
        rows[gene_key(g["id"])] = dict(chrom=g["chrom"], length=length, sites=len(i), reads=float(reads[i].sum()) if len(i) else 0., density=len(i)/kb,
                                      reads_per_kb=float(reads[i].sum())/kb if len(i) else 0.)
    return rows


def ortholog_groups(path):
    """{species: {gene: group}} from a TSV with columns group, species, gene (gene IDs as in gene_key)."""
    groups = {}
    with Path(path).open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            groups.setdefault(r["species"], {}).setdefault(gene_key(r["gene"]), set()).add(r["group"])
    # keep genes in exactly one group, and groups with exactly one gene per species
    out = {}
    for species, genes in groups.items():
        single = {g:next(iter(s)) for g, s in genes.items() if len(s) == 1}
        counts = {}
        for grp in single.values():counts[grp] = counts.get(grp, 0)+1
        out[species] = {g:grp for g, grp in single.items() if counts[grp] == 1}
    return out


def pair_type(a, b):
    if a.species != b.species:return "different species"
    if not a.strain or not b.strain:return "unknown"
    if a.strain != b.strain:return "different strain"
    if not a.pool or not b.pool:return "unknown"
    return "same pool (shared clones)" if a.pool == b.pool else "same strain, different pool"


def matched(a, b, ma, mb, orthologs, metric, covered_a=None, covered_b=None):
    """Value vectors for genes comparable between two libraries (on chromosomes covered by each library's insertions)."""
    if covered_a is not None:ma = {k:v for k, v in ma.items() if v["chrom"] in covered_a}
    if covered_b is not None:mb = {k:v for k, v in mb.items() if v["chrom"] in covered_b}
    if a.species == b.species:
        keys = sorted(set(ma) & set(mb)); return np.array([ma[k][metric] for k in keys]), np.array([mb[k][metric] for k in keys])
    if orthologs is None or a.species not in orthologs or b.species not in orthologs:return np.zeros(0), np.zeros(0)
    ga, gb = orthologs[a.species], orthologs[b.species]
    by_group_b = {grp:g for g, grp in gb.items() if g in mb}
    keys = [(g, by_group_b[grp]) for g, grp in ga.items() if g in ma and grp in by_group_b]
    return np.array([ma[x][metric] for x, _ in keys]), np.array([mb[y][metric] for _, y in keys])


def site_table(hits, chromosomes):
    """Sorted integer site keys (chromosome index << 32 | position) and reads; `chromosomes` is a shared name -> index map."""
    keys, reads = [], []
    for c, (ps, rs) in hits.items():
        code = chromosomes.setdefault(c, len(chromosomes))
        keys.append((np.int64(code) << 32) | ps.astype(np.int64)); reads.append(rs)
    keys = np.concatenate(keys) if keys else np.zeros(0, np.int64); reads = np.concatenate(reads) if reads else np.zeros(0)
    o = np.argsort(keys); return keys[o], reads[o]


def clone_evidence(ta, tb):
    """Site overlap (Jaccard), shared sites, and Spearman of log reads at shared sites, with an evidence label."""
    (ka, ra), (kb, rb) = ta, tb
    common, ia, ib = np.intersect1d(ka, kb, assume_unique=True, return_indices=True)
    union = len(ka)+len(kb)-len(common)
    rho = float(spearmanr(ra[ia], rb[ib])[0]) if len(common) >= 100 else None
    label = None if rho is None else "shared clones" if rho >= SHARED else "no shared clone abundance" if rho <= INDEPENDENT else "ambiguous"
    return dict(shared_sites=int(len(common)), site_jaccard=len(common)/union if union else None, shared_site_read_rho=rho, clone_evidence=label)


def conflict(pair, evidence):
    if evidence is None:return ""
    if pair == "same pool (shared clones)" and evidence == "no shared clone abundance":
        return "declared same pool, but read counts at shared sites do not correlate"
    if pair in ("same strain, different pool", "different strain") and evidence == "shared clones":
        return "declared independent libraries, but read counts at shared sites correlate as in split cultures"
    return ""


def correlate(samples, metric="density", orthologs=None, boot=500, seed=2026, min_reads=1, min_genes=100, gene_values=None):
    """Pairwise Spearman correlations with gene-bootstrap 95% CIs. Returns (pair rows, matrix, per-sample rows).
    A dict passed as gene_values receives each library's per-gene metrics (for gene-level plots)."""
    if metric not in ("density", "reads_per_kb"):raise ValueError("metric must be 'density' or 'reads_per_kb'")
    if len({s.name for s in samples}) != len(samples):raise ValueError("Sample names must be unique")
    annotation = {}
    for s in samples:
        if s.gff not in annotation:annotation[s.gff] = annotations(Path(s.gff).resolve())
    hits = {s.name:read_hits(s.hits, min_reads) for s in samples}
    chromosomes = {}; tables = {s.name:site_table(hits[s.name], chromosomes) for s in samples}
    metrics = {s.name:gene_metrics(annotation[s.gff], hits[s.name]) for s in samples}
    if gene_values is not None:gene_values.update(metrics)
    for s in samples:
        annotated = {g["chrom"] for g in annotation[s.gff]}; total = sum(len(p) for p, _ in hits[s.name].values())
        outside = sum(len(p) for c, (p, _) in hits[s.name].items() if c not in annotated)
        if total and outside/total > .1:
            raise ValueError(f"{s.name}: {outside/total:.0%} of insertion sites are on chromosomes absent from its GFF "
                             f"(hits {sorted(set(hits[s.name])-annotated)[:4]}; GFF {sorted(annotated)[:4]}); check that hits and GFF use one reference")
    info = [asdict(s) | dict(insertion_sites=int(sum(len(p) for p, _ in hits[s.name].values())),
                             chromosomes_with_insertions=len(hits[s.name]),
                             genes=len(metrics[s.name]), genes_with_insertions=sum(m["sites"] > 0 for m in metrics[s.name].values()))
            for s in samples]
    rng = np.random.default_rng(seed); rows = []; n = len(samples); R = np.eye(n)
    for i in range(n):
        for j in range(i+1, n):
            a, b = samples[i], samples[j]
            ca, cb = set(hits[a.name]), set(hits[b.name])
            if a.species == b.species:ca = cb = ca & cb
            x, y = matched(a, b, metrics[a.name], metrics[b.name], orthologs, metric, ca, cb)
            row = dict(sample_a=a.name, sample_b=b.name, pair_type=pair_type(a, b), genes=len(x),
                       condition_a=a.condition, condition_b=b.condition, species_a=a.species, species_b=b.species)
            if len(x) < min_genes:
                rows.append(row | dict(rho=None, ci_low=None, ci_high=None, note=f"fewer than {min_genes} comparable genes")); R[i, j] = R[j, i] = np.nan; continue
            rho = float(spearmanr(x, y)[0]); draws = []
            for _ in range(boot):
                k = rng.integers(0, len(x), len(x)); draws.append(spearmanr(x[k], y[k])[0])
            row |= dict(rho=rho, ci_low=float(np.nanquantile(draws, .025)) if draws else None, ci_high=float(np.nanquantile(draws, .975)) if draws else None, note="")
            if a.species == b.species:
                ev = clone_evidence(tables[a.name], tables[b.name]); row |= ev; row["flag"] = conflict(row["pair_type"], ev["clone_evidence"])
            rows.append(row); R[i, j] = R[j, i] = rho
    return rows, R, info


def order(R):
    if len(R) < 3 or not np.isfinite(R).all():return np.arange(len(R))
    D = 1-R; return leaves_list(linkage(D[np.triu_indices(len(R), 1)], "average", optimal_ordering=True))


def plot(R, labels, path, title):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import textwrap
    n = len(R); size = 2.5+.45*n
    fig, ax = plt.subplots(figsize=(size+2, size))
    finite = R[np.isfinite(R) & ~np.eye(n, dtype=bool)]; low = min(0., float(finite.min())) if finite.size else 0.
    im = ax.imshow(np.where(np.isfinite(R), R, np.nan), cmap="viridis", vmin=low, vmax=1)
    if n <= 25:
        for i in range(n):
            for j in range(n):
                if i != j and np.isfinite(R[i, j]):ax.text(j, i, f"{R[i, j]:.2f}", ha="center", va="center", fontsize=6, color="white" if R[i, j] < (low+1)/2 else "black")
    ax.set_xticks(range(n), labels, rotation=90, fontsize=7); ax.set_yticks(range(n), labels, fontsize=7)
    fig.colorbar(im, ax=ax, label="Spearman rho across genes", shrink=.7); ax.set_title(textwrap.fill(title, 70), fontsize=9)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def write_tsv(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)


def run(samples, output, metric="density", orthologs=None, boot=500, seed=2026, min_reads=1):
    output = Path(output)
    if output.exists() and any(output.iterdir()):raise ValueError("Choose a new output directory")
    groups = ortholog_groups(orthologs) if orthologs else None
    if len({s.species for s in samples}) > 1 and groups is None:raise ValueError("Libraries from several species need --orthologs")
    genes = {}
    rows, R, info = correlate(samples, metric, groups, boot, seed, min_reads, gene_values=genes)
    output.mkdir(parents=True, exist_ok=True)
    species = {s.name:s.species for s in samples}
    with gzip.open(output/"gene_values.tsv.gz", "wt", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["library", "species", "gene", "chrom", "length", "sites", "reads", "density", "reads_per_kb"])
        for name, per_gene in genes.items():
            for gene, m in per_gene.items():
                w.writerow([name, species[name], gene, m["chrom"], m["length"], m["sites"], round(m["reads"], 3), round(m["density"], 5), round(m["reads_per_kb"], 5)])
    o = order(R); names = [samples[k].name for k in o]
    labels = [f"{samples[k].name} | {'/'.join(x for x in (samples[k].strain, samples[k].pool, samples[k].condition) if x)}" for k in o]
    write_tsv(output/"samples.tsv", list(info[0]), info)
    write_tsv(output/"pairs.tsv", ["sample_a", "sample_b", "pair_type", "species_a", "species_b", "condition_a", "condition_b", "genes", "rho", "ci_low", "ci_high",
                                   "shared_sites", "site_jaccard", "shared_site_read_rho", "clone_evidence", "flag", "note"], rows)
    write_tsv(output/"correlation_matrix.tsv", ["sample"]+names, [dict(sample=a) | {b:(None if not np.isfinite(R[i, j]) else float(R[i, j])) for j, b in zip(o, names)} for i, a in zip(o, names)])
    plot(R[np.ix_(o, o)], labels, output/"correlation.png", f"Hermes libraries: Spearman correlation of gene insertion {'density' if metric == 'density' else 'reads per kb'}")
    summary = {}
    for t in PAIR_TYPES:
        vals = [r["rho"] for r in rows if r["pair_type"] == t and r["rho"] is not None]
        if vals:summary[t] = dict(pairs=len(vals), median_rho=float(np.median(vals)), min_rho=float(min(vals)), max_rho=float(max(vals)),
                                  median_site_jaccard=(float(np.median(js)) if (js := [r.get("site_jaccard") for r in rows if r["pair_type"] == t and r.get("site_jaccard") is not None]) else None),
                                  median_shared_site_read_rho=(float(np.median(cs)) if (cs := [r.get("shared_site_read_rho") for r in rows if r["pair_type"] == t and r.get("shared_site_read_rho") is not None]) else None))
    flags = [dict(sample_a=r["sample_a"], sample_b=r["sample_b"], flag=r["flag"], shared_site_read_rho=r["shared_site_read_rho"]) for r in rows if r.get("flag")]
    manifest = dict(status="success", revision=REVISION, metric=metric, orthologs=str(orthologs) if orthologs else None, bootstrap=boot, seed=seed,
                    min_reads=min_reads, samples=len(samples), by_pair_type=summary, flagged_pairs=flags,
                    clone_evidence_thresholds=dict(shared=SHARED, independent=INDEPENDENT, calibration="Gale 2020: same pool 0.89-0.99, independent pools 0.42-0.43, different strains 0.33"),
                    interpretation={"same pool (shared clones)":"shares clones: high correlation reflects library structure, not biological similarity",
                                    "same strain, different pool":"independent libraries: reproducibility of the insertion landscape",
                                    "different strain":"strain differences in essentiality/accessibility, plus library noise",
                                    "different species":"conservation of the insertion landscape over 1:1 orthologs"},
                    files=dict(samples="samples.tsv", pairs="pairs.tsv", matrix="correlation_matrix.tsv", genes="gene_values.tsv.gz", plot="correlation.png"))
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest
