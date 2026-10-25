"""Per-library QC table for a project: every QC metric in one table, flagged and commented, plus plot data.

Collects, for each included library, the read and mapping metrics (Hermes cleanup, MapFastq), the Tn-seq count
statistics computed from the hit file (tnseq_stats.statistics, over the project's empirical site catalogue), the
strand correlation, the SummaryTable gene-level metrics and the saved library diagnostics (MidLC, jackpots).

Flags (green / yellow / red; grey = not flagged) never use invented cut-offs:
- "rule of thumb": only where published Tn-seq QC guidance is scale-free and carries over to Hermes (skewness of the
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
    ("input_reads", "reads", "Reads sequenced", "low", "relative", "Reads entering Hermes cleanup (the junction read)."),
    ("retained_fraction", "reads", "Reads retained after cleanup", "low", "relative", "Fraction of reads kept after Hermes/adaptor trimming and the length filter."),
    ("adapter_fraction", "reads", "Reads with 3' adaptor", "high", "relative", "Fraction of reads with adaptor at the 3' end (short fragments or primer dimers)."),
    ("mapped_reads", "reads", "Mapped reads", "low", "relative", "Primary alignments that map to the reference."),
    ("percent_mapped", "reads", "Mapped (%)", "low", "relative", "Primary alignments that map to the reference."),
    ("percent_mapq", "reads", "Mapped with MAPQ above threshold (%)", "low", "relative", "Uniquely placed alignments; low values mean repeats or the wrong reference."),
    ("mean_mapq", "reads", "Mean MAPQ", "low", "relative", "Mean mapping quality of mapped primary alignments."),
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
    ("fr_corr", "insertions", "Strand correlation (W vs C)", None, None, "Correlation of read counts on the two strands at a position. For Hermes the two orientations are separate insertions, so it is shown, not flagged."),
    ("percent_features_hit", "genes", "Genes hit (%)", "low", "relative", "SummaryTable: features with at least one insertion."),
    ("mean_hits_per_feature", "genes", "Insertions per gene", "low", "relative", "SummaryTable: mean insertion sites per feature."),
    ("percent_hits_in_features", "genes", "Insertions in genes (%)", "both", "relative", "SummaryTable: share of insertion sites inside features."),
    ("neighborhood_spearman", "genes", "Insertion-neighbourhood correlation", "low", "relative", "SummaryTable: Spearman correlation between a gene's insertions and its neighbourhood's."),
    ("midlc_est", "saturation", "MidLC estimate", None, None, "Library diagnostics: reads at which half the discoverable sites are found."),
    ("depth_ratio", "saturation", "Depth / MidLC", None, None, "Reads at insertion sites divided by the MidLC estimate (above 1: sequenced past half-saturation)."),
    ("jackpot_frac_reads", "saturation", "Reads in jackpot sites", "high", "relative", "Library diagnostics: fraction of reads in the top jackpot sites."),
    ("n_jackpot_sites", "saturation", "Jackpot sites", "high", "relative", "Library diagnostics: number of sites flagged as jackpots."),
]
HINTS = {"skewness":BGC_HINT, "max_share":TTR_HINT, "top1pct_share":BGC_HINT, "jackpot_frac_reads":TTR_HINT,
         "pickands_tail_index":BGC_HINT}
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


def _hermes(path: Path) -> dict:
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
    rows, plots = [], {"histogram":[], "qq":[], "ranked":[], "lorenz":[]}
    for s in samples:
        name = str(s["sample"]); role = str(s.get("library_role") or s.get("condition") or s.get("guessed_condition") or "library").strip().lower()
        row = dict(sample=name, display_name=s.get("display_name") or name, role=role, pool=s.get("pool") or s.get("guessed_pool") or "",
                   background=s.get("background") or s.get("guessed_background") or "")
        row |= (_hermes(project/"mapfastq"/name/f"{name}.hermes_preprocessing.json") or _mapping_log(project/"mapfastq"/name/f"{name}_log.txt")) \
            | _mapping(project/"mapfastq"/name/f"{name}.mapping_stats.csv")
        if name in hits:
            row |= insertion_metrics(*hits[name], universe, genome_bp)
            data = plot_data(name, np.fromiter(hits[name][0].values(), dtype=float))
            for k in plots:plots[k] += data[k]
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
    output.mkdir(parents=True)
    keys = ["sample", "display_name", "role", "pool", "background"]+[m[0] for m in METRICS]+["max_site", "pti_valid_terms", "catalogue_sites", "occupied_sites", "stats_notes"]
    _write(output/"qc_table.tsv", keys, rows)
    _write(output/"qc_flags.tsv", list(flags[0]), flags)
    _write(output/"metrics.tsv", ["metric", "theme", "label", "bad_direction", "rule", "comment"],
           [dict(metric=k, theme=THEMES[t], label=l, bad_direction=d or "", rule=r or "descriptive", comment=c) for k, t, l, d, r, c in METRICS])
    (output/"plots").mkdir()
    for name, data in plots.items():
        _write(output/"plots"/f"{name}.tsv.gz", list(data[0]) if data else ["sample"], data)
    write_xlsx(output/"qc_table.xlsx", rows, flags)
    manifest = dict(status="success", project=config.get("project_id"), project_config=str(config_path), libraries=len(rows),
                    libraries_without_hits=missing, catalogue_sites=universe, genome_bp=genome_bp,
                    flag_counts={c:sum(f["flag"] == c for f in flags) for c in ("red", "yellow", "green", "grey")},
                    rules=dict(rule_of_thumb=RULE_OF_THUMB, relative=LEYS, yellow=ROBUST_YELLOW, red=ROBUST_RED, min_role_libraries=MIN_ROLE_LIBRARIES, min_relative_spread=MIN_RELATIVE_SPREAD,
                               comparison="within role (parents with parents, treated with treated)"),
                    created=datetime.now(timezone.utc).isoformat(),
                    files=dict(table="qc_table.tsv", flags="qc_flags.tsv", metrics="metrics.tsv", excel="qc_table.xlsx",
                               histogram="plots/histogram.tsv.gz", qq="plots/qq.tsv.gz", ranked="plots/ranked.tsv.gz",
                               lorenz="plots/lorenz.tsv.gz"))
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
