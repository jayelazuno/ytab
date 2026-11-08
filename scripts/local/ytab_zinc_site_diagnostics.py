#!/usr/bin/env python3
"""Label-free site- and gene-level agreement diagnostics for independent-pool designs (zinc), with Gale as reference.

Reads only normalization outputs (raw counts + factors), the reference and the GFF. Reads no fitness,
z-score, classifier, GO or method call tables. Writes output/essentiality/pooled_gumbel/zinc_diagnostics/.
"""
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.normalization.ttr import write_table

spec = importlib.util.spec_from_file_location("v2", ROOT/"scripts/local/ytab_resampling_benchmark_v2.py")
v2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(v2)
OUT = ROOT/"output/essentiality/pooled_gumbel/zinc_diagnostics"
SCRATCH = ROOT/"output/essentiality/pooled_gumbel/_unused_output"
GALE = dict(controls=["Gale_2020::SRR11561505", "Gale_2020::SRR11561504"], treatments=["Gale_2020::SRR11561503", "Gale_2020::SRR11561502"],
            normalized="output/normalization/gale_flu_ttr_validation_v1/methods/ttr", gff=v2.ZN["gff"], reference=v2.ZN["reference"])
OFFSETS = (0, 7)   # Hermes 8-bp target duplication: one insertion can be reported at p or p±7 depending on read strand


def dataset(name):
    if name == "zinc":
        d = ccg.load_design(ROOT/v2.ZN["normalizations"]["ttr"], "ttr", ROOT/v2.ZN["reference"], SCRATCH, v2.ZN["controls"], v2.ZN["treatments"],
                            min_controls=1, stratified=True)
        labels = [f"mock{p}" for p in range(1, 5)]+[f"Zn{p}" for p in range(1, 5)]; pools = [1, 2, 3, 4]*2
    else:
        d = ccg.load_design(ROOT/GALE["normalized"], "ttr", ROOT/GALE["reference"], SCRATCH, GALE["controls"], GALE["treatments"])
        labels = ["ctrl1", "ctrl2", "FLU1", "FLU2"]; pools = ["Cg1"]*4
    counts = d["counts"]; counts = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
    return d, counts.astype(np.int64), labels, pools


def partner_fraction(chrom, pos, a, b, delta):
    """Fraction of sites in mask a that have a site in mask b at offset ±delta (delta=0: same coordinate)."""
    hit = np.zeros(len(pos), dtype=bool)
    for c in np.unique(chrom[a]):
        on = chrom == c; target = np.sort(pos[on & b])
        if not len(target):continue
        for shift in {delta, -delta}:
            q = pos[on & a]+shift; k = np.clip(np.searchsorted(target, q), 0, len(target)-1)
            hit[np.flatnonzero(on & a)[target[k] == q]] = True
    return float(hit[a].mean()) if a.any() else None


def pair_rows(name, counts, factors, labels, pools, chrom, pos, gene_idx):
    present = counts > 0; rows = []
    gene = np.array([[counts[i, j].sum() for j in range(counts.shape[1])] for i in gene_idx], dtype=float)
    scaled_gene = gene*factors   # normalized gene sums
    for i in range(len(labels)):
        for j in range(i+1, len(labels)):
            both = present[:, i] & present[:, j]; union = present[:, i] | present[:, j]
            strong = counts[:, i] >= 10
            g = (gene[:, i] >= 10) & (gene[:, j] >= 10)
            ratio = np.log2(scaled_gene[g, i]/scaled_gene[g, j])
            kind = ("within-pool" if pools[i] == pools[j] else "between-pool")+": "+"-".join(sorted({labels[i].rstrip("0123456789"), labels[j].rstrip("0123456789")}))
            rows.append(dict(dataset=name, a=labels[i], b=labels[j], comparison=kind, site_jaccard=both.sum()/union.sum(),
                             strong_sites_a_present_b=partner_fraction(chrom, pos, strong, present[:, j], 0),
                             strong_sites_a_partner_at_7bp_b=partner_fraction(chrom, pos, strong, present[:, j], 7),
                             site_spearman_shared=float(spearmanr(counts[both, i], counts[both, j])[0]),
                             gene_spearman=float(spearmanr(scaled_gene[g, i], scaled_gene[g, j])[0]),
                             gene_log2ratio_robust_sd=float(1.4826*np.median(np.abs(ratio-np.median(ratio)))), genes_compared=int(g.sum())))
    return rows


def sample_rows(name, counts, labels):
    rows = []
    for j, label in enumerate(labels):
        c = counts[:, j]; nz = c[c > 0]; s = np.sort(nz)[::-1]
        rows.append(dict(dataset=name, sample=label, reads=int(c.sum()), sites=int(len(nz)), median_reads_per_site=float(np.median(nz)),
                         singleton_sites=float(np.mean(nz == 1)), top1pct_read_share=float(s[:max(1, len(s)//100)].sum()/s.sum())))
    return rows


def mock_noise_fits(d, counts, factors, pairs):
    """Site noise curve (selection-aware lognormal-NB) fitted on treatment-free pairs only."""
    out = {}; settings = ccg.Settings(); inverse = 1/factors
    for i, j in pairs:
        pair = counts[:, [i, j]]; keep = pair.sum(1) >= 1; e = inverse[[i, j]]/inverse[[i, j]].mean()
        fit = ccg.fit_background(pair[keep], e, settings)
        out[f"{i}-{j}"] = {k:fit[k] for k in ("m", "tau", "phi", "linear", "boundary", "sites")}
    return out


def clone_sharing():
    """Zinc: pairwise overlap of sites above read thresholds (1, 10, 100, 1000 reads, and >= 1 read per million)."""
    d, counts, labels, _ = dataset("zinc"); rows = []
    for name, mask in [(f">= {t} reads", counts >= t) for t in (1, 10, 100, 1000)]+[(">= 1 RPM", (counts >= counts.sum(0)/1e6) & (counts > 0))]:
        share = (counts*mask).sum(0)/counts.sum(0)
        for i in range(len(labels)):
            rows.append(dict(threshold=name, sample=labels[i], sites=int(mask[:, i].sum()), read_share=float(share[i]),
                             **{labels[j]:float((mask[:, i] & mask[:, j]).sum()/max((mask[:, i] | mask[:, j]).sum(), 1)) for j in range(len(labels))}))
    write_table(OUT/"clone_sharing.tsv", list(rows[0]), rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if "--clones-only" in sys.argv:clone_sharing(); return
    clone_sharing(); pairs, samples, fits = [], [], {}
    for name in ("zinc", "gale"):
        d, counts, labels, pools = dataset(name)
        chrom = np.array([c for c, _ in d["sites"]]); pos = np.array([p for _, p in d["sites"]], dtype=np.int64)
        genes = ccg.annotations(ROOT/v2.ZN["gff"]); idx = ccg.gene_sites(genes, d["sites"], ccg.Settings())
        pairs += pair_rows(name, counts, d["factors"], labels, pools, chrom, pos, idx); samples += sample_rows(name, counts, labels)
        mock_pairs = [(0, 1), (2, 3), (0, 2), (1, 3)] if name == "zinc" else [(0, 1)]
        fits[name] = {f"{labels[i]}-{labels[j]}":v for (i, j), v in zip(mock_pairs, mock_noise_fits(d, counts, d["factors"], mock_pairs).values())}
    write_table(OUT/"pairs.tsv", list(pairs[0]), pairs); write_table(OUT/"samples.tsv", list(samples[0]), samples)
    summary = {}
    for row in pairs:
        key = f"{row['dataset']} | {row['comparison']}"; summary.setdefault(key, []).append(row)
    summary = {k:{m:float(np.mean([r[m] for r in v])) for m in ("site_jaccard", "strong_sites_a_present_b", "strong_sites_a_partner_at_7bp_b",
               "site_spearman_shared", "gene_spearman", "gene_log2ratio_robust_sd")} | dict(pairs=len(v)) for k, v in summary.items()}
    result = dict(status="success", treatment_free_site_noise_fits=fits, comparison_means=summary,
                  blinding="normalization outputs only; no fitness, z-score, classifier, GO or call tables read")
    (OUT/"diagnostics.json").write_text(json.dumps(ccg._json_ready(result), indent=2)+"\n")
    print(json.dumps(ccg._json_ready(result), indent=1))


if __name__ == "__main__":main()
