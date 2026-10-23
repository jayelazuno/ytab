# Condition comparisons: treated vs control

YTAB compares treated libraries with their matched controls to find genes whose mutants are depleted (the gene helps
survive the treatment) or enriched (losing the gene helps). Three methods are offered. They answer slightly different
questions and were validated in a pre-registered benchmark.

| method | what it gives | normalization (default) | use it for |
|---|---|---|---|
| **Local-abundance z-score** | one z per gene, depleted/enriched calls by threshold; no p-values | BGC | ranking genes, comparison with published z-score screens |
| **Resampling** | site-restricted permutation test, p and BH q, both directions | BGC (trimmed-sum recentering) | hit lists with false-discovery control |
| **Control-conditioned Gumbel** | runs of depleted sites vs the controls' count model, p and q | TTR | depletion-only hit lists with strict error control |

## Input and design

A comparison design pairs each treated library with its control (background, pool). The CLI writes it with
`ytab comparisons`, and it can be edited in the app.

**Normalization:**
- **What it is:** any Site Normalization method.
- **How it enters the z-score:** the method's normalized site counts are summed per gene and expressed per million of
  the mean library total. The pseudocount (0.5) therefore means the same as with counts per million.
- **Default:** BGC (beta-geometric correction), which the benchmark supports (below).
- **Fallback:** total reads (CPM) needs no normalization run.

## Local-abundance z-score (Gale et al. 2020)

1. **Control noise model.** For each background, the two control libraries give A = mean abundance and M = log2 ratio
   (pseudocount 0.5). The running SD of M over 41 genes, ordered by A, is fitted as SD(A) = m1 + m2·A^m3. Low-count
   genes therefore get a wider noise band.
2. **Replicates averaged.** The replicate libraries of each condition are averaged, and each gene gets one
   log2FC = log2((treated + 0.5) / (control + 0.5)) and z = log2FC / SD(A). This is Gale's published procedure.
3. **Call.** |z| beyond the 99th percentile of the control–control |z| gives depleted or enriched.
4. **Top hits.** For labels, the heatmap and GO, the called genes with at least a 2-fold change are ranked by |z|
   (`top_hit_evidence.csv`).

The earlier per-pair rule ("consistent in all pairs" / "single pool") is available with `--combine pairs`, but made most
of the z-score's false calls in the benchmark.

## Resampling

A permutation test on the gene's normalized site counts, controls vs treated. Labels are permuted within each pool when
there are several. The test is two-sided, with Benjamini–Hochberg FDR. With BGC, the treated libraries are recentred
(trimmed sum) so that typical genes show no change. Adapted from TRANSIT's resampling (DeJesus et al. 2015).

## Control-conditioned Gumbel

Finds each gene's longest run of consecutive depleted sites in the treated libraries and compares it with the longest
run expected from the controls' count model. It uses the exact tail probability, with the Gumbel fit as a diagnostic.
It tests depletion only, needs two or more controls per background, and accepts scalar normalizations (TTR, total
reads, nonzero mean). Adapted from TRANSIT's Gumbel method (DeJesus et al. 2013). Details: `docs/control_conditioned_gumbel.md`.

## Validation

**Pre-registered simulation benchmark** (Gale-calibrated; 140 signal and 40 null experiments; amended before
unblinding). Full results: `output/treated_vs_control/benchmark_blind_2026-10-22/RESULT.md`.

- **Normalization:** BGC beats CPM and TTR for both the z-score (AUC +0.012 vs CPM) and resampling (+0.008 vs total
  reads).
- **Ranking:** the z-score with BGC ranks affected genes best (precision at K +0.086 over resampling).
- **Calls:** only resampling with BGC (FDR 0.063) and the Gumbel with TTR (FDR 0.010) pass the calibration gates. The
  z-score has no error control: 64 false calls per null experiment, FDR 0.20.

**Real data** (Gale fluconazole; `output/treated_vs_control/benchmark_phase2_2026-10-22/RESULT.md`):
- **Published hits:** all methods recover Gale's published hits. Of the 19 strict hits, the z-score calls 19,
  resampling 18 and the Gumbel 15. UPC2A ranks first with the z-score.
- **Swap null** (control + FLU vs control + FLU): the z-score and the Gumbel make no calls, and resampling with BGC
  makes 7.

## Outputs the app reads

Every run writes the same tables, whatever the normalization, library or species, and the command fails if one is
missing (`src/ytab/fitness/contract.py`):
- **z-score:**
  - `zscore_manifest.json`;
  - the per-gene summary (call, effect, z, abundances);
  - per-pair and combined statistics;
  - gene values;
  - library sizes;
  - the control–control z;
  - the noise-model fit;
  - top-hit evidence.
- **Resampling / Gumbel:** `manifest.json`, `genes.tsv`, and per group `sites.tsv.gz`.

## Commands

```sh
python scripts/local/ytab.py fitness --project-config P/config/project.yaml --normalization betageom --normalization-run P/normalization/RUN
python scripts/local/ytab.py treated-vs-control --project-config P/config/project.yaml --method resampling \
  --normalization-run P/normalization/RUN --normalization betageom --output P/treated_vs_control/RUN
```

## References

- DeJesus MA, Ioerger TR (2013). A hidden Markov model for identifying essential and growth-defect regions in bacterial
  genomes from transposon insertion sequencing data. BMC Bioinformatics 14:303.
- DeJesus MA, Zhang YJ, Sassetti CM, Rubin EJ, Sacchettini JC, Ioerger TR (2013). Bayesian analysis of gene
  essentiality based on sequencing of transposon insertion libraries. Bioinformatics 29:695–703 (the Gumbel method).
- DeJesus MA et al. (2015). TRANSIT. PLoS Comput Biol 11:e1004401.
- Gale AN et al. (2020). Identification of essential genes and fluconazole susceptibility genes in Candida glabrata by
  profiling Hermes transposon insertions. G3 10:3859–3870.
