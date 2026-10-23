# Control-conditioned Hermes Gumbel (experimental)

Status: **experimental run-based depletion evidence** (`ccg_v5_normalized_input`). It is meant as an
orthogonal method alongside the z-score fitness analysis. It is not a validated essentiality,
fitness or published Hermes method. It is the only Gumbel in YTAB: the earlier single-condition
Gumbel and its outputs were removed after this method was tested.

## Workflow position

1. **Normalize** with `ytab normalize`, choosing `ttr` (recommended default), `totreads` or `nzmean`.
2. **Run** `ytab control-gumbel` on that normalization output. The Gumbel runs no normalization of
   its own.

The count model needs raw integer counts, so a normalization enters as an **exposure** (a per-sample
depth multiplier on expected counts), not as rescaled counts. For observation *j*:
`exposure_j = (1/factor_j) / mean over controls of (1/factor)`. The runner reads the normalization's
`factors.tsv` and takes raw counts and coordinates from that normalization's source catalogue. It
refuses the input unless:
- `normalized_counts.npz == raw × factor` exactly;
- the catalogue's checksums still match those recorded at normalization;
- all selected tracks share one normalization group.

Distribution-reshaping methods (quantile, betageom, aBGC, zinfnb), `nonorm` and planned methods are
refused. They cannot be expressed as one depth multiplier per sample, or they were not tested with
this method.

| normalization | factor_j | exposure_j ∝ |
|---|---|---|
| `ttr` (default) | target / (density_j × trimmed positive-site mean_j) | density × trimmed mean (robust to a few extreme sites) |
| `totreads` | mean total / total_j | total insertion support |
| `nzmean` | mean positive-site mean / positive-site mean_j | mean support per occupied site |

## Why an extension of Gumbel

A presence/absence Gumbel (TRANSIT's Gumbel method; DeJesus et al. 2013) asks whether a gene has an unusually long run of **empty**
sites. It sees neither the matched control nor read abundance, so it cannot recognise
treatment-specific depletion while reads remain. In Gale, UPC2A keeps 7 of 68 sites occupied in FLU
yet loses about 97% of reads; the presence/absence version matched 0 of 20 published FLU-depleted
genes. This extension keeps the longest-run / extreme-value logic and changes what a run is made of:

| | Presence/absence Gumbel | Control-conditioned Gumbel |
|---|---|---|
| Opportunities | catalogue sites | sites supported by the combined controls |
| Site event | count < minread (empty) | combined treatment count improbably low given combined controls |
| Statistic | longest empty run | longest run of depleted sites |
| Null | theoretical Bernoulli Gumbel | no-effect run distribution from the control-conditioned count model |
| Output | E/EB/NE/U/S/NA | p, BH q, `depleted run` / `not significant` / `not tested` |

## Design: combined replicates

The control replicates are **summed into one control track**, which defines the null universe
(eligible sites) and the baseline. The treatment replicates are **summed into one treatment track**,
which is compared against it. Each condition's noise is estimated from its own replicates *before*
combining, using only how much they disagree with each other. A treatment effect that both
replicates share creates no disagreement, so it stays signal.

## Model

For site *i*, with control replicates *C<sub>ij</sub>*, treatment replicates *T<sub>ik</sub>* and
exposures *e* from the chosen normalization:

1. **Eligibility (controls only).** Combined control `ΣC_ij ≥ min_control_reads` (default 1).
   Optional `--eligibility detectable` also requires that a zero treatment count would be flagged.
2. **Per-condition count model.** `X ~ NB(e·λ_i)` with a shared two-term mean–variance curve
   `Var = μ(1 + a) + φμ²`. Here *a* is mean-proportional extra variance (bottleneck-like sampling of
   clones) and φ is the usual quadratic dispersion. `log λ_i ~ N(m, τ²)`.
   - Controls: m, τ, φ_C and a_C are fitted on the control replicate tracks by maximising the
     **zero-truncated** marginal likelihood `p(C_i)/P(ΣC_i ≥ min_control_reads)`. Sites enter because
     the controls saw them; ignoring that selection biases abundance upward (tested).
   - Treatment (`--treatment-dispersion replicates`, default): φ_T and a_T are fitted the same way on
     the treatment replicate tracks. `control` assumes control dispersion explicitly.
   - `--dispersion-model quadratic` fixes a = 0.

   Integrals use adaptive Gauss–Hermite quadrature around each posterior mode (32 nodes, checked
   against 64).
3. **Combined tracks.** A sum of independent replicates is approximated by an NB with the same mean,
   the same *a*, and dispersion `φ·Σe²/(Σe)²` (moment-matched; tested by simulation).
4. **Depletion encoding.** `F_i(t) = ∫ NB_CDF(t; E_T λ, φ_T', a_T) p(λ | ΣC_i) dλ`, and
   `D_i = 1[F_i(ΣT_i) ≤ α_site]`. α_site = 0.05 is an encoding setting, not gene significance.
   Equivalently `D_i = 1[ΣT_i ≤ t*_i]`, where `t*_i = −1` when even a zero would not be flagged: weak
   controls give uncertainty, not depletion. Under the null, site *i* is flagged with probability
   `π_i = F_i(t*_i)`.
5. **Statistic.** `R_g` = longest run of consecutive `D_i = 1` across the gene's eligible sites in
   genomic order (own CDS union, last 3 bases excluded, N/C trimming as in TRANSIT's Gumbel method). The last
   tied maximum wins. Span = `end − start + 1`, reported as a diagnostic only.
6. **Gene p-value (`--inference`).** Both modes use the same null: sites are conditionally independent
   given λ, with flag probabilities π_i.
   - `exact` (**default**): `p = P(R ≥ r)` computed exactly by a dynamic program over π_i.
     Deterministic and seed-free. A discrete Gumbel is still fitted to 1000 seeded simulations and
     checked; the result is reported in `gumbel_check` but does not gate.
   - `gumbel-mc`: a discrete Gumbel rounded/truncated on 0..n,
     `P(R=r) = [G(r+½) − G(r−½)] / [G(n+½) − G(−½)]`, is fitted by maximum likelihood to 20,000 seeded
     no-effect simulations (seeded per gene: same answer regardless of gene order). Then
     `p = [G(n+½) − G(r−½)] / [G(n+½) − G(−½)]`. The p is used only if the fit converged, which
     requires both checks below to pass (they allow for simulation error):
     - the 20,000 validation runs lie within a two-sample KS bound `c·√(1/n_val + 1/n_cal)`;
     - `|log(p_fit/p_exact)| ≤ log 2 + z·SE`, where SE is the delta-method Monte Carlo standard error
       of the fitted log tail. Below `fdr/family_size`, a disagreement is tolerated and the gene is
       labelled "magnitude unvalidated".

     Otherwise the status is `Gumbel null inadequate: …` and p = 1.
7. **Multiple testing.** The family is the genes with ≥ `min_sites` (3) eligible sites, fixed from
   control eligibility before any treatment result. BH at `fdr` = 0.05. Genes with no depleted or no
   detectable site enter with p = 1.

## Running

```bash
PYTHONPATH=src python scripts/local/ytab.py normalize --catalogue <catalogue> --method ttr \
  --group-by library_id --output <norm>          # optionally also --method totreads --method nzmean

PYTHONPATH=src python scripts/local/ytab.py control-gumbel \
  --normalized <norm>            # run directory (uses ttr unless --normalization totreads|nzmean), or <norm>/<method>
  --gff resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff \
  --reference resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.fna \
  --control Gale_2020::SRR11561505 --control Gale_2020::SRR11561504 \
  --treatment Gale_2020::SRR11561503 --treatment Gale_2020::SRR11561502 \
  --output <new directory>       # add --inference gumbel-mc for the Gumbel cross-check
```

All tracks must share project, `library_id`, pool and normalization group. Controls share one role
and treatments another. Runtime on the Gale matrix is about 10 min for `exact`; `gumbel-mc` takes
considerably longer. Outputs:
- `genes.tsv`: effect sizes, run, p_value, q, call, status, exact and fitted tails, `gumbel_check`.
- `sites.tsv.gz`: per-replicate and combined counts, predictive CDF, threshold, null flag
  probability, flag.
- `run_manifest.json`: normalization method, input, factors and hashes; exposures; parameters; both
  dispersion fits with integration and posterior-predictive checks; limitations.

Synthetic evidence: `scripts/local/ytab_control_conditioned_gumbel_synthetic.py --inference {exact,gumbel-mc}` →
`output/essentiality/control_conditioned_gumbel/synthetic_validation_v4_{exact,gumbel_mc}/`.
Tests: `tests/test_control_conditioned_gumbel.py`.

## Synthetic evidence

Each scenario has 8 experiments × 300 genes, with 2 control and 2 treatment tracks combined as in the
method, and known exposures. Generative settings are generic and not fitted to Gale. Values are
given as exact / gumbel-mc (20k), with 95% Clopper–Pearson intervals for exact.

| Scenario (no-effect unless stated) | unaffected p ≤ 0.05 | false / total BH calls (exact · mc) | site flags obs/exp | sensitivity |
|---|---|---|---|---|
| equal depth | 1.17% (0.78–1.68) / 1.17% | 0 / 0 · 0 / 0 | 0.95 | – |
| shallow treatment (0.3×) | 1.63% (1.16–2.22) / 1.63% | 0 / 0 · 0 / 0 | 1.00 | – |
| deep treatment (2–2.5×) | 1.38% (0.95–1.93) / 1.38% | 0 / 0 · 0 / 0 | 0.98 | – |
| low-read controls | 1.06% (0.68–1.55) / 1.06% | 0 / 0 · 0 / 0 | 0.99 | – |
| high dispersion (φ = 1) | 1.00% (0.64–1.49) / 1.00% | 0 / 0 · 0 / 0 | 0.93 | – |
| misspecified prior (5% jackpots) | 1.21% (0.81–1.73) / 1.21% | 0 / 0 · 0 / 0 | 0.99 | – |
| conditional (fixed controls) | 1.42% (0.98–1.97) / 1.46% | 0 / 0 · 0 / 0 | 1.02 | – |
| noisy treatment replicates (φ_T = 25 φ_C) | 1.38% (0.95–1.93) / 1.38% | 0 / 0 · 0 / 0 | 0.98 | – |
| same, **assuming control dispersion** | 49.9% / 49.8% | **934 / 934 · 931 / 931** | 6.17 | – |
| bottleneck treatment noise (a_T = 4) | 1.29% (0.88–1.83) / 1.29% | 0 / 0 · 0 / 0 | 1.02 | – |
| same, **quadratic-only model** | 2.46% (1.88–3.16) / 2.42% | 0 / 0 · 0 / 0 | 1.31 | – |
| complete loss, contiguous half of gene | 1.67% / 1.57% | 2 / 219 · 2 / 217 | 2.08 | 0.90 (0.86–0.94) / 0.90 |
| 4× depletion, contiguous | 0.93% / 0.93% | 0 / 164 · 0 / 164 | 1.82 | 0.68 (0.62–0.74) / 0.68 |
| 4× depletion, noisy treatment | 1.11% / 1.16% | 1 / 132 · 0 / 132 | 1.64 | 0.55 (0.48–0.61) / 0.55 |
| 4× depletion, bottleneck noise | 1.62% / 1.62% | 3 / 181 · 3 / 181 | 1.94 | 0.74 (0.68–0.80) / 0.74 |
| 2× depletion, contiguous | 1.25% / 1.25% | 0 / 30 · 0 / 32 | 1.28 | 0.13 (0.09–0.17) / 0.13 |
| complete loss, scattered (non-adjacent) | 1.20% / 1.20% | 0 / 14 · 0 / 15 | 1.81 | 0.06 (0.03–0.10) / 0.06 |

- **Calibration.** Under the model the test is conservative: about 1–1.6% of unaffected genes reach
  p ≤ 0.05 at nominal 5%, with at most 3 false BH calls per scenario (empirical FDR ≤ 2%). Site flags
  match expectation (obs/exp 0.93–1.02) in every correctly specified null. Across all 17 scenarios,
  7 genes failed the gumbel-mc adequacy check.
- **The two misspecified rows show why the design matters:**
  - Treating treatment noise as equal to control noise gives a 6× flag excess and over 900 false calls.
  - Omitting the bottleneck term leaves a 1.3× excess.
- **Sensitivity:**
  - Detection needs adjacent sites and roughly ≥ 4× depletion at α_site = 0.05.
  - For 4× contiguous depletion, sensitivity is 0.98–1.0 for genes with ≥ 20 sites and 0.40 for genes
    with 3–19 sites.
  - Scattered losses are mostly missed, by design of a run statistic; the z-score covers that case.

These numbers establish calibration only when the data follow the assumed model. They do not
validate genome-wide tiny tails.

## Gale FLU (combined controls 05+04 vs combined FLU 03+02)

See `output/essentiality/control_conditioned_gumbel/Gale_2020/INTERPRETATION.md`. Reproducibility
of the tested settings:

| | exact | gumbel-mc (20k) |
|---|---|---|
| between seeds (Jaccard) | 1.000 (identical p) | 0.97–0.98 |
| ttr / totreads / nzmean calls | 136 / 146 / 151 | 140–141 / 145–146 / 153–154 |
| strict published genes called (of 19) | 15 / 15 / 16 | 15 / 15 / 16 |

The earlier 1000-simulation Gumbel with a fixed ×2 adequacy gate was removed: its call sets
overlapped only 0.56 between seeds.

## Limitations

- Dispersion is a shared two-term curve per condition; site-specific dispersion is not modelled.
- Treatment dispersion from two replicates cannot separate replicate-specific effects from noise.
- Reproducible but non-gene-level (site/clone-level) treatment deviations also raise site flags. The
  run statistic needs adjacency, which limits but does not remove their influence.
- Under the null, sites are conditionally independent given λ. Spatially correlated noise (shared
  clone or positional effects) would make runs longer than the null allows. This is not modelled.
- Plug-in hyperparameters; their uncertainty is not propagated.
- The lognormal abundance prior is an assumption (posterior-predictive check reported).
- Scalar normalization exposures assume most sites are unaffected by treatment. The choice of scalar
  normalization shifts borderline calls (Gale: 3–11% of calls differ between pairs of methods).
- Overlapping genes share sites; BH assumes valid p-values.
