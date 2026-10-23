# Beta-geometric correction (BGC): optional Hermes adaptation

## Principle and data characteristics

BGC (`betageom`) is a nonlinear transformation for insertion-support distributions
dominated by a small number of high-count sites. Unlike a scalar multiplier, it
can change distribution shape. It is a candidate when a geometric target is a
reasonable approximation and disproportionate upper-tail counts are suspected
to be technical. It cannot determine whether high counts are technical artifacts
or real biological enrichment; the transformation can suppress genuine effects.

This is an experimental, user-selectable option, not a preference over TTR, CPM
or quantile. Evaluations below characterize behavior, not a method ranking.
More extensive validation is needed before recommendations. No app defaults,
existing fitness results or classifier workflows change.

Sources: [DeJesus and Ioerger, 2016](https://pubmed.ncbi.nlm.nih.gov/26932272/),
[TRANSIT documentation](https://transit.readthedocs.io/en/latest/method_normalization.html),
and the inspected TRANSIT `BetaGeomNorm` at commit
`50e23d6bec2605a24835fbd67f33691e330f3361` (`norm_tools.py`, lines 511–584).

## Model and parameter estimates

The positive-support model is:

```
p ~ Beta(alpha, beta)
X | p ~ Geometric(p), with X = 1, 2, ...
alpha = kappa*rho
beta  = kappa*(1-rho)
F_BG(x) = 1 - B(alpha, beta+floor(x)) / B(alpha,beta)
```

rho is the mean beta probability, not the inverse of the beta-geometric mean.
kappa is a concentration parameter; neither is viable-mutant abundance or a
Hermes motif preference parameter. No length/motif/bottleneck layer is fitted.

The inspected upstream fitting estimators are retained. For the sorted full
catalogue vector e_r, including its zeros, U sites and q_r=r/U (r=0,...,U-1):

```
rho = max(1 / trim_mean(e+1, 0.001), 0.0001)
A = sum(log(1-q_r)^2) / sum(e_r * log(1-q_r))
kappa = max((2*exp(A)-1)/(exp(A)+rho-1), 10)
```

These are upstream estimators, not a new maximum-likelihood or Bayesian fit.
The 0.1% trim is per tail and estimates rho only; sites/counts are not deleted.
The same chromosome/one-based coordinate catalogue and nonnegative raw support
are used as in the other methods, with strands summed at identical positions.
There is no Himar1 TA-site assumption.

**Catalogue dependence:** because zeros enter these estimates, adding all-zero
rows can change rho and kappa and thus the shape transformation—not just its
absolute scale. The full catalogue and source hashes are recorded. Testing a
10% zero-only expansion is a diagnostic of this sensitivity, not a proposed
alternative catalogue or a claim that those positions are biological targets.
No gene-only conditioning is silently substituted.

## Mapping and explicit final scaling

Each sample is fitted separately. Within an explicitly selected normalization
group, the target geometric mean is the equally weighted mean of samples'
raw positive-site means:

```
target_mean = mean_j(mean_i(raw positive counts in sample j))
p_target = 1 / target_mean
corrected_positive_count = Geometric_quantile(F_BG(raw_count), p_target)
```

Zeros remain zero. Equal counts receive equal corrected values. There is no
overall scalar factor for the nonlinear step. Whole-library counts are corrected
before gene aggregation; summing and then applying BGC is not equivalent.

The upstream default of final TTR scaling is supported explicitly: 5% per tail,
target 100 unless changed with `--trim`/`--target`. `--bgc-no-ttr` disables it.
This is part of the selected BGC method, not accidental double normalization.
CPM is never applied to its output. Fitted parameters, geometric target, engine,
post-TTR settings and factor are exported separately in `factors.tsv`.

## Numerical engines and differences from upstream

- `--bgc-cdf analytic` (default) computes the beta-geometric CDF directly using
  log-beta functions and expm1. It has no simulation noise and needs no seed.
- `--bgc-cdf simulation --bgc-seed 0 --bgc-draws 200000` uses the upstream-style
  beta/geometric sampling approximation. A local per-sample seed is derived from
  the specified seed and sorted raw counts so column order does not change the
  result. It does not change global RNG state. The generators differ from legacy
  TRANSIT; bitwise random-output parity is not claimed.
- Finite simulation can give ECDF=0 for a positive count. This adaptation keeps
  such a positive observation at a minimum corrected count of 1, rather than
  converting it to an unobserved insertion. Probability endpoints are bounded
  to avoid infinity; affected positive-site counts are reported. The analytic
  and simulation results therefore have explicit safeguards beyond upstream.
- Invalid/singular fits or simulation overflows fail clearly. The upstream's
  silent/printed geometric fallback is not used.

The retained parameter formulas are tested against the actual local TRANSIT
function when available. The exact CDF is independently checked against summed
beta-geometric probabilities and seeded mixtures; full stochastic output parity
with legacy TRANSIT is neither expected nor advertised.

## Commands

Example for Gale's four Cg1 cultures; substitute any compatible experiment's
matrix and new output paths. Zinc uses its four independent matched-pool groups.

```sh
python scripts/local/ytab.py normalize \
  --catalogue output/normalization/gale_flu_ttr_validation_v1/matrix \
  --output output/normalization/gale_flu_bgc_sites_v1 \
  --group-by library_id --method betageom --bgc-cdf analytic

python scripts/local/ytab.py normalize-genes \
  --validation output/normalization/gale_flu_ttr_validation_v1 \
  --output output/normalization/gale_flu_bgc_gene_validation_v1 \
  --existing-fitness output/projects/Gale_2020/treated_vs_parent/Cg1_FLU_vs_no_drug_raw_cpm \
  --site-normalization output/normalization/gale_flu_bgc_sites_v1

python scripts/local/ytab.py bgc-audit \
  --site-run output/normalization/gale_flu_bgc_sites_v1 \
  --output output/normalization/gale_flu_bgc_audit_v1 \
  --seed 0 --draws 200000
```

For zinc, corresponding paths start `zn_`; use `--compare-independent-pools`
for gene validation and its existing fitness result
`output/projects/Zn_toxicity_screen/treated_vs_parent/Zn_1_5mM_vs_mock`.
Cross-pool control differences remain library differences, not technical noise.

## Diagnostics and tests

Every run records rho, kappa, alpha, beta, A, input density, target geometric
mean, numerical safeguards and the fraction of positive-site support contributed
by the top 1% of occupied sites before/after correction. `bgc-audit` checks exact
analytic reproduction, fixed-seed simulation reproduction, unchanged source
hashes, zero patterns, monotonicity, analytic/simulation count differences and
catalogue sensitivity. Its CDF-discrepancy statistic is descriptive at observed
positive counts, not a fitted-model goodness-of-fit p-value.

60 focused tests cover these numerical/invariant properties, mixed-method CLI
output, nonlinear gene aggregation, existing CPM/raw parity and provenance.
Actual Gale and zinc audits passed. Their CDF-discrepancy statistics differ
substantially, underscoring that a valid implementation does not establish a good
biological model for every dataset. No fitted-model acceptance threshold or
automatic method selection is invented here.

Results and limitations: [Gale](../output/normalization/gale_flu_bgc_gene_validation_v1/INTERPRETATION.md)
and [zinc](../output/normalization/zn_bgc_gene_validation_v1/INTERPRETATION.md).
