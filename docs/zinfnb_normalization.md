# Zero-inflated negative-binomial scaling: mathematics and assumptions

`zinfnb` is an experimental, user-selectable scalar normalization. It fits each
sample's raw insertion-site counts independently, then divides by the fitted
negative-binomial (NB) component mean. Its availability is not a recommendation
over other methods. Existing fitness CPM, classifier MidLC and UI defaults are
unchanged; no normalization tab is introduced.

## Principle and assumptions

The model combines an NB count distribution with an additional zero component.
NB variability can exceed Poisson variability. This is a distributional model of
observed insertion support, not a mechanism for identifying essential genes,
missing mutants, PCR artifacts, or viable-mutant frequencies.

The fitting likelihood treats sites as independent observations from a common
within-sample distribution. Real genomic correlation, heterogeneous insertion
opportunity, PCR amplification and genuine selection can violate that
approximation. The further normalization assumption is that the fitted NB mean
is a suitable scale reference across compared samples. Convergence does not
establish either assumption. Scalar scaling cannot remove distributional skew,
restore mutants, or distinguish a genuine global shift from a depth difference.

All observed coding and intergenic coordinates enter the empirical catalogue;
zeros denote no observed support, not known inaccessible positions. There is no
extra TA-site, Hermes motif, gene-length or bottleneck weighting. Those features
are not separately fitted. Coordinate identity and strand merging follow the
same raw hit-file conventions as the other methods.

## Exact probability model

For nonnegative integer x, size n>0 and 0<p<1:

```
NB(x;n,p) = Gamma(x+n)/(Gamma(n)*Gamma(x+1)) * p^n * (1-p)^x
mu = n*(1-p)/p
Var_NB = mu + mu^2/n

P(X=0) = pi + (1-pi)*p^n
P(X=x>0) = (1-pi)*NB(x;n,p)
E[X] = (1-pi)*mu
```

Here pi is the fitted additional-zero mixture probability. It is not a measured
fraction of essential, inaccessible or non-insertion sites. NB itself also
generates zeros, so pi is not the observed zero fraction.

Let U be the full empirical catalogue size, z the number of zeros, and h_x the
histogram frequency of each observed positive count. We minimize:

```
L = -[z*log(pi+(1-pi)*p^n)
      + sum_(x>0) h_x*(log(1-pi)+log NB(x;n,p))] / U
```

Histogram compression is exactly equivalent to the full count-vector objective;
it is not binning into approximate count ranges or downsampling sites.

## Difference from the inspected TRANSIT source

The reference is `ZeroInflatedNBNorm` / `Fzinfnb` in
[TRANSIT norm_tools.py, commit 50e23d6](https://github.com/ioerger/transit/blob/50e23d6bec2605a24835fbd67f33691e330f3361/src/pytransit/norm_tools.py).
Its zero term uses `pi + NB(0)` without the multiplier `(1-pi)`, while positive
counts do include that multiplier. That is not a normalized ZINB probability
distribution. YTAB intentionally uses the proper mixture above; numerical parity
with that upstream objective is neither expected nor claimed. The upstream
inverse-NB-mean scaling interpretation is retained.

## Optimization and numerical parameters

We optimize log(mu) and log(n), with `p=n/(n+mu)`. For fixed mu,n, let q=z/U and
b=P_NB(0). The maximizing mixture weight is available analytically:

```
pi_hat = clip((q-b)/(1-b), 0, 1-1e-12)
```

This profiles pi rather than optimizing three parameters simultaneously.
Three deterministic starts use raw positive-count mean for mu and n in
`{0.2, 2, 20}`. L-BFGS-B minimizes mean negative log likelihood with
`maxiter=1000`, `ftol=1e-12`, `gtol=1e-7`.

Numerical safeguard bounds are log(mu) in [-18,25] and log(n) in [-14,18].
They are not biological priors or fitted library limits. The best converged,
finite-likelihood start is selected; counts of successful starts, optimizer
message and iterations are exported. All-start failure raises an error, with no
claim that the best of three local fits is a proven global optimum and no
silent fallback to another method. Estimates within 0.001 of a log-parameter
bound, or pi=0, receive a boundary warning. Pi=0 is the valid ordinary-NB boundary,
not an automatic failed fit. Empty samples and invalid raw inputs fail clearly.
No simulation, random seed, trimming or post-TTR step is used.

## Scaling, groups and gene summaries

```
factor_j = 1/mu_hat_j
normalized_site_ij = factor_j * raw_site_ij
normalized_gene_gj = factor_j * raw_gene_gj
```

It is the NB-component mean that is scaled to 1. The full mixture's fitted mean
after scaling is 1-pi; observed normalized totals and positive means need not
match across samples. This is not mixture-mean scaling or a group-target method.
Each sample is fitted separately, even in a one-sample group. Explicit group
metadata is retained for comparisons but does not influence these fits.

Gene validation verifies the site factors and equality of scaling before/after
the existing annotation-based aggregation. No gene-only refit or additional CPM
is applied. For shared-positive genes:

```
log2FC_zinfnb = log2FC_raw + log2(factor_treated/factor_control)
```

Consequently, fixed-set ranks, spread and correlations are unchanged by this
scalar method, although direction and fold-threshold counts may change. Zeros
remain zero; cross-method diagnostics report dropout/gain separately instead
of inventing a common pseudocount for incompatible count units.

## Diagnostics and catalogue sensitivity

Factors include fitted pi, mu, n, p, observed/fitted zero fractions, mean negative
log likelihood, convergence, successful starts, boundaries and maximum discrete
CDF discrepancy. The discrepancy compares both sides of observed count jumps;
it is descriptive, not a calibrated goodness-of-fit p-value.

Adding zero-only rows changes the likelihood and fitted pi. When the optimum
has interior pi, the positive-count distribution can constrain mu,n independently
of the added zeros, leaving factors nearly unchanged. At boundaries or weakly
identified fits, scale estimates can also change. Catalogue sensitivity must be
measured, not assumed absent or inevitable. Never infer saturation from the fit.

## Execution and outputs

```sh
python scripts/local/ytab.py normalize \
  --catalogue output/normalization/gale_flu_ttr_validation_v1/matrix \
  --output output/normalization/gale_flu_zinfnb_sites_v1 \
  --group-by library_id --method zinfnb

python scripts/local/ytab.py normalize-genes \
  --validation output/normalization/gale_flu_ttr_validation_v1 \
  --output output/normalization/gale_flu_zinfnb_gene_validation_v1 \
  --existing-fitness output/projects/Gale_2020/treated_vs_parent/Cg1_FLU_vs_no_drug_raw_cpm \
  --site-normalization output/normalization/gale_flu_zinfnb_sites_v1
```

Zinc uses analogous `zn_` paths and its existing fitness result, with explicit
`--compare-independent-pools` for descriptive within-background agreement.
Independent-pool mock differences are not technical control-noise replicates.
New output directories are required. Exports retain matrices, factor diagnostics,
method descriptions, source/code hashes, versions, gene metrics and plots.

Tests cover the proper mixture likelihood and total probability, synthetic
parameter recovery, deterministic fitting, preserved zeros/raw counts, grouping
independence, catalogue changes, invalid inputs, optimizer failure, provenance,
and site-to-gene aggregation. Actual-data evaluations are documented in
[Gale results](../output/normalization/gale_flu_zinfnb_gene_validation_v1/INTERPRETATION.md)
and [zinc results](../output/normalization/zn_zinfnb_gene_validation_v1/INTERPRETATION.md).
No new p-values, FDR, z-score fits, or significant fitness calls are generated.
