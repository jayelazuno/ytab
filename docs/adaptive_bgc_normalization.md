# Adaptive beta-geometric correction: principle and Hermes adaptation

`aBGC` is an optional experimental distribution correction, not a preferred
normalization. It searches fitting settings that make corrected positive counts
more geometric according to an explicit discrepancy score. It does not learn
which genes are correct, optimize biological agreement or estimate an FDR.
User selection should consider distribution assumptions, sample context,
candidate stability and treatment-effect changes. Recommendations are deferred
until broader testing of the method family is complete.

## Source and differences from ordinary BGC

The inspected reference is `AdaptiveBGCNorm` in TRANSIT's `norm_tools.py`,
commit `50e23d6bec2605a24835fbd67f33691e330f3361`, lines 320–412.
The implementation follows its positive-count parameter formulas, candidate
fraction grid and maximum-discrepancy score. It is not a claim of bitwise parity
with legacy random output or of independent validation in Hermes.

Unlike ordinary BGC, adaptive fitting excludes zeros, does not add 1 to counts,
and does not impose BGC's kappa>=10 floor. Its geometric correction target uses
each candidate's rho, not the group mean positive count. Each sample is fitted
independently; pool metadata and matched comparisons are retained. No motif,
gene-length, TA-site or viable-mutant-frequency model is introduced.

## Mathematics

Let z_r be the sorted positive counts in a sample, n their number, and
`q_r = (r-1)/n`, r=1,...,n. At least two occupied sites are required.

```
F = {0, 0.02, 0.04, ..., 0.30}         # trim fraction per tail
rho(f) = 1 / trim_mean(z, f)
A = sum(log(1-q_r)^2) / sum(z_r * log(1-q_r))
kappa(f) = (2*exp(A)-1)/(exp(A)+rho(f)-1)
alpha(f) = kappa(f)*rho(f)
beta(f)  = kappa(f)*(1-rho(f))
```

The same A is used for all candidates; rho and therefore kappa change with f.
Valid candidates require 0<rho<1 and finite positive kappa/alpha/beta. Invalid
fits or simulation failures are recorded with reasons and excluded. If no
candidate is valid, the run fails clearly. There is no fallback geometric model.

For each candidate, positive counts are mapped through the beta-geometric CDF
to a geometric target with that candidate's rho:

```
F_BG(x) = 1 - B(alpha, beta+x)/B(alpha,beta)
y_r(f) = Geom_quantile(F_BG(z_r), rho(f))
mean_y = mean(y_r(f))
g_r(f) = Geom_quantile(rank_probability_r, 1/mean_y)
score(f) = max_r((y_r(f)-g_r(f))^2/g_r(f))
```

The minimum-scoring valid candidate is selected; exact ties select the smaller
fraction. This is an optimum **within the supplied grid**, not a global optimum
or a statistical goodness-of-fit test. A boundary winner is reported as such
in the candidate table; no automatic grid expansion is performed.

## Tail convention is explicit

- `--abgc-tail endpoint` (default; `transit`, its earlier name, is still accepted) uses r/n for scoring quantiles. The final
  quantile is infinite; it is replaced with the largest finite quantile +100,
  retaining the inspected upstream rule. That constant is a numerical scoring
  convention, not a fitted biological parameter.
- `--abgc-tail midpoint` instead uses `(r-0.5)/n`, keeping every scoring quantile
  finite. This is a sensitivity option, explicitly different from upstream.

These choices can select different fractions. Both are implemented and tested;
neither is silently substituted or interpreted as proof of better biology.

## Numerical engine, boundaries and reproducibility

The default `--abgc-cdf analytic` uses the exact beta-geometric CDF. The optional
`--abgc-cdf simulation --abgc-seed 0 --abgc-draws 200000` uses the upstream-style
beta/geometric sampling approximation. Candidate/final RNG streams derive from
the sorted positive counts, base seed, candidate fraction and phase; they are
independent of column order and do not alter global RNG state. The final
simulation redraw mirrors the upstream search-then-redraw structure. Exact
legacy random-output parity is not claimed.

Zeros remain zero. A positive count remains positive even if finite simulation
gives ECDF=0. Endpoints are bounded and safeguard counts exported.

For the analytic engine, `F_BG(1)=rho` exactly. Floating log-beta subtraction can
otherwise place that percentile just above a geometric step, incorrectly mapping
one read to two. The implementation anchors that identity, snaps inverse
quantiles at floating-point integer boundaries, and enforces the theoretical
bound `corrected_count<=raw_count` before final scaling. This bound follows from
the beta mixture's survival function and Jensen's inequality. These numerical
rules are regression-tested and documented deviations from a naive legacy port.
Simulation can still cross discrete steps due to finite sampling; that is
reported in sensitivity diagnostics rather than hidden.

Optional final TTR is explicit and enabled by default, matching the upstream
workflow: 5% per tail and target 100. `--abgc-no-ttr` disables it. No CPM follows
the correction. Nonlinear site correction precedes aggregation with the existing
SummaryTable mapping. Adding zero-only catalogue rows leaves the positive fit
and unscaled map unchanged; final TTR's absolute scale changes proportionally
with the catalogue size.

## CLI and outputs

```sh
python scripts/local/ytab.py normalize \
  --catalogue output/normalization/gale_flu_ttr_validation_v1/matrix \
  --output output/normalization/gale_flu_abgc_sites_v1 \
  --group-by library_id --method aBGC \
  --abgc-cdf analytic --abgc-tail endpoint

python scripts/local/ytab.py normalize-genes \
  --validation output/normalization/gale_flu_ttr_validation_v1 \
  --output output/normalization/gale_flu_abgc_gene_validation_v1 \
  --existing-fitness output/projects/Gale_2020/treated_vs_parent/Cg1_FLU_vs_no_drug_raw_cpm \
  --site-normalization output/normalization/gale_flu_abgc_sites_v1

python scripts/local/ytab.py abgc-audit \
  --site-run output/normalization/gale_flu_abgc_sites_v1 \
  --output output/normalization/gale_flu_abgc_audit_v1 \
  --seed 0 --draws 200000
```

`--abgc-fractions` accepts an explicit comma-separated candidate grid. Fractions
must be distinct and within [0,0.5). The complete grid, engine, tail policy,
post-scaling settings and source/code hashes are recorded. Each sample has
all candidate rows in `aBGC/candidates.tsv`; status/reason, rho/kappa/alpha/beta/A,
score and selected flag are exported. `factors.tsv` records chosen parameters,
score, validity counts and post-TTR factor; the overall nonlinear factor is blank.

The audit reproduces analytic output, repeats fixed-seed simulation, tests a
10% zero-only catalogue expansion, and compares selected fractions under the
alternative tail rule. It exports sensitivity rather than inventing acceptance
cutoffs. In the actual data, candidate choices and low-count values can change
under finite simulation; the reported scores must not be overinterpreted.

70 focused normalization/CLI tests passed, covering estimators, manual tail-score
calculations, selection, invalid candidates, zero/tie/monotonic behavior,
one-read boundary precision, simulation order/reproducibility, zero-only row
invariance, candidate exports, audit and nonlinear gene aggregation. Gale and
zinc runs/audits passed, plots were inspected, and raw/CPM parity remained exact.

Results characterize behavior and assumptions, not a ranking:
[Gale](../output/normalization/gale_flu_abgc_gene_validation_v1/INTERPRETATION.md),
[zinc](../output/normalization/zn_abgc_gene_validation_v1/INTERPRETATION.md).
No Shiny UI, default method, accepted fitness result or classifier model changed.
