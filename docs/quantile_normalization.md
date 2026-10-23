# Experimental zero-preserving quantile normalization for Hermes

This is an optional positive-count distribution sensitivity analysis, not a
default replacement for CPM or TTR. No Shiny UI or existing fitness results are
changed. The method assumes comparable positive-count distributions are an
appropriate reference within the explicitly chosen normalization group; real
treatment-induced distribution changes may violate that assumption.

## Relationship to TRANSIT

The inspected TRANSIT checkout is commit
`50e23d6bec2605a24835fbd67f33691e330f3361`, specifically `QuantileNorm` in
`src/pytransit/norm_tools.py`. That implementation sorts all sites, averages the
sorted samples and uses a dense-rank lookup. It is not the tied rank-interval
averaging implemented here. TRANSIT's own
[normalization documentation](https://transit.readthedocs.io/en/latest/method_normalization.html)
warns that quantile normalization does not work well on Tn-seq data with many
zero-count TA-site ties. Its availability is not evidence of suitability for
Hermes or for a given experiment.

The Hermes adaptation deliberately differs: zeros remain zero and only observed
positive counts define the reference. No TA or motif-defined opportunity layer
is added. The adaptation has its own mathematical/invariant tests, **not a claim
of numerical parity with TRANSIT QuantileNorm**. Existing TTR/totreads/nzmean
parity checks remain separate.

## Mathematics and parameters

For each sample j, sort its n_j positive raw coordinate counts as z_(r),j.
For n_j > 1, its empirical quantile function Q_j(p) linearly interpolates between
the knots `((r-1)/(n_j-1), z_(r),j)`. A single positive count defines a constant
quantile function. Within an explicit group of J samples:

```
Q_group(p) = sum_j Q_j(p) / J
```

Each sample contributes equally, regardless of depth or occupied-site count.
For sample j, evaluate that reference at its positive-rank probabilities
`p_rj = (r-1)/(n_j-1)`, or p=0.5 for a singleton. If the original positive count
is tied over ranks a through b, assign the mean of the reference values over
that entire rank block to each tied site. Restore original coordinate order.
All original zeros remain zero. A one-sample group is the identity transform.

This operation is deterministic, nonnegative, monotonic among positive counts,
preserves ties and requires no random seed. No trimming, target=100 or TTR
post-scaling is applied. It has no single sample normalization factor. Unequal
occupied-site numbers, interpolation and ties mean that **complete sample
distributions need not become identical**. In particular, insertion density is
not forced to match.

For gene assignment matrix A and transformed site matrix X_quantile:

```
gene_counts_quantile = A @ X_quantile
```

Normalizing raw gene sums instead would generally give different results and
is not implemented. CPM is not applied after quantile normalization.

## Gale evaluation commands

The same existing four-sample assay matrix and coordinate catalogue were reused;
no alignment, raw-count generation or parent catalogue was repeated.

```sh
python scripts/local/ytab.py normalize \
  --catalogue output/normalization/gale_flu_ttr_validation_v1/matrix \
  --output output/normalization/gale_flu_quantile_sites_v1 \
  --group-by library_id --method quantile

python scripts/local/ytab.py normalize-genes \
  --validation output/normalization/gale_flu_ttr_validation_v1 \
  --output output/normalization/gale_flu_quantile_gene_validation_v1 \
  --existing-fitness output/projects/Gale_2020/treated_vs_parent/Cg1_FLU_vs_no_drug_raw_cpm \
  --site-normalization output/normalization/gale_flu_quantile_sites_v1
```

All four Cg1-derived cultures jointly estimate one reference distribution.
Consequently, improved apparent agreement would not constitute independent
validation: it can partly result from distribution matching itself. The design
does not represent four independent library inductions. No zinc quantile run
was performed in this step.

## Output and validation

Site output includes the normalized sparse matrix, method description, raw input
hashes, group membership and implementation hashes. Gene output includes the
raw/TTR/CPM/quantile matrices, comparison and agreement metrics, raw-count strata,
plots, and `effect_retention.tsv`. The retention table tracks descriptive
same-direction twofold CPM effects through quantile normalization; CPM is a
reference for sensitivity, not biological ground truth. No significant-hit calls,
p-values, FDR estimates or control-noise model refits are produced.

52 focused tests passed, including hand-calculated ties, unequal positive counts,
singleton behavior, group separation, sample-order equivariance, scaling
equivariance, reproducibility, zero/raw preservation, nonlinear aggregation order
and rejection of mismatched site provenance. On actual Gale data, exact site-zero
patterns, monotonicity and unchanged source hashes were independently verified.
Raw summary and existing CPM parity remained exact. The generated plot was
visually inspected.

See [the Gale results report](../output/normalization/gale_flu_quantile_gene_validation_v1/INTERPRETATION.md).
The initial evaluation does not establish superiority over CPM or TTR.
