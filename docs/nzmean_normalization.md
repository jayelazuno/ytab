# Nonzero-mean normalization: principle, mathematics and Hermes use

`nzmean` is a user-selectable scalar method that equalizes mean insertion support
among occupied coordinates within an explicitly chosen group. It is not a
preference over the other methods. Its core implementation was already present;
this step adds experiment/gene-level validation and dedicated guidance.

## Principle and suitable data characteristics

The intended reference is comparable average support per **observed positive
site**, not equal total reads, equal insertion density or equal gene recovery.
It is a candidate when that positive-site mean is an appropriate scale reference
and different occupied-site counts should not themselves be normalized away.
That is a scientific assumption, not a property established by the calculation.

The arithmetic mean is sensitive to very high counts. Technical outliers,
genuine mutant enrichment and broad biological abundance changes can all alter
the factor. Nzmean cannot distinguish them. It does not correct skew, impute
missing mutants, infer viable-mutant frequencies, or account separately for
Hermes motifs or gene length. Unlike TTR it does not trim; unlike whole-site
total-read scaling it does not equalize totals when occupied-site numbers differ.

## Exact formula

For raw nonnegative insertion-support count x_ij at coordinate i in sample j:

```
n_j = number of coordinates with x_ij > 0
R_j = sum_i x_ij
mu_j_positive = R_j / n_j

target_group = arithmetic mean of mu_j_positive across the group's samples
factor_j = target_group / mu_j_positive
normalized_x_ij = factor_j * x_ij
```

Each sample contributes equally to the group target, regardless of its depth or
occupied-site count. Therefore every normalized positive-site mean equals
`target_group`, but its normalized total is `n_j * target_group` and may differ
from another sample's total. A one-sample group has factor 1. Empty samples fail
clearly instead of receiving a fallback factor. No random seed, fitting, trimming
fraction or target=100 parameter applies to nzmean.

The formula follows TRANSIT `NZMeanNorm`, inspected at commit
`50e23d6bec2605a24835fbd67f33691e330f3361`, `norm_tools.py` lines 13–51. Numerical
parity with that actual class is covered by the shared-method tests when the
local source checkout is available. TRANSIT's `TAs_hit` variable is interpreted
as occupied observed Hermes coordinates here, not as a Himar1 target motif.

## Coordinate and grouping rules

Use the same raw, nonnegative canonical hit counts and reference identity as the
other methods. Counts at identical chromosome/one-based positions are summed
across strands. Both coding and intergenic insertion support enter the factor;
gene assignment happens separately. A zero means no observed support, not an
impossible insertion or an essential gene. Zero-only rows do not affect nzmean
factors or transformed positive values, although the catalogue remains important
for coordinate alignment and provenance.

Gale's four Cg1-derived assay cultures form one group. Zinc retains four separate
treated/mock groups, one per pool. No independent library catalogues are flattened
to calculate a shared factor. Within a group, changing its target by adding
samples changes the absolute scale; for a fixed pair the relative factor remains
`mu_control_positive / mu_treated_positive`.

For gene-assignment matrix A, nzmean commutes with aggregation:

```
normalized_gene_gj = sum_i A_gi * factor_j * raw_site_ij
                     = factor_j * raw_gene_gj
```

The factor is still estimated from **sites**, not from positive gene totals.
Gene validation checks the expected site factor and this aggregation identity.
Intergenic counts can influence the factor without being assigned to genes.
No CPM or other normalization is applied afterward.

For a shared-positive gene's treated/control log ratio:

```
log2FC_nzmean = log2FC_raw + log2(factor_treated/factor_control)
```

Consequently, scalar scaling shifts log-ratio centering but cannot change its
spread, ranks or correlation across a fixed shared-positive gene set. Direction
and fold-threshold classifications can change. Diagnostics exclude zero ratios
without inventing a cross-method pseudocount, and report dropout/gain separately.

## Commands and outputs

```sh
python scripts/local/ytab.py normalize \
  --catalogue output/normalization/gale_flu_ttr_validation_v1/matrix \
  --output output/normalization/gale_flu_nzmean_sites_v1 \
  --group-by library_id --method nzmean

python scripts/local/ytab.py normalize-genes \
  --validation output/normalization/gale_flu_ttr_validation_v1 \
  --output output/normalization/gale_flu_nzmean_gene_validation_v1 \
  --existing-fitness output/projects/Gale_2020/treated_vs_parent/Cg1_FLU_vs_no_drug_raw_cpm \
  --site-normalization output/normalization/gale_flu_nzmean_sites_v1
```

Zinc uses analogous `zn_` paths, its existing
`Zn_toxicity_screen/treated_vs_parent/Zn_1_5mM_vs_mock` fitness result, and
`--compare-independent-pools` for descriptive within-background agreement.
That option does not turn independent-pool mock differences into technical noise.

The site runner exports raw totals, occupied-site counts, positive means,
groups, factors, normalized matrix and source hashes. The gene runner exports
raw/TTR/CPM/nzmean matrices, count-stratified metrics, plots and provenance.
No fitting parameters or synthetic randomness are hidden in this method.

## Verification and scope

75 focused normalization/CLI tests passed. New tests cover hand calculations,
unequal densities/totals, outlier sensitivity, zero-only rows, pool separation,
single-sample identity, intergenic support and site-derived gene factors. Existing
TRANSIT parity, raw summaries and CPM parity remain satisfied. Actual Gale/zinc
outputs were checked for equal within-group positive means, unchanged zero
patterns and unchanged input hashes. Plots were inspected.

The evaluations characterize behavior, not method rankings or biological ground
truth. No p-values, FDR estimates, control-noise refits or new significant fitness
calls are produced. No existing app defaults, fitness results or classifier
workflows were changed.

Results: [Gale](../output/normalization/gale_flu_nzmean_gene_validation_v1/INTERPRETATION.md),
[zinc](../output/normalization/zn_nzmean_gene_validation_v1/INTERPRETATION.md).
