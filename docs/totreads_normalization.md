# Total-read normalization: principle, mathematics and Hermes use

`totreads` is a selectable scalar normalization, not a recommendation over other
methods. It equalizes total insertion-site support within an explicit group.
This step connects the existing shared implementation to gene-level validation
and evaluates Gale's fluconazole experiment and the zinc experiment.

## Principle and exact calculation

For raw nonnegative insertion-support counts x_ij:

```
R_j = sum_i x_ij
target_group = arithmetic mean of R_j across samples in the group
factor_j = target_group / R_j
normalized_x_ij = factor_j * x_ij
```

Every normalized sample total equals its group's target. Each sample contributes
equally to that target. A single-sample group has factor 1; a sample with no
support fails clearly. There is no trim fraction, fitting, random seed, or fixed
target of 100. This follows TRANSIT's total-read scaling, without interpreting
Hermes coordinates as Himar1 TA sites.

Counts are observed insertion support, including coding and intergenic sites,
not total raw FASTQ reads or measured viable-mutant frequencies. Identical
chromosome/one-based insertion positions are merged across strands consistently
with the other methods. Raw inputs and zeros remain unchanged. Adding rows that
are zero in every sample does not change the factors. No saturated target-site
universe, extra motif weights, or gene-length weights are required.

Gale's four Cg1-derived assay cultures share one group. Zinc's four independent
pools retain separate treated/mock groups. The algorithms use metadata, not
project-specific branches. The common empirical coordinate catalogue aligns
samples and preserves provenance; its unobserved entries are not known
biologically inaccessible sites.

## Interpretation and limitations

Totreads is a candidate when total observed insertion support is an appropriate
relative scale reference. It assumes differences in that total should be scaled
away. High-count outliers, strongly enriched mutants, global biological changes,
and composition differences can alter that reference. It cannot distinguish
these from sequencing-depth changes or restore absent mutants. Equal totals are
a mathematical consequence, not evidence of improved biological accuracy.

Unlike nzmean, totreads equalizes totals rather than occupied-site means. Unlike
TTR, it does not trim dominant counts or use a trimmed density-based estimate.
Unlike the existing fitness CPM workflow, its denominator includes all insertion
sites, not only gene-assigned support. It must not be followed by CPM again.

Site scaling commutes with the existing gene-assignment aggregation:

```
normalized_gene_gj = factor_j * raw_gene_gj
log2FC_totreads = log2FC_raw + log2(factor_treated/factor_control)
```

The gene runner verifies those identities using site-derived factors; it does
not estimate a new factor from gene totals. Intergenic support can affect the
factor without contributing to a protein-coding gene. For a fixed shared-positive
gene set, scaling shifts log-ratio centering, not spread, ranks, or correlation.
Fold-threshold and direction classifications can change. Zeros/dropouts are
reported separately, without adding a cross-method pseudocount.

## Execute

```sh
python scripts/local/ytab.py normalize \
  --catalogue output/normalization/gale_flu_ttr_validation_v1/matrix \
  --output output/normalization/gale_flu_totreads_sites_v1 \
  --group-by library_id --method totreads

python scripts/local/ytab.py normalize-genes \
  --validation output/normalization/gale_flu_ttr_validation_v1 \
  --output output/normalization/gale_flu_totreads_gene_validation_v1 \
  --existing-fitness output/projects/Gale_2020/treated_vs_parent/Cg1_FLU_vs_no_drug_raw_cpm \
  --site-normalization output/normalization/gale_flu_totreads_sites_v1
```

Zinc uses analogous `zn_` paths and its existing
`Zn_toxicity_screen/treated_vs_parent/Zn_1_5mM_vs_mock` result. Add
`--compare-independent-pools` for descriptive within-background effect agreement;
independent-pool mock differences are not technical control noise. Use new output
directories for reruns. Matrices, factors, source hashes, plots, and diagnostics
are exported independently of existing analysis results.

## Validation

79 focused normalization/CLI tests passed. Dedicated tests cover hand-calculated
factors, unequal occupied-site counts, zero-padding invariance, outlier
sensitivity, group isolation, single-sample identity, empty-sample rejection,
and intergenic support in site-derived gene factors. Existing shared-method
tests include TRANSIT parity when the optional source checkout is available.

Actual Gale and zinc runs preserved raw input hashes and zero patterns, achieved
equal within-group normalized site totals, reproduced raw gene summaries, and
passed existing CPM/default-effect parity checks. No existing fitness results,
classifier workflow, app defaults, or normalization UI were changed. No new
p-values, FDR estimates, z-score fits, or significant fitness calls were made.

Results: [Gale](../output/normalization/gale_flu_totreads_gene_validation_v1/INTERPRETATION.md),
[zinc](../output/normalization/zn_totreads_gene_validation_v1/INTERPRETATION.md).
