# Tn-seq statistics with Hermes QC

The output `tnseq_stats.tsv` retains TRANSIT's original column names/order and
tab-separated Excel-readable format. This is not an XLSX workbook. Hermes QC and
provenance columns are appended; the existing diagnostics are not replaced or
recomputed. No UI layout, normalization, classifier or fitness default changes.

## Mathematics and assumptions

For U empirical catalogue coordinates, n positive counts x and R=sum(x):

- `density=n/U`; this is empirical catalogue occupancy, not genome saturation.
- `mean_ct=R/U`, `NZmean=R/n`, `NZmedian=median(x)`.
- `max_ct=max(x)` and `total_cts=R` describe called-insertion support, not raw
  FASTQ yield or viable-mutant abundance.
- `skewness=m3/m2^(3/2)` and `kurtosis=m4/m2^2-3`, where mk is the mean kth
  central moment over positive counts. These retain SciPy/TRANSIT biased moment
  conventions; kurtosis is excess kurtosis. Constant/empty positives give blank
  undefined moments, not fabricated zeros.
- Pickands index is median(log2((y[m]-y[2m])/(y[2m]-y[4m]))) for m=10..99,
  descending counts y and upstream zero-based indices. Implicit zero counts
  follow positives. Only in-range terms with both positive spacings are used;
  `pti_valid_terms` reports the number. No valid terms means a blank estimate.

The inspected reference is TRANSIT `analysis/tnseq_stats.py` and
`tnseq_tools.get_data_stats` at commit 50e23d6bec2605a24835fbd67f33691e330f3361.
Unlike its formatted output, we retain precision, including fractional medians.
No Himar1 density, mean-count or tail-index acceptance thresholds are transferred
to Hermes. These are descriptive metrics, not significance tests or automatic
sample exclusions. Comparisons require consistent reference/catalogue identity;
zero-padding changes density and all-site mean. Positive moments are unchanged.

## Added Hermes fields

Existing summary columns are appended with `hermes_` prefixes: jackpot fraction,
top-site fraction/count, MidLC estimate and depth ratio, discovery-curve sampling
information, filtering/normalization settings and completion flags. Total support
and occupied-site counts already occur in the leading statistic fields.

All 16 saved +2/+7 base-pair fractions, background fractions and enrichment
ratios are copied where available. They retain the existing read-weighted,
coordinate-offset convention, not a new adjacent-TA or strand-aware motif test.
Saved discovery curves, motif tables, centromere profiles and TSS/TTS/tRNA
metaplots are linked by file paths. Curves remain curves rather than a fabricated
scalar. MidLC values are historical diagnostics, not validated new classifier
targets or saturated-library estimates.

Diagnostics join by explicit project/sample identity. Raw support totals and
occupied-site counts must match to reuse the saved metrics. Mismatched rows
retain their raw statistics but have saved metrics withheld and a clear status;
saved totals are included to expose the discrepancy. Historical input identity
is not revalidated by this join; that limitation is explicit in the status field.
Missing samples are labeled unavailable, not zero or failed. All reused files
and raw catalogue inputs are checksummed in the manifest. Original outputs stay
untouched. File links refer to the current filesystem and must accompany result
transfer if users want to inspect the source curves on another machine.

## CLI

```sh
python scripts/local/ytab.py tnseq-stats \
  --catalogue output/normalization/gale_our_parents_ttr_v1 \
  --diagnostics Gale_2020=output/projects/Gale_2020/library_diagnostics \
  --diagnostics H2O2_screen_v1=output/projects/H2O2_screen_v1/library_diagnostics \
  --diagnostics Zn_toxicity_screen=output/projects/Zn_toxicity_screen/library_diagnostics \
  --output output/qc/parent_tnseq_stats_v1
```

The generic runner accepts any existing raw catalogue and explicit diagnostic
directory mappings; there are no dataset-specific algorithm branches. It does
not parse strand-signed browser WIG as raw nonnegative support. Use the accepted
canonical hit-derived matrix. Future standalone WIG import should explicitly
handle coordinate/strand conventions rather than guess. New output directories
are required. The command can also be a step in the modular CLI master workflow.

Tests cover TRANSIT metric/PTI definitions, degenerate data, median precision,
original header order, saved-QC reuse, missing diagnostics and mismatch withholding.
