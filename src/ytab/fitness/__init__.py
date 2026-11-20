"""Comparative analyses: what a treatment does to a library, relative to its controls.

Every method here needs at least two conditions. They divide by what they measure, which is why their calls overlap
only in part, and why `agreement` reports a rank correlation beside every recovery figure:

  zscore                      the gene's change in abundance against the noise expected at that abundance
  resampling                  the same change, tested by permuting condition labels over the gene's sites
  control_conditioned_gumbel  the longest run of consecutive sites whose treated count is improbably low
  multicondition              permutation ANOVA and all-pairs resampling over several conditions
  genetic_interactions        whether a gene's treatment response depends on the background
  agreement                   how any two of the above agree, per comparison

Those that score a magnitude (zscore, resampling) are interchangeable apart from error control; the Gumbel,
which scores contiguity, answers a different question and will not agree call for call. None of
them classifies a gene as essential: single-condition essentiality lives in `ytab.essentiality`.
"""
