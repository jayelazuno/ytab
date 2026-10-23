# References

YTAB adapts several methods from TRANSIT. Where a YTAB method is adapted from TRANSIT, its documentation and code
cite the method and, for code, the source file at the inspected commit
(https://github.com/ioerger/transit, commit 50e23d6bec2605a24835fbd67f33691e330f3361). Each adaptation's
documentation says what was kept and what was changed for Hermes.

| YTAB method | Adapted from | Reference |
|---|---|---|
| TTR, totreads, nzmean, quantile, BGC, aBGC, ZINB normalization (`ytab normalize`) | TRANSIT normalization methods (`norm_tools.py`) | DeJesus et al. 2015 |
| Per-library Tn-seq statistics (`ytab tnseq-stats`, `ytab qc`) and their QC guidance and plots | TRANSIT `tnseq_stats` and quality control | DeJesus et al. 2015 |
| Site-restricted resampling (`ytab resampling`) | TRANSIT resampling | DeJesus et al. 2015 |
| Multi-condition tests and plots (`ytab multicondition`) | TRANSIT ANOVA, U-test, corrplot and heatmap | DeJesus et al. 2015 |
| Library correlation (`ytab library-correlation`) | TRANSIT corrplot | DeJesus et al. 2015 |
| Genetic interactions (`ytab genetic-interactions`) | TRANSIT GI | DeJesus et al. 2017 |
| Hermes window HMM (`ytab hermes-hmm`) | TRANSIT HMM | DeJesus & Ioerger 2013 |
| Control-conditioned Gumbel | TRANSIT Gumbel | DeJesus et al. 2013 (Bioinformatics) |

- DeJesus MA, Ambadipudi C, Baker R, Sassetti C, Ioerger TR. TRANSIT — a software tool for Himar1 TnSeq analysis.
  *PLoS Comput Biol* 2015;11(10):e1004401.
- DeJesus MA, Ioerger TR. A Hidden Markov Model for identifying essential and growth-defect regions in bacterial
  genomes from transposon insertion sequencing data. *BMC Bioinformatics* 2013;14:303.
- DeJesus MA, Zhang YJ, Sassetti CM, Rubin EJ, Sacchettini JC, Ioerger TR. Bayesian analysis of gene essentiality
  based on sequencing of transposon insertion libraries. *Bioinformatics* 2013;29(6):695–703.
- DeJesus MA, Nambi S, Smith CM, Baker RE, Sassetti CM, Ioerger TR. Statistical analysis of genetic interactions in
  Tn-Seq data. *Nucleic Acids Res* 2017;45(11):e93. doi:10.1093/nar/gkx128; PMID 28334803 (verified 2026-10-21).

Other methods cited in YTAB documentation:

- Leys C, Ley C, Klein O, Bernard P, Licata L. Detecting outliers: do not use standard deviation around the mean, use
  absolute deviation around the median. *J Exp Soc Psychol* 2013;49(4):764–766. (QC table, project-relative flags.)
- McCarthy DJ, Smyth GK. Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics*
  2009;25(6):765–771. (Genetic interactions, ROPE option.)
