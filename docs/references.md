# References

YTAB adapts several methods from TRANSIT. Where a YTAB method is adapted from TRANSIT, its documentation and code
cite the method and, for code, the source file at the inspected commit
(https://github.com/ioerger/transit, commit 50e23d6bec2605a24835fbd67f33691e330f3361). Each adaptation's
documentation says what was kept and what was changed for Hermes.

| YTAB method | Adapted from | Reference |
|---|---|---|
| TTR, totreads, nzmean, quantile, ZINB normalization (`ytab normalize`) | TRANSIT normalization methods (`norm_tools.py`) | DeJesus et al. 2015 |
| BGC and aBGC normalization (`ytab normalize`) | TRANSIT `BetaGeomNorm` and `AdaptiveBGCNorm` (`norm_tools.py`) | DeJesus & Ioerger 2016 |
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
- DeJesus MA, Ioerger TR. Normalization of transposon-mutant library sequencing datasets to improve identification of
  conditionally essential genes. *J Bioinform Comput Biol* 2016;14(3):1642004. doi:10.1142/S021972001642004X;
  PMID 26932272.
- DeJesus MA, Nambi S, Smith CM, Baker RE, Sassetti CM, Ioerger TR. Statistical analysis of genetic interactions in
  Tn-Seq data. *Nucleic Acids Res* 2017;45(11):e93. doi:10.1093/nar/gkx128; PMID 28334803 (verified 2026-10-21).

Methods adapted from the Hermes Tn-seq literature rather than from TRANSIT:

| YTAB method | Adapted from | Reference |
|---|---|---|
| Local-abundance z-score (`ytab fitness`) | the treated-vs-control z-score of the *C. glabrata* Hermes screens | Gale et al. 2020 |
| MidLC complexity, depth correction and jackpot fraction (`ytab library-diagnostics`, `ytab qc`) | the MidLC half-saturation measure and depth correction | Gale et al. 2020 |

Published datasets reprocessed by YTAB (sources, accessions and trimming: [Public datasets](public_datasets.md)):

- Gale AN, Sakhawala RM, Levitan A, Sharan R, Berman J, Timp W, Cunningham KW. Identification of essential genes and
  fluconazole susceptibility genes in *Candida glabrata* by profiling Hermes transposon insertions. *G3* 2020;10:3859–3870.
- Grech L, Jeffares DC, Sadée CY, Rodríguez-López M, Bitton DA, Hoti M, et al., Bähler J. Fitness landscape of the
  fission yeast genome. *Mol Biol Evol* 2019;36(8):1612–1623. doi:10.1093/molbev/msz113.
- Guo Y, Park JM, Cui B, Humes E, Gangadharan S, Hung S, FitzGerald PC, Hoe KL, Grewal SIS, Craig NL, Levin HL.
  Integration profiling of gene function with dense maps of transposon integration. *Genetics* 2013;195:599–609.
- Lee SY, Hung S, Esnault C, Pathak R, Johnson KR, Bankole O, Yamashita A, Zhang H, Levin HL. Dense transposon
  integration reveals essential cleavage and polyadenylation factors promote heterochromatin formation. *Cell Rep*
  2020;30(8):2686–2698. doi:10.1016/j.celrep.2020.01.094.
- Michel AH, Hatakeyama R, Kimmig P, Arter M, Peter M, Matos J, De Virgilio C, Kornmann B. Functional mapping of yeast
  genomes by saturated transposition. *eLife* 2017;6:e23570.
- Segal ES, Gritsenko V, Levitan A, Yadav B, Dror N, Steenwyk JL, et al., Berman J. Gene essentiality analyzed by in
  vivo transposon mutagenesis and machine learning in a stable haploid isolate of *Candida albicans*. *mBio*
  2018;9:e02048-18. doi:10.1128/mBio.02048-18.

Other methods cited in YTAB documentation:

- Leys C, Ley C, Klein O, Bernard P, Licata L. Detecting outliers: do not use standard deviation around the mean, use
  absolute deviation around the median. *J Exp Soc Psychol* 2013;49(4):764–766. (QC table, project-relative flags.)
- McCarthy DJ, Smyth GK. Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics*
  2009;25(6):765–771. (Genetic interactions, ROPE option.)
- Holm S. A simple sequentially rejective multiple test procedure. *Scand J Stat* 1979;6(2):65–70. (Genetic
  interactions, family-wise error option.)
- Efron B. Large-scale simultaneous hypothesis testing: the choice of a null hypothesis. *J Am Stat Assoc*
  2004;99(465):96–104. (Genetic interactions, optional empirical null.)
- Benjamini Y, Heller R. Screening for partial conjunction hypotheses. *Biometrics* 2008;64(4):1215–1222.
  PMID 18261164. (Pooled Gumbel, partial-conjunction null across pools.)
- Lund SP, Nettleton D, McCarthy DJ, Smyth GK. Detecting differential expression in RNA-sequence data using
  quasi-likelihood with shrunken dispersion estimates. *Stat Appl Genet Mol Biol* 2012;11(5):Article 8. PMID 23104842.
  (Pooled gene quasi-likelihood, empirical-Bayes dispersion shrinkage.)
