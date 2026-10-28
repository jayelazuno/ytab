# Patterson et al. 2018 - Yarrowia lipolytica Hermes (PRJNA356612)

Hermes saturation mutagenesis of the oleaginous yeast *Y. lipolytica* CLIB89 (W29), GEO GSE90993;
Patterson K et al. (2018) Metab Eng 48:184-196, doi:10.1016/j.ymben.2018.05.008.

Two independent trials (`pool` trial1, trial2). Each starts from the library before outgrowth (`gen0`, role parent)
and is grown to 20 and 80 generations in liquid YP with glucose or glycerol. A third pool (`facs`) was sorted by
BODIPY lipid signal: top 5% (`lipid_high`, treated) against bottom 95% (`lipid_low`, control).

Read structure: 100 bp single end, Hermes end `...TTGAAGTTCTCTG` at a fixed offset, then the genome; the
`splinkerette_pcr` design's automatic Hermes cleanup trims it (65% of the first 2,000 reads of SRR5084557 carry the
end; 65% of the trimmed reads map at MAPQ >= 20 to GCF_001761485.1).

Reference: `resources/species/lipolytica` (GCF_001761485.1 ASM176148v1, CLIB89(W29), annotation submitted by the same
laboratory; 6 chromosomes, 7,918 protein-coding genes).
