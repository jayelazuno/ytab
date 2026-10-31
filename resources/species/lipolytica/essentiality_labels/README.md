# Curated Y. lipolytica essentiality labels

Built by `scripts/local/ytab_curate_lipolytica_labels.py` from *S. cerevisiae* (SGD) viability through 1:1 orthologs
only. Rules are in the script docstring; sources and SHA-256 hashes are in `sources.json`; counts in `summary.json`.

| tier | evidence | essential | nonessential |
|---|---|---|---|
| silver | *S. cerevisiae* null-mutant viability (SGD), 1:1 similarity match | 404 | 761 |

8,570 annotated genes; 1,165 labelled, 7,405 unlabelled. **This is the weakest of the three new species' label sets,
in two distinct ways, and both belong in any sentence written about it.**

- **No gold tier.** No *Y. lipolytica* deletion collection or knockout set was obtained, so there is no direct
  evidence for either class. Every label is inferred from the *S. cerevisiae* ortholog.
- **Similarity, not curated orthology.** *C. neoformans* and *K. phaffii* both came with published ortholog tables
  (Billmyre S4 Table; Zhu Dataset S6). Patterson et al. 2018 has no PubMed Central copy, so no equivalent table was
  openly available. The orthologs here come instead from the CLIB89 RefSeq annotation's per-CDS notes — "similar to
  uniprot|&lt;acc&gt; Saccharomyces cerevisiae &lt;systematic name&gt;", inherited from the Génolevures annotation.
  Only "similar to" and "highly similar to" are used; "weakly similar to" (1,054 notes) and "some similarities with"
  (351) are excluded, and the pair must be 1:1 in both directions. Including the weak qualifiers would raise the
  count to 1,855 labels, and that was not done.
- **Coverage is low:** 1,165 of 8,570 genes. 3,774 CDS notes name an *S. cerevisiae* systematic name, 2,369 of them
  with a strong qualifier, of which 1,221 are 1:1 and 1,165 have a curated SGD viability.

## Patterson et al. 2018 is a comparator, never labels

Patterson et al. (2018) Metab Eng 48:184–196 call about 22% of genes essential from a Poisson analysis of Hermes
insertion density — computed on PRJNA356612, **the same 33 sequencing runs this project maps**. Scoring a YTAB model
against those calls would compare two analyses of one dataset. They are reported as agreement only, in the way
`gale2020_comparator` is for *C. glabrata*
(`resources/species/glabrata/essentiality_labels/README.md`).
