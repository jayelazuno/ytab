# Zaeske et al. 2026 - prototrophic C. glabrata Hermes pools (PRJNA1460454)

Zaeske AM et al. (2026) "Facile generation of Hermes insertion mutants in prototrophic *Candida glabrata* for use in
nutrient-limited environments", Microbiol Spectr 14(7):e0084826, doi:10.1128/spectrum.00848-26 (Cunningham
laboratory).

## Two enrichment schemes, and why it matters

- **Old scheme** (Gale et al. 2020; Nickels et al. 2024, 2025; and the deposit's `BG14` run): a *ura3* auxotroph
  carries a centromeric plasmid with the inducible transposase, NAT1 inside the transposon and a counter-selectable
  URA3. Transposants are enriched on **nourseothricin + 5-FOA**, which selects against cells keeping the plasmid.
  Every resulting mutant is *ura3*, which the paper shows is outcompeted by URA3+ prototrophs in the mouse gut.
- **New scheme** (this paper, the `BG2r`, `CBS138r` and `DPL1021r` runs): a recessive **rpl28-Q38E** cycloheximide
  resistance mutation goes into the chromosome of a **prototrophic** strain, and the launchpad plasmid
  (pCR-MET3-Hermes) carries wild-type **RPL28**, which is dominant and restores CHX sensitivity, in place of URA3.
  Transposants are enriched on **nourseothricin + cycloheximide**. The rpl28-Q38E marker itself is fitness-neutral
  in the mouse gut, and the scheme is portable across strains and clades.

## The cost of the new scheme: PDR1 and CDR1 become false essentials

PDR1 and CDR1 are required for cycloheximide resistance, so the CHX enrichment removes every mutant insert in them.
The authors state this ("all the strains lacked insertions in the PDR1 and CDR1 genes"). In YTAB the consequence is
direct: the label-free HMM calls both genes **essential** in every CHX-enriched library, and does not in any
FOA-enriched library of the same species.

| Library | Scheme | PDR1 | CDR1 |
|---|---|---|---|
| Zaeske BG2r + DPL1021r | CHX, prototroph | **ES** | **ES** |
| Zaeske CBS138r | CHX, prototroph | **ES** | **ES** |
| Gale 2020 BG14 | FOA, *ura3* | NE | GD |
| Nickels & Cunningham 2025 BG2 | FOA, *ura3* | NE | GD |
| Nickels et al. 2024 CBS138 | FOA, *ura3* | NE | NE |

Neither gene is essential: Nickels & Cunningham 2025 build and screen *pdr1* deletion pools, and this paper's own
*pdr1* mutant colonises mice. This is a selection artefact of the enrichment marker, of the same kind as the SATAY
donor-locus jackpot (R25), and it must be masked before essentiality is read from a CHX-enriched pool. It is a
recommendation for Table 2, keyed to the enrichment scheme rather than to the species.

## Runs

Five runs. `SRR38338547` is labelled BG14, the old-scheme *ura3* comparator, and duplicates SRR33002022 of Nickels &
Cunningham 2025; it is counted once under the earlier deposit and is not in the analysed set. The four analysed
libraries are the prototrophic ones: BG2r (two pools), DPL1021r and CBS138r. DPL1021 has no assembly of its own and
is mapped to the BG2 reference; CBS138r is in its own project against the CBS138 reference.

Measured density of the analysed libraries (R24): BG2r 16.5 and 13.1, CBS138r 6.9, DPL1021r 3.7 insertion sites per
coding kb, so only the BG2r pools approach the density the essentiality recommendation asks for.

The paper also reports that subtelomeric DNA is highly accessible to insertion in CBS138 and DPL1021 but not BG2,
tracking natural variation in the silencing gene SIR3.
