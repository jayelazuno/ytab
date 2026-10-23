# Identifying single Hermes mutants by nested arbitrary PCR

*Status: draft from the lab protocol `docs/20261016-hermes-mutant-screen.docx`, to be polished for the manuscript.
Items marked **[to confirm]** need a check at the bench.*

After a phenotypic screen, individual mutants are picked and their Hermes insertion site is found by two rounds of
nested arbitrary PCR (AP-PCR) followed by Sanger sequencing. This replaces the sonication and splinkerette-adaptor
ligation of the Tn-seq library preparation with arbitrary genomic priming. Only standard PCR and Sanger sequencing are
needed. The genomic flank in the read is then mapped by YTAB (**Screen Planning > Mutant Identification**) to give
the insertion coordinate and the disrupted gene.

Adapted from the arbitrary-PCR strategy of Das et al. (2005) and Martínez-García et al. (2014), and the nested AP-PCR
genotyping of individual transposon mutants of Alejaldre et al. (2026). It keeps the Hermes primers already validated
in the C. glabrata Tn-seq workflow.

## Primers

| primer | role | sequence (5′→3′) |
|---|---|---|
| HR-1 | Hermes-specific outer primer (PCR1) | AATACTTGCACTCAAAAGGCTTGACACC |
| ARB6 | degenerate arbitrary primer (PCR1) | GGCACGCGTCGACTAGTACNNNNNNNNNNACGCC |
| HR-2 (= "Sequencing_F1" in the protocol) | Hermes-specific nested primer (PCR2) and the Sanger sequencing primer | ACACTTGTGCTTATCTATGTGGC |
| ARB2 | nested arbitrary primer (PCR2); matches the 5′ tail of ARB6 | GGCACGCGTCGACTAGTAC |
Sanger sequencing uses **HR-2**, which the protocol also calls Sequencing_F1 (confirmed by the bench scientist,
2026-10-23). **[to confirm]** The protocol's draft lists a different sequence under Sequencing_F1
(GGCTTACGTTTGCCTGTGGCTTGTTGAAGTTCTCTG); remove or correct that line in the docx.

## Template

- **Genomic DNA** from each mutant, 20 ng per 50 µL PCR1. Alternatively use 5 µL of cell lysate.
- **Optimization round** (10-06-2026): DNA from 10 alkalinization mutants, at 23–80 ng/µL.

| mutant | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| ng/µL | 30 | 25 | 25 | 23 | 40 | 80 | 67 | 55 | 65 | 49 |

## PCR1: arbitrary priming (HR-1 + ARB6, OneTaq Master Mix, 50 µL)

| step | temperature | time | cycles |
|---|---|---|---|
| initial denaturation | 94 °C | 30 s | 1 |
| denaturation / **low-stringency annealing** / extension | 94 °C / **30 °C** / 68 °C | 30 s / 30 s / 90 s | 6 |
| denaturation / annealing / extension | 94 °C / 45 °C / 68 °C | 30 s / 30 s / 90 s | 30 |
| final extension | 68 °C | 5 min | 1 |
| hold | 4 °C | ∞ | |

Keep the 30 °C annealing: it lets the degenerate 3′ end of ARB6 prime in the genome. Do not raise it to the standard
OneTaq range.

## PCR2: nested amplification (HR-2 + ARB2, OneTaq, 50 µL; template 1 µL PCR1 or 40 ng purified PCR1)

| step | temperature | time | cycles |
|---|---|---|---|
| initial denaturation | 94 °C | 30 s | 1 |
| denaturation / annealing / extension | 94 °C / 52 °C / 68 °C | 30 s / 30 s / 90 s | 30 |
| final extension | 68 °C | 5 min | 1 |
| hold | 4 °C | ∞ | |

HR-2 is nested inside HR-1, and ARB2 anneals to the defined 5′ sequence that ARB6 introduced in PCR1.

## Product check and sequencing

1. **Gel:** run an aliquot of each PCR2 on an agarose gel. A smear or several bands can be normal for arbitrary PCR and
   does not mean failure. In the optimization round most mutants gave single bands (gel in the app: 100 bp ladder,
   M1–M9 and negative control).
2. **Sequence:** purify suitable products and Sanger-sequence them with HR-2.
3. **Analyse:** load the `.seq`/`.ab1` files, or a folder or `.tar` of them, into **Mutant Identification**.
   - The read is trimmed for quality and for the Hermes end.
   - The genomic flank is aligned to the reference.
   - The insertion is assigned to a gene, or called intergenic.
   - The result can optionally be linked to a project's screen results.

   Analysis method: `docs/mutant_identification.md`.

## Optimization result

The PCR adaptation worked. Mostly single bands were recovered and smearing was not a concern. The products were
purified and Sanger-sequenced, and their identification is the example in the app's Mutant Identification tab.

## References

- Das S, Noe JC, Paik S, Kitten T (2005). An improved arbitrary primed PCR method for rapid characterization of
  transposon insertion sites. J Microbiol Methods 63:89–94.
- Martínez-García E, Aparicio T, de Lorenzo V, Nikel PI (2014). New transposon tools tailored for metabolic engineering
  of Gram-negative microbial cell factories. Front Bioeng Biotechnol 2:46.
- Alejaldre L, Anhel AM, Grozinger L, Goñi-Moreno Á (2026). Automated workflow for genotyping individual transposon
  library variants. BMC Methods 3:3. PMID 41869086.
