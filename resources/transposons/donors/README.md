# Transposon donors: counting the reads that never left the donor

Induction does not always lose the donor. Reads that start at a transposon end and run into the donor's flanking
sequence are sequenced with the genomic junctions, and until now the pipeline dropped them without counting: the ones
whose flank is donor-only sequence (vector, foreign markers, terminators) fail to map, and the ones whose flank the
donor shares with the host genome (a promoter, marker or CEN/ARS fragment) are placed at that locus, where they look
like insertions. `ytab.mapping.donor_reads` counts both during mapping, and QC reports them under Reads & Mapping as
**Reads from the transposon donor (%)** and **Donor reads placed on the genome (%)**.

Two ways to count, both written to `<sample>.donor_reads.csv`:

1. **Donor sequence known** (`mapping.donor_fasta` in the project config; `scripts/local/ytab_attach_donors.py` sets
   it): the cleaned junction reads are aligned to the donor, and each read is donor-only, shared (donor and genome) or
   genomic.
2. **No donor sequence** (always computed): a donor junction is one fixed sequence, so donor-only reads share their
   first bases. Reads unmapped to the genome whose first 20 bases recur in at least 0.01% of reads are counted. An
   estimate that needs no sequence, for any transposon.

## The donors here

| folder | donor | source |
|---|---|---|
| `hermes_pCU-MET3-Hermes` | pCU-MET3-Hermes (10,679 bp) | **the real plasmid**, Cunningham lab (pLannotate-annotated GenBank) |
| `hermes_pHL2577_pombe` | pHL2577 (URA3, kanMX6) + pHL2574/pHL2578 | rebuilt from parts (Guo 2013; Lee 2020; Billmyre 2022) |
| `acds_pBK549_satay` | pBK549 / pBK257: ADE2::MiniDs + GAL1-Ac | rebuilt from parts (Michel 2017; Kingma 2025) |
| `acds_pBB51_neoformans` | pBB51 Ds-NEO minitransposon, pBB65 at URA5 | rebuilt from parts (Billmyre 2025) |
| `piggybac_pIMY015_phaffii` | pIMY015 + pIMY029 (PAOX1-PBase) | rebuilt from parts (Yao 2026) |

Only the Cunningham-lab plasmid is a real sequence; none of the others is deposited, so each is **rebuilt from its
published parts** with `scripts/local/ytab_build_donor.py` from a `recipe.yaml` (GenBank accessions in
`components/`, the transposon reference, and host-genome regions taken from the species' own reference). Each folder
keeps `recipe.yaml`, `donor.fasta` and `parts.json` (every part, its source and length).

## How good is a rebuilt donor?

Rebuilding pCU-MET3-Hermes from eight parts and comparing it with the real plasmid
(`hermes_pCU-MET3-Hermes/rebuild/validation.json`):

- **90.5%** of the real plasmid's 100-bp windows are covered by the parts;
- **100%** of the 60 windows within 150 bp of a transposon end are covered — and that is the sequence donor reads
  actually contain;
- on the same library, the rebuilt donor counts **6.13%** donor reads against **6.56%** for the real plasmid (93% of
  the count).

So a rebuilt donor measures the donor fraction to within about a tenth of its value, and the number should be read as
a close estimate rather than an exact count. Where a laboratory can supply the real plasmid, drop it in and rebuild
nothing.

## Measured on real libraries

Every donor was tested on a sample of reads from a library of its own study
(`scripts/local/ytab_donor_reads.py --fastq`, 50,000-100,000 reads from ENA):

| donor | library tested | donor reads | donor-only | placed on the genome |
|---|---|---|---|---|
| pCU-MET3-Hermes (real) | C. glabrata pool, this lab | 6.63% | 6.40% | 0.23% |
| pHL2577 (rebuilt) | Lee 2020, SRR8502440 | 0.79% | 0.79% | 0.00% |
| pBK549 (rebuilt) | Michel 2017 wild type 1, ERR1533147 | 13.36% | **8.62%** | **4.74%** |
| pBB51 (rebuilt) | Billmyre 2025, SRR29762409 | 1.96% | 1.96% | 0.00% |
| pIMY015 (rebuilt) | Yao 2026, SRR36343507 | 0.03% | 0.03% | 0.00% |

The SATAY library is the clearest case: 8.6% of its reads never left the donor and were discarded unseen, and a
further 4.7% are placed on the genome at the donor's own ADE2, URA3 and GAL1-10 fragments, where they are
indistinguishable from insertions. Those are the loci the field masks by hand; the count measures them instead.
pHL2577 is low because the S. pombe protocol selects against the donor before sequencing, and the K. phaffii number
is near zero - which may equally mean that the rebuilt pIMY015/pIMY029 misses the construct actually used, so treat a
near-zero from a rebuilt donor as "nothing found", not as proof of absence.

**Why there can be so many, and where they go.** Every read starts at the transposon junction by design, so the only
question is what follows it. When the donor was lost, the flank is genomic and the read maps; when it was not, the
flank is plasmid backbone, which is not in the genome, so the read fails to map and the pipeline discards it. It never
reaches a MAPQ filter and never becomes an insertion. In the eight libraries measured here, donor reads account for
56-95% of everything the aligner could not place (yH299-parent-pool4: 21.7% donor against 22.9% unmapped, i.e. 95%),
while the fraction that does map - the donor's own host fragments - stays at 0.1-0.3%. Aligning those reads back to
the plasmid shows them piled where that story predicts: 83% of them land within 250 bp of the masked Hermes end, the
sequence immediately outside the transposon.

**Reads that keep the transposon end.** The donor carries the transposon, so such a read would match any donor
whatever its flank. The transposon's own sequence is therefore masked inside the donor before indexing, and alignment
is local with at least 30 aligned bases, so only the flank can produce a hit. Lee 2020's reads start with the Hermes
end and still give 0.79%, not 100%, which is the check that this works.

## Projects without a donor sequence

The piggyBac studies (Chow 2024, Song 2026, Gao 2021) deliver the transposon from a **chromosomal donor at ARG4**, not
a plasmid, and Segal 2018 (C. albicans Ac/Ds) and Patterson 2018 (Y. lipolytica Hermes) do not describe a construct we
can rebuild. These keep the flank estimate, which needs no sequence.
