# Public Tn-seq datasets

YTAB reprocesses published transposon libraries from raw reads, so every species and transposon goes through the same
mapping, hit-file, SummaryTable and analysis steps. These datasets form the four-species essentiality benchmark
(C. glabrata, C. albicans, S. pombe, S. cerevisiae; Hermes and Ac/Ds) and supply condition, mutant and time-course
libraries for the condition-comparison, genetic-interaction and multi-condition methods.

Each dataset is described in `pipeline/hpc/datasets/<dataset>/`: `samples.csv` (one row per sequencing run, keyed by
its accession, with role, condition, background, pool, transposon and library design) and `dataset.env` (project ID
and species). The same `samples.csv` drives the download and becomes the project's sample sheet.

## Datasets

| dataset (project) | species, strain | transposon | runs | accession | reference |
|---|---|---|---|---|---|
| `gale_2020_glabrata` (`Gale_2020`) | C. glabrata BG2 | Hermes | 7: 3 parents, 2 no-drug controls, 2 fluconazole | SRA SRX8130913–SRX8130919 | Gale et al. 2020 |
| `segal_2018_albicans` (`Segal_2018_albicans`) | C. albicans YJB-T1081, stable haploid (SC5314) | Ac/Ds (Ds-NAT1) | 11 libraries (01–11) | BioProject PRJNA490565 (SRR7824838–SRR7824848) | Segal et al. 2018 |
| `pombe_hermes` (`Pombe_Hermes`) | S. pombe | Hermes | 26 (25 analysed) | see below | Guo et al. 2013; Grech et al. 2019; Lee et al. 2020 |
| `michel_2017_satay` (`Michel_2017_SATAY`) | S. cerevisiae BY4741/BY4742 ade2Δ | Ac/Ds (MiniDs, SATAY) | 7 | ArrayExpress E-MTAB-4885 (ENA PRJEB14804, ERR1533143–ERR1533149) | Michel et al. 2017 |

**S. pombe runs:**
- Guo et al. 2013 (SRR327340): the in vivo Hermes library; the one Segal et al. 2018 used. Its published insertion
  table is kept in `resources/species/pombe/tnseq_data/guo_2013_hermes/hermes_hits.csv`.
- SRR836165: Hermes transposition into naked DNA in vitro. It is not a fitness library and is left out of analyses
  (`include` False); it records target-sequence preference without chromatin or selection.
- Grech et al. 2019 (PRJEB27324, ERR2641199–ERR2641214): 16 paired-end high-density libraries; the authors also
  called constrained regions with a five-state HMM.
- Lee et al. 2020 (PRJNA517661, SRR8502433–SRR8502440): WT and otr1R::ura4 strains (two each) at generation 0
  (`control`) and 80 (`treated`), a time-course comparison.

**S. cerevisiae (SATAY) libraries** (Michel et al. 2017, Table 1):

| run | library | role in YTAB |
|---|---|---|
| ERR1533147 | Wild-type 1 | parent (pool WT1) |
| ERR1533148 | Wild-type 2 | parent (pool WT2) |
| ERR1533149 | Wild-type 2 + rapamycin | treated, condition rapamycin (pool WT2: grown from Wild-type 2) |
| ERR1533143 | dpl1Δ | treated, condition dpl1del |
| ERR1533144 | dpl1Δ psd2Δ | treated, condition dpl1del_psd2del |
| ERR1533146 | VPS13(D716H) | treated, condition VPS13_D716H |
| ERR1533145 | mmm1Δ VPS13(D716H) | treated, condition mmm1del_VPS13_D716H |

The mutants are `treated` so they pair with the wild-type libraries; dpl1Δ / dpl1Δ psd2Δ and VPS13(D716H) /
mmm1Δ VPS13(D716H) are genetic-interaction sets. Table 1 also lists a YEN1on library, which was not deposited.
Wild-type 1 and VPS13(D716H) were sequenced on two MiSeq runs each; ENA serves each as one run.

## Read structures and trimming

Each library design is a transposon profile entry (`resources/transposons/<transposon>/profile.json`); the sample
sheet's `transposon` and `library_design` columns select it. Sequences are written 5′→3′ as they appear in the read.

| transposon / design | read (junction mate) | 5′ trimming | 3′ trimming | datasets |
|---|---|---|---|---|
| hermes / `splinkerette_pcr` | Hermes end, genome | automatic Hermes cleanup (Hermes ends and Gale primers, `primers.json`) | adaptors in `primers.json` | Gale 2020 |
| hermes / `end_ligated_adapter` | Hermes end, genome, adapter, tag | `GTGGCTTGTTGAAGTTCTCTG` (end of Hermes; overlap 21; reads without it discarded) | `GTCCCTTAAGCGGAG`; `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` | Grech 2019, Lee 2020 |
| hermes / `terminal_10bp` | last 10 bp of the Hermes end, genome (50-bp reads) | `^AAGTTCTCTG` (anchored; reads without it discarded) | `GTCCCTTAAGCGGAG`; `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` | Guo 2013 |
| ac_ds / `ds_primer_tail` | 6-nt tag, Ds primer and end, genome | `GTATTTTACCGACCGTTACCGACCGTTTTCATCCCTA` (overlap 37; reads without it discarded) | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` | Segal 2018 |
| ac_ds / `satay_restriction` | genome from the MiniDs junction to a DpnII (GATC) or NlaIII (CATG) site, then circularised Ds | none | `GTATCGGTTTTCGATTACCGTATTTATCCCGTTCG` (after GATC); `TGTGATTTTACCGAACAAAAATACCGGTTCCCGTCCG` (after CATG); `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` | Michel 2017 |

Sources of the sequences: the Ds primer and end are those of the authors' pipeline (github.com/berman-lab/transposon-
pipeline, `MapFastq.py`); the Hermes end is the Hermes sequence (GenBank L34807.1) and the Gale sequencing primer;
the ligated adapter and the SATAY Ds sequences were read from the libraries themselves. All 3′ sequences of a design
are trimmed in one cutadapt pass, each read at the first one found. Reads are then mapped with Bowtie2 and filtered
at MAPQ 20, as for every YTAB library; the insertion is the first genomic base of the junction read.

## Validation

On read samples of each dataset (about 50,000–80,000 reads each), mapped with the designs above:

| library | reads with the expected sequence | aligned |
|---|---|---|
| Segal 2018 library 03 | 74% | 93% |
| Guo 2013 | 99% | 80% |
| Lee 2020 (10160B, T0) | 85% | 78% |
| Grech 2019 HL_D_2 / V3 | 99% | 89% / 44% (V3 is mostly donor-plasmid reads) |
| SATAY Wild-type 1 / + rapamycin | — | 73% / 51% |

- **Guo et al. 2013:** every YTAB insertion matches the authors' published position and strand.
- **SATAY:** YTAB insertions match the published wild-type track at genomic sites, with a constant 1-bp offset on the
  C strand (the track's convention).
- **Donor loci:** reads from the SATAY donor plasmid pile up on ADE2, URA3, CYC1 and GAL1-10, and the albicans Ds donor
  sits in the ADE2 promoter. These loci are masked in essentiality benchmarks.

## Downloading and processing (Argon)

```bash
bash pipeline/hpc/ytab_public_dataset.sh <dataset> download   # one SGE job; NCBI SRA, ENA fallback; resumable
bash pipeline/hpc/ytab_public_dataset.sh <dataset> status
bash pipeline/hpc/ytab_public_dataset.sh <dataset> submit --profile core --threads 8
```

The download job (`pipeline/hpc/ytab_download_public.txt`) resolves each accession and its read count through the
ENA portal API, then fetches each run from NCBI SRA (prefetch, vdb-validate, fasterq-dump, pigz; ERR runs are
mirrored there), accepting it only when the archive validates and the read count equals ENA's. A run that fails on
SRA is downloaded from ENA with MD5 checks (`SOURCE=ena` uses ENA for everything; it is slow from Argon). Both give
identical reads. FASTQs go to `/nfsscratch/jayelazuno/workspace/public_data/<dataset>/fastq/<study>/<run>/`
(`DATA_ROOT` changes the base). Paired runs keep only their two mate files. `submit` creates the project from the
manifest and submits reference preparation, mapping, hit files, SummaryTable, library diagnostics, report and export.
On any other machine, download the same accessions and run `python scripts/local/ytab.py run --profile core` on a
project initialized with the manifest as its sample sheet.

## References

- Gale AN, Sakhawala RM, Levitan A, Sharan R, Berman J, Timp W, Cunningham KW (2020). Identification of essential
  genes and fluconazole susceptibility genes in *Candida glabrata* by profiling Hermes transposon insertions. G3
  10:3859–3870.
- Grech L, Jeffares DC, Sadée CY, Rodríguez-López M, Bitton DA, Hoti M, et al., Bähler J (2019). Fitness landscape of
  the fission yeast genome. Mol Biol Evol 36:1612–1623. doi:10.1093/molbev/msz113.
- Guo Y, Park JM, Cui B, Humes E, Gangadharan S, Hung S, FitzGerald PC, Hoe KL, Grewal SIS, Craig NL, Levin HL
  (2013). Integration profiling of gene function with dense maps of transposon integration. Genetics 195:599–609.
- Lee SY, Hung S, Esnault C, Pathak R, Johnson KR, Bankole O, Yamashita A, Zhang H, Levin HL (2020). Dense transposon
  integration reveals essential cleavage and polyadenylation factors promote heterochromatin formation. Cell Rep
  30:2686–2698. doi:10.1016/j.celrep.2020.01.094.
- Michel AH, Hatakeyama R, Kimmig P, Arter M, Peter M, Matos J, De Virgilio C, Kornmann B (2017). Functional mapping
  of yeast genomes by saturated transposition. eLife 6:e23570.
- Segal ES, Gritsenko V, Levitan A, Yadav B, Dror N, Steenwyk JL, et al., Berman J (2018). Gene essentiality analyzed
  by in vivo transposon mutagenesis and machine learning in a stable haploid isolate of *Candida albicans*. mBio
  9:e02048-18.
