# Experimental designs of the public datasets

Which runs make which library, and which libraries each analysis uses, for every reprocessed study. Written by
`ytab datasets` from the sample sheets, `config/comparison_design.csv` and the analyses' manifests; search a run ID,
strain, drug or condition to find where it is used. From the command line: `ytab datasets --find TEXT`.

## Gale et al. 2020: Candida glabrata (Gale_2020)

BG14 (BG2 background); Hermes (NAT), splinkerette PCR. SRA SRX8130913–SRX8130919. 3 independent inductions (1 tube, 40 tubes, Cg1); 2 no-drug controls and 2 fluconazole cultures grown from Cg1.

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| SRR11561502 |  | treated | treated | BG14 | Cg1 |
| SRR11561503 |  | treated | treated | BG14 | Cg1 |
| SRR11561504 |  | control | control | BG14 | Cg1 |
| SRR11561505 |  | control | control | BG14 | Cg1 |
| SRR11561506 |  | parent | parent | BG14 | Cg3 |
| SRR11561507 |  | parent | parent | BG14 | Cg2 |
| SRR11561508 |  | parent | parent | BG14 | Cg1 |

**Analyses**

- **Essentiality (HMM)**, `essentiality_20261025_131741`
  - libraries: SRR11561508; SRR11561507; SRR11561506
  - `ytab essentiality --project-config output/projects/Gale_2020/config/project.yaml --libraries SRR11561508 SRR11561507 SRR11561506 --output output/projects/Gale_2020/essentiality/essentiality_20261025_131741`
- **Treated vs control (Resampling, betageom)**, `resampling_BGC`
  - BG14, treated: SRR11561505 (control); SRR11561504 (control); SRR11561503 (treated); SRR11561502 (treated)
  - `ytab treated-vs-control --project-config output/projects/Gale_2020/config/project.yaml --method resampling --normalization-run output/projects/Gale_2020/normalization/norm_20261022_all_methods --normalization betageom --output output/projects/Gale_2020/treated_vs_control/resampling_BGC`
- **Treated vs control (z-score, CPM)**, `Cg1_FLU_vs_no_drug_raw_cpm`
  - BG14, treated: SRR11561505 (control); SRR11561503 (treated); SRR11561504 (control); SRR11561502 (treated)
  - `ytab fitness --project-config output/projects/Gale_2020/config/project.yaml --analysis-id Cg1_FLU_vs_no_drug_raw_cpm --comparisons Cg1_FLU1_vs_control1,Cg1_FLU2_vs_control2`
- **Treated vs control (z-score, betageom)**, `treated_vs_parent_BGC_20Oct2026_2315`
  - BG14, treated: SRR11561505 (control); SRR11561503 (treated); SRR11561504 (control); SRR11561502 (treated)
  - `ytab fitness --project-config output/projects/Gale_2020/config/project.yaml --analysis-id treated_vs_parent_BGC_20Oct2026_2315 --comparisons Cg1_FLU1_vs_control1,Cg1_FLU2_vs_control2`
- **Treated vs control (z-score, betageom)**, `treated_vs_parent_BGC_22Oct2026_1048`
  - BG14, treated: SRR11561505 (control); SRR11561503 (treated); SRR11561504 (control); SRR11561502 (treated)
  - `ytab fitness --project-config output/projects/Gale_2020/config/project.yaml --analysis-id treated_vs_parent_BGC_22Oct2026_1048 --comparisons Cg1_FLU1_vs_control1,Cg1_FLU2_vs_control2`
- **Treated vs control (z-score, ttr)**, `treated_vs_parent_TTR_22Oct2026_1048`
  - BG14, treated: SRR11561505 (control); SRR11561503 (treated); SRR11561504 (control); SRR11561502 (treated)
  - `ytab fitness --project-config output/projects/Gale_2020/config/project.yaml --analysis-id treated_vs_parent_TTR_22Oct2026_1048 --comparisons Cg1_FLU1_vs_control1,Cg1_FLU2_vs_control2`
- **Multi-condition (betageom)**, `Cg1_parent_control_FLU_BGC`
  - parent, control, FLU: SRR11561508 (parent); SRR11561505 (control (reference)); SRR11561504 (control (reference)); SRR11561503 (FLU); SRR11561502 (FLU)
  - `ytab multicondition-project --project-config output/projects/Gale_2020/config/project.yaml --normalized None --normalization betageom --sample Gale_2020::SRR11561508=parent --sample Gale_2020::SRR11561505=control --sample Gale_2020::SRR11561504=control --sample Gale_2020::SRR11561503=FLU --sample Gale_2020::SRR11561502=FLU --reference-condition control --output output/projects/Gale_2020/multicondition/Cg1_parent_control_FLU_BGC`

## Segal et al. 2018: Candida albicans (Segal_2018_albicans)

YJB-T1081 (stable haploid of SC5314); Ac/Ds (Ds-NAT1), Ds primer tail. BioProject PRJNA490565 (SRR7824838–SRR7824848). 11 independent Ds libraries.

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| SRR7824842 |  | parent | parent | YJB-T1081 | Ca01 |
| SRR7824845 |  | parent | parent | YJB-T1081 | Ca02 |
| SRR7824843 |  | parent | parent | YJB-T1081 | Ca03 |
| SRR7824844 |  | parent | parent | YJB-T1081 | Ca04 |
| SRR7824839 |  | parent | parent | YJB-T1081 | Ca05 |
| SRR7824840 |  | parent | parent | YJB-T1081 | Ca06 |
| SRR7824841 |  | parent | parent | YJB-T1081 | Ca07 |
| SRR7824846 |  | parent | parent | YJB-T1081 | Ca08 |
| SRR7824847 |  | parent | parent | YJB-T1081 | Ca09 |
| SRR7824848 |  | parent | parent | YJB-T1081 | Ca10 |
| SRR7824838 |  | parent | parent | YJB-T1081 | Ca11 |

**Analyses**

- **Essentiality (HMM)**, `essentiality_20261025_131741`
  - libraries: SRR7824843; SRR7824841; SRR7824838
  - `ytab essentiality --project-config output/projects/Segal_2018_albicans/config/project.yaml --libraries SRR7824843 SRR7824841 SRR7824838 --output output/projects/Segal_2018_albicans/essentiality/essentiality_20261025_131741`

## Guo et al. 2013: Schizosaccharomyces pombe (Pombe_Hermes)

HL Hermes strain; Hermes, terminal 10 bp of the Hermes end. SRA SRR327340. 1 in vivo library (the library used by Segal et al. 2018).

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| SRR327340 |  | parent | parent | HL_Hermes | Guo2013 |

**Analyses:** parent libraries only (no condition analyses)


## Grech et al. 2019: Schizosaccharomyces pombe (Pombe_Hermes)

wild type; Hermes, end-ligated adapter (paired end). ENA PRJEB27324 (ERR2641199–ERR2641214). 8 high-density libraries, each sequenced twice.

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| ERR2641199 |  | parent | parent | Grech2019 | LG1_LG2 |
| ERR2641200 |  | parent | parent | Grech2019 | LG1_LG2 |
| ERR2641201 |  | parent | parent | Grech2019 | LG3_LG4 |
| ERR2641202 |  | parent | parent | Grech2019 | LG3_LG4 |
| ERR2641203 |  | parent | parent | Grech2019 | LG7_LG8 |
| ERR2641204 |  | parent | parent | Grech2019 | LG7_LG8 |
| ERR2641205 |  | parent | parent | Grech2019 | LG10_LG11 |
| ERR2641206 |  | parent | parent | Grech2019 | LG10_LG11 |
| ERR2641207 |  | parent | parent | Grech2019 | LG12_LG13 |
| ERR2641208 |  | parent | parent | Grech2019 | LG12_LG13 |
| ERR2641209 |  | parent | parent | Grech2019 | HL_B |
| ERR2641210 |  | parent | parent | Grech2019 | HL_B |
| ERR2641211 |  | parent | parent | Grech2019 | HL_D |
| ERR2641212 |  | parent | parent | Grech2019 | HL_D |
| ERR2641213 |  | parent | parent | Grech2019 | V2 |
| ERR2641214 |  | parent | parent | Grech2019 | V3 |

**Analyses**

- **Essentiality (HMM)**, `essentiality_20261025_131741`
  - libraries: ERR2641199; ERR2641200; ERR2641201; ERR2641202; ERR2641203; ERR2641204; ERR2641205; ERR2641206; ERR2641207; ERR2641208; ERR2641209; ERR2641210; ERR2641211; ERR2641212
  - `ytab essentiality --project-config output/projects/Pombe_Hermes/config/project.yaml --libraries ERR2641199 ERR2641200 ERR2641201 ERR2641202 ERR2641203 ERR2641204 ERR2641205 ERR2641206 ERR2641207 ERR2641208 ERR2641209 ERR2641210 ERR2641211 ERR2641212 --output output/projects/Pombe_Hermes/essentiality/essentiality_20261025_131741`

## Lee et al. 2020: Schizosaccharomyces pombe (Pombe_Hermes)

WT (YHL9451) and otr1R::ura4 (YHL10161); Hermes, end-ligated adapter. BioProject PRJNA517661 (SRR8502433–SRR8502440). 2 strains × 2 pools × generation 0 and 80 (FOA).

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| SRR8502438 |  | control | T0 | otr1R::ura4 | 10160B |
| SRR8502434 |  | treated | T80 | otr1R::ura4 | 10160B |
| SRR8502437 |  | control | T0 | otr1R::ura4 | 10161B |
| SRR8502433 |  | treated | T80 | otr1R::ura4 | 10161B |
| SRR8502440 |  | control | T0 | WT | 9451A |
| SRR8502436 |  | treated | T80 | WT | 9451A |
| SRR8502439 |  | control | T0 | WT | 9451C |
| SRR8502435 |  | treated | T80 | WT | 9451C |

**Analyses**

- **Treated vs control (Resampling, betageom)**, `T80_vs_T0_BGC`
  - otr1R::ura4, treated: SRR8502438 (control); SRR8502437 (control); SRR8502434 (treated); SRR8502433 (treated)
  - WT, treated: SRR8502440 (control); SRR8502439 (control); SRR8502436 (treated); SRR8502435 (treated)
  - `ytab treated-vs-control --project-config output/projects/Pombe_Hermes/config/project.yaml --method resampling --normalization-run output/projects/Pombe_Hermes/normalization/design_20261024 --normalization betageom --output output/projects/Pombe_Hermes/treated_vs_control/T80_vs_T0_BGC`
- **Treated vs control (z-score, CPM)**, `T80_vs_T0_WT`
  - WT, treated: SRR8502440 (control); SRR8502436 (treated); SRR8502439 (control); SRR8502435 (treated)
  - `ytab fitness --project-config output/projects/Pombe_Hermes/config/project.yaml --analysis-id T80_vs_T0_WT --comparisons WT_pool9451A_treated_vs_parent,WT_pool9451C_treated_vs_parent`
- **Treated vs control (z-score, CPM)**, `T80_vs_T0_otr1R`
  - otr1R::ura4, treated: SRR8502438 (control); SRR8502434 (treated); SRR8502437 (control); SRR8502433 (treated)
  - `ytab fitness --project-config output/projects/Pombe_Hermes/config/project.yaml --analysis-id T80_vs_T0_otr1R --comparisons otr1R_ura4_pool10160B_treated_vs_parent,otr1R_ura4_pool10161B_treated_vs_parent`
- **Multi-condition (betageom)**, `otr1R_vs_WT_T80_BGC`
  - WT, otr1R: SRR8502436 (WT (reference)); SRR8502435 (WT (reference)); SRR8502434 (otr1R); SRR8502433 (otr1R)
  - `ytab multicondition-project --project-config output/projects/Pombe_Hermes/config/project.yaml --normalized output/projects/Pombe_Hermes/normalization/between_strain_T80/experiment/methods --normalization betageom --sample Pombe_Hermes::SRR8502436=WT --sample Pombe_Hermes::SRR8502435=WT --sample Pombe_Hermes::SRR8502434=otr1R --sample Pombe_Hermes::SRR8502433=otr1R --reference-condition WT --output output/projects/Pombe_Hermes/multicondition/otr1R_vs_WT_T80_BGC`
- **Genetic interactions**, `T80_vs_T0_otr1R_vs_WT_BGC`
  - otr1R vs WT, strain x treatment: SRR8502440 (WT control); SRR8502439 (WT control); SRR8502438 (otr1R control); SRR8502437 (otr1R control); SRR8502436 (WT treated); SRR8502435 (WT treated); SRR8502434 (otr1R treated); SRR8502433 (otr1R treated)
  - `ytab genetic-interactions --normalized output/projects/Pombe_Hermes/normalization/conditions_v1/experiment/methods --normalization betageom --gff resources/species/pombe/reference_genome/schizosaccharomyces_pombe.chr.gff3 --reference resources/species/pombe/reference_genome/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa --reference-strain WT --query-strain otr1R --sample Pombe_Hermes::SRR8502440=WT:control --sample Pombe_Hermes::SRR8502439=WT:control --sample Pombe_Hermes::SRR8502438=otr1R:control --sample Pombe_Hermes::SRR8502437=otr1R:control --sample Pombe_Hermes::SRR8502436=WT:treated --sample Pombe_Hermes::SRR8502435=WT:treated --sample Pombe_Hermes::SRR8502434=otr1R:treated --sample Pombe_Hermes::SRR8502433=otr1R:treated --output output/projects/Pombe_Hermes/genetic_interactions/T80_vs_T0_otr1R_vs_WT_BGC`
- **Genetic interactions**, `T80_vs_T0_otr1R_vs_WT_TTR`
  - otr1R vs WT, strain x treatment: SRR8502440 (WT control); SRR8502439 (WT control); SRR8502438 (otr1R control); SRR8502437 (otr1R control); SRR8502436 (WT treated); SRR8502435 (WT treated); SRR8502434 (otr1R treated); SRR8502433 (otr1R treated)
  - `ytab genetic-interactions --normalized output/projects/Pombe_Hermes/normalization/conditions_v1/experiment/methods --normalization ttr --gff resources/species/pombe/reference_genome/schizosaccharomyces_pombe.chr.gff3 --reference resources/species/pombe/reference_genome/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa --reference-strain WT --query-strain otr1R --sample Pombe_Hermes::SRR8502440=WT:control --sample Pombe_Hermes::SRR8502439=WT:control --sample Pombe_Hermes::SRR8502438=otr1R:control --sample Pombe_Hermes::SRR8502437=otr1R:control --sample Pombe_Hermes::SRR8502436=WT:treated --sample Pombe_Hermes::SRR8502435=WT:treated --sample Pombe_Hermes::SRR8502434=otr1R:treated --sample Pombe_Hermes::SRR8502433=otr1R:treated --output output/projects/Pombe_Hermes/genetic_interactions/T80_vs_T0_otr1R_vs_WT_TTR`
- **Genetic interactions**, `T80_vs_T0_otr1R_vs_WT_total_reads`
  - otr1R vs WT, strain x treatment: SRR8502440 (WT control); SRR8502439 (WT control); SRR8502438 (otr1R control); SRR8502437 (otr1R control); SRR8502436 (WT treated); SRR8502435 (WT treated); SRR8502434 (otr1R treated); SRR8502433 (otr1R treated)
  - `ytab genetic-interactions --normalized output/projects/Pombe_Hermes/normalization/conditions_v1/experiment/methods --normalization totreads --gff resources/species/pombe/reference_genome/schizosaccharomyces_pombe.chr.gff3 --reference resources/species/pombe/reference_genome/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa --reference-strain WT --query-strain otr1R --sample Pombe_Hermes::SRR8502440=WT:control --sample Pombe_Hermes::SRR8502439=WT:control --sample Pombe_Hermes::SRR8502438=otr1R:control --sample Pombe_Hermes::SRR8502437=otr1R:control --sample Pombe_Hermes::SRR8502436=WT:treated --sample Pombe_Hermes::SRR8502435=WT:treated --sample Pombe_Hermes::SRR8502434=otr1R:treated --sample Pombe_Hermes::SRR8502433=otr1R:treated --output output/projects/Pombe_Hermes/genetic_interactions/T80_vs_T0_otr1R_vs_WT_total_reads`

## Michel et al. 2017: Saccharomyces cerevisiae (Michel_2017_SATAY)

BY4741 ade2Δ background; Ac/Ds (MiniDs), restriction digestion and circularization. ArrayExpress E-MTAB-4885; ENA PRJEB14804 (ERR1533143–ERR1533149). 2 wild-type libraries; rapamycin (grown from wild type 2); dpl1Δ, dpl1Δ psd2Δ, VPS13(D716H), mmm1Δ VPS13(D716H).

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| ERR1533147 |  | parent | parent | BY4741_ade2 | WT1 |
| ERR1533148 |  | parent | parent | BY4741_ade2 | WT2 |
| ERR1533149 |  | treated | rapamycin | BY4741_ade2 | WT2 |
| ERR1533143 |  | treated | dpl1del | BY4741_ade2 | dpl1del |
| ERR1533144 |  | treated | dpl1del_psd2del | BY4741_ade2 | dpl1del_psd2del |
| ERR1533146 |  | treated | VPS13_D716H | BY4741_ade2 | VPS13_D716H |
| ERR1533145 |  | treated | mmm1del_VPS13_D716H | BY4741_ade2 | mmm1del_VPS13_D716H |

**Analyses**

- **Essentiality (HMM)**, `essentiality_20261025_131741`
  - libraries: ERR1533147; ERR1533148
  - `ytab essentiality --project-config output/projects/Michel_2017_SATAY/config/project.yaml --libraries ERR1533147 ERR1533148 --output output/projects/Michel_2017_SATAY/essentiality/essentiality_20261025_131741`
- **Treated vs control (Resampling, betageom)**, `rapamycin_vs_WT2_BGC`
  - BY4741_ade2, treated: ERR1533148 (control); ERR1533149 (treated)
  - `ytab treated-vs-control --project-config output/projects/Michel_2017_SATAY/config/project.yaml --method resampling --normalization-run output/projects/Michel_2017_SATAY/normalization/conditions_v1 --normalization betageom --output output/projects/Michel_2017_SATAY/treated_vs_control/rapamycin_vs_WT2_BGC`
- **Treated vs control (Resampling, ttr)**, `rapamycin_vs_WT2_TTR`
  - BY4741_ade2, treated: ERR1533148 (control); ERR1533149 (treated)
  - `ytab treated-vs-control --project-config output/projects/Michel_2017_SATAY/config/project.yaml --method resampling --normalization-run output/projects/Michel_2017_SATAY/normalization/conditions_v1 --normalization ttr --output output/projects/Michel_2017_SATAY/treated_vs_control/rapamycin_vs_WT2_TTR`
- **Treated vs control (Resampling, totreads)**, `rapamycin_vs_WT2_total_reads`
  - BY4741_ade2, treated: ERR1533148 (control); ERR1533149 (treated)
  - `ytab treated-vs-control --project-config output/projects/Michel_2017_SATAY/config/project.yaml --method resampling --normalization-run output/projects/Michel_2017_SATAY/normalization/conditions_v1 --normalization totreads --output output/projects/Michel_2017_SATAY/treated_vs_control/rapamycin_vs_WT2_total_reads`
- **Multi-condition (betageom)**, `mutants_vs_WT_BGC`
  - WT, dpl1del, dpl1del_psd2del, VPS13_D716H, mmm1del_VPS13_D716H: ERR1533147 (WT (reference)); ERR1533148 (WT (reference)); ERR1533143 (dpl1del); ERR1533144 (dpl1del_psd2del); ERR1533146 (VPS13_D716H); ERR1533145 (mmm1del_VPS13_D716H)
  - `ytab multicondition-project --project-config output/projects/Michel_2017_SATAY/config/project.yaml --normalized output/projects/Michel_2017_SATAY/normalization/conditions_v1/experiment/methods --normalization betageom --sample Michel_2017_SATAY::ERR1533147=WT --sample Michel_2017_SATAY::ERR1533148=WT --sample Michel_2017_SATAY::ERR1533143=dpl1del --sample Michel_2017_SATAY::ERR1533144=dpl1del_psd2del --sample Michel_2017_SATAY::ERR1533146=VPS13_D716H --sample Michel_2017_SATAY::ERR1533145=mmm1del_VPS13_D716H --reference-condition WT --output output/projects/Michel_2017_SATAY/multicondition/mutants_vs_WT_BGC`
- **Multi-condition (ttr)**, `mutants_vs_WT_TTR`
  - WT, dpl1del, dpl1del_psd2del, VPS13_D716H, mmm1del_VPS13_D716H: ERR1533147 (WT (reference)); ERR1533148 (WT (reference)); ERR1533143 (dpl1del); ERR1533144 (dpl1del_psd2del); ERR1533146 (VPS13_D716H); ERR1533145 (mmm1del_VPS13_D716H)
  - `ytab multicondition-project --project-config output/projects/Michel_2017_SATAY/config/project.yaml --normalized output/projects/Michel_2017_SATAY/normalization/conditions_v1/experiment/methods --normalization ttr --sample Michel_2017_SATAY::ERR1533147=WT --sample Michel_2017_SATAY::ERR1533148=WT --sample Michel_2017_SATAY::ERR1533143=dpl1del --sample Michel_2017_SATAY::ERR1533144=dpl1del_psd2del --sample Michel_2017_SATAY::ERR1533146=VPS13_D716H --sample Michel_2017_SATAY::ERR1533145=mmm1del_VPS13_D716H --reference-condition WT --output output/projects/Michel_2017_SATAY/multicondition/mutants_vs_WT_TTR`
- **Multi-condition (totreads)**, `mutants_vs_WT_total_reads`
  - WT, dpl1del, dpl1del_psd2del, VPS13_D716H, mmm1del_VPS13_D716H: ERR1533147 (WT (reference)); ERR1533148 (WT (reference)); ERR1533143 (dpl1del); ERR1533144 (dpl1del_psd2del); ERR1533146 (VPS13_D716H); ERR1533145 (mmm1del_VPS13_D716H)
  - `ytab multicondition-project --project-config output/projects/Michel_2017_SATAY/config/project.yaml --normalized output/projects/Michel_2017_SATAY/normalization/conditions_v1/experiment/methods --normalization totreads --sample Michel_2017_SATAY::ERR1533147=WT --sample Michel_2017_SATAY::ERR1533148=WT --sample Michel_2017_SATAY::ERR1533143=dpl1del --sample Michel_2017_SATAY::ERR1533144=dpl1del_psd2del --sample Michel_2017_SATAY::ERR1533146=VPS13_D716H --sample Michel_2017_SATAY::ERR1533145=mmm1del_VPS13_D716H --reference-condition WT --output output/projects/Michel_2017_SATAY/multicondition/mutants_vs_WT_total_reads`

## Nickels et al. 2024: Candida glabrata (Nickels_2024_glabrata)

CBS138 ura3 derivatives (CBS138u, 2001u); Hermes (NAT), splinkerette PCR. BioProject PRJNA1133102. 4 parent pools; 8 micafungin cultures (2–32 ng/µL).

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| SRR29743993 |  | parent | parent | CBS138 | CBS138_Ind1 |
| SRR29743992 |  | parent | parent | CBS138 | CBS138_Ind2 |
| SRR29743986 |  | treated | micafungin_16ng | CBS138 | CBS138_Ind1 |
| SRR29743989 |  | treated | micafungin_2ng | CBS138 | CBS138_Ind1 |
| SRR29743988 |  | treated | micafungin_4ng | CBS138 | CBS138_Ind1 |
| SRR29743987 |  | treated | micafungin_8ng | CBS138 | CBS138_Ind1 |
| SRR29743985 |  | parent | parent | CBS138_2001 | CBS138_2001_Ind1 |
| SRR29743984 |  | parent | parent | CBS138_2001 | CBS138_2001_Ind2 |
| SRR29743991 |  | treated | micafungin_16ng | CBS138_2001 | CBS138_2001_Ind1 |
| SRR29743990 |  | treated | micafungin_32ng | CBS138_2001 | CBS138_2001_Ind1 |
| SRR29743983 |  | treated | micafungin_4ng | CBS138_2001 | CBS138_2001_Ind1 |
| SRR29743982 |  | treated | micafungin_8ng | CBS138_2001 | CBS138_2001_Ind1 |

**Analyses**

Rerun every condition analysis: `ytab run --workflow output/projects/Nickels_2024_glabrata/config/analyses_workflow.yaml`

- **Essentiality (HMM)**, `essentiality_CBS138_2001`
  - libraries: SRR29743984; SRR29743985
  - `ytab essentiality --project-config output/projects/Nickels_2024_glabrata/config/project.yaml --libraries SRR29743984 SRR29743985 SRR33001998 --output output/projects/Nickels_2024_glabrata/essentiality/essentiality_CBS138_2001`
- **Essentiality (HMM)**, `essentiality_CBS138_Nickels2024_CBS138`
  - libraries: SRR29743992; SRR29743993
  - `ytab essentiality --project-config output/projects/Nickels_2024_glabrata/config/project.yaml --libraries SRR29743992 SRR29743993 --output output/projects/Nickels_2024_glabrata/essentiality/essentiality_CBS138_Nickels2024_CBS138`
- **Treated vs control (Resampling, betageom)**, `treated_vs_control_BGC`
  - CBS138, micafungin_16ng: SRR29743993 (control); SRR29743986 (treated)
  - CBS138, micafungin_2ng: SRR29743993 (control); SRR29743989 (treated)
  - CBS138, micafungin_4ng: SRR29743993 (control); SRR29743988 (treated)
  - CBS138, micafungin_8ng: SRR29743993 (control); SRR29743987 (treated)
  - CBS138_2001, micafungin_16ng: SRR29743985 (control); SRR29743991 (treated)
  - CBS138_2001, micafungin_32ng: SRR29743985 (control); SRR29743990 (treated)
  - CBS138_2001, micafungin_4ng: SRR29743985 (control); SRR29743983 (treated)
  - CBS138_2001, micafungin_8ng: SRR29743985 (control); SRR29743982 (treated)
  - `ytab treated-vs-control --project-config output/projects/Nickels_2024_glabrata/config/project.yaml --method resampling --normalization-run output/projects/Nickels_2024_glabrata/normalization/conditions_v1 --normalization betageom --output output/projects/Nickels_2024_glabrata/treated_vs_control/treated_vs_control_BGC`
- **Multi-condition (betageom)**, `CBS138_2001_micafungin_series_BGC`
  - parent, micafungin_16ng, micafungin_32ng, micafungin_4ng, micafungin_8ng: SRR29743985 (parent (reference)); SRR29743991 (micafungin_16ng); SRR29743990 (micafungin_32ng); SRR29743983 (micafungin_4ng); SRR29743982 (micafungin_8ng)
  - `ytab multicondition-project --project-config output/projects/Nickels_2024_glabrata/config/project.yaml --normalized output/projects/Nickels_2024_glabrata/normalization/conditions_v1/experiment/methods --normalization betageom --sample Nickels_2024_glabrata::SRR29743985=parent --sample Nickels_2024_glabrata::SRR29743991=micafungin_16ng --sample Nickels_2024_glabrata::SRR29743990=micafungin_32ng --sample Nickels_2024_glabrata::SRR29743983=micafungin_4ng --sample Nickels_2024_glabrata::SRR29743982=micafungin_8ng --reference-condition parent --output output/projects/Nickels_2024_glabrata/multicondition/CBS138_2001_micafungin_series_BGC`
- **Multi-condition (betageom)**, `CBS138_micafungin_series_BGC`
  - parent, micafungin_2ng, micafungin_4ng, micafungin_8ng, micafungin_16ng: SRR29743993 (parent (reference)); SRR29743989 (micafungin_2ng); SRR29743988 (micafungin_4ng); SRR29743987 (micafungin_8ng); SRR29743986 (micafungin_16ng)
  - `ytab multicondition-project --project-config output/projects/Nickels_2024_glabrata/config/project.yaml --normalized output/projects/Nickels_2024_glabrata/normalization/conditions_v1/experiment/methods --normalization betageom --sample Nickels_2024_glabrata::SRR29743993=parent --sample Nickels_2024_glabrata::SRR29743989=micafungin_2ng --sample Nickels_2024_glabrata::SRR29743988=micafungin_4ng --sample Nickels_2024_glabrata::SRR29743987=micafungin_8ng --sample Nickels_2024_glabrata::SRR29743986=micafungin_16ng --reference-condition parent --output output/projects/Nickels_2024_glabrata/multicondition/CBS138_micafungin_series_BGC`

## Nickels & Cunningham 2025: Candida glabrata (Nickels_2025_glabrata)

BG2, BG2 PDR1 gain of function, BG2 pdr1Δ; CBS138 2001u; Hermes (NAT), splinkerette PCR. BioProject PRJNA1247003. 10 no-drug controls; 39 cultures in five echinocandins or ibrexafungerp (1–64 ng/µL).

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| SRR33001989 |  | control | anidulafungin_no_drug | BG2 | BG2_Ind3 |
| SRR33002022 |  | control | caspofungin_no_drug | BG2 | BG2_Ind3 |
| SRR33001993 |  | control | caspofungin_no_drug | BG2 | BG2_Ind3 |
| SRR33002027 |  | control | ibrexafungerp_no_drug | BG2 | BG2_Ind3 |
| SRR33002033 |  | control | micafungin_no_drug | BG2 | BG2_Ind3 |
| SRR33002031 |  | control | rezafungin_no_drug | BG2 | BG2_Ind3 |
| SRR33001985 |  | treated | anidulafungin_16ng | BG2 | BG2_Ind3 |
| SRR33001988 |  | treated | anidulafungin_2ng | BG2 | BG2_Ind3 |
| SRR33001987 |  | treated | anidulafungin_4ng | BG2 | BG2_Ind3 |
| SRR33001986 |  | treated | anidulafungin_8ng | BG2 | BG2_Ind3 |
| SRR33002016 |  | treated | caspofungin_16ng | BG2 | BG2_Ind3 |
| SRR33002020 |  | treated | caspofungin_1ng | BG2 | BG2_Ind3 |
| SRR33002019 |  | treated | caspofungin_2ng | BG2 | BG2_Ind3 |
| SRR33001992 |  | treated | caspofungin_2ng | BG2 | BG2_Ind3 |
| SRR33002018 |  | treated | caspofungin_4ng | BG2 | BG2_Ind3 |
| SRR33001991 |  | treated | caspofungin_4ng | BG2 | BG2_Ind3 |
| SRR33002017 |  | treated | caspofungin_8ng | BG2 | BG2_Ind3 |
| SRR33001990 |  | treated | caspofungin_8ng | BG2 | BG2_Ind3 |
| SRR33002025 |  | treated | ibrexafungerp_16ng | BG2 | BG2_Ind3 |
| SRR33002024 |  | treated | ibrexafungerp_32ng | BG2 | BG2_Ind3 |
| SRR33002023 |  | treated | ibrexafungerp_64ng | BG2 | BG2_Ind3 |
| SRR33002026 |  | treated | ibrexafungerp_8ng | BG2 | BG2_Ind3 |
| SRR33001999 |  | treated | micafungin_16ng | BG2 | BG2_Ind3 |
| SRR33002032 |  | treated | micafungin_2ng | BG2 | BG2_Ind3 |
| SRR33002021 |  | treated | micafungin_4ng | BG2 | BG2_Ind3 |
| SRR33002010 |  | treated | micafungin_8ng | BG2 | BG2_Ind3 |
| SRR33002030 |  | treated | rezafungin_2ng | BG2 | BG2_Ind3 |
| SRR33002029 |  | treated | rezafungin_4ng | BG2 | BG2_Ind3 |
| SRR33002028 |  | treated | rezafungin_8ng | BG2 | BG2_Ind3 |
| SRR33002004 |  | control | micafungin_no_drug | PDR1_GOF | PDR1_GOF_Ind1 |
| SRR33002000 |  | treated | micafungin_16ng | PDR1_GOF | PDR1_GOF_Ind1 |
| SRR33002003 |  | treated | micafungin_2ng | PDR1_GOF | PDR1_GOF_Ind1 |
| SRR33002002 |  | treated | micafungin_4ng | PDR1_GOF | PDR1_GOF_Ind1 |
| SRR33002001 |  | treated | micafungin_8ng | PDR1_GOF | PDR1_GOF_Ind1 |
| SRR33002009 |  | control | ibrexafungerp_no_drug | pdr1_KO | pdr1_KO_Ind1 |
| SRR33002015 |  | control | micafungin_no_drug | pdr1_KO | pdr1_KO_Ind1 |
| SRR33002007 |  | treated | ibrexafungerp_16ng | pdr1_KO | pdr1_KO_Ind1 |
| SRR33002006 |  | treated | ibrexafungerp_32ng | pdr1_KO | pdr1_KO_Ind1 |
| SRR33002005 |  | treated | ibrexafungerp_64ng | pdr1_KO | pdr1_KO_Ind1 |
| SRR33002008 |  | treated | ibrexafungerp_8ng | pdr1_KO | pdr1_KO_Ind1 |
| SRR33002011 |  | treated | micafungin_16ng | pdr1_KO | pdr1_KO_Ind1 |
| SRR33002014 |  | treated | micafungin_2ng | pdr1_KO | pdr1_KO_Ind1 |
| SRR33002013 |  | treated | micafungin_4ng | pdr1_KO | pdr1_KO_Ind1 |
| SRR33002012 |  | treated | micafungin_8ng | pdr1_KO | pdr1_KO_Ind1 |

**Analyses**

Rerun every condition analysis: `ytab run --workflow output/projects/Nickels_2025_glabrata/config/analyses_workflow.yaml`

- **Essentiality (HMM)**, `essentiality_H1_new_BG2`
  - libraries: SRR33001989; SRR33001993; SRR33002022; SRR33002027; SRR33002031; SRR33002033
  - `ytab essentiality --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --libraries SRR33001989 SRR33001993 SRR33002022 SRR33002027 SRR33002031 SRR33002033 SRR38338548 SRR38338549 --output output/projects/Nickels_2025_glabrata/essentiality/essentiality_H1_new_BG2`
- **Essentiality (HMM)**, `essentiality_glabrata_PDR1_GOF`
  - libraries: SRR33002004
  - `ytab essentiality --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --libraries SRR33002004 --output output/projects/Nickels_2025_glabrata/essentiality/essentiality_glabrata_PDR1_GOF`
- **Essentiality (HMM)**, `essentiality_glabrata_pdr1_KO`
  - libraries: SRR33002009; SRR33002015
  - `ytab essentiality --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --libraries SRR33002009 SRR33002015 --output output/projects/Nickels_2025_glabrata/essentiality/essentiality_glabrata_pdr1_KO`
- **Treated vs control (Resampling, betageom)**, `treated_vs_control_BGC`
  - BG2, anidulafungin_16ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33001985 (treated)
  - BG2, anidulafungin_2ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33001988 (treated)
  - BG2, anidulafungin_4ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33001987 (treated)
  - BG2, anidulafungin_8ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33001986 (treated)
  - BG2, caspofungin_16ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002016 (treated)
  - BG2, caspofungin_1ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002020 (treated)
  - BG2, caspofungin_2ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002019 (treated); SRR33001992 (treated)
  - BG2, caspofungin_4ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002018 (treated); SRR33001991 (treated)
  - BG2, caspofungin_8ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002017 (treated); SRR33001990 (treated)
  - BG2, ibrexafungerp_16ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002025 (treated)
  - BG2, ibrexafungerp_32ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002024 (treated)
  - BG2, ibrexafungerp_64ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002023 (treated)
  - BG2, ibrexafungerp_8ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002026 (treated)
  - BG2, micafungin_16ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33001999 (treated)
  - BG2, micafungin_2ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002032 (treated)
  - BG2, micafungin_4ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002021 (treated)
  - BG2, micafungin_8ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002010 (treated)
  - BG2, rezafungin_2ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002030 (treated)
  - BG2, rezafungin_4ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002029 (treated)
  - BG2, rezafungin_8ng: SRR33001989 (control); SRR33002022 (control); SRR33001993 (control); SRR33002027 (control); SRR33002033 (control); SRR33002031 (control); SRR33002028 (treated)
  - PDR1_GOF, micafungin_16ng: SRR33002004 (control); SRR33002000 (treated)
  - PDR1_GOF, micafungin_2ng: SRR33002004 (control); SRR33002003 (treated)
  - PDR1_GOF, micafungin_4ng: SRR33002004 (control); SRR33002002 (treated)
  - PDR1_GOF, micafungin_8ng: SRR33002004 (control); SRR33002001 (treated)
  - pdr1_KO, ibrexafungerp_16ng: SRR33002009 (control); SRR33002015 (control); SRR33002007 (treated)
  - pdr1_KO, ibrexafungerp_32ng: SRR33002009 (control); SRR33002015 (control); SRR33002006 (treated)
  - pdr1_KO, ibrexafungerp_64ng: SRR33002009 (control); SRR33002015 (control); SRR33002005 (treated)
  - pdr1_KO, ibrexafungerp_8ng: SRR33002009 (control); SRR33002015 (control); SRR33002008 (treated)
  - pdr1_KO, micafungin_16ng: SRR33002009 (control); SRR33002015 (control); SRR33002011 (treated)
  - pdr1_KO, micafungin_2ng: SRR33002009 (control); SRR33002015 (control); SRR33002014 (treated)
  - pdr1_KO, micafungin_4ng: SRR33002009 (control); SRR33002015 (control); SRR33002013 (treated)
  - pdr1_KO, micafungin_8ng: SRR33002009 (control); SRR33002015 (control); SRR33002012 (treated)
  - `ytab treated-vs-control --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --method resampling --normalization-run output/projects/Nickels_2025_glabrata/normalization/conditions_v1 --normalization betageom --output output/projects/Nickels_2025_glabrata/treated_vs_control/treated_vs_control_BGC`
- **Multi-condition (betageom)**, `BG2_anidulafungin_series_BGC`
  - anidulafungin_no_drug, anidulafungin_16ng, anidulafungin_2ng, anidulafungin_4ng, anidulafungin_8ng: SRR33001989 (anidulafungin_no_drug (reference)); SRR33001985 (anidulafungin_16ng); SRR33001988 (anidulafungin_2ng); SRR33001987 (anidulafungin_4ng); SRR33001986 (anidulafungin_8ng)
  - `ytab multicondition-project --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --normalized output/projects/Nickels_2025_glabrata/normalization/conditions_v1/experiment/methods --normalization betageom --sample Nickels_2025_glabrata::SRR33001989=anidulafungin_no_drug --sample Nickels_2025_glabrata::SRR33001985=anidulafungin_16ng --sample Nickels_2025_glabrata::SRR33001988=anidulafungin_2ng --sample Nickels_2025_glabrata::SRR33001987=anidulafungin_4ng --sample Nickels_2025_glabrata::SRR33001986=anidulafungin_8ng --reference-condition anidulafungin_no_drug --output output/projects/Nickels_2025_glabrata/multicondition/BG2_anidulafungin_series_BGC`
- **Multi-condition (betageom)**, `BG2_caspofungin_series_BGC`
  - caspofungin_no_drug, caspofungin_16ng, caspofungin_1ng, caspofungin_2ng, caspofungin_4ng, caspofungin_8ng: SRR33002022 (caspofungin_no_drug (reference)); SRR33001993 (caspofungin_no_drug (reference)); SRR33002016 (caspofungin_16ng); SRR33002020 (caspofungin_1ng); SRR33002019 (caspofungin_2ng); SRR33001992 (caspofungin_2ng); SRR33002018 (caspofungin_4ng); SRR33001991 (caspofungin_4ng); SRR33002017 (caspofungin_8ng); SRR33001990 (caspofungin_8ng)
  - `ytab multicondition-project --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --normalized output/projects/Nickels_2025_glabrata/normalization/conditions_v1/experiment/methods --normalization betageom --sample Nickels_2025_glabrata::SRR33002022=caspofungin_no_drug --sample Nickels_2025_glabrata::SRR33001993=caspofungin_no_drug --sample Nickels_2025_glabrata::SRR33002016=caspofungin_16ng --sample Nickels_2025_glabrata::SRR33002020=caspofungin_1ng --sample Nickels_2025_glabrata::SRR33002019=caspofungin_2ng --sample Nickels_2025_glabrata::SRR33001992=caspofungin_2ng --sample Nickels_2025_glabrata::SRR33002018=caspofungin_4ng --sample Nickels_2025_glabrata::SRR33001991=caspofungin_4ng --sample Nickels_2025_glabrata::SRR33002017=caspofungin_8ng --sample Nickels_2025_glabrata::SRR33001990=caspofungin_8ng --reference-condition caspofungin_no_drug --output output/projects/Nickels_2025_glabrata/multicondition/BG2_caspofungin_series_BGC`
- **Multi-condition (betageom)**, `BG2_ibrexafungerp_series_BGC`
  - ibrexafungerp_no_drug, ibrexafungerp_16ng, ibrexafungerp_32ng, ibrexafungerp_64ng, ibrexafungerp_8ng: SRR33002027 (ibrexafungerp_no_drug (reference)); SRR33002025 (ibrexafungerp_16ng); SRR33002024 (ibrexafungerp_32ng); SRR33002023 (ibrexafungerp_64ng); SRR33002026 (ibrexafungerp_8ng)
  - `ytab multicondition-project --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --normalized output/projects/Nickels_2025_glabrata/normalization/conditions_v1/experiment/methods --normalization betageom --sample Nickels_2025_glabrata::SRR33002027=ibrexafungerp_no_drug --sample Nickels_2025_glabrata::SRR33002025=ibrexafungerp_16ng --sample Nickels_2025_glabrata::SRR33002024=ibrexafungerp_32ng --sample Nickels_2025_glabrata::SRR33002023=ibrexafungerp_64ng --sample Nickels_2025_glabrata::SRR33002026=ibrexafungerp_8ng --reference-condition ibrexafungerp_no_drug --output output/projects/Nickels_2025_glabrata/multicondition/BG2_ibrexafungerp_series_BGC`
- **Multi-condition (betageom)**, `BG2_micafungin_series_BGC`
  - micafungin_no_drug, micafungin_16ng, micafungin_2ng, micafungin_4ng, micafungin_8ng: SRR33002033 (micafungin_no_drug (reference)); SRR33001999 (micafungin_16ng); SRR33002032 (micafungin_2ng); SRR33002021 (micafungin_4ng); SRR33002010 (micafungin_8ng)
  - `ytab multicondition-project --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --normalized output/projects/Nickels_2025_glabrata/normalization/conditions_v1/experiment/methods --normalization betageom --sample Nickels_2025_glabrata::SRR33002033=micafungin_no_drug --sample Nickels_2025_glabrata::SRR33001999=micafungin_16ng --sample Nickels_2025_glabrata::SRR33002032=micafungin_2ng --sample Nickels_2025_glabrata::SRR33002021=micafungin_4ng --sample Nickels_2025_glabrata::SRR33002010=micafungin_8ng --reference-condition micafungin_no_drug --output output/projects/Nickels_2025_glabrata/multicondition/BG2_micafungin_series_BGC`
- **Multi-condition (betageom)**, `BG2_rezafungin_series_BGC`
  - rezafungin_no_drug, rezafungin_2ng, rezafungin_4ng, rezafungin_8ng: SRR33002031 (rezafungin_no_drug (reference)); SRR33002030 (rezafungin_2ng); SRR33002029 (rezafungin_4ng); SRR33002028 (rezafungin_8ng)
  - `ytab multicondition-project --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --normalized output/projects/Nickels_2025_glabrata/normalization/conditions_v1/experiment/methods --normalization betageom --sample Nickels_2025_glabrata::SRR33002031=rezafungin_no_drug --sample Nickels_2025_glabrata::SRR33002030=rezafungin_2ng --sample Nickels_2025_glabrata::SRR33002029=rezafungin_4ng --sample Nickels_2025_glabrata::SRR33002028=rezafungin_8ng --reference-condition rezafungin_no_drug --output output/projects/Nickels_2025_glabrata/multicondition/BG2_rezafungin_series_BGC`
- **Multi-condition (betageom)**, `PDR1_GOF_micafungin_series_BGC`
  - micafungin_no_drug, micafungin_16ng, micafungin_2ng, micafungin_4ng, micafungin_8ng: SRR33002004 (micafungin_no_drug (reference)); SRR33002000 (micafungin_16ng); SRR33002003 (micafungin_2ng); SRR33002002 (micafungin_4ng); SRR33002001 (micafungin_8ng)
  - `ytab multicondition-project --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --normalized output/projects/Nickels_2025_glabrata/normalization/conditions_v1/experiment/methods --normalization betageom --sample Nickels_2025_glabrata::SRR33002004=micafungin_no_drug --sample Nickels_2025_glabrata::SRR33002000=micafungin_16ng --sample Nickels_2025_glabrata::SRR33002003=micafungin_2ng --sample Nickels_2025_glabrata::SRR33002002=micafungin_4ng --sample Nickels_2025_glabrata::SRR33002001=micafungin_8ng --reference-condition micafungin_no_drug --output output/projects/Nickels_2025_glabrata/multicondition/PDR1_GOF_micafungin_series_BGC`
- **Multi-condition (betageom)**, `pdr1_KO_ibrexafungerp_series_BGC`
  - ibrexafungerp_no_drug, ibrexafungerp_16ng, ibrexafungerp_32ng, ibrexafungerp_64ng, ibrexafungerp_8ng: SRR33002009 (ibrexafungerp_no_drug (reference)); SRR33002007 (ibrexafungerp_16ng); SRR33002006 (ibrexafungerp_32ng); SRR33002005 (ibrexafungerp_64ng); SRR33002008 (ibrexafungerp_8ng)
  - `ytab multicondition-project --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --normalized output/projects/Nickels_2025_glabrata/normalization/conditions_v1/experiment/methods --normalization betageom --sample Nickels_2025_glabrata::SRR33002009=ibrexafungerp_no_drug --sample Nickels_2025_glabrata::SRR33002007=ibrexafungerp_16ng --sample Nickels_2025_glabrata::SRR33002006=ibrexafungerp_32ng --sample Nickels_2025_glabrata::SRR33002005=ibrexafungerp_64ng --sample Nickels_2025_glabrata::SRR33002008=ibrexafungerp_8ng --reference-condition ibrexafungerp_no_drug --output output/projects/Nickels_2025_glabrata/multicondition/pdr1_KO_ibrexafungerp_series_BGC`
- **Multi-condition (betageom)**, `pdr1_KO_micafungin_series_BGC`
  - micafungin_no_drug, micafungin_16ng, micafungin_2ng, micafungin_4ng, micafungin_8ng: SRR33002015 (micafungin_no_drug (reference)); SRR33002011 (micafungin_16ng); SRR33002014 (micafungin_2ng); SRR33002013 (micafungin_4ng); SRR33002012 (micafungin_8ng)
  - `ytab multicondition-project --project-config output/projects/Nickels_2025_glabrata/config/project.yaml --normalized output/projects/Nickels_2025_glabrata/normalization/conditions_v1/experiment/methods --normalization betageom --sample Nickels_2025_glabrata::SRR33002015=micafungin_no_drug --sample Nickels_2025_glabrata::SRR33002011=micafungin_16ng --sample Nickels_2025_glabrata::SRR33002014=micafungin_2ng --sample Nickels_2025_glabrata::SRR33002013=micafungin_4ng --sample Nickels_2025_glabrata::SRR33002012=micafungin_8ng --reference-condition micafungin_no_drug --output output/projects/Nickels_2025_glabrata/multicondition/pdr1_KO_micafungin_series_BGC`

## Nickels & Cunningham 2025: Candida glabrata (Nickels_2025_glabrata_CBS138)

BG2, BG2 PDR1 gain of function, BG2 pdr1Δ; CBS138 2001u; Hermes (NAT), splinkerette PCR. BioProject PRJNA1247003. 10 no-drug controls; 39 cultures in five echinocandins or ibrexafungerp (1–64 ng/µL).

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| SRR33001998 |  | control | caspofungin_no_drug | CBS138_2001 | CBS138_2001_Ind1 |
| SRR33001994 |  | treated | caspofungin_16ng | CBS138_2001 | CBS138_2001_Ind1 |
| SRR33001997 |  | treated | caspofungin_2ng | CBS138_2001 | CBS138_2001_Ind1 |
| SRR33001996 |  | treated | caspofungin_4ng | CBS138_2001 | CBS138_2001_Ind1 |
| SRR33001995 |  | treated | caspofungin_8ng | CBS138_2001 | CBS138_2001_Ind1 |

**Analyses**

Rerun every condition analysis: `ytab run --workflow output/projects/Nickels_2025_glabrata_CBS138/config/analyses_workflow.yaml`

- **Essentiality (HMM)**, `essentiality_CBS138_2001`
  - libraries: SRR33001998
  - `ytab essentiality --project-config output/projects/Nickels_2025_glabrata_CBS138/config/project.yaml --libraries SRR29743984 SRR29743985 SRR33001998 --output output/projects/Nickels_2025_glabrata_CBS138/essentiality/essentiality_CBS138_2001`
- **Treated vs control (Resampling, betageom)**, `treated_vs_control_BGC`
  - CBS138_2001, caspofungin_16ng: SRR33001998 (control); SRR33001994 (treated)
  - CBS138_2001, caspofungin_2ng: SRR33001998 (control); SRR33001997 (treated)
  - CBS138_2001, caspofungin_4ng: SRR33001998 (control); SRR33001996 (treated)
  - CBS138_2001, caspofungin_8ng: SRR33001998 (control); SRR33001995 (treated)
  - `ytab treated-vs-control --project-config output/projects/Nickels_2025_glabrata_CBS138/config/project.yaml --method resampling --normalization-run output/projects/Nickels_2025_glabrata_CBS138/normalization/conditions_v1 --normalization betageom --output output/projects/Nickels_2025_glabrata_CBS138/treated_vs_control/treated_vs_control_BGC`
- **Multi-condition (betageom)**, `CBS138_2001_caspofungin_series_BGC`
  - caspofungin_no_drug, caspofungin_16ng, caspofungin_2ng, caspofungin_4ng, caspofungin_8ng: SRR33001998 (caspofungin_no_drug (reference)); SRR33001994 (caspofungin_16ng); SRR33001997 (caspofungin_2ng); SRR33001996 (caspofungin_4ng); SRR33001995 (caspofungin_8ng)
  - `ytab multicondition-project --project-config output/projects/Nickels_2025_glabrata_CBS138/config/project.yaml --normalized output/projects/Nickels_2025_glabrata_CBS138/normalization/conditions_v1/experiment/methods --normalization betageom --sample Nickels_2025_glabrata_CBS138::SRR33001998=caspofungin_no_drug --sample Nickels_2025_glabrata_CBS138::SRR33001994=caspofungin_16ng --sample Nickels_2025_glabrata_CBS138::SRR33001997=caspofungin_2ng --sample Nickels_2025_glabrata_CBS138::SRR33001996=caspofungin_4ng --sample Nickels_2025_glabrata_CBS138::SRR33001995=caspofungin_8ng --reference-condition caspofungin_no_drug --output output/projects/Nickels_2025_glabrata_CBS138/multicondition/CBS138_2001_caspofungin_series_BGC`

## Zaeske et al. 2026: Candida glabrata (Zaeske_2026_glabrata)

BG14, BG2r, DPL1021r, CBS138r; Hermes (NAT), splinkerette PCR. BioProject PRJNA1460454. 5 parent pools (one deposit duplicates a Nickels & Cunningham run and is counted once).

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| SRR38338547 (not analysed) |  | parent | parent | BG14 | BG14_pool1 |
| SRR38338549 |  | parent | parent | BG2r | BG2r_pool1 |
| SRR38338548 |  | parent | parent | BG2r | BG2r_pool2 |
| SRR38338545 |  | parent | parent | DPL1021r | DPL1021r_pool1 |

**Analyses**

- **Essentiality (HMM)**, `essentiality_H1_new_BG2`
  - libraries: SRR38338548; SRR38338549
  - `ytab essentiality --project-config output/projects/Zaeske_2026_glabrata/config/project.yaml --libraries SRR33001989 SRR33001993 SRR33002022 SRR33002027 SRR33002031 SRR33002033 SRR38338548 SRR38338549 --output output/projects/Zaeske_2026_glabrata/essentiality/essentiality_H1_new_BG2`
- **Essentiality (HMM)**, `essentiality_glabrata_DPL1021r`
  - libraries: SRR38338545
  - `ytab essentiality --project-config output/projects/Zaeske_2026_glabrata/config/project.yaml --libraries SRR38338545 --output output/projects/Zaeske_2026_glabrata/essentiality/essentiality_glabrata_DPL1021r`

## Zaeske et al. 2026: Candida glabrata (Zaeske_2026_glabrata_CBS138)

BG14, BG2r, DPL1021r, CBS138r; Hermes (NAT), splinkerette PCR. BioProject PRJNA1460454. 5 parent pools (one deposit duplicates a Nickels & Cunningham run and is counted once).

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| SRR38338546 |  | parent | parent | CBS138r | CBS138r_pool1 |

**Analyses**

- **Essentiality (HMM)**, `essentiality_CBS138r_Zaeske`
  - libraries: SRR38338546
  - `ytab essentiality --project-config output/projects/Zaeske_2026_glabrata_CBS138/config/project.yaml --libraries SRR38338546 --output output/projects/Zaeske_2026_glabrata_CBS138/essentiality/essentiality_CBS138r_Zaeske`

## Billmyre et al. 2022: Schizosaccharomyces pombe (Billmyre_2022_pombe)

wild type; Hermes, left-end adapter with UMI. BioProject PRJNA758956. glycerol, mitotic, short and long meiotic outgrowth; 2 libraries each, each sequenced twice.

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| 3786_glycerol_libSRR15669607 | SRR15669607 SRR15669608 | parent | glycerol | 3786 | 3786 |
| 3786_glycerol_libSRR15669609 | SRR15669609 SRR15669610 | parent | glycerol | 3786 | 3786 |
| 3786_mitotic_libSRR15669611 | SRR15669611 SRR15669612 | control | mitotic | 3786 | 3786 |
| 3786_mitotic_libSRR15669619 | SRR15669619 SRR15669620 | control | mitotic | 3786 | 3786 |
| 3786_meiotic_longOutgrowth_libSRR15669613 | SRR15669613 SRR15669614 | treated | meiotic_longOutgrowth | 3786 | 3786 |
| 3786_meiotic_longOutgrowth_libSRR15669615 | SRR15669615 SRR15669616 | treated | meiotic_longOutgrowth | 3786 | 3786 |
| 3786_meiotic_shortOutgrowth_libSRR15669605 | SRR15669605 SRR15669606 | treated | meiotic_shortOutgrowth | 3786 | 3786 |
| 3786_meiotic_shortOutgrowth_libSRR15669617 | SRR15669617 SRR15669618 | treated | meiotic_shortOutgrowth | 3786 | 3786 |

**Analyses**

Rerun every condition analysis: `ytab run --workflow output/projects/Billmyre_2022_pombe/config/analyses_workflow.yaml`

- **Essentiality (HMM)**, `essentiality_pombe_Billmyre`
  - libraries: 3786_glycerol_libSRR15669607; 3786_glycerol_libSRR15669609; 3786_mitotic_libSRR15669611; 3786_mitotic_libSRR15669619
  - `ytab essentiality --project-config output/projects/Billmyre_2022_pombe/config/project.yaml --libraries 3786_glycerol_libSRR15669607 3786_glycerol_libSRR15669609 3786_mitotic_libSRR15669611 3786_mitotic_libSRR15669619 --output output/projects/Billmyre_2022_pombe/essentiality/essentiality_pombe_Billmyre`
- **Treated vs control (Resampling, betageom)**, `treated_vs_control_BGC`
  - 3786, meiotic_longOutgrowth: 3786_mitotic_libSRR15669611 (control); 3786_mitotic_libSRR15669619 (control); 3786_meiotic_longOutgrowth_libSRR15669613 (treated); 3786_meiotic_longOutgrowth_libSRR15669615 (treated)
  - 3786, meiotic_shortOutgrowth: 3786_mitotic_libSRR15669611 (control); 3786_mitotic_libSRR15669619 (control); 3786_meiotic_shortOutgrowth_libSRR15669605 (treated); 3786_meiotic_shortOutgrowth_libSRR15669617 (treated)
  - `ytab treated-vs-control --project-config output/projects/Billmyre_2022_pombe/config/project.yaml --method resampling --normalization-run output/projects/Billmyre_2022_pombe/normalization/conditions_v1 --normalization betageom --output output/projects/Billmyre_2022_pombe/treated_vs_control/treated_vs_control_BGC`
- **Treated vs control (z-score, CPM)**, `meiotic_longOutgrowth_vs_mitotic`
  - 3786, meiotic_longOutgrowth: 3786_mitotic_libSRR15669611 (control); 3786_meiotic_longOutgrowth_libSRR15669613 (treated); 3786_mitotic_libSRR15669619 (control); 3786_meiotic_longOutgrowth_libSRR15669615 (treated)
  - `ytab fitness --project-config output/projects/Billmyre_2022_pombe/config/project.yaml --analysis-id meiotic_longOutgrowth_vs_mitotic --comparisons 3786_meiotic_longOutgrowth_libSRR15669613_vs_3786_mitotic_libSRR15669611,3786_meiotic_longOutgrowth_libSRR15669613_vs_3786_mitotic_libSRR15669619,3786_meiotic_longOutgrowth_libSRR15669615_vs_3786_mitotic_libSRR15669611,3786_meiotic_longOutgrowth_libSRR15669615_vs_3786_mitotic_libSRR15669619`
- **Treated vs control (z-score, CPM)**, `meiotic_shortOutgrowth_vs_mitotic`
  - 3786, meiotic_shortOutgrowth: 3786_mitotic_libSRR15669611 (control); 3786_meiotic_shortOutgrowth_libSRR15669605 (treated); 3786_mitotic_libSRR15669619 (control); 3786_meiotic_shortOutgrowth_libSRR15669617 (treated)
  - `ytab fitness --project-config output/projects/Billmyre_2022_pombe/config/project.yaml --analysis-id meiotic_shortOutgrowth_vs_mitotic --comparisons 3786_meiotic_shortOutgrowth_libSRR15669605_vs_3786_mitotic_libSRR15669611,3786_meiotic_shortOutgrowth_libSRR15669605_vs_3786_mitotic_libSRR15669619,3786_meiotic_shortOutgrowth_libSRR15669617_vs_3786_mitotic_libSRR15669611,3786_meiotic_shortOutgrowth_libSRR15669617_vs_3786_mitotic_libSRR15669619`
- **Multi-condition (betageom)**, `3786_conditions_BGC`
  - mitotic, meiotic_longOutgrowth, meiotic_shortOutgrowth: 3786_mitotic_libSRR15669611 (mitotic (reference)); 3786_mitotic_libSRR15669619 (mitotic (reference)); 3786_meiotic_longOutgrowth_libSRR15669613 (meiotic_longOutgrowth); 3786_meiotic_longOutgrowth_libSRR15669615 (meiotic_longOutgrowth); 3786_meiotic_shortOutgrowth_libSRR15669605 (meiotic_shortOutgrowth); 3786_meiotic_shortOutgrowth_libSRR15669617 (meiotic_shortOutgrowth)
  - `ytab multicondition-project --project-config output/projects/Billmyre_2022_pombe/config/project.yaml --normalized output/projects/Billmyre_2022_pombe/normalization/conditions_v1/experiment/methods --normalization betageom --sample Billmyre_2022_pombe::3786_mitotic_libSRR15669611=mitotic --sample Billmyre_2022_pombe::3786_mitotic_libSRR15669619=mitotic --sample Billmyre_2022_pombe::3786_meiotic_longOutgrowth_libSRR15669613=meiotic_longOutgrowth --sample Billmyre_2022_pombe::3786_meiotic_longOutgrowth_libSRR15669615=meiotic_longOutgrowth --sample Billmyre_2022_pombe::3786_meiotic_shortOutgrowth_libSRR15669605=meiotic_shortOutgrowth --sample Billmyre_2022_pombe::3786_meiotic_shortOutgrowth_libSRR15669617=meiotic_shortOutgrowth --reference-condition mitotic --output output/projects/Billmyre_2022_pombe/multicondition/3786_conditions_BGC`

## Kingma et al. 2025: Saccharomyces cerevisiae (Kingma_2025_SATAY)

W303 background: wild type (yWT01) and nrp1Δ (yLIC136); Ac/Ds (MiniDs), restriction digestion and circularization; Ds-end adapter. ENA PRJEB80469 (ERR13705597–ERR13705640). 11 wild-type libraries; 2 nrp1Δ libraries, each 2 digests × 2 technical replicates × 4 lanes (32 runs); 1 run of an undescribed strain left out.

| library | runs | role | condition | strain | pool |
|---|---|---|---|---|---|
| nrp1D_1_nrp1_deletion_libERR13705597 | ERR13705597 ERR13705598 ERR13705599 ERR13705600 ERR13705601 ERR13705602 ERR13705603 ERR13705604 ERR13705605 ERR13705606 ERR13705607 ERR13705608 ERR13705609 ERR13705610 ERR13705611 ERR13705612 | parent | nrp1_deletion | yLIC136_nrp1D | nrp1D_1 |
| nrp1D_2_nrp1_deletion_libERR13705613 | ERR13705617 ERR13705618 ERR13705619 ERR13705620 ERR13705613 ERR13705614 ERR13705615 ERR13705616 ERR13705625 ERR13705626 ERR13705627 ERR13705628 ERR13705621 ERR13705622 ERR13705623 ERR13705624 | parent | nrp1_deletion | yLIC136_nrp1D | nrp1D_2 |
| ERR13705640 (not analysed) |  | parent | unannotated | yLIC137 | yLIC137_7 |
| ERR13705629 |  | parent | wild_type | yWT01 | FD11 |
| ERR13705630 |  | parent | wild_type | yWT01 | FD12 |
| ERR13705631 |  | parent | wild_type | yWT01 | FD7 |
| ERR13705632 |  | parent | wild_type | yWT01 | FD7 |
| ERR13705633 |  | parent | wild_type | yWT01 | FD7 |
| ERR13705634 |  | parent | wild_type | yWT01 | FD7 |
| ERR13705635 |  | parent | wild_type | yWT01 | FD7 |
| ERR13705636 |  | parent | wild_type | yWT01 | FD7 |
| ERR13705637 |  | parent | wild_type | yWT01 | FD9 |
| ERR13705638 |  | parent | wild_type | yWT01 | FD9 |
| ERR13705639 |  | parent | wild_type | yWT01 | FD9 |

**Analyses**

Rerun every condition analysis: `ytab run --workflow output/projects/Kingma_2025_SATAY/config/analyses_workflow.yaml`

- **Essentiality (HMM)**, `essentiality_cerevisiae_Kingma_nrp1`
  - libraries: nrp1D_1_nrp1_deletion_libERR13705597; nrp1D_2_nrp1_deletion_libERR13705613
  - `ytab essentiality --project-config output/projects/Kingma_2025_SATAY/config/project.yaml --libraries nrp1D_1_nrp1_deletion_libERR13705597 nrp1D_2_nrp1_deletion_libERR13705613 --output output/projects/Kingma_2025_SATAY/essentiality/essentiality_cerevisiae_Kingma_nrp1`
- **Essentiality (HMM)**, `essentiality_cerevisiae_Kingma_wild_type`
  - libraries: ERR13705629; ERR13705630; ERR13705631; ERR13705632; ERR13705633; ERR13705634; ERR13705635; ERR13705636; ERR13705637; ERR13705638; ERR13705639
  - `ytab essentiality --project-config output/projects/Kingma_2025_SATAY/config/project.yaml --libraries ERR13705629 ERR13705630 ERR13705631 ERR13705632 ERR13705633 ERR13705634 ERR13705635 ERR13705636 ERR13705637 ERR13705638 ERR13705639 --output output/projects/Kingma_2025_SATAY/essentiality/essentiality_cerevisiae_Kingma_wild_type`
- **Multi-condition (betageom)**, `nrp1_deletion_vs_wild_type_BGC`
  - nrp1_deletion, wild_type: ERR13705629 (wild_type (reference)); ERR13705630 (wild_type (reference)); ERR13705631 (wild_type (reference)); ERR13705632 (wild_type (reference)); ERR13705633 (wild_type (reference)); ERR13705634 (wild_type (reference)); ERR13705635 (wild_type (reference)); ERR13705636 (wild_type (reference)); ERR13705637 (wild_type (reference)); ERR13705638 (wild_type (reference)); ERR13705639 (wild_type (reference)); nrp1D_1_nrp1_deletion_libERR13705597 (nrp1_deletion); nrp1D_2_nrp1_deletion_libERR13705613 (nrp1_deletion)
  - `ytab multicondition-project --project-config output/projects/Kingma_2025_SATAY/config/project.yaml --normalized output/projects/Kingma_2025_SATAY/normalization/conditions_v1/experiment/methods --normalization betageom --sample Kingma_2025_SATAY::ERR13705629=wild_type --sample Kingma_2025_SATAY::ERR13705630=wild_type --sample Kingma_2025_SATAY::ERR13705631=wild_type --sample Kingma_2025_SATAY::ERR13705632=wild_type --sample Kingma_2025_SATAY::ERR13705633=wild_type --sample Kingma_2025_SATAY::ERR13705634=wild_type --sample Kingma_2025_SATAY::ERR13705635=wild_type --sample Kingma_2025_SATAY::ERR13705636=wild_type --sample Kingma_2025_SATAY::ERR13705637=wild_type --sample Kingma_2025_SATAY::ERR13705638=wild_type --sample Kingma_2025_SATAY::ERR13705639=wild_type --sample Kingma_2025_SATAY::nrp1D_1_nrp1_deletion_libERR13705597=nrp1_deletion --sample Kingma_2025_SATAY::nrp1D_2_nrp1_deletion_libERR13705613=nrp1_deletion --reference-condition wild_type --output output/projects/Kingma_2025_SATAY/multicondition/nrp1_deletion_vs_wild_type_BGC`
