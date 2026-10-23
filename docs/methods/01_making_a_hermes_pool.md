# Making a Hermes Tn-seq pool in yeast

*This page is the editable version of the lab protocol (converted from `docs/20251016-Hermes-Tn-seq-Induction .docx`);
edit it here and the app shows the change. Day-by-day notes and QC values
are to be added from the current induction (started 2026-10-19/20). Items marked **[to confirm]** need a check by the
bench scientist before publication.*

This guide follows one induction from transformation to frozen pool. Double every quantity for two inductions in
parallel. The aims are a pool that is saturated (insertions in as many genes as possible) and even (no single early
insertion taking over), with a record of the QC values YTAB needs downstream.

## Principle

The plasmid pCU-MET3-HERMES carries:
- the **Hermes transposase** under the **MET3 promoter**, repressed by cysteine and methionine, so transposition is
  switched on by growing cells without Cys and Met;
- a **URA3** marker, which keeps selection on the plasmid and lets 5-FOA select against it later;
- the **Hermes transposon** carrying a **NAT** resistance cassette.

Cells that have both inserted the transposon in the genome and lost the plasmid grow on **SCD + NAT + 5-FOA**.
- NAT selects for the transposon.
- 5-FOA kills cells that still have URA3, i.e. still have the plasmid.

Induction is split across **40 tubes** so that a rare early transposition, which then multiplies, can take over only
one fortieth of the pool.

## Materials (one induction)

| item | amount |
|---|---|
| YPD plates | 9 (+12 optional) |
| SCD−Ura plates | 10 (+12 optional) |
| SCD + NAT + 5-FOA plates | 8 (+12 optional) |
| YPD liquid | 15–50 mL |
| SCD−cmu (−Cys −Met −Ura) | 110 mL |
| 1× PBS pH 7.4 / 2× PBS | ~250 mL / 30 mL |
| SCD (complete) | 2.4 L |
| autoclaved green-cap glass tubes (disposable) | 41 |
| autoclaved red-cap glass tubes (reusable) | 2–5 |
| autoclaved 250 mL glass bottle with lid | 1 |
| autoclaved 2 L Erlenmeyer flasks | 4 |
| autoclaved 450 mL centrifuge jars | 4 |
| autoclaved 1.5 mL tubes | 24 (+24 optional) |
| sterile 50 mL / 15 mL conical tubes | 2 / 12 |
| 5-FOA | 2.4 g |
| NAT (clonNAT) | 300 mg |
| 30% glycerol, filter-sterilized | 30 mL |
| TE-LiOAc / PEG-LiTE | 1.1 mL / 800 µL (make 10 mL each, fresh) |
| 100× methionine | 100 µL |

## Media and reagents

| reagent | total | made as | sterilize | used |
|---|---|---|---|---|
| SCD | 2.4 L | 4 × 600 mL | filter | days 3–6 |
| SCD−cmu | 110 mL | 1 × 110 mL (make ~1 L) | filter | day 0 inoculum |
| PBS (1×) | ~250 mL | 1 bottle | filter | washes, dilutions |
| 2× PBS | 30 mL | 1 × 30 mL | filter | freezing |
| 30% glycerol | 30 mL | 1 × 30 mL | filter | freezing |
| YPD | 15–50 mL | 1 bottle | autoclave | transformation |
| 5-FOA, 120 mg/mL in DMSO | 2.4 g in 20 mL | 4 × 5 mL | none needed | plates, days 4–6 |
| NAT, 100 mg/mL in water | 300 mg in 3 mL | 4 × 600 µL | filter | plates, days 4–6 |
| TE-LiOAc | 1.1 mL | 10 mL fresh | filter | transformation |
| PEG-LiTE | 800 µL | 10 mL fresh | filter | transformation |
| 100× methionine | 100 µL | 1 mL | filter | control tube |

One 5-FOA aliquot (600 mg) and one NAT aliquot (60 mg) drug one 600 mL bottle of SCD (final 1 mg/mL 5-FOA, 100 µg/mL
NAT). Three pairs go to the enrichment days and one to the selective plates.

**Recipes:**
- **SCD (per L):** 6.8 g YNB without amino acids (with ammonium sulfate), 2.0 g complete amino-acid supplement (CSM),
  960 mL water. Autoclave, cool to 55 °C, then add 40 mL dextrose (50% stock **[to confirm]**). Make 2.4 L as one batch
  and split it into 4 × 600 mL, so every enrichment day gets identical medium.
- **SCD−cmu (per L):** as SCD, with a −Cys −Met −Ura supplement. Residual methionine is the most common cause of a weak
  induction.
- **100× methionine:** 2 mg/mL L-methionine (100× a 20 mg/L SC recipe), filtered and stored at −20 °C.
- **5-FOA:**
  - Dissolve 2.4 g in DMSO to 20 mL; a 37 °C bath or 20 min on a nutator helps it dissolve. Dissolve it completely
    before aliquoting.
  - DMSO freezes at about 18 °C: thaw for 30 min.
  - Filtering is not needed; if you filter, use a PTFE or nylon membrane.
- **NAT:** 300 mg in 3 mL water, filtered, 4 × 600 µL. NAT decays even at 4 °C, so add it to the medium on the day of
  use.
- **TE-LiOAc (fresh):** 1 mL 10× TE (100 mM Tris-HCl, 10 mM EDTA, pH 7.5), 1 mL 10× LiOAc (1 M, pH 7.5), 8 mL water.
- **PEG-LiTE (fresh):** 8 mL 50% PEG 3350, 1 mL 10× TE, 1 mL 10× LiOAc. Pipette the PEG slowly and mix thoroughly by
  inversion.
- **Freezing medium:** 30 mL 2× PBS + 30 mL 30% glycerol = 60 mL of 15% glycerol in 1× PBS.

**Plates:**
- **Base media:**
  - YPD (9);
  - SCD−Ura (10);
  - SCD + NAT + 5-FOA (8). This one is on complete SCD: cells need uracil available for 5-FOA to report plasmid loss.
- **Adding the drugs:** cool the agar to 55 °C before adding them; NAT loses activity when hot. Add the 5-FOA, then the
  NAT, and pour promptly.
- **Storage:** keep plates at 4 °C in the dark and use them within a few weeks.

## Prep checklist

Make ahead:
- [ ] SCD 4 × 600 mL
- [ ] SCD−cmu
- [ ] PBS, 2× PBS
- [ ] 30% glycerol
- [ ] 5-FOA 4 aliquots
- [ ] NAT 4 aliquots
- [ ] 100× methionine
- [ ] 10× TE and 10× LiOAc stocks
- [ ] 50% PEG 3350
- [ ] plates: 9 YPD, 10 SCD−Ura, 8 SCD + NAT + 5-FOA

Autoclave ahead:
- [ ] 41 green-cap tubes
- [ ] 2–5 red-cap tubes
- [ ] 250 mL bottle
- [ ] 4 × 2 L flasks
- [ ] 4 centrifuge jars
- [ ] inoculation sticks
- [ ] 24 × 1.5 mL tubes

**The one rule worth repeating:** drug the SCD the morning you use it. A day-old bottle of NAT medium lets untransposed
cells survive enrichment.

## Step 1 — Transform pCU-MET3-HERMES into the parent strain (ura3⁻)

1. **Streak** the ura3⁻ parent from the −80 °C stock onto YPD and grow 2 days at 30 °C.
2. **Grow to early log:** inoculate one colony into YPD overnight. Dilute in the morning and grow to OD600 0.3–0.4
   (about 3.5 h). Use 5–10 mL per transformation.
3. **Pellet and wash:**
   - Spin at 3000 × g, 5 min.
   - Meanwhile boil the salmon-sperm DNA at 99 °C for 10 min.
   - Resuspend in 500 µL water, move to a 1.5 mL tube and spin 30 s.
   - Wash in 1 mL TE-LiOAc.
4. **Transform:**
   - Resuspend in 90 µL water with 150–300 ng plasmid and 10 µL salmon-sperm DNA.
   - Add 900 µL PEG-LiTE and mix; C. glabrata tolerates short vortex pulses.
   - Incubate 30–60 min at 30 °C, mixing every 15 min.
   - Add 100 µL DMSO and heat-shock at 42 °C for 20 min.
5. **Plate:** spin, wash in 1 mL PBS, resuspend in 200 µL PBS and plate on SCD−Ura. Incubate 2 days at 30 °C.
   **[to confirm]** The protocol also mentions 37 °C for C. glabrata.
6. **Restreak** 4–6 colonies on SCD−Ura for single colonies and grow 2 more days. Induce soon; fresh cells work best.

**QC 1:**
- The untransformed parent does not grow on SCD−Ura (negative control), while transformants do.
- Record the strain, plasmid batch and transformant colonies used.

## Step 2 — Induce the transposon (days 1–3)

1. **Inoculate:** put 110 mL SCD−cmu into the 250 mL bottle, add one large transformant colony and swirl hard, without
   shaking up and down.
2. **Fill the tubes:** dispense 2.5 mL into each of the 40 green-cap tubes with a 5 mL serological pipette.
3. **Set up the control:** put the rest into tube 41 and note its volume. Add 1/100 volume of 100× methionine (e.g. 70 µL
   to 7 mL); this represses the transposase.
4. **Induce:** shake at 30 °C, 225 rpm, for 3 days. Transposition is most efficient in stationary phase.

**Day 3, pool:**

5. **Pool the tubes:** collect all 40 tubes into a centrifuge jar.
6. **Plate (CFU checkpoint 1):**
   - Mix 100 µL of the pool into 900 µL PBS and make seven 1:10 serial dilutions.
   - Plate them on YPD, SCD−Ura and SCD + NAT + 5-FOA (Table 1).
   - Do the same for control tube 41.
7. **Rinse the tubes:** rinse each tube with 500 µL PBS into the jar and spin at 3000 × g for 15 min.
8. **Recover:**
   - Resuspend the pellet in 100 mL SCD and move it to a 2 L flask.
   - Rinse the jar with 100 mL SCD twice more and bring the flask to 600 mL.
   - Shake 24 h at 30 °C. This recovery in SCD lets cells express NAT resistance and use up any remaining Ura3 protein
     before selection.

## Step 3 — Enrich for transposed cells (days 4–5)

1. **Plate (CFU checkpoint):** take 100 µL for the dilution series.
2. **Pellet:** spin the 600 mL at 3000 × g for 15 min, rinsing the flask 3 × 25 mL PBS into the jars.
3. **Prepare the selective medium:**
   - Thaw the 5-FOA (~30 min) and NAT (~5 min).
   - Vortex, spin down and add both to 600 mL SCD.
4. **Resuspend:** decant into 10% bleach. Resuspend the pellets in SCD + NAT + 5-FOA in 3 × 100 mL washes into a fresh
   2 L flask.
5. **Grow:** shake 24 h at 30 °C.
6. **Clean up:** bleach the glassware for at least 15 min.
7. **Repeat** steps 1–6 on day 5.

## Step 4 — Regrow and freeze the pool (days 6–7)

1. **Day 6:**
   - Take 100 µL for dilutions (CFU checkpoint).
   - Spin 60 mL (a 1:10 transfer) in two 50 mL tubes at 3000 × g, 5 min.
   - Resuspend in SCD + NAT + 5-FOA and wash 3 × 30 mL into a 2 L flask.
   - Shake 24 h at 30 °C.
2. **Day 7:**
   - Take 100 µL for dilutions (final CFU checkpoint).
   - Pellet all 600 mL, washing the flask 3 × with PBS.
   - Resuspend in 30 mL 2× PBS + 30 mL 30% glycerol.
3. **Freeze:**
   - Aliquot into 12 labelled 15 mL tubes: strain, induction number, initials, date.
   - Freeze at −80 °C.
   - Do not freeze-thaw more than 3–4 times, and mark the tube at every thaw.

## QC during induction: CFU worksheet

Cells/mL = CFU / 0.2 × dilution factor (200 µL plated). Example: 87 colonies at 10⁻⁶ gives 87 / 0.2 × 10⁶ =
4.35 × 10⁸ cells/mL.

| plate | measures |
|---|---|
| YPD | total viable cells |
| SCD−Ura | cells still carrying the plasmid |
| SCD + NAT + 5-FOA | cells with a genomic insertion that have lost the plasmid, i.e. pool members |

From these:
- **Transposition frequency** = (SCD + NAT + 5-FOA) / YPD. Compare the induced pool with the methionine control (tube
  41), which shows background (leaky transposase, spontaneous 5-FOA resistance).
- **Enrichment** is working when the SCD + NAT + 5-FOA / YPD ratio rises towards 1 over days 4–6, and SCD−Ura counts
  fall.
- **Complexity:** the number of independent transposition events in the pool is about the SCD + NAT + 5-FOA CFU/mL on
  day 3, times the culture volume. This bounds how many distinct insertion sites sequencing can find. Compare it with
  YTAB's Screen Planning target (about 3× the represented genes for 95% gene coverage).

| day | sample | YPD CFU/mL | SCD−Ura CFU/mL | SCD + NAT + 5-FOA CFU/mL | ratio NAT/FOA : YPD | notes |
|---|---|---|---|---|---|---|
| 3 | induced pool | | | | | |
| 3 | Met control (41) | | | | | |
| 4 | | | | | | |
| 5 | | | | | | |
| 6 | | | | | | |
| 7 | final pool | | | | | |

Table 1 (dilutions to plate per medium and day) is to be added from the docx **[to confirm]**.

## Induction log (current induction)

- [x] Day 0 (2026-10-19/20 **[to confirm]**): streaked yH001 (ura3⁻) and the control. The negative control on SCD−Ura
  confirms the strain is ura3⁻.
- [ ] Day 1: transformation / inoculation
- [ ] Day 2
- [ ] Day 3: pool 40 tubes, CFU checkpoint 1, SCD recovery
- [ ] Day 4: enrichment 1 (SCD + NAT + 5-FOA), CFU
- [ ] Day 5: enrichment 2, CFU
- [ ] Day 6: regrowth, CFU
- [ ] Day 7: final CFU, freeze 12 aliquots

## After the pool: sequencing and YTAB

YTAB takes over from the frozen pool:
1. Sequence parent libraries, and treated libraries for screens.
2. **Run** (mapping, hit files, SummaryTables).
3. **Quality Control** (library QC table, saturation, sequence bias).
4. **Normalization** (BGC by default).
5. **Essentiality** (HMM).
6. **Condition Comparisons** (z-score for ranking; resampling or Gumbel for error-controlled hits).
7. **Screen Planning** (colonies needed) and **Mutant Identification** (single mutants by arbitrary PCR).

The methods pages in the app describe each step.
