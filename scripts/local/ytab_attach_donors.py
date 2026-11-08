#!/usr/bin/env python3
"""Attach the transposon donor of each project (mapping.donor_fasta in config/project.yaml), so remapping counts the
reads that come from the donor rather than from an insertion (ytab.mapping.donor_reads).

A donor belongs to a study, not to a transposon: the same Hermes is delivered on pCU-MET3-Hermes in the Cunningham-lab
C. glabrata pools and on pHL2577 in the S. pombe libraries. DONORS maps each project to the donor built under
resources/transposons/donors/; projects whose donor is unknown keep the recurring-unmapped-flank estimate, which needs
no sequence.

  python scripts/local/ytab_attach_donors.py [--write]      (default: show what would change)
"""
from __future__ import annotations
import argparse
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[2]; DONORS_DIR = Path("resources/transposons/donors")
DONORS = {
    # Cunningham-lab Hermes pools: the lab's own pCU-MET3-Hermes (real sequence)
    "Gale_2020":"hermes_pCU-MET3-Hermes", "Nickels_2024_glabrata":"hermes_pCU-MET3-Hermes",
    "Nickels_2025_glabrata":"hermes_pCU-MET3-Hermes", "Nickels_2025_glabrata_CBS138":"hermes_pCU-MET3-Hermes",
    "Nickels_2026_glabrata_FKS2":"hermes_pCU-MET3-Hermes", "Nickels_2026_glabrata_tolerance":"hermes_pCU-MET3-Hermes",
    "Zaeske_2026_glabrata":"hermes_pCU-MET3-Hermes", "Zaeske_2026_glabrata_CBS138":"hermes_pCU-MET3-Hermes",
    "H2O2_screen_v1":"hermes_pCU-MET3-Hermes", "Zn_toxicity_screen":"hermes_pCU-MET3-Hermes",
    "local_glabrata_smoke_v1":"hermes_pCU-MET3-Hermes",
    # S. pombe Hermes: pHL2577 (rebuilt from parts)
    "Pombe_Hermes":"hermes_pHL2577_pombe", "Billmyre_2022_pombe":"hermes_pHL2577_pombe",
    # Ac/Ds
    "Michel_2017_SATAY":"acds_pBK549_satay", "Kingma_2025_SATAY":"acds_pBK549_satay",
    "Billmyre_2025_neoformans":"acds_pBB51_neoformans",
    # piggyBac
    "Yao_2026_phaffii":"piggybac_pIMY015_phaffii",
    # No donor sequence yet (chromosomal piggyBac donor at ARG4; Hermes donor not described):
    #   Chow_2024_glabrata, Song_2026_auris, Gao_2021_auris, Segal_2018_albicans, Patterson_2018_lipolytica
}


def main(write: bool):
    rows = []
    for project, donor in sorted(DONORS.items()):
        config = ROOT/"output/projects"/project/"config/project.yaml"
        fasta = DONORS_DIR/donor/"donor.fasta"
        if not config.is_file() or not (ROOT/fasta).is_file():
            rows.append((project, donor, "missing")); continue
        text = config.read_text(); data = yaml.safe_load(text)
        current = (data.get("mapping") or {}).get("donor_fasta")
        if current == str(fasta):rows.append((project, donor, "already set")); continue
        if write:
            data.setdefault("mapping", {})["donor_fasta"] = str(fasta)
            config.write_text(yaml.safe_dump(data, sort_keys=False))
        rows.append((project, donor, "set" if write else "would set"))
    for project, donor, state in rows:print(f"{project:34} {donor:28} {state}")
    print(f"\n{len(rows)} projects; donors under {DONORS_DIR}. Projects not listed keep the flank estimate.")


if __name__ == "__main__":
    a = argparse.ArgumentParser(); a.add_argument("--write", action="store_true"); main(a.parse_args().write)
