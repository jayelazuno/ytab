"""Transposon profiles: how reads of each transposon and library design are trimmed before mapping.

resources/transposons/<transposon>/profile.json lists the transposon's library designs, each with the MapFastq
settings it needs (hermes_cleanup, trim_5p_seq, trim_5p_overlap, trim_5p_times, discard_untrimmed, trim_3p_seq). A project selects
one with `transposon` and `library_design` in its config's `mapping` block or as sample-sheet columns; explicit
settings there override the design's. Without a selection the transposon is Hermes with its default design.
"""
from __future__ import annotations

import json
from pathlib import Path

DEFAULT_TRANSPOSON = "hermes"
MAPPING_KEYS = ("hermes_cleanup", "trim_5p_seq", "trim_5p_overlap", "trim_5p_times", "discard_untrimmed", "trim_3p_seq", "donor_fasta")


def profile_path(transposon: str, root: Path) -> Path:
    return Path(root)/"resources"/"transposons"/transposon/"profile.json"


def available_transposons(root: Path) -> list[str]:
    return sorted(p.parent.name for p in (Path(root)/"resources"/"transposons").glob("*/profile.json"))


def target_site(transposon: str, root: Path) -> dict | None:
    """{motif, duplication_bp} when the transposon inserts only at a fixed motif, otherwise None.

    piggyBac declares TTAA with a 4 bp duplication; the hAT elements (Hermes, Ac/Ds) declare nothing, because they
    insert at nearly any position. A recorded position is on target when the motif starts there or `duplication_bp-1`
    bases before it, the two ends of the one duplicated site.
    """
    path = profile_path(transposon, root)
    if not path.is_file():return None
    site = json.loads(path.read_text()).get("target_site") or {}
    motif = str(site.get("motif") or "").upper()
    return dict(motif=motif, duplication_bp=int(site.get("duplication_bp") or len(motif))) if motif else None


def on_target(sequence: str, position: int, site: dict) -> bool:
    """Is a 1-based recorded hit position one of the two ends of an occupied target site in `sequence`?"""
    motif, offset = site["motif"], site["duplication_bp"]-1
    # Both starts are clamped at 0: a negative slice start would wrap to the end of the chromosome and compare
    # against the wrong sequence for the first few bases.
    return sequence[position-1:position-1+len(motif)] == motif or \
           (position-1-offset >= 0 and sequence[position-1-offset:position-1-offset+len(motif)] == motif)


def resolve_design(settings: dict, root: Path) -> dict:
    """{transposon, library_design, profile, mapping}: the design's settings under any explicit ones in `settings`."""
    transposon = str(settings.get("transposon") or DEFAULT_TRANSPOSON).strip()
    path = profile_path(transposon, root)
    if not path.is_file():
        raise ValueError(f"Unknown transposon '{transposon}'. Available: {', '.join(available_transposons(root)) or 'none'}")
    profile = json.loads(path.read_text())
    design = str(settings.get("library_design") or profile["default_design"]).strip()
    if design not in profile["designs"]:
        raise ValueError(f"Unknown library design '{design}' for {transposon}. Available: {', '.join(profile['designs'])}")
    mapping = dict(profile["designs"][design].get("mapping") or {})
    for key in MAPPING_KEYS:
        if settings.get(key) not in (None, ""):
            mapping[key] = settings[key]
    if isinstance(mapping.get("trim_3p_seq"), str):
        mapping["trim_3p_seq"] = [s for s in mapping["trim_3p_seq"].split(";") if s]
    return dict(transposon=transposon, library_design=design, profile=path, mapping=mapping)
