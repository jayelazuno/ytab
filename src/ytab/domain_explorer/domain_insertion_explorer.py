"""YTAB-native gene/domain insertion explorer.

This module adapts the useful plotting concepts from the legacy
``DomainFigures.py`` guide into current YTAB project structure.  It does not
import from ``codex/`` and it does not run or modify scientific pipeline
calculations.  Inputs are existing project annotations plus existing
CreateHitFile/combined-hit outputs.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


REPO_ROOT = Path(__file__).resolve().parents[3]
IDENTIFIER_FIELDS = (
    "gene_id",
    "ID",
    "Name",
    "gene",
    "locus_tag",
    "standard_name",
    "common_name",
    "Alias",
    "aliases",
    "product",
    "systematic_name",
    "feature_name",
)
_PROJECT_CONTEXT_CACHE: dict[str, dict[str, Any]] = {}
_GENE_LOOKUP_CACHE: dict[str, dict[str, Any]] = {}


@dataclass(frozen=True)
class GeneRecord:
    gene_id: str
    display_name: str
    standard_name: str
    common_name: str
    systematic_name: str
    chromosome: str
    start: int
    end: int
    strand: str
    product: str
    aliases: tuple[str, ...]
    exons: tuple[tuple[int, int], ...]
    domains: tuple[tuple[int, int], ...]
    domain_source: str

    def as_candidate(self) -> dict[str, Any]:
        return {
            "gene_id": self.gene_id,
            "display_name": self.display_name,
            "standard_name": self.standard_name,
            "common_name": self.common_name,
            "systematic_name": self.systematic_name,
            "chromosome": self.chromosome,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "product": self.product,
        }

    def as_dict(self) -> dict[str, Any]:
        data = self.as_candidate()
        data.update(
            {
                "aliases": list(self.aliases),
                "exons": [list(x) for x in self.exons],
                "domains": [list(x) for x in self.domains],
                "domain_source": self.domain_source,
            }
        )
        return data


def _repo_path(path: str | Path | None, base: Path | None = None) -> Path | None:
    if path is None or str(path) == "":
        return None
    p = Path(path)
    if p.is_absolute():
        return p
    return ((base or REPO_ROOT) / p).resolve()


def _read_project_config(project_config: str | Path) -> dict[str, Any]:
    path = _repo_path(project_config)
    if path is None or not path.is_file():
        raise FileNotFoundError(f"Project config not found: {project_config}")
    with path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    cfg["_project_config_path"] = str(path)
    cfg["_repo_root"] = str(REPO_ROOT)
    cfg["_project_root"] = str(_repo_path(cfg.get("output_project_dir"), REPO_ROOT) or path.parents[1])
    return cfg


def _load_reference_resolved(cfg: dict[str, Any]) -> dict[str, Any]:
    project_root = Path(cfg["_project_root"])
    resolved = project_root / "config" / "reference_resolved.json"
    if resolved.is_file():
        with resolved.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        return data
    return dict(cfg.get("reference") or {})


def _safe_text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _uniq(values: list[str] | tuple[str, ...] | set[str]) -> tuple[str, ...]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        text = _safe_text(value)
        if not text or text in seen or text.lower() in {"none", "nan"}:
            continue
        seen.add(text)
        out.append(text)
    return tuple(out)


def _ranges_from_rangeset(value: Any) -> tuple[tuple[int, int], ...]:
    try:
        return tuple((int(a), int(b)) for a, b in value)
    except Exception:
        return tuple()


def _feature_to_record(feature: Any) -> GeneRecord:
    aliases = set(getattr(feature, "all_names", set()) or set())
    for attr in ("primary_name", "standard_name", "common_name", "name", "feature_name", "protein_id"):
        aliases.add(_safe_text(getattr(feature, attr, "")))
    gene_id = _safe_text(getattr(feature, "standard_name", "")) or _safe_text(getattr(feature, "primary_name", ""))
    common = _safe_text(getattr(feature, "common_name", ""))
    display = common or _safe_text(getattr(feature, "name", "")) or gene_id
    domains = _ranges_from_rangeset(getattr(feature, "domains", ()))
    strand = _safe_text(getattr(feature, "strand", ""))
    strand = {"W": "+", "C": "-"}.get(strand, strand)
    return GeneRecord(
        gene_id=gene_id,
        display_name=display,
        standard_name=gene_id,
        common_name=common,
        systematic_name=gene_id,
        chromosome=_safe_text(getattr(feature, "chromosome", "")),
        start=int(getattr(feature, "start")),
        end=int(getattr(feature, "stop")),
        strand=strand,
        product=_safe_text(getattr(feature, "description", "")),
        aliases=_uniq(tuple(aliases)),
        exons=_ranges_from_rangeset(getattr(feature, "exons", ())) or ((int(getattr(feature, "start")), int(getattr(feature, "stop"))),),
        domains=domains,
        domain_source="InterProScan/core feature DB" if domains else "",
    )


def _load_core_feature_records(reference: dict[str, Any], species: str) -> list[GeneRecord]:
    try:
        from ytab.core import GenomicFeatures
    except Exception:
        return []

    gff = _repo_path(reference.get("gff"))
    # This viewer needs coordinates/features/domains, not chromosome sequence.
    # Avoid FASTA parsing here so Shiny searches/plots stay responsive.
    fasta = None
    species_lower = species.lower()
    try:
        if "glab" in species_lower and gff and gff.is_file():
            return []
        elif "alb" in species_lower and reference.get("feature_table"):
            db = GenomicFeatures.AlbicansFeatureDB(
                str(_repo_path(reference.get("feature_table"))),
                str(fasta) if fasta else None,
                str(gff) if gff and gff.is_file() else None,
                str(_repo_path(reference.get("domain_table"))) if reference.get("domain_table") else None,
            )
        elif "pomb" in species_lower and gff and gff.is_file():
            db = GenomicFeatures.PombeFeatureDB(str(gff), str(fasta) if fasta else None)
        else:
            return []
        return [_feature_to_record(feature) for feature in db.get_all_features() if getattr(feature, "start", -1) > 0]
    except Exception:
        return []


def _parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
        else:
            continue
        result[key.strip()] = value.strip().strip('"')
    return result


def _load_gff_gene_records(reference: dict[str, Any]) -> list[GeneRecord]:
    gff = _repo_path(reference.get("gff")) or _repo_path(reference.get("gtf"))
    if not gff or not gff.is_file():
        return []
    records: list[GeneRecord] = []
    with gff.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in {"gene", "pseudogene"}:
                continue
            attrs = _parse_attributes(parts[8])
            gene_id = attrs.get("locus_tag") or attrs.get("gene_id") or attrs.get("ID") or attrs.get("Name")
            if not gene_id:
                continue
            name = attrs.get("gene") or attrs.get("Name") or gene_id
            product = attrs.get("product") or attrs.get("description") or ""
            aliases = _uniq([gene_id, name, attrs.get("ID", ""), attrs.get("Alias", "")])
            records.append(
                GeneRecord(
                    gene_id=gene_id,
                    display_name=name,
                    standard_name=gene_id,
                    common_name=name if name != gene_id else "",
                    systematic_name=gene_id,
                    chromosome=parts[0],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    strand=parts[6],
                    product=product,
                    aliases=aliases,
                    exons=((int(parts[3]), int(parts[4])),),
                    domains=tuple(),
                    domain_source="",
                )
            )
    return records


def _load_feature_table_records(reference: dict[str, Any]) -> list[GeneRecord]:
    table = _repo_path(reference.get("feature_table"))
    if not table or not table.is_file():
        return []
    records: list[GeneRecord] = []
    with table.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("!")), delimiter="\t")
        for row in reader:
            feature = _safe_text(row.get("# feature") or row.get("feature"))
            if feature and feature != "gene":
                continue
            gene_id = _safe_text(row.get("locus_tag") or row.get("GeneID") or row.get("name"))
            if not gene_id:
                continue
            start = _safe_text(row.get("start"))
            end = _safe_text(row.get("end"))
            chrom = _safe_text(row.get("genomic_accession") or row.get("chromosome"))
            if not start or not end or not chrom:
                continue
            symbol = _safe_text(row.get("symbol"))
            name = _safe_text(row.get("name"))
            product = _safe_text(row.get("product_accession") or row.get("attributes"))
            aliases = _uniq([gene_id, symbol, name, row.get("GeneID", "")])
            records.append(
                GeneRecord(
                    gene_id=gene_id,
                    display_name=symbol or name or gene_id,
                    standard_name=gene_id,
                    common_name=symbol,
                    systematic_name=gene_id,
                    chromosome=chrom,
                    start=min(int(start), int(end)),
                    end=max(int(start), int(end)),
                    strand={"plus": "+", "minus": "-", "+": "+", "-": "-"}.get(_safe_text(row.get("strand")), _safe_text(row.get("strand"))),
                    product=product,
                    aliases=aliases,
                    exons=((min(int(start), int(end)), max(int(start), int(end))),),
                    domains=tuple(),
                    domain_source="",
                )
            )
    return records


def load_project_gene_context(project_config: str | Path) -> dict[str, Any]:
    """Load project/reference context and all available gene records."""
    cfg = _read_project_config(project_config)
    cache_key = cfg["_project_config_path"]
    if cache_key in _PROJECT_CONTEXT_CACHE:
        return _PROJECT_CONTEXT_CACHE[cache_key]
    reference = _load_reference_resolved(cfg)
    species = _safe_text(cfg.get("species") or reference.get("species"))
    records = _load_core_feature_records(reference, species)
    source = "ytab.core.GenomicFeatures"
    if not records:
        records = _load_gff_gene_records(reference)
        source = "GFF annotation"
    if not records:
        records = _load_feature_table_records(reference)
        source = "feature table"
    if not records:
        raise ValueError("No gene annotation records were found for this project reference.")
    context = {
        "project_config": cfg["_project_config_path"],
        "project_id": cfg.get("project_id") or Path(cfg["_project_root"]).name,
        "project_root": cfg["_project_root"],
        "species": species,
        "reference": reference,
        "annotation_source": source,
        "genes": [record.as_dict() for record in records],
        "_records": records,
    }
    _PROJECT_CONTEXT_CACHE[cache_key] = context
    return context


def _glabrata_annotation_lookup_path() -> Path:
    return REPO_ROOT / "resources" / "comparative" / "orthology" / "20260610_Cgla_CAGL_to_Scer_annotation_lookup.csv"


def _add_lookup_index(index: dict[str, list[int]], key: str, ix: int) -> None:
    key = _safe_text(key)
    if key:
        index.setdefault(key, []).append(ix)


def _add_glabrata_annotation_aliases(records: list[GeneRecord], exact: dict[str, list[int]], lower: dict[str, list[int]]) -> None:
    path = _glabrata_annotation_lookup_path()
    if not path.is_file():
        return
    by_gwk: dict[str, int] = {}
    for ix, record in enumerate(records):
        for name in _uniq(record.aliases + (record.gene_id, record.standard_name, record.systematic_name)):
            if name.startswith("GWK60_"):
                by_gwk[name.lower()] = ix
    if not by_gwk:
        return
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            gwk = _safe_text(row.get("gwk60_id_clean"))
            ix = by_gwk.get(gwk.lower())
            if ix is None:
                continue
            cagl_id = _safe_text(row.get("cagl_id"))
            scer_gene_name = _safe_text(row.get("scer_gene_name"))
            cgla_common_name = _safe_text(row.get("cgla_common_name"))
            display_name = scer_gene_name or cagl_id or records[ix].display_name
            if display_name and display_name != records[ix].display_name:
                records[ix] = replace(
                    records[ix],
                    display_name=display_name,
                    standard_name=cagl_id or records[ix].standard_name,
                    common_name=cgla_common_name or records[ix].common_name,
                    aliases=_uniq(records[ix].aliases + (cagl_id, scer_gene_name, cgla_common_name,
                                                          _safe_text(row.get("qng_id")),
                                                          _safe_text(row.get("scer_gene_id")))),
                )
            for field in ("cagl_id", "qng_id", "cgla_common_name", "scer_gene_id", "scer_gene_name"):
                value = _safe_text(row.get(field))
                if not value:
                    continue
                _add_lookup_index(exact, value, ix)
                _add_lookup_index(lower, value.lower(), ix)


def build_gene_lookup(project_config: str | Path) -> dict[str, Any]:
    """Build a searchable gene lookup from current project annotation."""
    context = load_project_gene_context(project_config)
    if context["project_config"] in _GENE_LOOKUP_CACHE:
        return _GENE_LOOKUP_CACHE[context["project_config"]]
    records: list[GeneRecord] = context["_records"]
    exact: dict[str, list[int]] = {}
    lower: dict[str, list[int]] = {}
    for i, record in enumerate(records):
        names = set(record.aliases) | {
            record.gene_id,
            record.display_name,
            record.standard_name,
            record.common_name,
            record.systematic_name,
        }
        for name in _uniq(names):
            exact.setdefault(name, []).append(i)
            lower.setdefault(name.lower(), []).append(i)
    if _safe_text(context.get("species")).lower() == "glabrata":
        _add_glabrata_annotation_aliases(records, exact, lower)
    lookup = {
        "project_id": context["project_id"],
        "species": context["species"],
        "annotation_source": context["annotation_source"],
        "records": records,
        "exact": exact,
        "lower": lower,
    }
    _GENE_LOOKUP_CACHE[context["project_config"]] = lookup
    return lookup


def _candidate_rows(records: list[GeneRecord], indexes: list[int], match_type: str) -> list[dict[str, Any]]:
    seen: set[int] = set()
    out: list[dict[str, Any]] = []
    for ix in indexes:
        if ix in seen:
            continue
        seen.add(ix)
        row = records[ix].as_candidate()
        row["match_type"] = match_type
        out.append(row)
    return out


def query_gene(project_config: str | Path, query: str) -> list[dict[str, Any]]:
    """Query genes by exact, case-insensitive exact, then partial match."""
    query = _safe_text(query)
    if not query:
        return []
    lookup = build_gene_lookup(project_config)
    records: list[GeneRecord] = lookup["records"]
    if query in lookup["exact"]:
        return _candidate_rows(records, lookup["exact"][query], "exact")
    lowered = query.lower()
    if lowered in lookup["lower"]:
        return _candidate_rows(records, lookup["lower"][lowered], "case_insensitive_exact")
    matches: list[int] = []
    for ix, record in enumerate(records):
        haystack = " ".join(_uniq(record.aliases + (record.product, record.display_name))).lower()
        if lowered in haystack:
            matches.append(ix)
    return _candidate_rows(records, matches, "partial")


def resolve_gene(project_config: str | Path, query: str) -> dict[str, Any]:
    """Resolve a query to one gene or return candidates for ambiguous/no-match cases."""
    candidates = query_gene(project_config, query)
    if not candidates:
        return {"status": "no_match", "query": query, "candidates": []}
    exact = [row for row in candidates if row.get("match_type") in {"exact", "case_insensitive_exact"}]
    if len(exact) == 1:
        return {"status": "resolved", "query": query, "gene": exact[0], "candidates": candidates}
    if len(candidates) == 1:
        return {"status": "resolved", "query": query, "gene": candidates[0], "candidates": candidates}
    return {"status": "ambiguous", "query": query, "candidates": candidates}


def _sample_rows(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    rows = cfg.get("samples")
    if isinstance(rows, list) and rows:
        return [dict(row) for row in rows]
    sheet = _repo_path(cfg.get("sample_sheet"))
    if not sheet or not sheet.is_file():
        sheet = Path(cfg["_project_root"]) / "config" / "sample_sheet.csv"
    if not sheet.is_file():
        return []
    with sheet.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _infer_track_role(row: dict[str, Any], sample: str) -> str:
    values = [
        _safe_text(row.get("library_role")),
        _safe_text(row.get("fitness_role")),
        _safe_text(row.get("control_or_treated")),
        _safe_text(row.get("condition")),
        _safe_text(row.get("treatment")),
        _safe_text(row.get("guessed_condition")),
        sample,
    ]
    text = " ".join(values).lower()
    if re.search(r"treated|h2o2|zn-treated|1[_ .-]?5\s*mm|1_5mm|1\.5mm", text):
        return "treated"
    if re.search(r"mock|parent|control", text):
        return "parent"
    return ""


def _infer_track_pool(row: dict[str, Any], sample: str) -> str:
    for field in ("pool", "pool_id", "replicate_id", "guessed_pool"):
        value = _safe_text(row.get(field))
        if value:
            match = re.search(r"(\d+)", value)
            if match:
                return match.group(1)
    match = re.search(r"pool[_ -]?(\d+)", sample, flags=re.IGNORECASE)
    return match.group(1) if match else ""


def _track_sort_key(track: dict[str, Any]) -> tuple[int, int, int, str]:
    role_rank = {"parent": 0, "treated": 1}.get(_safe_text(track.get("role")).lower(), 2)
    pool = _safe_text(track.get("pool"))
    pool_rank = int(pool) if pool.isdigit() else 9999
    return (role_rank, pool_rank, int(track.get("sample_order", 9999)), str(track.get("sample", "")))


def _order_tracks_for_display(tracks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if any(_safe_text(track.get("role")) for track in tracks):
        return sorted(tracks, key=_track_sort_key)
    return sorted(tracks, key=lambda track: int(track.get("sample_order", 9999)))


def _matched_pair_pools(tracks: list[dict[str, Any]]) -> list[str]:
    pools = sorted(
        {
            _safe_text(track.get("pool"))
            for track in tracks
            if _safe_text(track.get("pool"))
        },
        key=lambda value: (int(value) if value.isdigit() else 9999, value),
    )
    return [
        pool
        for pool in pools
        if any(track.get("role") == "parent" and _safe_text(track.get("pool")) == pool for track in tracks)
        and any(track.get("role") == "treated" and _safe_text(track.get("pool")) == pool for track in tracks)
    ]


def _matched_pair_tracks(tracks: list[dict[str, Any]], pool: str | None = None) -> list[dict[str, Any]]:
    pools = [pool] if pool else _matched_pair_pools(tracks)
    selected: list[dict[str, Any]] = []
    for current_pool in pools:
        parents = [
            track
            for track in tracks
            if track.get("role") == "parent" and _safe_text(track.get("pool")) == current_pool
        ]
        treated = [
            track
            for track in tracks
            if track.get("role") == "treated" and _safe_text(track.get("pool")) == current_pool
        ]
        selected.extend(_order_tracks_for_display(parents))
        selected.extend(_order_tracks_for_display(treated))
    return selected


def resolve_track_preset(
    project_config: str | Path,
    track_source: str = "raw",
    track_preset: str = "all",
    samples: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Resolve a display preset into existing tracks without renaming samples."""
    preset = (track_preset or "all").lower()
    tracks = list_insertion_tracks(project_config, track_source)
    if preset == "custom":
        selected = [sample for sample in (samples or []) if sample]
        if not selected or "all" in {sample.lower() for sample in selected}:
            return _order_tracks_for_display(tracks)
        selected_tracks: list[dict[str, Any]] = []
        for sample in selected:
            for track in tracks:
                if track["sample"] == sample or track["track_name"] == sample:
                    selected_tracks.append(track)
                    break
        if not selected_tracks:
            raise ValueError("Custom track preset did not match any available tracks.")
        return selected_tracks
    raw_like = [track for track in tracks if track.get("track_source") == "raw"]
    if preset == "all":
        return _order_tracks_for_display(raw_like if raw_like else tracks)
    if preset in {"parents", "treated"}:
        role = "parent" if preset == "parents" else "treated"
        selected_tracks = [track for track in raw_like if track.get("role") == role]
        if not selected_tracks:
            raise ValueError(
                f"Could not resolve {preset} because parent/treated metadata were not available for this project. "
                "Use --samples or --track-preset custom."
            )
        return _order_tracks_for_display(selected_tracks)
    if preset in {"matched", "matched_pairs", "pairs"}:
        selected_tracks = _matched_pair_tracks(raw_like)
        if not selected_tracks:
            raise ValueError(
                "Could not resolve matched pairs because parent/treated pool metadata were not available for this project. "
                "Use --samples or --track-preset custom."
            )
        return selected_tracks
    match = re.fullmatch(r"pool([0-9]+)_pair", preset)
    if match:
        pool = match.group(1)
        selected_tracks = _matched_pair_tracks(raw_like, pool=pool)
        if not selected_tracks:
            raise ValueError(
                f"Could not resolve {preset} because parent/treated pool metadata were not available for this project. "
                "Use --samples or --track-preset custom."
            )
        return selected_tracks
    raise ValueError(f"Unknown track preset: {track_preset}")


def list_insertion_tracks(project_config: str | Path, track_source: str = "raw") -> list[dict[str, Any]]:
    """List available raw and/or combined insertion hit tracks."""
    cfg = _read_project_config(project_config)
    project_root = Path(cfg["_project_root"])
    selected_source = (track_source or "raw").lower()
    tracks: list[dict[str, Any]] = []
    if selected_source in {"raw", "all"}:
        for sample_order, row in enumerate(_sample_rows(cfg)):
            sample = _safe_text(row.get("sample"))
            hit_file = _repo_path(row.get("hit_file"))
            if sample and hit_file and hit_file.is_file():
                tracks.append(
                    {
                        "sample": sample,
                        "track_name": sample,
                        "track_source": "raw",
                        "source_file": str(hit_file),
                        "role": _infer_track_role(row, sample),
                        "pool": _infer_track_pool(row, sample),
                        "sample_order": sample_order,
                    }
                )
        if not tracks:
            for sample_order, hit_file in enumerate(sorted((project_root / "create_hit_file").glob("*/*_hits.txt"))):
                sample = hit_file.parent.name
                row = {"sample": sample}
                tracks.append(
                    {
                        "sample": sample,
                        "track_name": sample,
                        "track_source": "raw",
                        "source_file": str(hit_file),
                        "role": _infer_track_role(row, sample),
                        "pool": _infer_track_pool(row, sample),
                        "sample_order": sample_order,
                    }
                )
    if selected_source in {"combined", "combined_parent", "all"}:
        combined_root = project_root / "combined_hits"
        for hit_file in sorted(combined_root.glob("*/combined_parent_hits*.txt")):
            target = hit_file.parent.name
            tracks.append(
                {
                    "sample": f"combined_parent_{target}",
                    "track_name": f"Combined parent {target}",
                    "track_source": "combined_parent",
                    "source_file": str(hit_file),
                    "role": "parent",
                    "pool": "",
                    "sample_order": 9999,
                }
            )
    return tracks


def _read_hits_for_track(track: dict[str, Any]) -> list[dict[str, Any]]:
    path = Path(track["source_file"])
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = _safe_text(row.get("Chromosome") or row.get("chromosome") or row.get("chrom"))
            pos = _safe_text(row.get("Hit position") or row.get("hit_pos") or row.get("position"))
            count = _safe_text(row.get("Hit count") or row.get("hit_count") or row.get("read_count"))
            if not chrom or not pos:
                continue
            try:
                position = int(float(pos))
                read_count = int(float(count)) if count else 1
            except ValueError:
                continue
            rows.append(
                {
                    "sample": track["sample"],
                    "track_name": track["track_name"],
                    "track_source": track["track_source"],
                    "chromosome": chrom,
                    "position": position,
                    "read_count": read_count,
                    "strand": _safe_text(row.get("Source") or row.get("strand")),
                    "source_file": track["source_file"],
                }
            )
    return rows


def load_insertion_tracks(
    project_config: str | Path,
    samples: list[str] | None = None,
    track_source: str = "raw",
    track_preset: str = "custom",
) -> list[dict[str, Any]]:
    """Load insertion hit records from selected existing project tracks."""
    tracks = resolve_track_preset(project_config, track_source, track_preset, samples)
    loaded: list[dict[str, Any]] = []
    for track in tracks:
        track = dict(track)
        track["records"] = _read_hits_for_track(track)
        loaded.append(track)
    return loaded


def _record_by_gene_id(project_config: str | Path, gene_id: str) -> GeneRecord | None:
    lookup = build_gene_lookup(project_config)
    for record in lookup["records"]:
        if gene_id == record.gene_id or gene_id in record.aliases:
            return record
    return None


def collect_gene_insertions(
    project_config: str | Path,
    gene_id: str,
    samples: list[str] | None = None,
    flank_bp: int = 1000,
    track_source: str = "raw",
    track_preset: str = "custom",
    show_domains: bool = True,
) -> dict[str, Any]:
    """Collect gene model, optional domains, tracks, and overlapping insertions."""
    context = load_project_gene_context(project_config)
    record = _record_by_gene_id(project_config, gene_id)
    if record is None:
        resolved = resolve_gene(project_config, gene_id)
        if resolved.get("status") == "resolved":
            record = _record_by_gene_id(project_config, resolved["gene"]["gene_id"])
    if record is None:
        raise ValueError(f"Gene could not be resolved: {gene_id}")
    flank_bp = max(0, int(flank_bp))
    start = max(1, int(record.start) - flank_bp)
    end = int(record.end) + flank_bp
    gene_start = int(record.start)
    gene_end = int(record.end)
    tracks = load_insertion_tracks(
        project_config,
        samples=samples,
        track_source=track_source,
        track_preset=track_preset,
    )
    selected_rows: list[dict[str, Any]] = []
    track_summaries: list[dict[str, Any]] = []
    for track in tracks:
        overlapping: list[dict[str, Any]] = []
        for row in track["records"]:
            if row["chromosome"] != record.chromosome or not (start <= row["position"] <= end):
                continue
            enriched = dict(row)
            enriched["displayed_region_start"] = start
            enriched["displayed_region_end"] = end
            enriched["gene_start"] = gene_start
            enriched["gene_end"] = gene_end
            enriched["in_display_region"] = True
            enriched["in_gene"] = gene_start <= int(row["position"]) <= gene_end
            if enriched["in_gene"]:
                enriched["region_class"] = "gene_body"
            elif int(row["position"]) < gene_start:
                enriched["region_class"] = "left_flank"
            else:
                enriched["region_class"] = "right_flank"
            overlapping.append(enriched)
        gene_body = [row for row in overlapping if row["in_gene"]]
        flank_context = [row for row in overlapping if not row["in_gene"]]
        selected_rows.extend(overlapping)
        track_summaries.append(
            {
                "sample": track["sample"],
                "track_name": track["track_name"],
                "track_source": track["track_source"],
                "source_file": track["source_file"],
                "role": track.get("role", ""),
                "pool": track.get("pool", ""),
                "sample_order": track.get("sample_order", 9999),
                "insertions_in_region": len(overlapping),
                "insertions_in_gene": len(gene_body),
                "insertions_in_flank": len(flank_context),
                "reads_in_region": int(sum(row["read_count"] for row in overlapping)),
                "reads_in_gene": int(sum(row["read_count"] for row in gene_body)),
            }
        )
    domains = list(record.domains) if show_domains else []
    return {
        "project_id": context["project_id"],
        "project_config": context["project_config"],
        "project_root": context["project_root"],
        "species": context["species"],
        "annotation_source": context["annotation_source"],
        "gene": record.as_dict(),
        "region": {"chromosome": record.chromosome, "start": start, "end": end, "flank_bp": flank_bp},
        "count_region": {"start": gene_start, "end": gene_end},
        "track_source": track_source,
        "track_preset": track_preset,
        "selected_samples": samples or ["all"],
        "tracks": track_summaries,
        "insertions": selected_rows,
        "domains": [list(x) for x in domains],
        "domain_source": record.domain_source if domains else "",
        "domains_available": bool(domains),
        "domain_message": "Domain annotations available." if domains else "No domain annotation available for this gene.",
    }


def _format_track_label(track_name: str, site_count: int, label_mode: str = "full", show_site_counts: bool = True) -> str:
    """Return a display-only track label without changing sample IDs."""
    label = str(track_name)
    if label_mode == "compact" and len(label) > 34:
        label = label[:16] + "…" + label[-15:]
    if show_site_counts:
        label = f"{label}  n={site_count}"
    return label


def _gene_title(gene: dict[str, Any]) -> str:
    display = str(gene.get("display_name") or "").strip()
    gene_id = str(gene.get("gene_id") or "").strip()
    cagl_id = str(gene.get("standard_name") or "").strip()
    if cagl_id and not cagl_id.upper().startswith("CAGL"):
        cagl_id = ""
    if display and cagl_id and display != cagl_id:
        return f"{display} ({cagl_id})"
    if display and gene_id and display != gene_id:
        return f"{display} ({gene_id})"
    return display or gene_id or "Selected gene"


def draw_gene_domain_insertion_figure(
    payload: dict[str, Any],
    output_png: str | Path,
    width_px: int = 1800,
    dpi: int = 150,
    label_mode: str = "full",
    show_site_counts: bool = True,
) -> str:
    """Draw a stacked insertion/gene/domain PNG with matplotlib Agg."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    from matplotlib.patches import FancyArrow, Rectangle

    output_png = Path(output_png)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    gene = payload["gene"]
    region = payload["region"]
    tracks = payload["tracks"]
    insertions = payload["insertions"]
    width_px = max(900, int(width_px))
    dpi = max(96, int(dpi))
    label_mode = label_mode if label_mode in {"full", "compact"} else "full"
    n_tracks = max(1, len(tracks))
    fig_width = width_px / float(dpi)
    fig_height = max(5.5, 1.0 + 0.45 * n_tracks + 1.8)
    track_ratio = max(2.0, n_tracks * 0.45)
    fig, (ax_tracks, ax_gene) = plt.subplots(
        2,
        1,
        figsize=(fig_width, fig_height),
        dpi=dpi,
        sharex=True,
        gridspec_kw={"height_ratios": [track_ratio, 1.2], "hspace": 0.12},
    )

    region_start = int(region["start"])
    region_end = int(region["end"])
    gene_start = int(gene["start"])
    gene_end = int(gene["end"])
    ax_tracks.set_xlim(region_start, region_end)
    ax_tracks.set_ylim(0.4, n_tracks + 0.6)
    ax_gene.set_ylim(0, 1)

    strand = str(gene.get("strand") or "unknown")
    subtitle = (
        f"{region['chromosome']}:{region_start:,}–{region_end:,}"
        f" | strand {strand}"
        f" | flank {int(region.get('flank_bp', 0)):,} bp"
    )
    fig.suptitle(_gene_title(gene), fontsize=15, fontweight="bold", y=0.985)
    fig.text(0.5, 0.945, subtitle, ha="center", va="top", fontsize=10.5, color="#555555")

    for ax in (ax_tracks, ax_gene):
        ax.set_facecolor("white")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#777777")
        ax.spines["bottom"].set_color("#777777")
        ax.spines["left"].set_linewidth(1.0)
        ax.spines["bottom"].set_linewidth(1.0)
        ax.grid(axis="x", color="#d0d0d0", alpha=0.25, linewidth=0.9)
        ax.tick_params(axis="x", labelsize=9.5, width=1.0)
    ax_tracks.tick_params(axis="y", labelsize=8.5, length=0)
    ax_gene.tick_params(axis="y", left=False, labelleft=False)

    by_track: dict[str, list[dict[str, Any]]] = {}
    for row in insertions:
        by_track.setdefault(row["sample"], []).append(row)
    track_labels: list[str] = []
    y_positions: list[float] = []
    for i, track in enumerate(tracks):
        y = n_tracks - i
        y_positions.append(y)
        ax_tracks.hlines(y, region_start, region_end, color="#d8d8d8", linewidth=1.0)
        rows = by_track.get(track["sample"], [])
        track_labels.append(
            _format_track_label(
                track["track_name"],
                int(track.get("insertions_in_gene", 0)),
                label_mode,
                show_site_counts,
            )
        )
        for row in rows:
            read_count = max(0, int(row.get("read_count") or 0))
            linewidth = min(3.0, max(1.0, math.log10(read_count + 1))) if read_count else 1.2
            in_gene = bool(row.get("in_gene"))
            ax_tracks.vlines(
                int(row["position"]),
                y - 0.25,
                y + 0.25,
                color="#303030" if in_gene else "#9a9a9a",
                linewidth=linewidth,
                alpha=0.95 if in_gene else 0.45,
            )
    ax_tracks.set_yticks(y_positions)
    ax_tracks.set_yticklabels(track_labels)
    ax_tracks.set_ylabel("Insertion tracks", fontsize=11, fontweight="bold")

    gene_y = 0.40
    gene_h = 0.25
    gene_mid = gene_y + gene_h / 2
    ax_gene.hlines(gene_mid, gene_start, gene_end, color="#1d1d1d", linewidth=1.0, alpha=0.65)
    ax_gene.add_patch(
        Rectangle(
            (gene_start, gene_y),
            gene_end - gene_start,
            gene_h,
            facecolor="#1f77b4",
            edgecolor="#1d1d1d",
            linewidth=1.3,
        )
    )
    for exon_start, exon_end in gene.get("exons") or []:
        ax_gene.add_patch(
            Rectangle(
                (int(exon_start), gene_y - 0.06),
                int(exon_end) - int(exon_start),
                gene_h + 0.12,
                facecolor="#0b5ea8",
                edgecolor="#1d1d1d",
                linewidth=1.0,
            )
        )

    if payload.get("domains_available"):
        for dom_start, dom_end in payload.get("domains", []):
            ax_gene.add_patch(
                Rectangle(
                    (int(dom_start), gene_y + 0.055),
                    int(dom_end) - int(dom_start),
                    gene_h - 0.11,
                    facecolor="#111111",
                    edgecolor="#111111",
                    linewidth=0.9,
                )
            )
        domain_caption = f"Domains shown from: {payload.get('domain_source') or 'annotation'}"
    else:
        domain_caption = "No domain annotation available for this gene."

    if payload.get("show_direction", True):
        direction = 1 if gene.get("strand") != "-" else -1
        arrow_len = max(1, min((gene_end - gene_start) * 0.22, (region_end - region_start) * 0.08))
        arrow_start = gene_end - arrow_len if direction == 1 else gene_start + arrow_len
        ax_gene.add_patch(
            FancyArrow(
                arrow_start,
                gene_mid,
                arrow_len * direction,
                0,
                width=0.035,
                head_width=0.12,
                head_length=max(1, arrow_len * 0.28),
                color="#222222",
                length_includes_head=True,
            )
        )

    ax_gene.text(
        (gene_start + gene_end) / 2,
        gene_y + gene_h + 0.12,
        _gene_title(gene),
        fontsize=10.5,
        fontweight="bold",
        ha="center",
        va="bottom",
        color="#1f1f1f",
    )
    ax_gene.text(
        0.01,
        -0.42,
        domain_caption,
        transform=ax_gene.transAxes,
        fontsize=9.5,
        color="#666666",
        ha="left",
        va="top",
        clip_on=False,
    )
    ax_gene.set_xlabel(f"{region['chromosome']} coordinate (bp)", fontsize=11, fontweight="bold")
    ax_gene.xaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f"{int(value):,}"))
    fig.subplots_adjust(left=0.30, right=0.98, top=0.90, bottom=0.16)
    fig.savefig(output_png, bbox_inches="tight")
    plt.close(fig)
    return str(output_png)


def _safe_slug(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return value.strip("._") or "gene"


def _hash_params(data: dict[str, Any]) -> str:
    payload = json.dumps(data, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


def _write_insertions_table(rows: list[dict[str, Any]], table_path: Path) -> None:
    table_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample",
        "track_name",
        "track_source",
        "chromosome",
        "position",
        "read_count",
        "strand",
        "source_file",
        "displayed_region_start",
        "displayed_region_end",
        "gene_start",
        "gene_end",
        "in_display_region",
        "in_gene",
        "region_class",
    ]
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def run_gene_domain_explorer(
    project_config: str | Path,
    query: str,
    samples: list[str] | None = None,
    track_source: str = "raw",
    track_preset: str = "all",
    flank_bp: int = 1000,
    show_domains: bool = True,
    show_direction: bool = True,
    width_px: int = 1800,
    dpi: int = 150,
    label_mode: str = "full",
    show_site_counts: bool = True,
    output: str | Path | None = None,
    force: bool = False,
) -> dict[str, Any]:
    """Resolve a gene, generate/reuse PNG/table/manifest, and return paths."""
    resolved = resolve_gene(project_config, query)
    if resolved["status"] != "resolved":
        return resolved
    gene_id = resolved["gene"]["gene_id"]
    payload = collect_gene_insertions(
        project_config,
        gene_id,
        samples=samples,
        flank_bp=flank_bp,
        track_source=track_source,
        track_preset=track_preset,
        show_domains=show_domains,
    )
    payload["query"] = query
    payload["show_direction"] = bool(show_direction)
    requested_samples = samples or ["all"]
    selected_samples = [track["sample"] for track in payload["tracks"]]
    cache_params = {
        "project_config": str(_repo_path(project_config)),
        "query": query,
        "gene_id": gene_id,
        "samples": requested_samples,
        "track_source": track_source,
        "track_preset": track_preset,
        "flank_bp": int(flank_bp),
        "show_domains": bool(show_domains),
        "show_direction": bool(show_direction),
        "width_px": int(width_px),
        "dpi": int(dpi),
        "label_mode": label_mode,
        "show_site_counts": bool(show_site_counts),
    }
    key = _hash_params(cache_params)
    project_root = Path(payload["project_root"])
    base = project_root / "gene_domain_explorer"
    slug = _safe_slug(gene_id)
    figure_path = Path(output) if output else base / "figures" / f"{slug}.{track_source}.{key}.png"
    if not figure_path.is_absolute():
        figure_path = (REPO_ROOT / figure_path).resolve()
    table_path = base / "tables" / f"{slug}.{track_source}.{key}.insertions.csv"
    manifest_path = base / "manifests" / f"{slug}.{track_source}.{key}.json"
    for directory in (figure_path.parent, table_path.parent, manifest_path.parent):
        directory.mkdir(parents=True, exist_ok=True)
    cached = figure_path.is_file() and table_path.is_file() and manifest_path.is_file() and not force
    if not cached:
        draw_gene_domain_insertion_figure(
            payload,
            figure_path,
            width_px=width_px,
            dpi=dpi,
            label_mode=label_mode,
            show_site_counts=show_site_counts,
        )
        _write_insertions_table(payload["insertions"], table_path)
    manifest = {
        "project_id": payload["project_id"],
        "project_config": str(_repo_path(project_config)),
        "species": payload["species"],
        "query": query,
        "resolved_gene_id": gene_id,
        "resolved_gene_name": payload["gene"]["display_name"],
        "chromosome": payload["gene"]["chromosome"],
        "start": payload["gene"]["start"],
        "end": payload["gene"]["end"],
        "flank_bp": int(flank_bp),
        "track_source": track_source,
        "track_preset": track_preset,
        "selected_samples": selected_samples,
        "requested_samples": requested_samples,
        "displayed_region_start": payload["region"]["start"],
        "displayed_region_end": payload["region"]["end"],
        "gene_start": payload["gene"]["start"],
        "gene_end": payload["gene"]["end"],
        "count_region": payload["count_region"],
        "count_definition": "insertion sites inside gene_start..gene_end only; flank insertions are displayed but not counted",
        "domain_source": payload["domain_source"],
        "domains_available": bool(payload["domains_available"]),
        "figure_path": str(figure_path),
        "table_path": str(table_path),
        "manifest_path": str(manifest_path),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "status": "cached" if cached else "success",
        "domain_message": payload["domain_message"],
        "track_count": len(payload["tracks"]),
        "insertion_count": len(payload["insertions"]),
        "track_summaries": payload["tracks"],
        "width_px": int(width_px),
        "dpi": int(dpi),
        "label_mode": label_mode,
        "show_site_counts": bool(show_site_counts),
    }
    if not cached:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    else:
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["status"] = "cached"
        except Exception:
            pass
    return {"status": manifest.get("status", "success"), "gene": payload["gene"], "payload": payload, "manifest": manifest}
