"""Discover species reference resources and report preparation readiness."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass
class ReferenceInfo:
    species: str
    species_dir: Path
    reference_dir: Path
    fasta: Path | None = None
    feature_table: Path | None = None
    gff: Path | None = None
    gtf: Path | None = None
    gffutils_db: Path | None = None
    bowtie2_index_prefix: Path | None = None
    bowtie2_index_complete: bool = False
    bowtie2_index_type: str | None = None
    centromere_bed: Path | None = None
    trna_bed: Path | None = None
    orthology_file: Path | None = None
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {key: str(value) if isinstance(value, Path) else value for key, value in asdict(self).items()}


def _repo_root(repo_root: Path | None) -> Path:
    return (repo_root or Path(__file__).resolve().parents[3]).expanduser().resolve()


def list_available_species(resources_dir: Path) -> list[str]:
    resources_dir = Path(resources_dir)
    if not resources_dir.is_dir():
        return []
    reference_patterns = (
        "*.fna", "*.fa", "*.fasta", "*.gff", "*.gff3", "*.gtf",
        "*feature_table*.txt", "*chromosomal_feature*.tab",
    )
    available = []
    for species_dir in resources_dir.iterdir():
        reference_dir = species_dir / "reference_genome"
        if (species_dir.is_dir() and not species_dir.name.startswith(".")
                and reference_dir.is_dir()
                and any(_glob_many(reference_dir, reference_patterns))):
            available.append(species_dir.name)
    return sorted(available)


def reference_is_selectable(info: ReferenceInfo) -> bool:
    return info.reference_dir.is_dir() and bool(
        info.fasta or info.gff or info.gtf or info.feature_table
    )


def reference_can_be_prepared(info: ReferenceInfo) -> bool:
    return bool(info.fasta and (info.gff or info.gtf or info.feature_table))


def reference_is_runnable(info: ReferenceInfo) -> bool:
    return reference_can_be_prepared(info) and info.bowtie2_index_complete


def find_bowtie2_index_prefix(reference_dir: Path) -> tuple[Path | None, bool, str | None]:
    reference_dir = Path(reference_dir)
    for extension in ("bt2", "bt2l"):
        for first in sorted(reference_dir.glob(f"*.1.{extension}")):
            prefix = Path(str(first)[: -len(f".1.{extension}")])
            required = [
                Path(f"{prefix}.1.{extension}"), Path(f"{prefix}.2.{extension}"),
                Path(f"{prefix}.3.{extension}"), Path(f"{prefix}.4.{extension}"),
                Path(f"{prefix}.rev.1.{extension}"), Path(f"{prefix}.rev.2.{extension}"),
            ]
            if all(path.is_file() for path in required):
                return prefix, True, extension
    return None, False, None


def _select(candidates: list[Path], label: str, species: str, warnings: list[str]) -> Path | None:
    unique = sorted(set(candidates), key=lambda path: (species.lower() not in path.name.lower(), len(path.name), path.name.lower()))
    if not unique:
        return None
    if len(unique) > 1:
        warnings.append(f"Multiple {label} candidates found; selected {unique[0].name}.")
    return unique[0]


def _select_centromere_bed(reference_dir: Path, species: str, warnings: list[str]) -> Path | None:
    candidates = sorted(set(_glob_many(reference_dir, ("centromeres*.bed",))), key=lambda path: path.name.lower())
    if not candidates:
        return None
    preferred = [path for path in candidates if "ncbi" in path.name.lower()]
    selected = (preferred or candidates)[0]
    if len(candidates) > 1:
        warnings.append(f"Multiple centromere BED candidates found; selected {selected.name}.")
    return selected


def _glob_many(directory: Path, patterns: tuple[str, ...]) -> list[Path]:
    return [path for pattern in patterns for path in directory.glob(pattern) if path.is_file()]


def resolve_reference(species: str, repo_root: Path | None = None) -> ReferenceInfo:
    root = _repo_root(repo_root)
    species_dir = root / "resources" / "species" / species
    if not species_dir.is_dir():
        available = ", ".join(list_available_species(root / "resources" / "species")) or "none"
        raise ValueError(f"Unknown species '{species}'. Available species: {available}")
    reference_dir = species_dir / "reference_genome"
    warnings: list[str] = []
    if not reference_dir.is_dir():
        warnings.append(f"Reference genome directory not found: {reference_dir}")

    fasta = _select(_glob_many(reference_dir, ("*.fna", "*.fa", "*.fasta")), "FASTA", species, warnings)
    gff = _select(_glob_many(reference_dir, ("*.gff", "*.gff3")), "GFF", species, warnings)
    gtf = _select(_glob_many(reference_dir, ("*.gtf",)), "GTF", species, warnings)
    feature = _select(_glob_many(reference_dir, ("*feature_table*.txt", "*chromosomal_feature*.tab")), "feature table", species, warnings)
    database = _select(_glob_many(reference_dir, ("*.gffutils_db.sqlite",)), "gffutils DB", species, warnings)
    centromere = _select_centromere_bed(reference_dir, species, warnings)
    trna = _select(_glob_many(reference_dir, ("*tRNA*.bed", "tRNAs.bed")), "tRNA BED", species, warnings)
    orthology = _select(_glob_many(species_dir, ("*orthologs*.txt",)), "orthology file", species, warnings)
    prefix, complete, index_type = find_bowtie2_index_prefix(reference_dir)

    if fasta is None:
        warnings.append("Reference FASTA missing. Cannot build Bowtie2 index.")
    if feature is None and gff is None and gtf is None:
        warnings.append("No GFF/GTF/feature table found. Feature annotation will not be available.")
    if complete:
        warnings.append("Bowtie2 index found; indexing skipped.")
    elif fasta is not None:
        warnings.append("No complete Bowtie2 index found. Run reference preparation before alignment.")

    return ReferenceInfo(
        species=species, species_dir=species_dir, reference_dir=reference_dir,
        fasta=fasta, feature_table=feature, gff=gff, gtf=gtf,
        gffutils_db=database, bowtie2_index_prefix=prefix,
        bowtie2_index_complete=complete, bowtie2_index_type=index_type,
        centromere_bed=centromere, trna_bed=trna,
        orthology_file=orthology, warnings=warnings,
    )
