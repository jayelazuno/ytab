"""Project-scoped representation inventory for the analytical planner.

Uses the existing raw-coordinate/CDS adapter, never the legacy cross-project audit.
No simulations, classifier filtering or changes to Poisson equations are made.
"""
import csv
import hashlib
import json
import tempfile
from pathlib import Path

from ytab.pipeline.project_config import load_project_config
from .library import is_parent_sample, load_library


def nonprotein_loci(gff):
    """Honor explicit gene biotypes/pseudogenes, as in the accepted audit."""
    excluded = set()
    with Path(gff).open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "pseudogene"}:
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            biotype = attrs.get("gene_biotype", attrs.get("biotype", ""))
            if (fields[2] == "pseudogene" or attrs.get("pseudo", "").lower() == "true" or
                    (biotype and biotype != "protein_coding")):
                excluded.add(attrs.get("locus_tag") or attrs.get("gene_id"))
    return excluded


def parent_rows(rows):
    included = [r for r in rows if str(r.get("include", "true")).lower() in {"true", "yes", "1"}]
    eligible = [r for r in included if is_parent_sample(r)]
    def explicit_parent(row):
        for key in ("library_role", "control_or_treated", "condition", "guessed_condition"):
            value = str(row.get(key, "")).strip().lower()
            if value:
                return value in {"parent", "parental", "parent library"}
        return False
    parents = [r for r in eligible if explicit_parent(r)]
    return parents or eligible


def project_summary(project_config):
    config_path = Path(project_config).resolve()
    cfg = load_project_config(config_path)
    root = Path(cfg["_repo_root"])
    resolve = lambda value: Path(value) if Path(value).is_absolute() else root / value
    project = resolve(cfg["output_project_dir"])
    sheet = resolve(cfg["sample_sheet"])
    with sheet.open() as handle:
        rows = parent_rows(list(csv.DictReader(handle)))
    if not rows:
        raise ValueError("No included parent/control samples are defined for this project.")
    resolved = project / "config/reference_resolved.json"
    reference = json.loads(resolved.read_text()) if resolved.exists() else cfg["reference"]
    inputs = [config_path, sheet, resolve(reference["gff"]), Path(__file__), Path(__file__).with_name("library.py")]
    if resolved.exists():
        inputs.append(resolved)
    for row in rows:
        hit = resolve(row.get("hit_file") or ".")
        if not hit.is_file():
            hit = project / "create_hit_file" / row["sample"] / (row["sample"] + "_hits.txt")
        if not hit.is_file():
            raise ValueError(f"Parent insertion catalogue is missing for {row['sample']}; run CreateHitFile for this project.")
        inputs.append(hit)
    # Identity and file revisions are project-scoped, including annotation/code.
    signature = hashlib.sha256(json.dumps([
        (str(p.resolve()), p.stat().st_size, p.stat().st_mtime_ns) for p in inputs
    ]).encode()).hexdigest()
    cache = project / "screen_coverage" / "project_summary" / (signature + ".json")
    if cache.is_file():
        return json.loads(cache.read_text())
    library = load_library(config_path, [r["sample"] for r in rows], universe="all")
    excluded = nonprotein_loci(resolve(reference["gff"]))
    counts = [g["informative_sites"] for g in library["genes"] if g["gene_id"] not in excluded]
    out = dict(project_config=str(config_path), project_id=cfg["project_id"],
        samples=[r["sample"] for r in rows], annotated=len(counts),
        represented=sum(n > 0 for n in counts), ge2=sum(n >= 2 for n in counts),
        ge3=sum(n >= 3 for n in counts), unique_sites=library["unique_sites"],
        coding_sites=sum(counts), sources=library["sources"],
        definition=library["informative_definition"], signature=signature)
    cache.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=cache.parent, suffix=".tmp", delete=False) as handle:
        json.dump(out, handle, allow_nan=False)
        temporary = Path(handle.name)
    temporary.replace(cache)
    return out
