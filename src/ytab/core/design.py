"""Reading a YTAB `normalize` run into a validated treated-vs-control design.

`load_design` is the single entry point every comparative method uses, so they all read normalized input,
check the catalogue against its checksums and select observations identically. `load_normalization` locates one
method inside a normalization run. The normalization name tuples say which methods are expressible as scalar
count-model exposures (NORMALIZATIONS) and which reshape the counts instead (NONLINEAR_NORMALIZATIONS)."""

import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from ytab.normalization.ttr import sha256

NORMALIZATIONS = ("ttr", "totreads", "nzmean")   # scalar methods whose factors are valid count-model exposures
DEFAULT_NORMALIZATION = "ttr"

NONLINEAR_NORMALIZATIONS = ("betageom", "aBGC")   # resampling only: count-reshaping, no scalar factor


def load_normalization(normalized, normalization=None, nonlinear=False):
    """Locate one scalar method inside a YTAB `normalize` output (run directory or its method folder).

    nonlinear=True (resampling only) also accepts the beta-geometric corrections, which reshape counts.
    """
    normalized = Path(normalized).resolve()
    if (normalized/"method.json").is_file():
        folder, root = normalized, normalized.parent
    else:
        root = normalized
        folder = root/(normalization or DEFAULT_NORMALIZATION)
    if not (root/"run_manifest.json").is_file() or not (folder/"method.json").is_file():
        raise ValueError(f"Expected a YTAB normalization output with {folder.name}/method.json and run_manifest.json: {normalized}")
    method = json.loads((folder/"method.json").read_text())["id"]
    if normalization and method != normalization:
        raise ValueError(f"Normalized input is {method!r}, not the requested {normalization!r}")
    if method not in NORMALIZATIONS+(NONLINEAR_NORMALIZATIONS if nonlinear else ()):
        raise ValueError(f"{method!r} is not supported here; use one of {', '.join(NORMALIZATIONS)} "
                         "(scalar factors only: distribution-reshaping methods cannot be expressed as count-model exposures)")
    manifest = json.loads((root/"run_manifest.json").read_text())
    if manifest.get("status") != "success" or method not in manifest.get("methods", []):
        raise ValueError(f"Normalization run is incomplete or does not list {method!r}")
    with (root/"observations.tsv").open() as handle:observations = [r["observation_id"] for r in csv.DictReader(handle, delimiter="\t")]
    with (folder/"factors.tsv").open() as handle:factors = list(csv.DictReader(handle, delimiter="\t"))
    if [r["observation_id"] for r in factors] != observations or len(set(observations)) != len(observations):
        raise ValueError("Normalization factors and observations disagree")
    return dict(method=method, root=root, folder=folder, manifest=manifest, observations=observations,
                factors={r["observation_id"]:float(r["factor"]) if r["factor"] != "" else None for r in factors},
                groups={r["observation_id"]:r["normalization_group"] for r in factors})

def load_design(normalized, normalization, reference, output, controls, treatments, min_controls=2, min_treatments=1, stratified=False,
                nonlinear=False, check_roles=True, independent=False):
    """Validated controls/treatments from a YTAB scalar normalization output and its source catalogue.

    stratified=True (resampling only) allows several independent libraries, each with >=1 control and
    >=1 treatment sharing one normalization group; comparisons are then made within library.

    Shared by the control-conditioned Gumbel and resampling so both read normalized input identically.
    Returns raw integer counts and normalized counts (selected columns), factors, sites and provenance.
    nonlinear=True (resampling only) accepts betageom/aBGC output: factors are None, and the normalized
    counts are checked to keep every raw zero and the within-sample order of raw counts.
    check_roles=False (multi-condition tests only): conditions are user-assigned, so catalogue roles are not compared.
    independent=True (unstratified resampling only): the samples may come from different, independently made libraries
    (e.g. one library per strain, as in SATAY mutant screens); labels are permuted among samples within each site, as in
    TRANSIT's resampling. Differences in library composition then add to the noise, so the test is for library-level designs.
    """
    norm = load_normalization(normalized, normalization, nonlinear)
    catalogue = Path(norm["manifest"].get("catalogue", "")).resolve()
    reference, output = Path(reference).resolve(), Path(output).resolve()
    if any(output == x or x in output.parents or output in x.parents for x in (catalogue, norm["root"])) or (output.exists() and any(output.iterdir())):
        raise ValueError("Choose a new output directory, separate from the normalization and its catalogue")
    controls = list(controls); treatments = [treatments] if isinstance(treatments, str) else list(treatments)
    selected = controls+treatments
    if len(controls) < min_controls or len(treatments) < min_treatments or len(set(selected)) != len(selected):
        raise ValueError(f"Select >={min_controls} unique control observations and distinct treatment observations explicitly")
    if not (catalogue/"run_manifest.json").is_file():
        raise ValueError(f"Normalization source catalogue not found: {catalogue}")
    stale = [name for name, digest in norm["manifest"].get("source_hashes", {}).items() if sha256(catalogue/name) != digest]
    if stale or not norm["manifest"].get("source_hashes"):
        raise ValueError(f"Catalogue changed since normalization (checksum mismatch: {stale}); renormalize")
    source = json.loads((catalogue/"run_manifest.json").read_text())
    if source.get("status") != "success" or sha256(reference) != source.get("reference_sha256") or norm["manifest"].get("reference_sha256") != source.get("reference_sha256"):
        raise ValueError("Successful catalogue and checksum-identical reference required")
    with (catalogue/"observations.tsv").open() as handle:metadata = list(csv.DictReader(handle, delimiter="\t"))
    ids = [r["observation_id"] for r in metadata]
    if len(set(ids)) != len(ids):raise ValueError("Duplicate observation IDs in catalogue")
    missing = [x for x in selected if x not in ids or x not in norm["factors"]]
    if missing:raise ValueError(f"Observations missing from the catalogue or normalization: {missing}")
    chosen = [metadata[ids.index(x)] for x in selected]
    for field in ("project",) if stratified or independent else ("project", "library_id", "pool"):
        if len({r.get(field, "") for r in chosen}) != 1 or not chosen[0].get(field):
            raise ValueError(f"Controls and treatments must share one explicit {field}")
    if stratified:
        libraries = {r.get("library_id", "") for r in chosen}
        for library in libraries:
            members = [x for x, r in zip(selected, chosen) if r.get("library_id", "") == library]
            if not library or not set(members) & set(controls) or not set(members) & set(treatments):
                raise ValueError(f"Each library needs an explicit library_id with >=1 control and >=1 treatment: {library!r}")
            if len({norm["groups"][x] for x in members}) != 1:
                raise ValueError(f"Samples of library {library!r} must come from one normalization group")
    elif len({norm["groups"][x] for x in selected}) != 1:
        raise ValueError("Controls and treatments must come from one normalization group")
    control_roles = {r.get("role", "") for r in chosen[:len(controls)]}
    treatment_roles = {r.get("role", "") for r in chosen[len(controls):]}
    if check_roles and (len(control_roles) != 1 or len(treatment_roles) != 1 or control_roles == treatment_roles):
        raise ValueError("Controls must share one role, treatments another, and the roles must differ")
    raw = sparse.load_npz(catalogue/"raw_counts.npz").tocsc(); raw.sum_duplicates()
    if (raw.shape[1] != len(ids) or raw.shape[0] != source.get("universe_sites") or np.any(~np.isfinite(raw.data))
            or np.any(raw.data < 0) or np.any(raw.data != np.floor(raw.data))):
        raise ValueError("Raw integer matrix dimensions/counts disagree with catalogue")
    with gzip.open(catalogue/"sites.tsv.gz", "rt") as handle:
        sites = [(r["chromosome"], int(r["position_1based"])) for r in csv.DictReader(handle, delimiter="\t")]
    if len(sites) != raw.shape[0] or len(set(sites)) != len(sites) or any(p < 1 for _, p in sites):
        raise ValueError("Coordinate catalogue duplicates or dimensions invalid")
    counts = raw[:, [ids.index(x) for x in selected]].toarray().astype(np.int64)
    scaled = sparse.load_npz(norm["folder"]/"normalized_counts.npz").tocsc()
    if scaled.shape[0] != raw.shape[0]:raise ValueError("Normalized matrix and catalogue universe disagree")
    scaled = scaled[:, [norm["observations"].index(x) for x in selected]].toarray()
    if norm["method"] in NONLINEAR_NORMALIZATIONS:
        factors = None
        for j in range(counts.shape[1]):
            order = np.lexsort((scaled[:, j], counts[:, j]))
            if (np.any((scaled[:, j] > 0) != (counts[:, j] > 0)) or np.any(~np.isfinite(scaled[:, j]))
                    or np.any(np.diff(scaled[order, j]) < -1e-9*max(1., scaled[:, j].max()))):
                raise ValueError("Normalized counts do not preserve the raw zeros and count order of this catalogue")
    else:
        factors = np.array([norm["factors"][x] for x in selected])
        if np.any(~np.isfinite(factors)) or np.any(factors <= 0):
            raise ValueError("Normalization factors must be finite and positive")
        if not np.allclose(scaled, counts*factors, rtol=1e-9, atol=1e-9):
            raise ValueError("Normalized counts are not raw catalogue counts times the reported factors; input is not a valid scalar normalization")
    return dict(norm=norm, catalogue=catalogue, output=output, controls=controls, treatments=treatments, selected=selected,
                libraries=[r.get("library_id", "") for r in chosen],
                chosen=chosen, counts=counts, scaled=scaled, factors=factors, sites=sites)
