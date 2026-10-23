"""Control-conditioned Hermes Gumbel for independent pools (biological replicates).

Each pool is an independent library with its own clones, so sites are never matched across pools.
Within each pool the mock and treated cultures share clones: the pool's noise curve
(Var = mu(1+a) + phi mu^2, lognormal abundance prior) is fitted on its mock/treated pair with the
selection-aware marginal likelihood (treatment effects affect a minority of sites; any residual
effect inflates the fitted noise, which is conservative). A site is flagged when the treated count is
improbable given the pool's mock count; each gene's longest run of flagged sites is scored per pool
by the exact conditional tail. Replication across pools is tested with a partial-conjunction
hypothesis (Benjamini & Heller 2008): H0 = fewer than u of the gene's n testable pools have a run
effect; p is Fisher's combination of the n-u+1 largest pool p-values. BH controls the FDR over genes.
Input is a YTAB scalar normalization (ttr, totreads or nzmean); no normalization is computed here.
"""
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2

from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.normalization.ttr import sha256, write_table

MODEL = "Hermes control-conditioned Gumbel for independent pools"
REVISION = "pooled_ccg_v1"


@dataclass(frozen=True)
class Settings:
    direction: str = "depleted"         # depleted | enriched
    min_pools: int | None = None        # u; None = majority of the gene's testable pools, floor(n/2)+1
    min_testable_pools: int = 2         # genes testable in fewer pools cannot show replication
    alpha_site: float = .05
    min_control_reads: int = 1
    min_sites: int = 3                  # per pool
    nodes: int = 32
    fdr: float = .05
    seed: int = 2026
    nterm: float = 0.
    cterm: float = 0.
    ignore_stop: bool = True
    dispersion_model: str = "linear-quadratic"

    def validate(self):
        if (self.direction not in {"depleted", "enriched"} or (self.min_pools is not None and (not isinstance(self.min_pools, int) or self.min_pools < 1)) or
                not isinstance(self.min_testable_pools, int) or self.min_testable_pools < 1 or not 0 < self.alpha_site < .5 or
                not isinstance(self.min_sites, int) or self.min_sites < 1 or not 0 < self.fdr < 1 or
                self.dispersion_model not in {"linear-quadratic", "quadratic"} or not 0 <= self.nterm < 100 or not 0 <= self.cterm < 100):
            raise ValueError("Invalid direction, min_pools, min_testable_pools, alpha_site, min_sites, FDR, dispersion model or terminal percentages")

    def ccg(self):
        return ccg.Settings(alpha_site=self.alpha_site, min_control_reads=self.min_control_reads, min_sites=self.min_sites, nodes=self.nodes,
                            seed=self.seed, nterm=self.nterm, cterm=self.cterm, ignore_stop=self.ignore_stop, dispersion_model=self.dispersion_model)


def required_pools(n, settings):
    return min(n, settings.min_pools) if settings.min_pools else n//2+1


def partial_conjunction(pvalues, u):
    """Fisher partial-conjunction p for 'at least u of n effects' from independent p-values."""
    p = np.sort(np.clip(np.asarray(pvalues, dtype=float), 1e-300, 1))[u-1:]
    return float(chi2.sf(-2*np.log(p).sum(), 2*len(p)))


def enrichment_thresholds(locations, log_weights, treatment_exposure, dispersion, alpha):
    """Smallest t with P(T >= t) <= alpha, and the null flag probability P(T >= t*)."""
    n = len(locations); args = (treatment_exposure, dispersion)
    upper = lambda t, u, w:1-ccg.predictive_cdf(t-1, u, w, *args)
    lo = np.zeros(n, dtype=np.int64); hi = np.maximum(1, np.ceil(2*(np.exp(log_weights)*treatment_exposure*np.exp(locations)).sum(1))).astype(np.int64)
    while True:
        high = upper(hi, locations, log_weights) > alpha
        if not high.any():break
        lo[high] = hi[high]; hi[high] *= 2
    while np.any(hi-lo > 1):
        active = np.flatnonzero(hi-lo > 1); middle = (lo[active]+hi[active])//2
        high = upper(middle, locations[active], log_weights[active]) > alpha
        lo[active[high]] = middle[high]; hi[active[~high]] = middle[~high]
    return hi, upper(hi, locations, log_weights)


def pool_model(mock, treated, exposure_mock, exposure_treated, settings):
    """Fit one pool's noise curve on its mock/treated pair; per-site thresholds and null flag probabilities."""
    base = settings.ccg(); pair = np.column_stack([mock, treated]); e = np.array([exposure_mock, exposure_treated])
    supported = pair.sum(1) >= settings.min_control_reads
    fit = ccg.fit_background(pair[supported], e, base)
    dispersion = (fit["phi"], fit["linear"])
    eligible = mock >= settings.min_control_reads
    threshold = np.full(len(mock), -1, dtype=np.int64); probability = np.zeros(len(mock))
    index = np.flatnonzero(eligible)
    if len(index):
        profiles, inverse = np.unique(mock[index], return_inverse=True); inverse = np.ravel(inverse)
        u, w, _ = ccg.site_posterior(profiles[:, None], np.array([exposure_mock]), fit["m"], fit["tau"], dispersion, settings.nodes)
        if settings.direction == "depleted":t, p = ccg.depletion_thresholds(u, w, exposure_treated, dispersion, settings.alpha_site)
        else:t, p = enrichment_thresholds(u, w, exposure_treated, dispersion, settings.alpha_site)
        threshold[index], probability[index] = t[inverse], p[inverse]
    flagged = (treated <= threshold) if settings.direction == "depleted" else (treated >= threshold) & (threshold >= 0)
    flagged &= eligible & (probability > 0)
    return dict(fit=fit, eligible=eligible, threshold=threshold, probability=probability, flagged=flagged)


def analyze(mocks, treated, exposures, genes, gene_indices, positions, settings=None, pool_names=None):
    """mocks, treated: sites x P raw integer counts (one column per pool); exposures: P x 2 (mock, treated)."""
    settings = settings or Settings(); settings.validate()
    mocks = ccg._integer_counts(mocks, "Mock counts"); treated = ccg._integer_counts(treated, "Treated counts")
    exposures = np.asarray(exposures, dtype=float)
    if mocks.ndim != 2 or mocks.shape != treated.shape or mocks.shape[1] < 2 or exposures.shape != (mocks.shape[1], 2) or np.any(exposures <= 0):
        raise ValueError("Need sites x >=2 pools of mock and treated counts and one positive (mock, treated) exposure pair per pool")
    if len(genes) != len(gene_indices) or len({g["id"] for g in genes}) != len(genes):
        raise ValueError("Gene identities and site lists must align and be unique")
    P = mocks.shape[1]; names = pool_names or [f"pool{j+1}" for j in range(P)]
    models = [pool_model(mocks[:, j], treated[:, j], *exposures[j], settings) for j in range(P)]
    rows = []
    for gene, idx in zip(genes, gene_indices):
        idx = np.asarray(idx, dtype=np.int64)
        row = dict(ORF=gene["id"], Name=gene["name"], Description=gene["description"], chromosome=gene["chrom"], strand=gene["strand"])
        pool_p, details = [], []
        for name, model in zip(names, models):
            sites = idx[model["eligible"][idx]]
            if len(sites) < settings.min_sites:details.append(f"{name}:n={len(sites)},untested"); continue
            r, _, _ = ccg.longest_run(model["flagged"][sites])
            tail = ccg.exact_run_tails(model["probability"][sites], max(r, 1)+1)
            p = float(tail[r-1]) if r else 1.
            pool_p.append(p); details.append(f"{name}:n={len(sites)},flagged={int(model['flagged'][sites].sum())},r={r},p={p:.3g}")
        n = len(pool_p); row |= dict(testable_pools=n, pools_p_le_05=int(sum(p <= .05 for p in pool_p)), pool_detail=";".join(details),
                                     pool_p_values=";".join(f"{p:.6g}" for p in pool_p))
        if n < settings.min_testable_pools:
            rows.append(row | dict(status="too few testable pools", in_family=False)); continue
        u = required_pools(n, settings)
        rows.append(row | dict(required_pools=u, p_value=partial_conjunction(pool_p, u), status="tested", in_family=True))
    family = [r for r in rows if r["in_family"]]
    if family:
        for row, q in zip(family, ccg.benjamini_hochberg([r["p_value"] for r in family])):row["q_bh"] = float(q)
    label = f"{settings.direction} run (replicated)"
    for row in rows:
        row["call"] = "not tested" if not row["in_family"] else label if row["q_bh"] <= settings.fdr else "not significant"
    diagnostics = dict(direction=settings.direction, pools=names, family_size=len(family),
        pool_fits={name:{k:m["fit"][k] for k in ("m", "tau", "phi", "linear", "boundary", "sites")} | dict(
            eligible_sites=int(m["eligible"].sum()), flagged_sites=int(m["flagged"].sum()), expected_flagged=float(m["probability"][m["eligible"]].sum()))
            for name, m in zip(names, models)},
        calls={c:sum(r["call"] == c for r in rows) for c in (label, "not significant", "not tested")})
    return rows, diagnostics


def run(normalized, gff, reference, output, controls, treatments, settings=None, normalization=None):
    """File-based runner: one or more mocks and treated cultures per independent library (replicates within a library are summed)."""
    settings = settings or Settings(); settings.validate()
    design = ccg.load_design(normalized, normalization, reference, output, controls, treatments, min_controls=2, min_treatments=2, stratified=True)
    J = len(design["controls"]); libraries = list(dict.fromkeys(design["libraries"]))
    if len(libraries) < 2:raise ValueError("Independent-pool analysis needs >=2 libraries; use control-gumbel for replicates within one library")
    counts, inverse = design["counts"], 1/design["factors"]
    mocks, treated, exposures = [], [], []
    for library in libraries:
        c = [j for j in range(J) if design["libraries"][j] == library]; t = [j for j in range(J, len(design["selected"])) if design["libraries"][j] == library]
        mocks.append(counts[:, c].sum(1)); treated.append(counts[:, t].sum(1))
        exposures.append((inverse[c].sum()/inverse[c].sum(), inverse[t].sum()/inverse[c].sum()))   # relative to the pool's mock
    genes = ccg.annotations(Path(gff).resolve()); sites = design["sites"]
    positions = np.array([p for _, p in sites], dtype=np.int64)
    rows, diagnostics = analyze(np.column_stack(mocks), np.column_stack(treated), exposures, genes, ccg.gene_sites(genes, sites, settings), positions,
                                settings, libraries)
    output = design["output"]; output.mkdir(parents=True, exist_ok=True)
    leading = ["ORF", "Name", "Description", "testable_pools", "required_pools", "pools_p_le_05", "p_value", "q_bh", "call", "status", "pool_detail"]
    write_table(output/"genes.tsv", leading+list(dict.fromkeys(k for r in rows for k in r if k not in leading)), rows)
    norm = design["norm"]
    manifest = dict(status="success", analysis=MODEL, algorithm_revision=REVISION, parameters=asdict(settings), libraries=libraries,
        controls=design["chosen"][:J], treatments=design["chosen"][J:], exposures={lib:dict(mock=e[0], treated=e[1]) for lib, e in zip(libraries, exposures)},
        normalization=dict(method=norm["method"], input=str(norm["folder"]), catalogue=str(design["catalogue"]),
            factors={x:float(f) for x, f in zip(design["selected"], design["factors"])},
            hashes={name:sha256(norm["folder"]/name) for name in ("normalized_counts.npz", "factors.tsv", "method.json")}
                   | dict(run_manifest=sha256(norm["root"]/"run_manifest.json")),
            check="normalized_counts == raw catalogue counts * factor (verified)"),
        model=dict(noise="per pool: lognormal-NB fitted on the pool's mock/treated pair, Var = mu(1+a) + phi mu^2, selection-aware",
            site_flag=f"{settings.direction}: treated count beyond the alpha_site tail of its predictive given the pool's mock",
            pool_test="longest run of flagged sites in the gene within the pool; exact conditional tail",
            replication="Fisher partial conjunction: at least u of the gene's n testable pools (u = floor(n/2)+1 unless min_pools set)",
            multiple_testing="BH over genes with >= min_testable_pools testable pools"),
        diagnostics=diagnostics, implementation_sha256=sha256(__file__), ccg_implementation_sha256=sha256(Path(ccg.__file__)),
        limitations=["Sites are never matched across pools; replication is at gene level",
            "Noise is fitted on mock/treated pairs, so genuine effects inflate it (conservative) and unequal mock/treated noise is not separated",
            "Pools are assumed independent; libraries from one parental strain may share background effects",
            "Scalar normalization exposures assume most sites are unaffected by treatment",
            "Sites are conditionally independent given abundance within a pool"])
    (output/"run_manifest.json").write_text(json.dumps(ccg._json_ready(manifest), indent=2, allow_nan=False)+"\n")
    return manifest
