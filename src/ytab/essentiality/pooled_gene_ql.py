"""Gene-level negative-binomial quasi-likelihood test for independent pools (biological replicates).

For designs where each pool is an independent library and its mock and treated cultures do not share
site-level structure (independent outgrowth/bottlenecks), site-level tests have no power. Reads are
summed over each gene's sites (N/C-terminal trimming) and modelled jointly over pools. The count is
either the number of distinct clones (sites above a reads-per-million threshold; default, for bottlenecked
cultures where per-clone read counts are dominated by amplification) or the summed reads:

    y[g, pool, arm] ~ NB(mu, phi_g),  log mu = log exposure[pool, arm] + alpha[g, pool] + beta[g] * treated

alpha (one intercept per pool) absorbs each pool's representation of the gene, so pools need not be
comparable to each other; beta is the treatment log fold change shared by pools. phi follows an
abundance trend of Cox-Reid adjusted-profile-likelihood estimates. Gene-specific extra variation,
including disagreement between pools, enters a quasi-likelihood dispersion (residual deviance /
(P - 1)) that is shrunk by empirical Bayes towards an abundance trend (Lund et al. 2012; edgeR
glmQLFit/glmQLFTest; limma squeezeVar). The one-sided test uses the signed root of the QL F
statistic on d0 + P - 1 degrees of freedom; BH controls the FDR. A gene is called when q <= fdr and
its per-pool log2 fold change points in the tested direction in at least u = floor(P/2)+1 pools.
Input is a YTAB scalar normalization (ttr, totreads or nzmean), used as exposures; no normalization is
computed here.
"""
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.special import digamma, gammaln, polygamma
from scipy.stats import t as student_t

from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.normalization.ttr import sha256, write_table

MODEL = "Hermes gene-level NB quasi-likelihood test for independent pools"
REVISION = "pooled_gene_ql_v1"
LOG_PHI_GRID = np.linspace(np.log(1e-4), np.log(1e2), 121)


@dataclass(frozen=True)
class Settings:
    direction: str = "depleted"     # depleted | enriched
    min_pools: int | None = None    # u for the direction-consistency rule; None = floor(P/2)+1
    counts: str = "clones"          # clones: distinct sites with >= min_rpm reads per million; reads: summed raw reads
    min_rpm: float = 1.             # clone threshold, reads per million of the culture's total (clones only)
    min_mock_total: int = 5         # treatment-blind filter: mock count summed over pools >= this ...
                                    # ... and mock count > 0 in >= floor(P/2)+1 pools
    fdr: float = .05
    trend_bins: int = 20
    max_iterations: int = 50
    nterm: float = 0.
    cterm: float = 0.
    ignore_stop: bool = True

    def validate(self):
        if (self.direction not in {"depleted", "enriched"} or (self.min_pools is not None and (not isinstance(self.min_pools, int) or self.min_pools < 1)) or
                self.counts not in {"clones", "reads"} or not self.min_rpm > 0 or
                not isinstance(self.min_mock_total, int) or self.min_mock_total < 1 or not 0 < self.fdr < 1 or
                not isinstance(self.trend_bins, int) or self.trend_bins < 1 or not 0 <= self.nterm < 100 or not 0 <= self.cterm < 100):
            raise ValueError("Invalid direction, min_pools, counts, min_rpm, min_mock_total, FDR, trend bins or terminal percentages")


def design_matrix(P, treatment=True):
    """Columns: one intercept per pool, then the treatment indicator. Rows: mock1..mockP, treated1..treatedP."""
    X = np.vstack([np.eye(P), np.eye(P)])
    return np.column_stack([X, np.r_[np.zeros(P), np.ones(P)]]) if treatment else X


def nb_loglik(y, mu, phi):
    r = 1/phi
    return gammaln(y+r)-gammaln(r)-gammaln(y+1)+r*np.log(r/(r+mu))+y*np.log(np.maximum(mu, 1e-300)/(r+mu))


def nb_deviance(y, mu, phi):
    r = 1/phi; mu = np.maximum(mu, 1e-300)
    term = np.where(y > 0, y*np.log(np.where(y > 0, y, 1)/mu), 0.)
    return 2*(term-(y+r)*np.log((1+phi*y)/(1+phi*mu)))


def fit_glm(y, offset, X, phi, iterations=50):
    """Batched NB IRLS. y: G x N, offset: N, X: N x K, phi: G (or scalar). Returns beta (G x K), mu (G x N), XtWX (G x K x K)."""
    G, N = y.shape; phi = np.broadcast_to(np.asarray(phi, dtype=float), (G,))[:, None]
    pinv = np.linalg.pinv(X)
    beta = (np.log((y+.5)/np.exp(offset)))@pinv.T
    for _ in range(iterations):
        eta = np.clip(beta@X.T+offset, -50, 50); mu = np.exp(eta)
        w = mu/(1+phi*mu); z = eta-offset+(y-mu)/mu
        XtWX = np.einsum("nk,gn,nl->gkl", X, w, X)+1e-10*np.eye(X.shape[1])
        new = np.linalg.solve(XtWX, np.einsum("nk,gn->gk", X, w*z)[..., None])[..., 0]
        new = np.clip(new, -40, 40)
        if np.max(np.abs(new-beta)) < 1e-8:beta = new; break
        beta = new
    eta = np.clip(beta@X.T+offset, -50, 50); mu = np.exp(eta); w = mu/(1+phi*mu)
    return beta, mu, np.einsum("nk,gn,nl->gkl", X, w, X)+1e-10*np.eye(X.shape[1])


def cox_reid_phi(y, offset, X, iterations=50):
    """Gene-wise Cox-Reid adjusted profile likelihood maximiser on LOG_PHI_GRID (parabolic refinement)."""
    apl = np.empty((len(y), len(LOG_PHI_GRID)))
    for k, log_phi in enumerate(LOG_PHI_GRID):
        phi = np.exp(log_phi); _, mu, XtWX = fit_glm(y, offset, X, phi, iterations)
        apl[:, k] = nb_loglik(y, mu, phi).sum(1)-.5*np.linalg.slogdet(XtWX)[1]
    best = np.argmax(apl, axis=1); inner = (best > 0) & (best < len(LOG_PHI_GRID)-1)
    log_phi = LOG_PHI_GRID[best].copy(); step = LOG_PHI_GRID[1]-LOG_PHI_GRID[0]
    g = np.flatnonzero(inner); a, b, c = apl[g, best[g]-1], apl[g, best[g]], apl[g, best[g]+1]
    denominator = a-2*b+c; shift = np.where(denominator < 0, .5*(a-c)/np.where(denominator < 0, denominator, -1), 0.)
    log_phi[g] += np.clip(shift, -.5, .5)*step
    return log_phi, ~inner


def binned_trend(x, values, bins):
    """Median of values in equal-count bins of x, linearly interpolated at x (constant beyond the end bins)."""
    order = np.argsort(x, kind="stable"); groups = np.array_split(order, max(1, min(bins, len(x)//30)))
    centres = np.array([np.median(x[g]) for g in groups]); medians = np.array([np.median(values[g]) for g in groups])
    return np.interp(x, centres, medians) if len(groups) > 1 else np.full(len(x), medians[0])


def trigamma_inverse(x):
    """Newton solve of trigamma(y) = x (limma::trigammaInverse)."""
    x = np.asarray(x, dtype=float); y = .5+1/x
    for _ in range(50):
        tri = polygamma(1, y); delta = tri*(1-tri/x)/polygamma(2, y); y = y+delta
        if np.all(-delta/y < 1e-8):break
    return y


def squeeze_variances(s2, df, covariate):
    """limma fitFDist with a quadratic abundance trend: prior df d0 and trended prior variance s0^2; posterior variances."""
    s2 = np.maximum(s2, 1e-8*np.median(s2[s2 > 0]) if np.any(s2 > 0) else 1e-8)
    e = np.log(s2)-digamma(df/2)+np.log(df/2)
    basis = np.column_stack([np.ones_like(covariate), covariate, covariate**2])
    coef, *_ = np.linalg.lstsq(basis, e, rcond=None); fitted = basis@coef
    excess = np.sum((e-fitted)**2)/max(len(e)-basis.shape[1], 1)-polygamma(1, df/2)
    if excess > 0:
        d0 = float(2*trigamma_inverse(excess)); s0 = np.exp(fitted+digamma(d0/2)-np.log(d0/2))
        return (d0*s0+df*s2)/(d0+df), d0, s0
    s0 = np.exp(fitted); return s0, np.inf, s0


def analyze(mocks, treated, exposures, genes, gene_indices, settings=None, pool_names=None):
    """mocks, treated: sites x P raw integer counts over the whole catalogue (one column per pool).

    exposures: P x 2 (mock, treated) normalization exposures, any common scale; used for counts="reads".
    For counts="clones" the exposure of a culture is its total number of clones (computed here from the
    same raw counts), because a scalar read normalization does not describe clone detection.
    """
    settings = settings or Settings(); settings.validate()
    mocks = ccg._integer_counts(mocks, "Mock counts"); treated = ccg._integer_counts(treated, "Treated counts")
    exposures = np.asarray(exposures, dtype=float)
    if mocks.ndim != 2 or mocks.shape != treated.shape or mocks.shape[1] < 2 or exposures.shape != (mocks.shape[1], 2) or np.any(exposures <= 0):
        raise ValueError("Need sites x >=2 pools of mock and treated counts and one positive (mock, treated) exposure pair per pool")
    if len(genes) != len(gene_indices) or len({g["id"] for g in genes}) != len(genes):
        raise ValueError("Gene identities and site lists must align and be unique")
    P = mocks.shape[1]; names = pool_names or [f"pool{j+1}" for j in range(P)]
    counts = np.column_stack([mocks, treated])
    if settings.counts == "clones":
        counts = (counts >= settings.min_rpm*counts.sum(0)/1e6) & (counts > 0); e = counts.sum(0).astype(float)
        if np.any(e == 0):raise ValueError("A culture has no clones above min_rpm")
    else:e = np.r_[exposures[:, 0], exposures[:, 1]]
    y = np.array([counts[np.asarray(i, dtype=np.int64)].sum(0) for i in gene_indices], dtype=float).reshape(len(genes), 2*P)
    sites = np.array([len(i) for i in gene_indices])
    e = e/np.exp(np.mean(np.log(e))); offset = np.log(e)
    u = min(P, settings.min_pools) if settings.min_pools else P//2+1
    tested = (sites > 0) & (y[:, :P].sum(1) >= settings.min_mock_total) & ((y[:, :P] > 0).sum(1) >= u)
    rows = []
    for g, gene in enumerate(genes):
        rows.append(dict(ORF=gene["id"], Name=gene["name"], Description=gene["description"], chromosome=gene["chrom"], strand=gene["strand"],
                         sites=int(sites[g]), mock_reads=";".join(str(int(v)) for v in y[g, :P]), treated_reads=";".join(str(int(v)) for v in y[g, P:]),
                         in_family=bool(tested[g]), status="tested" if tested[g] else "too little mock signal (filter) or no sites"))
    family = np.flatnonzero(tested); diagnostics = dict(direction=settings.direction, pools=names, family_size=int(len(family)))
    if not len(family):
        for row in rows:row["call"] = "not tested"
        return rows, diagnostics
    Y = y[family]; X1, X0 = design_matrix(P), design_matrix(P, treatment=False); df = P-1
    abundance = np.log2(np.mean(Y/e, axis=1)+.5)
    log_phi, at_bound = cox_reid_phi(Y, offset, X1, settings.max_iterations)
    phi = np.exp(binned_trend(abundance, log_phi, settings.trend_bins))
    beta1, mu1, _ = fit_glm(Y, offset, X1, phi, settings.max_iterations); _, mu0, _ = fit_glm(Y, offset, X0, phi, settings.max_iterations)
    deviance1 = nb_deviance(Y, mu1, phi[:, None]).sum(1); deviance0 = nb_deviance(Y, mu0, phi[:, None]).sum(1)
    lr = np.maximum(deviance0-deviance1, 0.)
    s2 = deviance1/df; s2_post, d0, s0 = squeeze_variances(s2, df, abundance)
    total_df = min(d0+df, 1e6); statistic = np.sign(beta1[:, -1])*np.sqrt(lr/s2_post)
    p = student_t.cdf(statistic, total_df) if settings.direction == "depleted" else student_t.sf(statistic, total_df)
    q = ccg.benjamini_hochberg(p)
    pool_lfc = np.log2((Y[:, P:]+.5)/e[P:])-np.log2((Y[:, :P]+.5)/e[:P])
    agree = (pool_lfc < 0).sum(1) if settings.direction == "depleted" else (pool_lfc > 0).sum(1)
    label = f"{settings.direction} (replicated)"
    for k, g in enumerate(family):
        rows[g] |= dict(log2FC=float(beta1[k, -1]/np.log(2)), pool_log2FC=";".join(f"{v:.4g}" for v in pool_lfc[k]), pools_in_direction=int(agree[k]),
                        required_pools=u, log2_abundance=float(abundance[k]), nb_phi_genewise=float(np.exp(log_phi[k])), nb_phi_trend=float(phi[k]),
                        ql_dispersion=float(s2[k]), ql_dispersion_posterior=float(s2_post[k]), statistic=float(statistic[k]),
                        p_value=float(p[k]), q_bh=float(q[k]))
    for row in rows:
        row["call"] = ("not tested" if not row["in_family"] else
                       label if row["q_bh"] <= settings.fdr and row["pools_in_direction"] >= row["required_pools"] else "not significant")
    diagnostics |= dict(counts=settings.counts, exposures=e.tolist(), prior_df=d0, residual_df=df, required_pools=u, nb_phi_trend_range=[float(phi.min()), float(phi.max())],
                        genewise_phi_at_grid_bound=int(at_bound.sum()), calls={c:sum(r["call"] == c for r in rows) for c in (label, "not significant", "not tested")})
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
        mocks.append(counts[:, c].sum(1)); treated.append(counts[:, t].sum(1)); exposures.append((inverse[c].sum(), inverse[t].sum()))
    genes = ccg.annotations(Path(gff).resolve()); sites = design["sites"]
    rows, diagnostics = analyze(np.column_stack(mocks), np.column_stack(treated), exposures, genes, ccg.gene_sites(genes, sites, settings), settings, libraries)
    output = design["output"]; output.mkdir(parents=True, exist_ok=True)
    leading = ["ORF", "Name", "Description", "log2FC", "pool_log2FC", "pools_in_direction", "required_pools", "p_value", "q_bh", "call", "status"]
    write_table(output/"genes.tsv", leading+list(dict.fromkeys(k for r in rows for k in r if k not in leading)), rows)
    norm = design["norm"]
    manifest = dict(status="success", analysis=MODEL, algorithm_revision=REVISION, parameters=asdict(settings), libraries=libraries,
        controls=design["chosen"][:J], treatments=design["chosen"][J:], exposures={lib:dict(mock=e[0], treated=e[1]) for lib, e in zip(libraries, exposures)},
        normalization=dict(method=norm["method"], input=str(norm["folder"]), catalogue=str(design["catalogue"]),
            factors={x:float(f) for x, f in zip(design["selected"], design["factors"])},
            hashes={name:sha256(norm["folder"]/name) for name in ("normalized_counts.npz", "factors.tsv", "method.json")}
                   | dict(run_manifest=sha256(norm["root"]/"run_manifest.json")),
            check="normalized_counts == raw catalogue counts * factor (verified)"),
        model=dict(counts=("clones: distinct catalogue sites with >= min_rpm reads per million in the gene; exposure = the culture's clone total"
                           if settings.counts == "clones" else "raw reads summed over the gene's catalogue sites; exposure = normalization")+" (N/C-terminal trimming)",
            glm="NB GLM per gene: pool intercepts + shared treatment log fold change; offsets = log normalization exposures",
            nb_dispersion="Cox-Reid APL gene-wise estimates, binned-median abundance trend",
            test="quasi-likelihood F (residual deviance / (P-1), EB-squeezed to a quadratic abundance trend); one-sided signed root, t on d0+P-1 df",
            replication=f"call = BH q <= fdr and per-pool log2FC in the tested direction in >= u pools (u = floor(P/2)+1 unless min_pools set)",
            filter="mock count summed over pools >= min_mock_total and > 0 in >= u pools (treatment-blind)"),
        diagnostics=diagnostics, implementation_sha256=sha256(__file__),
        limitations=["Gene-level: within-gene patterns (domains, partial loss) are not modelled",
            "Pools are assumed independent; libraries from one parental strain may share background effects",
            "Scalar normalization exposures assume most genes are unaffected by treatment",
            "With P pools the residual df is P-1; power relies on empirical-Bayes sharing across genes"])
    (output/"run_manifest.json").write_text(json.dumps(ccg._json_ready(manifest), indent=2, allow_nan=False)+"\n")
    return manifest
