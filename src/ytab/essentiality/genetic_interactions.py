"""Hermes genetic interactions, adapted from TRANSIT's GI method (DeJesus et al. 2017): does a gene's treatment response
differ between two strains?

Design (as in TRANSIT's GI): strain A (reference, e.g. wild type) and strain B (query, e.g. a mutant), each
with condition 1 (control) and condition 2 (treated). The interaction for gene g is
    delta_g = log2FC_B(g) - log2FC_A(g).
TRANSIT pools site x replicate counts and samples a normal posterior of each group mean; on Hermes the
clone-abundance variation between sites swamps the posterior, and it is nearly powerless
(output/essentiality/gi/fit_check/). This adaptation keeps TRANSIT's quantity and interaction classes but
estimates each strain's response site by site, inside the pool that was split into control and treated:

- Normalization: BGC (or any YTAB normalization) per library; then each treated culture is rescaled to its
  paired control by the iterative-median factor (multicondition.iterative_median_factor), so every
  replicate response is centred and depth offsets between pairs do not inflate the replicate variance.
- Replicate responses: within each library, control and treated samples are paired in order (equal numbers
  required). For pair r, R_r(g) = mean over the gene's sites eligible in that library of
  log2((treated + c)/(control + c)), c = pseudocount. A pair is one control and one treated culture of one pool,
  so founder abundance cancels site by site.
- Strain response: R(g) = mean_r R_r(g); s^2(g) = var_r R_r(g) (>= 2 pairs per strain).
- Variance of a pair's response, from two independent estimates per gene and strain:
    s^2_g: variance of R_r over the strain's pairs (m - 1 degrees of freedom);
    v_g: a site-based estimate. For pairs of one library (cultures split from one pool share sites): the
      pooled variance over sites of each site's log ratio around its mean over the library's pairs, divided by
      the number of sites (n(m - 1) df). For pairs from different libraries (different sites): the variance over
      sites of the log ratios divided by the sites (n - 1 df);
    combined: sigma^2_g = ((m-1) s^2_g + df_v (kappa v_g + tau^2))/((m-1) + df_v), df = (m-1) + df_v, where
      E[s^2_g] = kappa v_g + tau^2 is fitted across genes by least squares (kappa > 0, tau^2 >= 0): kappa calibrates
      the site-based scale to the replicate scale, tau^2 is the pair-level floor (e.g. each pair's recentering error
      shifts all sites alike, which site-based variances miss). Small genes keep few degrees of freedom.
- Test: t = (R_B - R_A) / sqrt(s~^2_A/m_A + s~^2_B/m_B), Welch-Satterthwaite df; two-sided p.
- Empirical null (Efron 2004; optional, off by default): t is converted to z (same two-sided p); the null centre is the median z
  and the null scale max(1, 1.4826 MAD(z)) over tested genes, assuming most genes do not interact. p-values are
  recomputed against N(centre, scale^2), then BH. The scale is never below 1, so the correction cannot make the test
  more liberal. Model p-values (p_model) are reported as well.
- Interval (in place of TRANSIT's HDI bounds of delta): delta_lower/upper = delta -/+ t_{df, (1+confidence)/2} * se, the
  Welch t interval matching the test (model scale).
- ROPE (region of practical equivalence, from TRANSIT's GI, where it is +-0.5), optional --rope tau > 0: tests H0 |delta| <= tau
  against |delta| > tau (McCarthy & Smyth 2009, TREAT):
      p = P(T_df > (|delta| - tau)/se) + P(T_df < (-|delta| - tau)/se),
  the largest rejection probability over the null interval (attained at its boundary), so the test is valid at
  level alpha for every |delta| <= tau. tau = 0 gives the ordinary two-sided p.
- Multiple testing (TRANSIT's GI offers BFDR or FWER): Benjamini-Hochberg FDR (default) or Holm (1979) family-wise
  error control (--correction holm): adjusted p_(i) = max_{j<=i} min(1, (n-j+1) p_(j)); valid under any dependence
  when p-values are exact. Status (output/essentiality/gi/options_check/): interval and ROPE validated; Holm NOT
  validated (FWER 8.3% vs 5%: family-wise control needs p-values exact far in the tail, where the t approximation
  is least accurate). BH remains the default.
- Call: adjusted p <= FDR/alpha (and |delta| >= min_delta, default 0) -> interaction class (as defined in TRANSIT's GI): Aggravating (delta < 0),
  Alleviating (delta >= 0 and |log2FC_B| < |log2FC_A|), Suppressive (delta >= 0 and |log2FC_B| > |log2FC_A|).
- Clone-sharing check per library: Spearman of log reads at sites detected in both members of each pair; the
  site-paired response assumes control and treated share clones (Gale split cultures: 0.89-0.99).
"""
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy import stats

from ytab.essentiality import multicondition as mc
from ytab.essentiality.control_conditioned_gumbel import _json_ready, annotations, benjamini_hochberg, gene_sites, load_design
from ytab.normalization.ttr import sha256, write_table

MODEL = "Hermes site-paired genetic interactions"
REVISION = "hermes_gi_v1"
CLASSES = ("Aggravating", "Alleviating", "Suppressive")


@dataclass(frozen=True)
class Settings:
    pseudocount: float = 5.          # on the BGC/TTR scale (target mean 100 per site)
    fdr: float = .05
    min_delta: float = 0.            # optional |delta| floor for calls; 0 = test only
    min_sites: int = 1
    recenter: str = "iterative-median"
    clone_warning: float = .5        # warn when control/treated read correlation at shared sites is below this
    empirical_null: bool = False     # optional Efron empirical null (centre and inflate-only scale) on the z-scores
    rope: float = 0.                 # ROPE: test |delta| > rope (TREAT); 0 = test delta != 0
    correction: str = "bh"           # bh (FDR) | holm (FWER)
    confidence: float = .95          # level of the interval for delta
    nterm: float = 0.
    cterm: float = 0.
    ignore_stop: bool = True

    def validate(self):
        if (not self.pseudocount > 0 or not 0 < self.fdr < 1 or self.min_delta < 0 or self.min_sites < 1 or
                self.recenter not in {"none", "iterative-median", "trimmed-sum"} or self.rope < 0 or self.correction not in {"bh", "holm"} or
                not 0 < self.confidence < 1):
            raise ValueError("Invalid GI settings")
        if self.empirical_null and self.rope > 0:raise ValueError("--rope and --empirical-null cannot be combined")


def holm(p):
    """Holm (1979) step-down adjusted p-values (family-wise error rate)."""
    p = np.asarray(p, float); n = len(p); order = np.argsort(p)
    adj = np.minimum(1., np.maximum.accumulate((n-np.arange(n))*p[order]))
    out = np.empty(n); out[order] = adj; return out


def treat_p(delta, se, df, tau):
    """TREAT p-value for H0 |delta| <= tau (McCarthy & Smyth 2009); equals the two-sided t p-value when tau = 0."""
    a = np.abs(delta)
    return stats.t.sf((a-tau)/se, df)+stats.t.cdf((-a-tau)/se, df)


def classify(delta, lfc_b, lfc_a):
    """Interaction classes as defined in TRANSIT's GI (classify_interaction; B = query strain, A = reference strain)."""
    if delta < 0:return "Aggravating"
    if abs(lfc_b) < abs(lfc_a):return "Alleviating"
    if abs(lfc_b) > abs(lfc_a):return "Suppressive"
    return "N/A"


def pairs_by_library(strains, conditions, libraries):
    """{strain: [(library, control column, treated column), ...]} pairing control and treated in order within each library."""
    out = {}
    for lib in dict.fromkeys(libraries):
        cols = [j for j, l in enumerate(libraries) if l == lib]; strain = {strains[j] for j in cols}
        if len(strain) != 1:raise ValueError(f"Library {lib!r} mixes strains {sorted(strain)}; a library is one strain's pool")
        c = [j for j in cols if conditions[j] == "control"]; t = [j for j in cols if conditions[j] == "treated"]
        if not c or len(c) != len(t):
            raise ValueError(f"Library {lib!r} needs equal numbers of control and treated samples (has {len(c)} and {len(t)}); they are paired in order")
        out.setdefault(strain.pop(), []).extend((lib, a, b) for a, b in zip(c, t))
    return out


def clone_sharing(raw, pairs):
    """Spearman of log reads at sites detected in both samples of each pair."""
    out = []
    for lib, c, t in pairs:
        both = (raw[:, c] > 0) & (raw[:, t] > 0)
        out.append(dict(library=lib, shared_sites=int(both.sum()),
                        read_rho=float(stats.spearmanr(raw[both, c], raw[both, t])[0]) if both.sum() >= 100 else None))
    return out


def analyze(values, strains, conditions, libraries, genes, gene_indices, reference="A", query="B", settings=None, raw=None):
    """values: sites x samples normalized counts; strains (two labels), conditions ('control'/'treated') and
    libraries: one label per column. Returns (rows, diagnostics)."""
    settings = settings or Settings(); settings.validate()
    values = np.array(values, float); strains, conditions, libraries = list(strains), list(conditions), list(libraries)
    if not {reference, query} == set(strains):raise ValueError(f"Samples must belong to strains {reference!r} and {query!r}")
    if not set(conditions) <= {"control", "treated"}:raise ValueError("Conditions must be 'control' or 'treated'")
    pairs = pairs_by_library(strains, conditions, libraries)
    for s in (reference, query):
        if len(pairs.get(s, [])) < 2:raise ValueError(f"Strain {s!r} needs >= 2 control/treated pairs (replicate cultures or pools)")
    eligible = {lib:values[:, [j for j, l in enumerate(libraries) if l == lib]].sum(1) > 0 for lib in dict.fromkeys(libraries)}
    factors = {}
    for s in (reference, query):                              # recentre each treated culture to its paired control
        for lib, c, t in pairs[s]:
            rows_lib = np.flatnonzero(eligible[lib])
            idx_lib = [np.searchsorted(rows_lib, np.asarray(i, dtype=np.int64)[eligible[lib][np.asarray(i, dtype=np.int64)]]) for i in gene_indices]
            f = (mc.iterative_median_factor(values[rows_lib][:, [c]], values[rows_lib][:, [t]], idx_lib) if settings.recenter == "iterative-median" else
                 mc.trimmed_sum_factor(values[rows_lib][:, [c]], values[rows_lib][:, [t]], idx_lib, .2) if settings.recenter == "trimmed-sum" else 1.)
            values[:, t] *= f; factors[f"{s}:{lib}:{t}"] = f
    resp, n_sites, site_var, site_df = {}, {}, {}, {}
    pc = settings.pseudocount
    for s in (reference, query):
        R = np.full((len(gene_indices), len(pairs[s])), np.nan); N = np.zeros_like(R)
        by_lib = {}
        for r, (lib, c, t) in enumerate(pairs[s]):by_lib.setdefault(lib, []).append((r, c, t))
        V = np.full((len(gene_indices), len(by_lib)), np.nan); W = np.zeros_like(V); F = np.zeros_like(V)
        for b, (lib, members) in enumerate(by_lib.items()):
            for g, idx in enumerate(gene_indices):
                idx = np.asarray(idx, dtype=np.int64); i = idx[eligible[lib][idx]]
                if not len(i):continue
                D = np.column_stack([np.log2((values[i, t]+pc)/(values[i, c]+pc)) for _, c, t in members])   # sites x pairs
                for k, (r, _, _) in enumerate(members):R[g, r] = D[:, k].mean(); N[g, r] = len(i)
                if D.shape[1] > 1:                            # shared sites: replicate deviations at each site
                    V[g, b] = ((D-D.mean(1, keepdims=True))**2).sum()/(len(i)*(D.shape[1]-1))/len(i); F[g, b] = len(i)*(D.shape[1]-1)
                elif len(i) > 1:V[g, b] = D[:, 0].var(ddof=1)/len(i); F[g, b] = len(i)-1
                W[g, b] = D.shape[1]
        resp[s], n_sites[s] = R, N
        with np.errstate(invalid="ignore"):
            site_var[s] = (np.where(np.isfinite(V), V, 0)*W).sum(1)/np.maximum((np.isfinite(V)*W).sum(1), 1)
            site_df[s] = F.sum(1)
    m = {s:resp[s].shape[1] for s in resp}
    ok = np.ones(len(gene_indices), bool)
    for s in resp:ok &= np.isfinite(resp[s]).all(1) & (n_sites[s].min(1) >= settings.min_sites)
    mean = {s:np.nanmean(resp[s], 1) for s in resp}; var, df_s, kappa = {}, {}, {}
    for s in resp:
        s2 = np.nanvar(resp[s], 1, ddof=1); dv = np.where(ok, site_df[s], 0.)
        use = ok & (dv > 0)
        if use.sum() >= 10:
            k_, t_ = np.linalg.lstsq(np.column_stack([site_var[s][use], np.ones(use.sum())]), s2[use], rcond=None)[0]
            if t_ < 0:k_, t_ = float(np.sum(s2[use])/np.sum(site_var[s][use])), 0.
            kappa[s] = (float(max(k_, 1e-6)), float(t_))
        else:kappa[s] = (1., 0.)
        sigma2 = np.where(dv > 0, ((m[s]-1)*s2+dv*(kappa[s][0]*site_var[s]+kappa[s][1]))/((m[s]-1)+dv), s2)
        var[s] = sigma2/m[s]; df_s[s] = (m[s]-1)+dv
    se2 = var[reference]+var[query]
    with np.errstate(invalid="ignore", divide="ignore"):
        df = se2**2/(var[reference]**2/df_s[reference]+var[query]**2/df_s[query])
    delta = mean[query]-mean[reference]
    with np.errstate(invalid="ignore", divide="ignore"):
        t = delta/np.sqrt(se2); p_model = np.where(ok, treat_p(delta, np.sqrt(se2), df, settings.rope), np.nan)
        half = stats.t.ppf((1+settings.confidence)/2, df)*np.sqrt(se2)
    z = np.where(ok, np.sign(t)*stats.norm.isf(np.clip(p_model, 1e-300, 1)/2), np.nan)
    centre, scale = 0., 1.
    if settings.empirical_null and ok.sum() >= 50:
        centre = float(np.median(z[ok])); scale = float(max(1., 1.4826*np.median(np.abs(z[ok]-centre))))
    p = np.where(ok, 2*stats.norm.sf(np.abs(z-centre)/scale), np.nan) if settings.empirical_null else p_model
    q = np.full(len(p), np.nan); q[ok] = benjamini_hochberg(p[ok]) if settings.correction == "bh" else holm(p[ok])
    rows = []
    for g, gene in enumerate(genes):
        row = dict(ORF=gene["id"], Name=gene["name"], sites_reference=float(n_sites[reference][g].mean()), sites_query=float(n_sites[query][g].mean()))
        if not ok[g]:rows.append(row | dict(status="not tested (too few eligible sites in some pair)", call="not tested")); continue
        called = q[g] <= settings.fdr and abs(delta[g]) >= settings.min_delta
        rows.append(row | dict(log2FC_reference=float(mean[reference][g]), log2FC_query=float(mean[query][g]), delta_log2FC=float(delta[g]),
                               delta_lower=float(delta[g]-half[g]), delta_upper=float(delta[g]+half[g]),
                               se=float(np.sqrt(se2[g])), t=float(t[g]), df=float(df[g]), p_model=float(p_model[g]), p_value=float(p[g]),
                               p_adjusted=float(q[g]), **({"q_bh":float(q[g])} if settings.correction == "bh" else {}),
                               interaction=classify(delta[g], mean[query][g], mean[reference][g]),
                               call=classify(delta[g], mean[query][g], mean[reference][g]) if called else "No interaction", status="tested"))
    sharing = {s:clone_sharing(raw, pairs[s]) for s in pairs} if raw is not None else None
    weak = [x for s in (sharing or {}) for x in sharing[s] if x["read_rho"] is not None and x["read_rho"] < settings.clone_warning]
    diagnostics = dict(model=MODEL, revision=REVISION, settings=asdict(settings), reference_strain=reference, query_strain=query,
                       pairs={s:[dict(library=l, control=c, treated=t_) for l, c, t_ in pairs[s]] for s in pairs},
                       recenter_factors=factors, site_variance_scale={s:dict(kappa=k[0], tau2=k[1]) for s, k in kappa.items()}, tested=int(ok.sum()), empirical_null=dict(centre=centre, scale=scale),
                       calls={c:sum(r["call"] == c for r in rows) for c in CLASSES+("No interaction", "not tested")},
                       clone_sharing=sharing,
                       warning=(f"In {len(weak)} control/treated pair(s) read counts at shared sites correlate below {settings.clone_warning} "
                                "(split cultures of one pool: 0.89-0.99). Control and treated may not share clones; the site-paired "
                                "response is then noisier, so power is reduced (replicate variance keeps the test calibrated)") if weak else None)
    return rows, diagnostics


def run(normalized, gff, reference_fasta, output, samples, reference="A", query="B", settings=None, normalization=None):
    """samples: (observation_id, strain, condition) triples; libraries come from the normalization catalogue."""
    settings = settings or Settings(); settings.validate()
    samples = [(str(o), str(s), str(c)) for o, s, c in samples]
    controls = [o for o, _, c in samples if c == "control"]; treated = [o for o, _, c in samples if c == "treated"]
    if len(controls)+len(treated) != len(samples):raise ValueError("Each sample's condition must be 'control' or 'treated'")
    design = load_design(normalized, normalization, reference_fasta, output, controls, treated, min_controls=1, stratified=True,
                         nonlinear=True, check_roles=False)
    meta = {o:(s, c) for o, s, c in samples}; selected = design["selected"]
    strains = [meta[o][0] for o in selected]; conditions = [meta[o][1] for o in selected]
    genes = annotations(Path(gff).resolve()); idx = gene_sites(genes, design["sites"], settings)
    rows, diag = analyze(design["scaled"].copy(), strains, conditions, design["libraries"], genes, idx, reference, query, settings, raw=design["counts"])
    output = design["output"]; output.mkdir(parents=True, exist_ok=True)
    lead = ["ORF", "Name", "log2FC_reference", "log2FC_query", "delta_log2FC", "delta_lower", "delta_upper", "se", "t", "df", "p_value", "p_adjusted",
            "q_bh", "call", "interaction", "status"]
    write_table(output/"genes.tsv", lead+list(dict.fromkeys(k for r in rows for k in r if k not in lead)), rows)
    norm = design["norm"]
    manifest = dict(status="success", analysis=MODEL, algorithm_revision=REVISION, parameters=asdict(settings),
                    samples=[dict(observation_id=o, strain=meta[o][0], condition=meta[o][1], library=l) for o, l in zip(selected, design["libraries"])],
                    normalization=dict(method=norm["method"], input=str(norm["folder"]), catalogue=str(design["catalogue"]),
                                       hashes={n:sha256(norm["folder"]/n) for n in ("normalized_counts.npz", "factors.tsv", "method.json")}),
                    diagnostics=diag, implementation_sha256=sha256(__file__),
                    limitations=["Needs >= 2 control/treated pairs per strain (replicate cultures of one pool, or several pools)",
                                 "Control and treated of a pair must come from one pool (shared clones); otherwise power drops (see clone_sharing)",
                                 "Sites within a gene are averaged as units; the test's variance comes from replicate pairs",
                                 "Recentering assumes most genes do not respond to treatment in each strain"])
    (output/"run_manifest.json").write_text(json.dumps(_json_ready(manifest), indent=2, allow_nan=False)+"\n")
    return manifest
