"""Hermes window HMM for gene essentiality, adapted from TRANSIT's HMM (DeJesus & Ioerger 2013) and trained by Baum-Welch.

TRANSIT's HMM (DeJesus & Ioerger 2013) walks TA sites with four states (ES, GD, NE, GA), geometric
emissions whose means are fixed fractions of the sample's mean read count, and transitions set from
the insertion density; nothing is fitted. Hermes has no fixed site universe (weak nTnnnnAn consensus,
chromatin-biased), so this module:
- observes fixed genomic windows (default 100 bp): y[w, l] = distinct insertion sites in window w of
  library l (several independent libraries are emitted jointly, conditionally independent given the state);
- masks windows that reads cannot map to uniquely (k-mer uniqueness of the reference; missing emission);
- models transposon accessibility with window covariates (sequence preference, AT content, fraction outside CDS); the
  sequence preference is the Hermes nTnnnnAn count (v2) or, in v3, a target-site model learned from the libraries'
  own insertion sites (learn_target_model), so any transposon can be modelled:
      y[w, l] | state k, u_w ~ Poisson(exp(a_l + x_w.b) * r_k * u_w),  u_w ~ Gamma(mean 1, variance psi),
      shared by the libraries (window accessibility common to independent libraries), giving a
      negative-multinomial joint emission (an NB for one library);   r_ES < r_GD < r_NE = 1 < r_GA;
- fits transitions, initial probabilities and all emission parameters by Baum-Welch (EM), deterministically;
- reports per-gene posterior state probabilities (window posteriors weighted by overlap with the trimmed
  CDS) and, as the essentiality score, P_ES_run: the exact posterior probability that the trimmed CDS
  contains a run of >= min_es_windows consecutive ES windows (mappable windows count toward the run; an
  unmappable ES window neither extends nor breaks it), from a dynamic program over the posterior Markov
  chain. ES call when P_ES_run >= es_call_threshold; otherwise the largest mean posterior of GD/NE/GA
  (v2; v1 called the largest mean posterior, which misses genes with an essential domain).
Forward-backward uses prefix products of scaled 4x4 matrices (exact up to floating point).
"""
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import digamma, gammaln

MODEL = "Hermes window HMM"
REVISION = "hermes_hmm_v2"
REVISION_LEARNED = "hermes_hmm_v3"        # v2 with the learned target-site covariate in place of the motif count
TARGET_OFFSETS = (-10, 17)               # target-site model: offsets around the hit coordinate, read outward from the transposon
TARGET_FLANK = 100                       # local background: bases within this distance of the insertion sites
STATES = ("ES", "GD", "NE", "GA")
COMPLEMENT = np.array([3, 2, 1, 0, 4], dtype=np.uint8)   # A C G T N -> T G C A N


@dataclass(frozen=True)
class Settings:
    window: int = 100
    min_reads: int = 1                  # a site counts when it has >= min_reads reads
    kmer: int = 50                      # mappability k-mer length
    unmappable_fraction: float = .5     # mask a window when >= this fraction of its k-mer starts are not unique
    max_iterations: int = 300
    tolerance: float = 1e-8             # relative log-likelihood change
    mstep_iterations: int = 25          # L-BFGS iterations per M-step (generalized EM)
    init_r: tuple = (.02, .3, 1., 3.)
    init_stay: float = .9
    nterm: float = 0.                   # gene trimming (% of CDS length)
    cterm: float = 10.
    min_gene_mappable: float = .5       # genes with less mappable trimmed CDS are reported but not called
    min_es_windows: int = 3             # v2: essential run length (windows) for P_ES_run
    es_call_threshold: float = .5

    def validate(self):
        if (not isinstance(self.window, int) or self.window < 10 or self.min_reads < 1 or not 16 <= self.kmer <= 64 or
                not 0 < self.unmappable_fraction <= 1 or self.max_iterations < 1 or not self.tolerance > 0 or
                len(self.init_r) != 4 or not all(np.diff(self.init_r) > 0) or self.init_r[2] != 1 or not 0 < self.init_stay < 1 or
                not 0 <= self.nterm < 100 or not 0 <= self.cterm < 100 or self.nterm+self.cterm >= 100 or not 0 <= self.min_gene_mappable <= 1 or
                not isinstance(self.min_es_windows, int) or self.min_es_windows < 1 or not 0 < self.es_call_threshold < 1):
            raise ValueError("Invalid HMM settings")


# ---------------------------------------------------------------- genome features

def encode(sequence):
    table = np.full(256, 4, dtype=np.uint8)
    for i, c in enumerate(b"ACGT"):table[c] = i; table[c+32] = i
    return table[np.frombuffer(sequence.encode(), dtype=np.uint8)]


def read_fasta(path):
    seqs, name, parts = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:seqs[name] = "".join(parts)
            name, parts = line[1:].split()[0], []
        else:parts.append(line.strip())
    if name:seqs[name] = "".join(parts)
    return seqs


def _pack(codes, k):
    """Exact 2-bit packing of every k-mer start (k <= 64) into two uint64 words; invalid (N-containing) starts flagged."""
    n = len(codes)-k+1
    if n <= 0:return np.zeros((0, 2), dtype=np.uint64), np.zeros(0, dtype=bool)
    words = np.zeros((n, 2), dtype=np.uint64); valid = np.ones(n, dtype=bool)
    for j in range(k):
        c = codes[j:j+n]; valid &= c < 4; w = 0 if j < 32 else 1
        words[:, w] = (words[:, w] << np.uint64(2)) | (c & 3).astype(np.uint64)
    return words, valid


def unique_kmer_starts(sequences, k):
    """Per chromosome, bool array over positions: the k-mer starting there occurs once among both strands of the genome."""
    packed = {}
    for name, codes in sequences.items():
        f, fv = _pack(codes, k); r, rv = _pack(COMPLEMENT[codes[::-1]], k)
        packed[name] = (f, fv, r, rv)
    allwords = np.concatenate([np.concatenate([f[fv], r[rv]]) for f, fv, r, rv in packed.values()])
    view = np.ascontiguousarray(allwords).view(np.dtype((np.void, 16))).ravel()
    keys, counts = np.unique(view, return_counts=True)
    out = {}
    for name, (f, fv, r, rv) in packed.items():
        u = np.zeros(len(sequences[name]), dtype=bool)
        fview = np.ascontiguousarray(f).view(np.dtype((np.void, 16))).ravel()
        idx = np.searchsorted(keys, fview); idx = np.clip(idx, 0, len(keys)-1)
        u[:len(f)] = fv & (keys[idx] == fview) & (counts[idx] == 1)
        out[name] = u
    return out


def learn_target_model(sequences, sites, offsets=TARGET_OFFSETS, flank=TARGET_FLANK):
    """Order-0 position weight matrix of the insertion target against the local background (label-free).

    sites: {chromosome: (positions 1-based, strands "W"/"C")}, each distinct site once. A W site at p is read as genome
    p+j; a C site as the complement of genome p-j (the sequence read outward from the transposon). The background Q is
    the oriented base composition within +-flank bp of the sites, so regional composition does not enter the motif.
    Returns offsets, P (L x 4), Q (4), logodds (L x 5, N = 0), information (bits per offset) and the site count."""
    js = np.arange(offsets[0], offsets[1]+1); counts = np.zeros((len(js), 4)); bg = np.zeros(4); n = 0
    for name, (pos, strand) in sites.items():
        if name not in sequences:continue
        codes = sequences[name]; L = len(codes); i0 = np.asarray(pos, dtype=np.int64)-1
        rev = np.asarray(strand) == "C"; n += len(i0)
        onehot = np.zeros((L+1, 4)); onehot[np.arange(L)[codes < 4]+1, codes[codes < 4]] = 1; cs = np.cumsum(onehot, 0)
        lo, hi = np.clip(i0-flank, 0, L), np.clip(i0+flank+1, 0, L)
        local = cs[hi]-cs[lo]; bg += local[~rev].sum(0)+local[rev][:, ::-1].sum(0)   # complement of A C G T is T G C A
        for k, j in enumerate(js):
            idx = np.where(rev, i0-j, i0+j); ok = (idx >= 0) & (idx < L); b = codes[idx[ok]]
            b = np.where(rev[ok], COMPLEMENT[b], b); b = b[b < 4]
            counts[k] += np.bincount(b, minlength=4)[:4]
    P = (counts+1)/(counts.sum(1, keepdims=True)+4); Q = bg/bg.sum()
    logodds = np.zeros((len(js), 5)); logodds[:, :4] = np.log(P/Q)
    return dict(offsets=js, P=P, Q=Q, logodds=logodds, information=(P*np.log2(P/Q)).sum(1), sites=int(n), flank=flank)


def target_propensity(codes, model):
    """Per position: mean over both strands of exp(target-model score) for an insertion at that coordinate."""
    js, lo = model["offsets"], model["logodds"]; L = len(codes); pad = max(abs(js[0]), abs(js[-1]))
    padded = np.concatenate([np.full(pad, 4, np.uint8), codes, np.full(pad, 4, np.uint8)]); comp = COMPLEMENT[padded]
    forward = np.zeros(L); reverse = np.zeros(L)
    for k, j in enumerate(js):
        forward += lo[k][padded[pad+j:pad+j+L]]                 # W: genome i+j
        reverse += lo[k][comp[pad-j:pad-j+L]]                   # C: complement of genome i-j
    return (np.exp(forward)+np.exp(reverse))/2


def window_features(sequences, cds_intervals, settings, target=None):
    """Per chromosome: window starts, covariates (log1p motif count, AT fraction, fraction outside CDS; with a learned
    target model also log mean target propensity) and mappability mask."""
    unique = unique_kmer_starts(sequences, settings.kmer); W = settings.window; out = {}
    for name, codes in sequences.items():
        n = len(codes); nw = -(-n//W); pad = nw*W-n
        motif = np.zeros(n, dtype=float); motif[:n-7] = (codes[1:n-6] == 3) & (codes[6:n-1] == 0)   # nTnnnnAn (palindromic consensus)
        at = ((codes == 0) | (codes == 3)).astype(float); valid = (codes < 4).astype(float)
        incds = np.zeros(n+1, dtype=float)
        for a, b in cds_intervals.get(name, []):incds[a-1] += 1; incds[b] -= 1
        incds = (np.cumsum(incds)[:n] > 0).astype(float)
        nonunique = (~unique[name]).astype(float)
        f = lambda x:np.concatenate([x, np.zeros(pad)]).reshape(nw, W).sum(1)
        length = f(np.ones(n))
        out[name] = dict(start=np.arange(nw)*W+1, length=length, motif=f(motif), at=f(at*valid)/np.maximum(f(valid), 1),
                         outside_cds=1-f(incds)/length, masked=f(nonunique)/length >= settings.unmappable_fraction)
        if target is not None:out[name]["target"] = np.log(f(target_propensity(codes, target))/length)
    return out


def window_counts(sites, starts, window, n_windows):
    """sites: dict chrom -> array of positions (1-based) -> counts per window."""
    y = np.zeros(n_windows, dtype=float); idx = (np.asarray(sites, dtype=np.int64)-1)//window
    np.add.at(y, idx[(idx >= 0) & (idx < n_windows)], 1.); return y


# ---------------------------------------------------------------- forward-backward via scaled prefix products

def _prefix(M):
    """P[t] = M[0] @ ... @ M[t], each rescaled to max 1; returns (P, log scales). Hillis-Steele scan."""
    P = M.copy(); logs = np.zeros(len(M)); s = P.max((1, 2)); P /= s[:, None, None]; logs += np.log(s)
    offset = 1
    while offset < len(P):
        new = P.copy(); new[offset:] = P[:-offset] @ P[offset:]
        nl = logs.copy(); nl[offset:] = logs[:-offset]+logs[offset:]
        s = new[offset:].max((1, 2)); new[offset:] /= s[:, None, None]; nl[offset:] += np.log(s)
        P, logs, offset = new, nl, offset*2
    return P, logs


def forward_backward(logB, A, pi, full=False):
    """logB: T x K emission log-probabilities (0 rows = missing). Returns gamma (T x K), xi sum (K x K), log-likelihood
    (and, with full=True, the per-window scaled emissions B and normalized backward variables beta)."""
    T, K = logB.shape; c = logB.max(1); B = np.exp(logB-c[:, None])
    v0 = pi*B[0]
    if T == 1:
        z = v0.sum(); out = ((v0/z)[None], np.zeros((K, K)), float(np.log(z)+c.sum()))
        return out+(B, np.ones((1, K))/K) if full else out
    M = A[None]*B[1:, None, :]                        # M[t-1] = A diag(B_t), t = 1..T-1
    P, lp = _prefix(M)
    alpha = np.vstack([v0[None], np.einsum("k,tkj->tj", v0, P)]); alpha /= alpha.sum(1, keepdims=True)
    end = v0 @ P[-1]; loglik = float(np.log(end.sum())+lp[-1]+c.sum())
    S, _ = _prefix(np.transpose(M[::-1], (0, 2, 1)))   # S[t] = (M_{T-1} ... M_{T-1-t})^T
    beta = np.vstack([S[::-1].sum(1), np.ones((1, K))]); beta /= beta.sum(1, keepdims=True)
    gamma = alpha*beta; gamma /= gamma.sum(1, keepdims=True)
    xi = alpha[:-1, :, None]*A[None]*(B[1:]*beta[1:])[:, None, :]; xi /= xi.sum((1, 2), keepdims=True)
    return (gamma, xi.sum(0), loglik, B, beta) if full else (gamma, xi.sum(0), loglik)


def run_probability(gamma0, A, B, beta, kept, m, state=0):
    """P(a run of >= m consecutive `state` windows | data) over a window span, from the posterior Markov chain.

    gamma0: posterior at the first window; B, beta: scaled emissions and normalized backward variables for the
    span (rows 1.. give the posterior transitions P(j | i) ∝ A_ij B_j beta_j); kept: mappable windows.
    A `state` window that is not kept neither extends nor breaks a run."""
    K = len(gamma0); p = np.zeros((K, m)); done = 0.
    for j in range(K):
        c = (1 if kept[0] else 0) if j == state else 0
        if c >= m:done += gamma0[j]
        else:p[j, c] += gamma0[j]
    for t in range(1, len(kept)):
        T = A*(B[t]*beta[t])[None, :]; T /= T.sum(1, keepdims=True)
        q = np.zeros((K, m)); q[:, 0] = p.sum(1)@T; q[state, 0] = 0.
        mass = p.T@T[:, state]                         # by previous run count c, moving into `state`
        for c in range(m):
            c2 = c+1 if kept[t] else c
            if c2 >= m:done += mass[c]
            else:q[state, c2] += mass[c]
        p = q
    return float(done)


# ---------------------------------------------------------------- emissions

def nb_logpmf(y, mu, phi):
    n = 1/phi
    return gammaln(y+n)-gammaln(n)-gammaln(y+1)+y*np.log(mu/(mu+n))+n*np.log(n/(mu+n))


def nm_logpmf(Y, MU, psi):
    """Negative multinomial: Poisson(MU[:, l] u) with shared u ~ Gamma(mean 1, variance psi); rows of Y and MU are windows."""
    nu = 1/psi; S = Y.sum(1); M = MU.sum(1)
    return (gammaln(nu+S)-gammaln(nu)-gammaln(Y+1).sum(1)+(Y*np.log(MU)).sum(1)-(nu+S)*np.log(nu+M)+nu*np.log(nu))


def unpack(theta, L, P):
    u = theta[:3]; log_gd = -np.exp(u[0]); log_es = log_gd-np.exp(u[1]); log_ga = np.exp(u[2])
    logr = np.array([log_es, log_gd, 0., log_ga])
    return logr, theta[3:3+L], theta[3+L:3+L+P], theta[3+L+P]


def pack(logr, a, b, logpsi):
    return np.concatenate([[np.log(-logr[1]), np.log(logr[1]-logr[0]), np.log(logr[3])], a, b, [logpsi]])


def emission_loglik(Y, X, keep, theta):
    """T x K log-probabilities; Y: T x L counts, X: T x P covariates, keep: T bool (False = masked)."""
    L, P = Y.shape[1], X.shape[1]; logr, a, b, logpsi = unpack(theta, L, P)
    base = a[None, :]+(X@b)[:, None]                   # T x L
    out = np.zeros((len(Y), 4))
    for k in range(4):out[:, k] = nm_logpmf(Y, np.exp(base+logr[k]), np.exp(logpsi))
    out[~keep] = 0.
    return out


def count_histograms(Y):
    """Distinct window totals S and their frequencies, and sum of lgamma(y+1) (the state-independent terms)."""
    v, c = np.unique(Y.sum(1), return_counts=True)
    return v, c, float(gammaln(Y+1).sum())


def _objective(theta, Y, X, G, hist):
    """Negative expected complete-data emission log-likelihood and gradient (G: T x K posteriors on kept windows).

    lgamma(nu+S)-lgamma(nu)-sum lgamma(y+1) does not depend on the state and enters with total weight 1 per window,
    so it is evaluated once per distinct window total S (hist)."""
    L, P = Y.shape[1], X.shape[1]; logr, a, b, logpsi = unpack(theta, L, P)
    nu = np.exp(-logpsi); base = a[None, :]+(X@b)[:, None]; S = Y.sum(1)
    v, c, lgy = hist
    value = (c*(gammaln(nu+v)-gammaln(nu))).sum()-lgy; g_nu = (c*(digamma(nu+v)-digamma(nu))).sum()
    g_logr = np.zeros(4); g_a = np.zeros(L); g_b = np.zeros(P)
    for k in range(4):
        logmu = base+logr[k]; mu = np.exp(logmu); M = mu.sum(1); w = G[:, k]; nm = nu+M
        value += (w*((Y*logmu).sum(1)-(nu+S)*np.log(nm)+nu*np.log(nu))).sum()
        d = w[:, None]*(Y-((nu+S)/nm)[:, None]*mu)      # d/dlog mu
        g_logr[k] = d.sum(); g_a += d.sum(0); g_b += X.T@d.sum(1)
        g_nu += (w*(-np.log(nm)-(nu+S)/nm+np.log(nu)+1)).sum()
    u = theta[:3]; e0, e1, e2 = np.exp(u)               # chain rule for the ordered r parametrisation
    g_u = np.array([-(g_logr[1]+g_logr[0])*e0, -g_logr[0]*e1, g_logr[3]*e2])
    grad = np.concatenate([g_u, g_a, g_b, [-nu*g_nu]])
    return -value, -grad


# ---------------------------------------------------------------- fitting

def fit(chromosomes, settings=None, theta=None, A=None, pi=None):
    """chromosomes: list of dicts with Y (T x L), X (T x P), keep (T). Baum-Welch; returns parameters, posteriors, trace."""
    settings = settings or Settings(); settings.validate()
    Y = np.vstack([c["Y"] for c in chromosomes]); X = np.vstack([c["X"] for c in chromosomes]); keep = np.concatenate([c["keep"] for c in chromosomes])
    L, P = Y.shape[1], X.shape[1]
    if theta is None:
        a0 = np.log(np.maximum(Y[keep].mean(0), 1e-3)); theta = pack(np.log(np.array(settings.init_r, dtype=float)), a0, np.zeros(P), np.log(.5))
    if A is None:A = np.full((4, 4), (1-settings.init_stay)/3); np.fill_diagonal(A, settings.init_stay)
    if pi is None:pi = np.full(4, .25)
    trace = []; previous = None; hist = count_histograms(Y[keep])
    for iteration in range(settings.max_iterations):
        logB = emission_loglik(Y, X, keep, theta); offsets = np.cumsum([0]+[len(c["Y"]) for c in chromosomes])
        gammas, xis, starts, loglik, chains = [], np.zeros((4, 4)), np.zeros(4), 0., []
        for i in range(len(chromosomes)):
            g, x, ll, Bs, bt = forward_backward(logB[offsets[i]:offsets[i+1]], A, pi, full=True)
            gammas.append(g); xis += x; starts += g[0]; loglik += ll; chains.append((Bs, bt))
        A_used = A
        trace.append(loglik)
        if previous is not None and abs(loglik-previous) <= settings.tolerance*abs(previous):break
        previous = loglik
        G = np.vstack(gammas)
        A = xis/xis.sum(1, keepdims=True); pi = starts/starts.sum()
        # generalized EM: a bounded number of quasi-Newton steps from the current parameters (Q never decreases)
        res = minimize(_objective, theta, args=(Y[keep], X[keep], G[keep], hist), jac=True, method="L-BFGS-B",
                       options=dict(maxiter=settings.mstep_iterations, ftol=1e-13, gtol=1e-9))
        if res.fun <= _objective(theta, Y[keep], X[keep], G[keep], hist)[0]:theta = res.x
    G = np.vstack(gammas)
    logr, a, b, logpsi = unpack(theta, L, P)
    return dict(theta=theta, A=A_used, pi=pi, chains=chains, r=np.exp(logr), a=a, b=b, psi=float(np.exp(logpsi)), loglik=trace[-1], trace=trace,
                iterations=len(trace), converged=len(trace) < settings.max_iterations, gamma=G)


def viterbi(logB, A, pi):
    T, K = logB.shape; logA = np.log(A); delta = np.log(pi)+logB[0]; back = np.zeros((T, K), dtype=np.int8)
    for t in range(1, T):
        s = delta[:, None]+logA; back[t] = s.argmax(0); delta = s.max(0)+logB[t]
    path = np.zeros(T, dtype=np.int8); path[-1] = delta.argmax()
    for t in range(T-1, 0, -1):path[t-1] = back[t, path[t]]
    return path


# ---------------------------------------------------------------- genes

def gene_windows(gene, window, settings):
    """(window index, overlap bp) of the gene's trimmed CDS (own CDS union; N/C trimming in transcript orientation)."""
    pos = np.concatenate([np.arange(a, b+1) for a, b in sorted(gene["blocks"])]); length = len(pos)
    o = np.arange(length) if gene["strand"] != "-" else length-1-np.arange(length)
    keep = (o >= length*settings.nterm/100) & (o < length*(1-settings.cterm/100))
    w, n = np.unique((pos[keep]-1)//window, return_counts=True)
    return list(zip(w.tolist(), n.tolist()))


def gene_table(genes, chrom_index, gamma_by_chrom, keep_by_chrom, settings, chain_by_chrom=None, A=None):
    """Per gene: mean window posteriors, P_ES_run (v2, when chain_by_chrom and A are given) and the call."""
    rows = []
    for gene in genes:
        ws = gene_windows(gene, settings.window, settings) if gene["chrom"] in chrom_index else []
        row = dict(ORF=gene["id"], Name=gene["name"], chromosome=gene["chrom"], strand=gene["strand"], windows=len(ws))
        if not ws:rows.append(row | dict(call="N/A", status="no windows")); continue
        idx = np.array([w for w, _ in ws]); ov = np.array([o for _, o in ws], dtype=float)
        kp = keep_by_chrom[gene["chrom"]][idx]; mappable = float(ov[kp].sum()/ov.sum())
        row["mappable_fraction"] = mappable
        if not kp.any():rows.append(row | dict(call="N/A", status="unmappable")); continue
        gam = gamma_by_chrom[gene["chrom"]]
        g = (gam[idx[kp]]*ov[kp, None]).sum(0)/ov[kp].sum()
        row |= {f"P_{s}":float(p) for s, p in zip(STATES, g)}
        called = mappable >= settings.min_gene_mappable
        if chain_by_chrom is not None:
            span = np.arange(idx.min(), idx.max()+1); inside = np.isin(span, idx[kp])   # only the gene's mappable windows count
            Bs, beta = chain_by_chrom[gene["chrom"]]
            m = max(1, min(settings.min_es_windows, int(kp.sum())))
            row["P_ES_run"] = run_probability(gam[span[0]], A, Bs[span], beta[span], inside, m)
            row["es_run_windows"] = m
            call = "ES" if row["P_ES_run"] >= settings.es_call_threshold else STATES[1+int(g[1:].argmax())]
        else:
            call = STATES[int(g.argmax())]
        row |= dict(call=call if called else "uncertain", status="called" if called else "mostly unmappable")
        rows.append(row)
    return rows


def _json_ready(x):
    if isinstance(x, dict):return {k:_json_ready(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):return [_json_ready(v) for v in x]
    if isinstance(x, np.ndarray):return _json_ready(x.tolist())
    if isinstance(x, (np.floating, np.integer)):return x.item()
    return x
