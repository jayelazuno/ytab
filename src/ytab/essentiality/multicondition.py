"""Hermes multi-condition tests: site-restricted permutation ANOVA and all-pairs resampling.

Adapts the purpose of TRANSIT's ANOVA (does a gene respond in any of K conditions?) and U-test (compare
two conditions), both from DeJesus et al. 2015, to Hermes designs, where one library (pool) is split into several conditions so every
insertion site is the same clone in every condition. TRANSIT's versions pool site x replicate counts as
independent observations; with Hermes clone-abundance variation that loses nearly all power (see
output/essentiality/transit_tests/fit_check/RESULT.md). Here, as in resampling v3:

- Permutation units are (site, library): under the null that the gene does not respond to condition,
  the condition labels of one library's samples are exchangeable at each site, independently across units.
- anova(): statistic T = sum_k N_k (m_k - m)^2, the between-condition sum of squares of the condition
  means m_k over eligible site x sample normalized values (N_k values in condition k, m the grand mean).
  N_k and m are permutation-invariant, so T* depends only on the condition sums S_k. For K = 2,
  T = (N_1 N_2 / N) (m_2 - m_1)^2, so the test equals two-sided resampling v3. p-values are Monte Carlo
  over the exact permutation distribution (seeded per gene), in stages: fresh draws at each stage, and a
  gene moves on only while p <= the stage threshold c. Each stage's p = (b+1)/(B+1) is valid, and so is
  the staged p: for alpha >= c, staged p <= alpha implies the first stage's p <= alpha; for alpha < c,
  it requires the next stage's p (fresh draws, independent of the first stage given the data) <= alpha.
- Each non-reference condition is recentered to the reference by one scalar factor. v2 (default):
  iterative median, the median gene log ratio (genes with pooled mean >= the median), re-estimated on
  genes within 2.5 MAD of it until stable. v1: resampling v3's trimmed-sum factor.
- all_pairs(): resampling v3 for every pair of conditions (within library), one BH family over all pairs.
  v2 (default): every condition is first recentered to the reference as above, and resampling v3 runs
  without its own recentering. v1: resampling v3's per-pair trimmed-sum recentering.

Status: v1 (output/essentiality/multicondition/benchmark/RESULT.md): all_pairs passed; anova failed with
opposite effects in two conditions (8.0% type I error, 13.7% FDR), caused by trimmed-sum miscentering
under enrichment. v2 passed its pre-registered evaluation (output/essentiality/multicondition/benchmark_v2/
RESULT.md): use anova for "responds to any condition" unless its replicate-noise warning fires; then use all_pairs.
"""
from dataclasses import asdict, dataclass
from itertools import combinations
import hashlib
import json
from math import factorial
from pathlib import Path

import numpy as np

from ytab.essentiality import resampling
from ytab.essentiality.control_conditioned_gumbel import _json_ready, annotations, benjamini_hochberg, gene_sites, load_design
from ytab.normalization.ttr import sha256, write_table

MODEL = "Hermes site-restricted permutation ANOVA"
REVISION = "hermes_anova_v2"


@dataclass(frozen=True)
class Settings:
    reference: str = ""                 # condition used for recentering and log2FC (default: the first condition)
    stages: tuple = (999, 9999, 99999)  # Monte Carlo draws per stage (fresh draws each stage)
    thresholds: tuple = (.01, .001)     # continue to the next stage while p <= threshold
    max_enumerate: int = 5000           # enumerate a library's distinct label assignments up to this many
    min_sites: int = 1
    pseudocount: float = 1.
    fdr: float = .05
    seed: int = 2026
    recenter: str = "iterative-median"  # none | iterative-median (v2) | trimmed-sum (v1)
    recenter_trim: float = .2
    nterm: float = 0.                   # gene site selection (ytab gene_sites conventions)
    cterm: float = 0.
    ignore_stop: bool = True
    asymmetry_warning: float = 2.

    def validate(self):
        if (len(self.stages) < 1 or len(self.thresholds) != len(self.stages)-1 or any(int(b) < 99 for b in self.stages) or
                any(not 0 < t < 1 for t in self.thresholds) or self.recenter not in {"none", "iterative-median", "trimmed-sum"} or not 0 <= self.recenter_trim < .5 or
                not isinstance(self.min_sites, int) or self.min_sites < 1 or not self.pseudocount >= 0 or not 0 < self.fdr < 1 or
                not isinstance(self.seed, int) or self.seed < 0 or self.max_enumerate < 1 or not self.asymmetry_warning >= 1):
            raise ValueError("Invalid multi-condition settings")


def trimmed_sum_factor(reference, other, gene_indices, trim):
    """Resampling v3 trimmed-sum recentering: factor that rescales `other` (sites x samples) to `reference`."""
    a_all, b_all = reference.mean(1), other.mean(1); ratios, pooled = [], []
    for idx in gene_indices:
        idx = np.asarray(idx, dtype=np.int64)
        if len(idx):ratios.append((a_all[idx].mean(), b_all[idx].mean())); pooled.append(a_all[idx].mean()+b_all[idx].mean())
    ratios, pooled = np.array(ratios), np.array(pooled); keep = (pooled >= np.median(pooled)) & (ratios[:, 0] > 0) & (ratios[:, 1] > 0)
    if keep.sum() < 10:return 1.
    log_ratio = np.log2(ratios[keep, 1]/ratios[keep, 0])
    lo, hi = np.quantile(log_ratio, [trim, 1-trim]); core = (log_ratio >= lo) & (log_ratio <= hi)
    return float(ratios[keep][core, 0].sum()/ratios[keep][core, 1].sum())


def iterative_median_factor(reference, other, gene_indices, mad_k=2.5, iterations=20):
    """Factor that rescales `other` to `reference` (sites x samples): 2^-m, with m the median log2 ratio of
    gene means over genes with pooled mean >= the median, re-estimated on genes within mad_k MAD of m."""
    a_all, b_all = reference.mean(1), other.mean(1); a, b = [], []
    for idx in gene_indices:
        idx = np.asarray(idx, dtype=np.int64)
        if len(idx):a.append(a_all[idx].mean()); b.append(b_all[idx].mean())
    a, b = np.array(a), np.array(b); keep = (a+b >= np.median(a+b)) & (a > 0) & (b > 0)
    if keep.sum() < 10:return 1.
    r = np.log2(b[keep]/a[keep]); m = float(np.median(r))
    for _ in range(iterations):
        mad = 1.4826*np.median(np.abs(r-m)); core = np.abs(r-m) <= mad_k*mad
        if core.sum() < 10:break
        new = float(np.median(r[core]))
        if abs(new-m) < 1e-9:break
        m = new
    return float(2**-m)


def recenter(values, code, reference, gene_indices, method, trim=.2):
    """Rescale each non-reference condition's columns to the reference; returns (values, factors by condition code)."""
    values = np.array(values, dtype=float); factors = {}
    for k in np.unique(code):
        if k == reference:factors[int(k)] = 1.; continue
        ref, other = values[:, code == reference], values[:, code == k]
        f = (iterative_median_factor(ref, other, gene_indices) if method == "iterative-median" else
             trimmed_sum_factor(ref, other, gene_indices, trim) if method == "trimmed-sum" else 1.)
        values[:, code == k] *= f; factors[int(k)] = f
    return values, factors


def assignments(labels, limit):
    """Distinct arrangements of a label multiset (A x n), or None when there are more than `limit`."""
    labels = np.asarray(labels); values, counts = np.unique(labels, return_counts=True)
    total = factorial(len(labels))
    for c in counts:total //= factorial(int(c))
    if total > limit:return None
    out, current, remaining = [], [], dict(zip(values.tolist(), counts.tolist()))
    def walk():
        if len(current) == len(labels):out.append(list(current)); return
        for v in values.tolist():
            if remaining[v]:
                remaining[v] -= 1; current.append(v); walk(); current.pop(); remaining[v] += 1
    walk()
    return np.array(out)


def _gene_rng(seed, gene_id, stage):
    key = int.from_bytes(hashlib.sha256(str(gene_id).encode()).digest()[:8], "big")
    return np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(key, 11, stage)))


def _statistic(S, N):
    """T = sum_k N_k (m_k - m)^2 from condition sums S (... x K) and counts N (K); conditions with N_k = 0 are skipped."""
    use = N > 0; total = S[..., use].sum(-1); n = N[use].sum()
    return (S[..., use]**2/N[use]).sum(-1)-total**2/n


def _draws(rng, blocks, K, B, chunk=4_000_000):
    """B permutation draws of the condition sums (B x K) over independent (site, library) units."""
    out = np.zeros((B, K))
    for kind, data in blocks:
        u = len(data[0]); rows = max(1, chunk//max(1, u*K*(1 if kind == "enumerated" else len(data[1]))))
        for begin in range(0, B, rows):
            b = min(rows, B-begin)
            if kind == "enumerated":   # data: (unit x assignment x K sums,)
                sums = data[0]; pick = rng.integers(0, sums.shape[1], size=(b, u))
                out[begin:begin+b] += sums[np.arange(u)[None, :], pick].sum(1)
            else:                      # data: (unit x sample values, labels)
                values, labels = data; perm = np.argsort(rng.random((b, u, len(labels))), axis=2)
                shuffled = np.take_along_axis(np.broadcast_to(values, (b,)+values.shape), perm, axis=2)
                out[begin:begin+b] += np.einsum("bun,nk->bk", shuffled, np.eye(K)[labels])
    return out


def anova(values, conditions, genes, gene_indices, libraries=None, settings=None):
    """values: sites x samples normalized counts; conditions / libraries: one label per sample column.
    Returns (rows, diagnostics); rows follow `genes`."""
    settings = settings or Settings(); settings.validate()
    values = np.asarray(values, dtype=float); conditions = list(conditions)
    libraries = ["library"]*values.shape[1] if libraries is None else list(libraries)
    if values.ndim != 2 or len(conditions) != values.shape[1] or len(libraries) != values.shape[1]:
        raise ValueError("Need sites x samples values with one condition and one library label per sample")
    if np.any(~np.isfinite(values)) or np.any(values < 0):raise ValueError("Normalized counts must be finite and nonnegative")
    if len(genes) != len(gene_indices) or len({g["id"] for g in genes}) != len(genes):raise ValueError("Gene identities and site lists must align and be unique")
    names = list(dict.fromkeys(conditions)); K = len(names)
    reference = settings.reference or names[0]
    if K < 2 or reference not in names:raise ValueError("Need >= 2 conditions and a reference among them")
    code = np.array([names.index(c) for c in conditions])
    strata = []
    for library in dict.fromkeys(libraries):
        cols = [j for j, x in enumerate(libraries) if x == library]
        if len(set(code[cols])) < 2:raise ValueError(f"Library {library!r} needs samples from >= 2 conditions")
        strata.append(dict(library=library, cols=np.array(cols), labels=code[cols]))
    values, by_code = recenter(values, code, names.index(reference), gene_indices, settings.recenter, settings.recenter_trim)
    factors = {names[k]:f for k, f in by_code.items()}
    for st in strata:
        st["eligible"] = values[:, st["cols"]].sum(1) > 0
        st["arrangements"] = assignments(st["labels"], settings.max_enumerate)
    rows = []
    for gene, idx in zip(genes, gene_indices):
        idx = np.asarray(idx, dtype=np.int64); blocks, S, N, n = [], np.zeros(K), np.zeros(K), 0
        per_condition = [[] for _ in range(K)]
        for st in strata:
            i = idx[st["eligible"][idx]]
            if not len(i):continue
            V = values[np.ix_(i, st["cols"])]; n += len(i)
            for k in range(K):
                m = st["labels"] == k; S[k] += V[:, m].sum(); N[k] += len(i)*m.sum(); per_condition[k].append(V[:, m].ravel())
            if st["arrangements"] is not None:
                A = st["arrangements"]; onehot = (A[:, :, None] == np.arange(K)).astype(float)   # A x n x K
                blocks.append(("enumerated", (np.einsum("un,ank->uak", V, onehot),)))
            else:blocks.append(("sampled", (V, st["labels"])))
        row = dict(ORF=gene["id"], Name=gene["name"], chromosome=gene.get("chrom"), n_sites=n)
        if n < settings.min_sites:
            rows.append(row | dict(status="too few eligible sites" if n else "no eligible sites", in_family=False)); continue
        means = {name:float(np.concatenate(per_condition[k]).mean()) if per_condition[k] else None for k, name in enumerate(names)}
        ref_mean = means[reference]
        for name in names:
            row[f"mean_{name}"] = means[name]
            if name != reference:
                row[f"log2FC_{name}"] = (float(np.log2((means[name]+settings.pseudocount)/(ref_mean+settings.pseudocount)))
                                         if means[name] is not None and ref_mean is not None and (settings.pseudocount > 0 or min(means[name], ref_mean) > 0) else None)
        t_obs = float(_statistic(S, N)); tolerance = 1e-9*max(1., abs(t_obs))
        for stage, B in enumerate(settings.stages):
            draws = _draws(_gene_rng(settings.seed, gene["id"], stage), blocks, K, int(B))
            p = float((np.sum(_statistic(draws, N) >= t_obs-tolerance)+1)/(B+1))
            if stage == len(settings.thresholds) or p > settings.thresholds[stage]:break
        lfc = {name:row.get(f"log2FC_{name}") for name in names if name != reference}
        top = max((x for x in lfc if lfc[x] is not None), key=lambda x:abs(lfc[x]), default=None)
        rows.append(row | dict(statistic=t_obs, p_value=p, permutations=int(B), top_condition=top, status="tested", in_family=True))
    family = [r for r in rows if r["in_family"]]
    for row, q in zip(family, benjamini_hochberg([r["p_value"] for r in family]) if family else []):row["q_bh"] = float(q)
    for row in rows:row["call"] = "not tested" if not row["in_family"] else "varies" if row["q_bh"] <= settings.fdr else "not significant"
    spreads = None
    if len(strata) == 1:
        spreads = {name:resampling.replicate_spread(values[:, code == k]) for k, name in enumerate(names)}
    measured = [s for s in (spreads or {}).values() if s]
    asymmetry = float(max(measured)/min(measured)) if len(measured) >= 2 else None
    diagnostics = dict(model=MODEL, revision=REVISION, settings=asdict(settings), conditions=names, reference=reference,
                       samples_per_condition={name:int((code == k).sum()) for k, name in enumerate(names)},
                       libraries=[st["library"] for st in strata], recenter_factors=factors,
                       enumerated_libraries=[st["library"] for st in strata if st["arrangements"] is not None],
                       replicate_spread=spreads, noise_asymmetry=asymmetry,
                       warning=(f"replicate noise differs {asymmetry:.1f}x between conditions; p-values may be mildly anti-conservative"
                                if asymmetry and asymmetry > settings.asymmetry_warning else None),
                       family_size=len(family), calls={c:sum(r["call"] == c for r in rows) for c in ("varies", "not significant", "not tested")})
    return rows, diagnostics


def all_pairs(values, conditions, genes, gene_indices, positions, libraries=None, settings=None, fdr=.05,
              recenter_method="iterative-median", reference=""):
    """Resampling v3 for every pair of conditions (first of the pair = control), within library.
    recenter_method "iterative-median" (v2): all conditions recentered to `reference` (default: the first
    condition) once, then resampling without recentering; "v3" (v1): resampling's own per-pair recentering;
    "none": no recentering.
    Returns (long rows with a 'pair' column, per-pair diagnostics); q_all_pairs is BH over all pairs' tests."""
    values = np.asarray(values, dtype=float); conditions = list(conditions)
    libraries = ["library"]*values.shape[1] if libraries is None else list(libraries)
    names = list(dict.fromkeys(conditions)); rows, diagnostics = [], {}
    if recenter_method == "v3":settings = settings or resampling.Settings(recenter="trimmed-sum")
    elif recenter_method == "none":settings = settings or resampling.Settings(recenter="none")
    elif recenter_method == "iterative-median":
        code = np.array([names.index(c) for c in conditions])
        values, by_code = recenter(values, code, names.index(reference or names[0]), gene_indices, "iterative-median")
        diagnostics["recenter_factors"] = {names[k]:f for k, f in by_code.items()}
        settings = settings or resampling.Settings(recenter="none")
        if settings.recenter != "none":raise ValueError("iterative-median all_pairs needs resampling settings with recenter='none'")
    else:raise ValueError("recenter_method must be 'iterative-median', 'v3' or 'none'")
    for a, b in combinations(names, 2):
        libs = [x for x in dict.fromkeys(libraries) if any(c == a and l == x for c, l in zip(conditions, libraries))
                and any(c == b and l == x for c, l in zip(conditions, libraries))]
        if not libs:diagnostics[f"{b} vs {a}"] = "no library contains both conditions"; continue
        ca = [j for j, (c, l) in enumerate(zip(conditions, libraries)) if c == a and l in libs]
        cb = [j for j, (c, l) in enumerate(zip(conditions, libraries)) if c == b and l in libs]
        stratified = len(libs) > 1
        r, d = resampling.analyze(values[:, ca], values[:, cb], genes, gene_indices, positions, settings,
                                  [libraries[j] for j in ca] if stratified else None, [libraries[j] for j in cb] if stratified else None)
        diagnostics[f"{b} vs {a}"] = d
        rows += [x | dict(pair=f"{b} vs {a}", control=a, treatment=b) for x in r]
    family = [r for r in rows if r.get("in_family")]
    for row, q in zip(family, benjamini_hochberg([r["p_value"] for r in family]) if family else []):row["q_all_pairs"] = float(q)
    for row in rows:
        row["call_all_pairs"] = ("not tested" if not row.get("in_family") else "not significant" if row["q_all_pairs"] > fdr else
                                 "depleted" if row["delta_mean"] < 0 else "enriched" if row["delta_mean"] > 0 else "not significant")
    return rows, diagnostics


def gene_sample_sums(values, gene_indices, libraries=None):
    """Per gene and sample: sum of normalized values over the gene's sites eligible in that sample's library,
    and the number of those sites (genes x samples each). Any grouping's site x sample mean is sum(S)/sum(N)."""
    values = np.asarray(values, dtype=float)
    libraries = ["library"]*values.shape[1] if libraries is None else list(libraries)
    S = np.zeros((len(gene_indices), values.shape[1])); N = np.zeros_like(S)
    for library in dict.fromkeys(libraries):
        cols = np.array([j for j, x in enumerate(libraries) if x == library]); eligible = values[:, cols].sum(1) > 0
        for g, idx in enumerate(gene_indices):
            idx = np.asarray(idx, dtype=np.int64); i = idx[eligible[idx]]
            if len(i):S[g, cols] = values[np.ix_(i, cols)].sum(0); N[g, cols] = len(i)
    return S, N


def run(normalized, gff, reference, output, samples, settings=None, normalization=None, stratified=False, tests=("anova", "all-pairs"),
        independent=False):
    """File-based runner on a YTAB normalization output (BGC recommended; ttr, totreads, nzmean and aBGC accepted).

    samples: (observation_id, condition) pairs; settings.reference names the reference condition (default:
    the first condition listed). stratified=True: several libraries (pools), each split into the reference
    and >= 1 other condition; labels are permuted only within library. independent=True: each condition may be its own,
    independently made library (strain screens such as SATAY mutants); labels are permuted among all samples within each
    site (unstratified). Writes anova_genes.tsv and/or
    pairs_genes.tsv and run_manifest.json to a new output directory.
    """
    settings = settings or Settings(); settings.validate()
    tests = tuple(tests)
    if not tests or any(t not in ("anova", "all-pairs") for t in tests):raise ValueError("tests must be 'anova' and/or 'all-pairs'")
    samples = [(str(o), str(c)) for o, c in samples]
    if len({o for o, _ in samples}) != len(samples):raise ValueError("Each observation may be listed once")
    names = list(dict.fromkeys(c for _, c in samples)); ref = settings.reference or names[0]
    if len(names) < 2 or ref not in names:raise ValueError("Assign samples to >= 2 conditions, including the reference condition")
    condition = dict(samples)
    controls = [o for o, c in samples if c == ref]; treatments = [o for o, c in samples if c != ref]
    if stratified and independent:raise ValueError("Choose stratified (pools split into conditions) or independent libraries, not both")
    design = load_design(normalized, normalization, reference, output, controls, treatments, min_controls=1, stratified=stratified,
                         nonlinear=True, check_roles=False, independent=independent)
    output, sites, values = design["output"], design["sites"], design["scaled"]
    conditions = [condition[x] for x in design["selected"]]; libraries = design["libraries"] if stratified else None
    settings = Settings(**(asdict(settings) | dict(reference=ref)))
    genes = annotations(Path(gff).resolve()); idx = gene_sites(genes, sites, settings)
    positions = np.array([p for _, p in sites], dtype=np.int64)
    output.mkdir(parents=True, exist_ok=True); diagnostics, files = {}, {}
    code = np.array([names.index(c) for c in conditions])
    centred, _ = recenter(values, code, names.index(ref), idx, settings.recenter, settings.recenter_trim)
    S, N = gene_sample_sums(centred, idx, design["libraries"] if stratified else None)
    write_table(output/"gene_sample_sums.tsv", ["ORF", "Name", "observation_id", "condition", "library", "eligible_sites", "sum"],
                [dict(ORF=g["id"], Name=g["name"], observation_id=x, condition=condition[x], library=l, eligible_sites=int(N[i, j]), sum=float(S[i, j]))
                 for i, g in enumerate(genes) for j, (x, l) in enumerate(zip(design["selected"], design["libraries"]))])
    files["gene_sample_sums"] = "gene_sample_sums.tsv"
    if "anova" in tests:
        rows, diagnostics["anova"] = anova(values, conditions, genes, idx, libraries, settings)
        lead = ["ORF", "Name", "n_sites", "statistic", "p_value", "q_bh", "call", "top_condition", "permutations"]
        write_table(output/"anova_genes.tsv", lead+list(dict.fromkeys(k for r in rows for k in r if k not in lead)), rows)
        files["anova"] = "anova_genes.tsv"
    if "all-pairs" in tests:
        method = {"iterative-median":"iterative-median", "trimmed-sum":"v3", "none":"none"}[settings.recenter]
        pair_settings = resampling.Settings(recenter="trimmed-sum" if method == "v3" else "none", fdr=settings.fdr, min_sites=settings.min_sites,
                                            pseudocount=settings.pseudocount, nterm=settings.nterm, cterm=settings.cterm, ignore_stop=settings.ignore_stop)
        rows, pdiag = all_pairs(values, conditions, genes, idx, positions, libraries, pair_settings, settings.fdr, method, ref)
        lead = ["pair", "control", "treatment", "ORF", "Name", "n_sites", "log2FC", "delta_mean", "p_value", "q_bh", "q_all_pairs", "call", "call_all_pairs"]
        write_table(output/"pairs_genes.tsv", lead+list(dict.fromkeys(k for r in rows for k in r if k not in lead)), rows)
        diagnostics["all_pairs"] = dict(recenter_factors=pdiag.get("recenter_factors")) | {
            pair:(d if isinstance(d, str) else dict(calls=d["calls"], family_size=d["family_size"], treatment_rescale=d["treatment_rescale"],
                                                    replicate_noise=d["replicate_noise"]))
            for pair, d in pdiag.items() if pair != "recenter_factors"}
        files["all_pairs"] = "pairs_genes.tsv"
    warning = (diagnostics.get("anova") or {}).get("warning")
    norm = design["norm"]
    manifest = dict(status="success", analysis=MODEL, algorithm_revision=REVISION, parameters=asdict(settings), tests=list(tests), files=files,
        samples=[dict(observation_id=x, condition=condition[x], library=l) for x, l in zip(design["selected"], design["libraries"])],
        conditions=names, reference_condition=ref, stratified=stratified, independent_libraries=independent,
        normalization=dict(method=norm["method"], input=str(norm["folder"]), catalogue=str(design["catalogue"]),
            hashes={name:sha256(norm["folder"]/name) for name in ("normalized_counts.npz", "factors.tsv", "method.json")}
                   | dict(run_manifest=sha256(norm["root"]/"run_manifest.json"))),
        recommendation=("Replicate noise differs between conditions: ANOVA p-values are anti-conservative here (about 10% instead of 5% "
                        "in the pre-registered evaluation). Use the all-pairs results." if warning else
                        "ANOVA answers 'does the gene respond to any condition'; all-pairs answers 'which conditions differ, and in which direction'."),
        validation="output/essentiality/multicondition/benchmark_v2/RESULT.md (pre-registered; passed)",
        diagnostics=diagnostics, implementation_sha256=sha256(__file__),
        limitations=["Unequal replicate noise between conditions makes ANOVA p-values anti-conservative (warning issued); all-pairs stays near nominal",
                     "Sites are treated as independent within a gene; sample-level effects shared across sites are not modelled",
                     "Recentering assumes most genes do not respond to each condition",
                     "Stratified designs need every library to contain the reference condition and >= 1 other condition",
                     "ANOVA p-values are Monte Carlo (seeded, staged; minimum about 1e-5)"])
    (output/"run_manifest.json").write_text(json.dumps(_json_ready(manifest), indent=2, allow_nan=False)+"\n")
    return manifest
