"""Pool-stratified insertion-positive mutant planning from audited catalogues.

Each draw is an insertion-positive mutant pick from the pool's full observed
genomic insertion-coordinate catalogue. Coding coordinates can recover genes;
noncoding/unassigned coordinates still consume picks. Coordinates identify
alleles globally, while probabilities are normalized within each independent
pool. No essentiality classifier is read.
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import yaml
from scipy.optimize import minimize

from .library import is_parent_sample

LABELS = {"equal": "Hermes Empirical Monte Carlo — Equal Insertion",
          "abundance": "Hermes Empirical Monte Carlo — Observed-Abundance Sensitivity"}
STRATEGIES = {"equal": "Equal physical-pool allocation baseline",
    "complexity": "Complexity-normalized allocation",
    "balanced": "Balanced saturation-aware allocation",
    "insertion_proportional": "Total-site-proportional allocation",
    "site_proportional": "Total-site-proportional allocation",
    "coding_proportional": "Coding-diversity-proportional allocation",
    "optimized": "Maximize immediate ≥1-gene recovery",
    "manual": "Manual pool counts"}


def pool_identity(row, library_namespace):
    """Identify an independently constructed pool, not a sequencing replicate.

    Explicit physical pool metadata takes precedence; otherwise the pool
    number is namespaced by library and species. Background/replicate labels
    remain provenance and do not split the same numbered pool.
    """
    number = str(row.get("pool_id") or row.get("pool") or "").strip()
    provenance = str(row.get("parent_lineage") or row.get("background") or "").strip()
    species = str(row.get("species") or "").strip()
    explicit = str(row.get("physical_pool_id") or "").strip()
    if not species or not number:
        raise ValueError("Pool identity needs species and pool number")
    namespace = str(row.get("library_id") or library_namespace)
    return "::".join((namespace, species, explicit or ("pool" + number)))


def load_pool_catalogue(audit_dir, repo_root, selected_sources=None, selected_pools=None):
    audit_dir, repo_root = Path(audit_dir), Path(repo_root)
    audit_path = audit_dir / "parent_union_audit.json"
    audit = json.loads(audit_path.read_text())
    inventory_path = audit_dir / "parent_union_gene_inventory.csv"
    site_path = audit_dir / "parent_union_site_provenance.csv"
    with inventory_path.open() as handle:
        genes = [g for g in csv.DictReader(handle) if g["biotype"] == "protein_coding"]
    if len(genes) != audit["gene_annotation"]["protein_coding_genes"]:
        raise ValueError("Protein-coding inventory disagrees with accepted audit")
    all_genes = {g["gene_id"]: g for g in genes}
    reference = Path(audit["reference"])
    reference_hash = hashlib.sha256(reference.read_bytes()).hexdigest()
    allowed = {f"{r['project']}::{r['sample']}" for r in audit["samples"]}
    chosen = allowed if selected_sources is None else set(selected_sources)
    if not chosen or not chosen <= allowed:
        raise ValueError("Select eligible sources in the audited parent catalogue")
    metadata, files, sheets = {}, [audit_path, inventory_path, site_path], {}
    def resolve(path):
        p = Path(path)
        return p if p.is_absolute() else repo_root / p
    for observation in audit["samples"]:
        source = observation["project"] + "::" + observation["sample"]
        if source not in chosen:
            continue
        project = observation["project"]
        if project not in sheets:
            config = repo_root / "output/projects" / project / "config/project.yaml"
            cfg = yaml.safe_load(config.read_text())
            sheet = resolve(cfg["sample_sheet"])
            with sheet.open() as handle:
                sheets[project] = {r["sample"]: r for r in csv.DictReader(handle)}
            files.extend([config, sheet])
            ref = resolve(cfg["reference"]["gff"])
            if hashlib.sha256(ref.read_bytes()).hexdigest() != reference_hash:
                raise ValueError("Parent observations must use the audited mapping reference")
        row = sheets[project].get(observation["sample"])
        if row is None or not is_parent_sample(row) or str(row.get("include", "true")).lower() not in {"true", "yes", "1"}:
            raise ValueError("Treated, excluded or missing source in parent catalogue: " + source)
        if str(row.get("pool_id") or row.get("pool")) != str(observation["pool"]):
            raise ValueError("Physical pool metadata changed since accepted audit")
        identity = pool_identity(row, reference_hash[:12])
        metadata[source] = dict(source=source, project=project, sample=row["sample"],
            pool_id=identity, pool_number=str(row.get("pool_id") or row["pool"]),
            parent_replicate_provenance=row.get("background", ""),
            species=row["species"])
    pool_ids = sorted({r["pool_id"] for r in metadata.values()})
    if selected_pools is not None:
        if not set(selected_pools) <= set(pool_ids) or not selected_pools:
            raise ValueError("Unknown or empty physical-pool selection")
        pool_ids = [p for p in pool_ids if p in selected_pools]
        metadata = {s:r for s,r in metadata.items() if r["pool_id"] in pool_ids}
    sources = sorted(metadata)
    coordinates, gene_ids, reads = [], [], []
    seen = set()
    with site_path.open() as handle:
        for row in csv.DictReader(handle):
            coordinate = (row["chromosome"], int(row["position"]))
            if coordinate in seen:
                raise ValueError("Duplicate global coordinate in accepted audit")
            seen.add(coordinate)
            support = json.loads(row["source_read_support_json"])
            if not set(support) <= allowed:
                raise ValueError("Unaudited source found in coordinate provenance")
            values = [float(support.get(s, 0)) for s in sources]
            if not any(values):
                continue
            coordinates.append(coordinate)
            gene_ids.append(row["unambiguous_CDS_locus"] if row["unambiguous_CDS_locus"] in all_genes else "")
            reads.append(values)
    if len(seen) != audit["union"]["unique_sites"]:
        raise ValueError("Coordinate count disagrees with audit")
    # Selecting fewer pools/observations must not silently redefine the target
    # universe. Retain every coding-represented gene in the accepted global audit;
    # inaccessible genes have p=0 and lower the selected-pool recovery ceiling.
    represented = sorted(g for g,row in all_genes.items() if int(row["union_distinct_CDS_sites"])>0)
    index = {g:i for i,g in enumerate(represented)}
    gids = np.array([index.get(g, -1) for g in gene_ids], dtype=int)
    reads = np.asarray(reads, dtype=float)
    if np.any(~np.isfinite(reads)) or np.any(reads < 0) or np.any(reads.sum(axis=0) <= 0):
        raise ValueError("Invalid sequencing support")
    presence, probabilities = pool_distributions(reads, [metadata[s]["pool_id"] for s in sources], pool_ids)
    n = len(represented)
    coding = gids >= 0
    gene_presence = np.zeros((n, len(pool_ids)), dtype=bool)
    site_counts = np.zeros((n, len(pool_ids)), dtype=int)
    gene_probabilities = {}
    for model, values in probabilities.items():
        gene_probabilities[model] = np.stack([np.bincount(gids[coding], weights=values[coding,j], minlength=n)
                                            for j in range(len(pool_ids))], axis=1)
    for j in range(len(pool_ids)):
        site_counts[:,j] = np.bincount(gids[coding & presence[:,j]], minlength=n)
        gene_presence[:,j] = site_counts[:,j] > 0
    inventory = []
    for j,pool in enumerate(pool_ids):
        obs = [metadata[s] for s in sources if metadata[s]["pool_id"] == pool]
        inventory.append(dict(pool_id=pool, label="Pool " + obs[0]["pool_number"],
            provenance_confidence="Investigator-designated independently constructed insertion pool",
            parent_replicate_provenance=obs[0]["parent_replicate_provenance"],
            observations=[x["source"] for x in obs], observation_count=len(obs),
            sites=int(presence[:,j].sum()), coding_sites=int((coding & presence[:,j]).sum()),
            noncoding_sites=int((~coding & presence[:,j]).sum()),
            coding_fraction=float((coding & presence[:,j]).sum()/presence[:,j].sum()),
            genes=int(gene_presence[:,j].sum()), unique_genes=int((gene_presence[:,j] & (gene_presence.sum(axis=1)==1)).sum()),
            unique_coding_alleles=int((coding & presence[:,j] & (presence.sum(axis=1)==1)).sum()),
            shared_coding_alleles=int((coding & presence[:,j] & (presence.sum(axis=1)>1)).sum()),
            ge2=int((site_counts[:,j]>=2).sum()), ge3=int((site_counts[:,j]>=3).sum()),
            accessible_fraction=float(gene_presence[:,j].mean())))
    overlaps = []
    for i in range(len(pool_ids)):
        for j in range(i+1,len(pool_ids)):
            sites = presence[:,i] & presence[:,j]; genes_both = gene_presence[:,i] & gene_presence[:,j]
            overlaps.append(dict(pool_a=pool_ids[i], pool_b=pool_ids[j], shared_sites=int(sites.sum()),
                site_jaccard=float(sites.sum()/(presence[:,i]|presence[:,j]).sum()),
                shared_coding_sites=int((sites & coding).sum()), shared_genes=int(genes_both.sum()),
                gene_jaccard=float(genes_both.sum()/(gene_presence[:,i]|gene_presence[:,j]).sum())))
    global_counts = np.array([int(all_genes[g]["union_distinct_CDS_sites"]) for g in represented])
    summary = dict(annotated=len(all_genes), represented=n, sites=len(seen), coding_sites=int(global_counts.sum()),
        ge2=int((global_counts>=2).sum()), ge3=int((global_counts>=3).sum()), missing=len(all_genes)-n,
        reference_representation=n/len(all_genes),selected_sites=len(gids),selected_coding_sites=int(coding.sum()),
        accessible_genes=int(np.any(gene_presence,axis=1).sum()),selected_pool_ceiling=float(np.any(gene_presence,axis=1).mean()))
    if chosen == allowed and selected_pools is None:
        expected = audit["union"]
        for field,key in [("sites","unique_sites"),("coding_sites","unambiguous_protein_coding_sites"),
                          ("represented","protein_genes_with_coding_insertion"),("ge2","protein_genes_with_ge2_sites"),("ge3","protein_genes_with_ge3_sites")]:
            if summary[field] != expected[key]:
                raise ValueError("Pool-aware catalogue changed accepted global count: " + field)
    return dict(summary=summary, pools=inventory, overlaps=overlaps, observations=[metadata[s] for s in sources],
        pool_ids=pool_ids, gene_ids=represented, gene_pool_presence=gene_presence,
        gene_lengths=np.array([int(all_genes[g]["coding_length"]) for g in represented],dtype=int),
        gene_coordinate_counts=site_counts,
        coordinates=coordinates, allele_genes=gids, presence=presence, probabilities=probabilities,
        gene_probabilities=gene_probabilities, observation_reads=reads,
        sources=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files])


def pool_distributions(reads, observation_pools, pool_ids=None):
    """Build pick probabilities over each pool's full observed catalogue.

    Parent observations are normalized individually before pooling. Noncoding
    and unassigned insertion coordinates remain in each pool's sampling
    universe.
    """
    reads = np.asarray(reads, dtype=float)
    if reads.ndim != 2 or reads.shape[1] != len(observation_pools) or np.any(~np.isfinite(reads)) or np.any(reads<0) or np.any(reads.sum(axis=0)<=0):
        raise ValueError("Invalid parent observations")
    pools = sorted(set(observation_pools)) if pool_ids is None else pool_ids
    q = reads / reads.sum(axis=0)
    membership = np.stack([np.any(reads[:,np.array(observation_pools)==p]>0,axis=1) for p in pools],axis=1)
    counts = membership.sum(axis=0)
    if np.any(counts == 0):
        raise ValueError("Each selected pool must contain an observed insertion coordinate")
    equal = membership / counts
    abundance = np.stack([q[:,np.array(observation_pools)==p].mean(axis=1) for p in pools],axis=1)
    abundance /= abundance.sum(axis=0)
    return membership, dict(equal=equal, abundance=abundance)


def _probabilities(p):
    p = np.asarray(p, dtype=float)
    if p.ndim != 2 or not p.shape[0] or not p.shape[1] or np.any(~np.isfinite(p)) or np.any(p<0) or np.any(p>1) or np.any(p.sum(axis=0)>1+1e-10):
        raise ValueError("Invalid per-pool gene probabilities")
    return p


def expected_recovery(p, allocation):
    """Exact independent-pool binomial marginals; NOT a Poisson model."""
    p = _probabilities(p); k = np.asarray(allocation, dtype=float)
    if k.shape != (p.shape[1],) or np.any(~np.isfinite(k)) or np.any(k<0):
        raise ValueError("Invalid allocation")
    # Treat p=1 separately to avoid 0 * -inf for unsampled pools.
    log_miss = np.log1p(-np.minimum(p, 1-np.finfo(float).eps)) @ k
    recovery = -np.expm1(log_miss)
    recovery[np.any((p==1) & (k>0)[None,:],axis=1)] = 1
    return recovery


def integer_allocation(total, weights, minimum=None, maximum=None):
    if int(total)!=total or total<0:
        raise ValueError("Total colony budget must be a nonnegative integer")
    weights=np.asarray(weights,dtype=float); n=len(weights)
    low=np.zeros(n,dtype=int) if minimum is None else np.asarray(minimum,dtype=int)
    high=np.full(n,int(total),dtype=int) if maximum is None else np.asarray(maximum,dtype=int)
    if low.shape!=(n,) or high.shape!=(n,) or np.any(low<0) or np.any(high<low) or low.sum()>total or high.sum()<total:
        raise ValueError("Infeasible pool constraints")
    if np.any(~np.isfinite(weights)) or np.any(weights<0) or not weights.sum():
        raise ValueError("Invalid allocation weights")
    # Continuous bounded apportionment, then largest-remainder integer rounding.
    # Equal allocation differs by at most one colony when K is not divisible by n.
    desired=low.astype(float);remaining=float(total-desired.sum())
    while remaining>1e-7:
        room=high-desired;active=room>1e-9;w=np.where(active,weights,0.)
        if not w.sum():w=active.astype(float)
        shares=remaining*w/w.sum();desired+=np.minimum(room,shares)
        remaining=max(0.,float(total-desired.sum()))
    out=np.minimum(high,np.floor(desired+1e-8).astype(int));remaining=int(total-out.sum())
    while remaining>0:
        order=np.argsort(-np.where(out<high,desired-out,-np.inf),kind="stable")
        for j in order:
            if out[j]<high[j] and remaining:out[j]+=1;remaining-=1
    return out


def optimize_allocation(p,total,site_counts=None,coding_counts=None,minimum=None,maximum=None):
    """Concave expected-gene objective; continuous solution + integer exchanges.

    Keep baseline candidates as a numerical safety gate. No assertion of a
    globally optimal integer solution: return relaxation bound and gap.
    """
    p=_probabilities(p); n=p.shape[1]
    weights={"equal":np.ones(n),"insertion_proportional":np.ones(n) if site_counts is None else site_counts,
             "coding_proportional":np.ones(n) if coding_counts is None else coding_counts}
    baseline={name:integer_allocation(total,w,minimum,maximum) for name,w in weights.items()}
    score=lambda k:float(expected_recovery(p,k).sum())
    if total==0:return dict(allocation=baseline["equal"],expected=0.,baselines=baseline,relaxation_gap=0.,converged=True)
    low=np.zeros(n) if minimum is None else np.asarray(minimum,dtype=float)
    high=np.full(n,total) if maximum is None else np.asarray(maximum,dtype=float)
    a=np.log1p(-np.minimum(p,1-np.finfo(float).eps))
    def objective(x):
        miss=np.exp(a@(x*total))
        return float(miss.mean()),total*(a.T@miss)/len(p)
    best=max(baseline.values(),key=score)
    fit=minimize(objective,best/total,jac=True,method="SLSQP",bounds=list(zip(low/total,high/total)),
        constraints=[dict(type="eq",fun=lambda x:x.sum()-1,jac=lambda x:np.ones(n))],options=dict(ftol=1e-12,maxiter=300))
    candidate=np.maximum(low,np.minimum(high,np.floor(np.maximum(0,fit.x)*total))).astype(int)
    while candidate.sum()<total:
        gains=[score(candidate+np.eye(n,dtype=int)[j])-score(candidate) if candidate[j]<high[j] else -np.inf for j in range(n)]
        candidate[np.argmax(gains)]+=1
    while candidate.sum()>total:
        losses=[score(candidate)-score(candidate-np.eye(n,dtype=int)[j]) if candidate[j]>low[j] else np.inf for j in range(n)]
        candidate[np.argmin(losses)]-=1
    if score(candidate)>score(best):best=candidate
    # Discrete one-colony exchange refinement, not arbitrary weighted objectives.
    for _ in range(100):
        current=score(best);move=None;gain=1e-9
        for i in range(n):
            for j in range(n):
                if i==j or best[i]<=low[i] or best[j]>=high[j]:continue
                trial=best.copy();trial[i]-=1;trial[j]+=1;delta=score(trial)-current
                if delta>gain:move=trial;gain=delta
        if move is None:break
        best=move
    # A supporting hyperplane gives a valid global upper bound on the concave
    # objective, including when the solver reports a numerical warning.
    miss=1-expected_recovery(p,best);gradient=-(a.T@miss)
    bound_alloc=low.copy();remaining=total-low.sum()
    for j in np.argsort(-gradient):
        extra=min(remaining,high[j]-low[j]);bound_alloc[j]+=extra;remaining-=extra
    gap=max(0.,float(gradient@(bound_alloc-best)))
    return dict(allocation=best,expected=score(best),baselines=baseline,relaxation_gap=gap,converged=bool(fit.success))


def target_allocation(p,target,site_counts,coding_counts,minimum=None,maximum=None,cap=10**10):
    """Integer-budget search on optimized exact expectation, for MC bracketing."""
    if not 0<target<1:raise ValueError("Target must be between 0 and 1")
    if np.mean(np.any(p>0,axis=1))<target:return None
    low=0 if minimum is None else int(np.sum(minimum))
    high=max(1,low,len(p))
    if maximum is not None:cap=min(cap,int(np.sum(maximum)))
    def fit(k):return optimize_allocation(p,k,site_counts,coding_counts,minimum,maximum)
    high=min(high,cap)
    while fit(high)["expected"]/len(p)<target and high<cap:high=min(cap,high*2)
    if fit(high)["expected"]/len(p)<target:return None
    while low<high:
        mid=(low+high)//2
        if fit(mid)["expected"]/len(p)>=target:high=mid
        else:low=mid+1
    result=fit(low);result["k"]=low
    return result


def equal_allocation_target(p,target,minimum=None,maximum=None,cap=10**10):
    """First exact-mean crossing under equal allocation, not optimization.

    Equal allocation's deterministic remainder gives nested per-pool depths.
    Exact expectations bracket a separately reported Monte Carlo mean crossing.
    """
    p=_probabilities(p)
    if not 0<target<1:raise ValueError("Target must be between 0 and 1")
    if np.any(p>0,axis=1).mean()<target:return None
    low=0 if minimum is None else int(np.sum(minimum))
    if maximum is not None:cap=min(cap,int(np.sum(maximum)))
    if low>cap:return None
    def allocation(k):return integer_allocation(k,np.ones(p.shape[1]),minimum,maximum)
    def mean(k):return expected_recovery(p,allocation(k)).mean()
    high=min(cap,max(1,low,len(p)))
    while mean(high)<target and high<cap:high=min(cap,2*high)
    if mean(high)<target:return None
    while low<high:
        mid=(low+high)//2
        if mean(mid)>=target:high=mid
        else:low=mid+1
    return dict(k=low,allocation=allocation(low),expected=float(mean(low)*len(p)))


def complexity_allocation(total,site_counts,minimum=None,maximum=None):
    """Allocate K in proportion to full observed genomic insertion complexity."""
    weights=np.asarray(site_counts,dtype=float)
    if weights.ndim!=1 or len(weights)==0 or np.any(weights<=0):
        raise ValueError("Each selected pool must have positive observed insertion complexity")
    return integer_allocation(total,weights,minimum,maximum)


def pool_allocation(strategy,total,pool_site_counts,p,site_probabilities=None,
                    pool_gene_presence=None,minimum=None,maximum=None):
    """Resolve one of the documented allocation strategies for arbitrary J."""
    n_pool=len(pool_site_counts)
    if strategy=="equal":return integer_allocation(total,np.ones(n_pool),minimum,maximum)
    if strategy=="complexity":return complexity_allocation(total,pool_site_counts,minimum,maximum)
    if strategy=="balanced":
        if pool_gene_presence is None:raise ValueError("Balanced allocation requires pool-specific represented genes")
        return balanced_saturation_allocation(total,p,pool_gene_presence,minimum,maximum)
    if strategy=="optimized":return optimize_allocation(p,total,pool_site_counts,None,minimum,maximum)["allocation"]
    if strategy=="coding_proportional":
        raise ValueError("Coding-proportional allocation requires a separate coding-site weight vector")
    if strategy=="manual":raise ValueError("Manual allocation is supplied directly")
    raise ValueError("Unknown allocation strategy: "+str(strategy))


def pool_gene_saturation(p,k,pool_gene_presence=None):
    """Expected within-pool represented-gene fraction at depth k.

    Denominator is each pool's own represented-gene repertoire N_j. Input
    probabilities are per insertion-positive pick from the full pool catalogue;
    remaining probability mass represents noncoding/unassigned coordinates.
    """
    p=_probabilities(p)
    depths=np.asarray(k,dtype=float)
    if p.ndim!=2 or np.any(~np.isfinite(depths)) or np.any(depths<0) or (depths.ndim>0 and depths.shape!=(p.shape[1],)):
        raise ValueError("Invalid pool saturation input")
    mask=(p>0) if pool_gene_presence is None else np.asarray(pool_gene_presence,dtype=bool)
    if mask.shape!=p.shape:raise ValueError("Pool presence mask must match gene probabilities")
    fractions=[]
    for j in range(p.shape[1]):
        eligible=mask[:,j];nj=int(eligible.sum())
        if nj==0:fractions.append(0.);continue
        values=p[eligible,j]
        depth=float(depths if depths.ndim==0 else depths[j])
        fractions.append(0. if depth==0 else float(np.mean(-np.expm1(int(depth)*np.log1p(-np.minimum(values,1.))))))
    return np.asarray(fractions)


def balanced_saturation_allocation(total,p,pool_gene_presence,minimum=None,maximum=None,
                                   complexity_weights=None):
    """Maximize the minimum pool saturation, with deterministic residual use.

    Maximizes min_j S_j(K_j) under integer K_j and sum K_j=K. Every selected
    pool receives at least one colony. After the maximum achievable minimum
    saturation is reached, any residual colonies that cannot raise that
    minimum are apportioned by largest remainder using observed complexity
    weights (or equally if none are supplied). The primary objective is
    max-min saturation, not global gene recovery.
    """
    p=_probabilities(p);j=p.shape[1]
    if int(total)!=total or total<j:raise ValueError("Balanced allocation requires at least one colony per selected pool")
    lower=np.ones(j,dtype=int) if minimum is None else np.maximum(1,np.asarray(minimum,dtype=int))
    upper=np.full(j,total,dtype=int) if maximum is None else np.asarray(maximum,dtype=int)
    if lower.shape!=(j,) or upper.shape!=(j,) or lower.sum()>total or upper.sum()<total or np.any(lower>upper):
        raise ValueError("Balanced allocation constraints are infeasible")
    presence=np.asarray(pool_gene_presence,dtype=bool)
    if presence.shape!=p.shape:raise ValueError("Pool gene repertoire mask must match probabilities")
    log_miss=[np.log1p(-np.minimum(p[presence[:,pool],pool],1.)) for pool in range(j)]
    def saturation(pool,depth):
        return 0. if depth<=0 else float(np.mean(-np.expm1(int(depth)*log_miss[pool])))
    ceiling=np.array([saturation(pool,10**10) for pool in range(j)])
    def at_least(pool,target):
        if saturation(pool,lower[pool])>=target:return int(lower[pool])
        lo=int(lower[pool]);hi=int(min(total,upper[pool]))
        if saturation(pool,hi)<target:return hi+1
        while lo<hi:
            mid=(lo+hi)//2
            if saturation(pool,mid)>=target:hi=mid
            else:lo=mid+1
        return lo
    lo,hi=0.,float(np.min(ceiling))
    # S_j is a fraction in [0,1]; 24 bisection steps resolve the common
    # threshold to <6.0e-8, comfortably below practical coverage precision.
    for _ in range(24):
        target=(lo+hi)/2
        needed=np.array([at_least(pool,target) for pool in range(j)],dtype=int)
        if np.any(needed>upper) or needed.sum()>total:hi=target
        else:lo=target
    allocation=np.array([at_least(pool,lo) for pool in range(j)],dtype=int)
    allocation=np.maximum(allocation,lower)
    if allocation.sum()>total:
        raise RuntimeError("Balanced saturation solver exceeded colony budget")
    residual=int(total-allocation.sum())
    if residual:
        weights=np.ones(j) if complexity_weights is None else np.asarray(complexity_weights,dtype=float)
        allocation+=integer_allocation(residual,weights,maximum=upper-allocation)
    return allocation


def saturation_grid(max_k,extra=()):
    """Dense deterministic early-to-plateau colony-depth grid."""
    max_k=max(1,int(max_k))
    geometric=np.geomspace(1,max_k,16).astype(int)
    linear=np.linspace(0,max_k,17).astype(int)
    return np.unique(np.r_[0,geometric,linear,np.asarray(list(extra),dtype=int)])


def pool_distinct_allele_expectations(gids,site_probabilities,k_values,n_genes):
    """Exact marginal expected genes hit by >=1/2/3 distinct coordinates.

    Probability mass outside a gene is retained as q. For two or three
    represented coordinates, occupancy probabilities use inclusion-exclusion.
    This is pool-specific and is not a colony-count redundancy model.
    """
    gids=np.asarray(gids,dtype=int);probs=np.asarray(site_probabilities,dtype=float)
    ks=np.asarray(k_values,dtype=int)
    if probs.ndim!=1 or len(gids)!=len(probs) or np.any(probs<0) or probs.sum()>1+1e-10 or np.any(ks<0):
        raise ValueError("Invalid pool coordinate probabilities")
    output=np.zeros((len(ks),3),dtype=float)
    def pow_mass(base):
        base=np.asarray(base,dtype=float)
        return np.where(ks[:,None]==0,1.,np.power(np.clip(base[None,:],0,1),ks[:,None]))
    for gene in np.unique(gids[(gids>=0)&(probs>0)]):
        ps=probs[(gids==gene)&(probs>0)]
        q=max(0.,1-float(ps.sum()))
        p0=np.power(q,ks,dtype=float)
        if np.allclose(ps,ps[0],rtol=1e-11,atol=1e-14):
            # Equal-insertion pools have equal per-coordinate probabilities;
            # closed occupancy expressions avoid a quadratic site-pair loop.
            m=len(ps);site=float(ps[0]);q1=q+site;q2=q+2*site
            p1=m*(np.power(q1,ks)-p0)
            p2=(m*(m-1)/2)*(np.power(q2,ks)-2*np.power(q1,ks)+p0) if m>=2 else np.zeros(len(ks))
        else:
            q1=q+ps
            p1=(np.power(q1[None,:],ks[:,None])-p0[:,None]).sum(axis=1)
            ii,jj=np.triu_indices(len(ps),1)
            if len(ii):
                q2=q+ps[ii]+ps[jj]
                pair=(np.power(q2[None,:],ks[:,None])-np.power((q+ps[ii])[None,:],ks[:,None])-
                      np.power((q+ps[jj])[None,:],ks[:,None])+p0[:,None])
                p2=pair.sum(axis=1)
            else:p2=np.zeros(len(ks))
        ge1=1-p0
        ge2=np.maximum(0.,ge1-p1)
        ge3=np.maximum(0.,ge2-p2)
        output+=np.column_stack((ge1,ge2,ge3))
    return np.minimum(np.maximum(output,0.),n_genes)


def expected_coordinate_recovery(site_probabilities,allocation):
    """Expected unique global coordinates and pool-derived coordinate occurrences."""
    p=np.asarray(site_probabilities,dtype=float);k=np.asarray(allocation,dtype=int)
    if p.ndim!=2 or k.shape!=(p.shape[1],) or np.any(k<0) or np.any(p<0) or np.any(p.sum(axis=0)>1+1e-10):
        raise ValueError("Invalid coordinate/pool allocation")
    if not k.any():return dict(unique_global_coordinates=0.,pool_derived_coordinate_occurrences=0.)
    misses=np.ones(p.shape[0]);occurrences=0.
    for j,depth in enumerate(k):
        if depth:
            q=p[:,j]
            misses*=np.exp(depth*np.log1p(-np.minimum(q,1.)))
            occurrences+=float((-np.expm1(depth*np.log1p(-np.minimum(q,1.)))).sum())
    return dict(unique_global_coordinates=float((1-misses).sum()),
                pool_derived_coordinate_occurrences=occurrences)


def recovery_denominators(row,represented_genes,reference_genes):
    """Explicit units alongside legacy fraction fields; never swap universes."""
    out=dict(row,represented_denominator=int(represented_genes),reference_denominator=int(reference_genes))
    for field in ("represented","lower_genes","upper_genes","exact_expected_genes"):
        value=row.get(field)
        if value is not None:
            out[field+"_union_percent"]=100*value/represented_genes
            out[field+"_reference_percent"]=100*value/reference_genes
    return out


def summarize_counts(counts,n):
    values=np.asarray(counts,dtype=float)
    lo,hi=np.quantile(values,[.025,.975])
    return dict(mean=float(values.mean()/n),represented=float(values.mean()),missed=float(n-values.mean()),
        sd_genes=float(values.std(ddof=1)),sd=float(values.std(ddof=1)/n),lower=float(lo/n),upper=float(hi/n),
        lower_genes=float(lo),upper_genes=float(hi))


def monte_carlo_pools(p,allocations,iterations=1000,seed=2026):
    """Exact pool-stratified multinomial gene simulation (lossless aggregation).

    Probability rows are represented genes under full-catalogue draws; missing
    probability mass is insertion-positive noncoding/unassigned picks. Common
    pool-specific draw streams support arbitrary allocations.
    """
    p=_probabilities(p); allocations=np.asarray(allocations)
    if allocations.ndim!=2 or allocations.shape[1]!=p.shape[1] or np.any(allocations<0) or np.any(allocations!=np.floor(allocations)) or iterations<2:
        raise ValueError("Invalid simulation allocation/iterations")
    allocations=allocations.astype(np.int64);n=len(p);counts=np.empty((len(allocations),iterations),dtype=int)
    rngs=[np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(p.shape[1])]
    for start in range(0,iterations,64):
        size=min(64,iterations-start);recovered=np.zeros((len(allocations),size,n),dtype=bool)
        for j,rng in enumerate(rngs):
            probs=np.r_[p[:,j],max(0.,1-p[:,j].sum())];probs/=probs.sum()
            observed=np.zeros((size,n),dtype=bool);previous=0
            for depth in sorted(set(allocations[:,j])):
                if depth>previous:observed|=rng.multinomial(int(depth-previous),probs,size=size)[:,:n]>0
                recovered[allocations[:,j]==depth]|=observed[None,:,:];previous=depth
        counts[:,start:start+size]=recovered.sum(axis=2)
    return [summarize_counts(x,n) for x in counts]


def monte_carlo_pool_alleles(gids,site_p,allocations,n,iterations=1000,seed=2026):
    """Actual coordinate draws per pool; globally deduplicate allele identities.

    Site probabilities span each full pool catalogue. Noncoding/unassigned mass
    is sampled too (equivalently, coding picks are thinned by their total mass)
    but cannot recover a gene. Shared coordinate identities are deduplicated.
    """
    gids=np.asarray(gids,dtype=int);site_p=np.asarray(site_p,dtype=float);allocations=np.asarray(allocations)
    if site_p.ndim!=2 or len(gids)!=len(site_p) or np.any(site_p<0) or np.any(~np.isfinite(site_p)) or not np.allclose(site_p.sum(axis=0),1) or allocations.ndim!=2 or allocations.shape[1]!=site_p.shape[1] or np.any(allocations<0) or np.any(allocations!=np.floor(allocations)) or iterations<2:
        raise ValueError("Invalid coordinate simulation")
    coding=gids>=0;genes=gids[coding];probs=site_p[coding];mass=probs.sum(axis=0)
    cdfs=[np.cumsum(probs[:,j]/mass[j]) if mass[j] else None for j in range(site_p.shape[1])]
    values=np.zeros((len(allocations),iterations,3),dtype=int)
    rngs=[np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(site_p.shape[1])]
    depths=[sorted(set(allocations[:,j].astype(int))) for j in range(site_p.shape[1])]
    for rep in range(iterations):
        recovered=np.zeros((len(allocations),len(genes)),dtype=bool)
        for j,rng in enumerate(rngs):
            seen=np.zeros(len(genes),dtype=bool);previous=0
            for depth in depths[j]:
                useful=rng.binomial(int(depth-previous),float(mass[j]))
                while useful:
                    block=min(useful,100000)
                    draw=np.searchsorted(cdfs[j],rng.random(block),side="right")
                    seen[np.minimum(draw,len(genes)-1)]=True;useful-=block
                recovered[allocations[:,j]==depth]|=seen[None,:];previous=depth
        for row,seen in enumerate(recovered):
            count=np.bincount(genes[seen],minlength=n)
            values[row,rep]=[(count>=r).sum() for r in (1,2,3)]
    return [[dict(r=r+1,**summarize_counts(values[i,:,r],n)) for r in range(3)] for i in range(len(allocations))]
