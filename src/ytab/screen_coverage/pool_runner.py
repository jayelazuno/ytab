"""Reproducible pool-aware scientific checkpoint and Shiny worker payload."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from .pool_coverage import (LABELS, STRATEGIES, load_pool_catalogue, expected_recovery,
    optimize_allocation, equal_allocation_target, integer_allocation, recovery_denominators,
    complexity_allocation, balanced_saturation_allocation, pool_gene_saturation,
    pool_distinct_allele_expectations, expected_coordinate_recovery, saturation_grid,
    monte_carlo_pools, monte_carlo_pool_alleles)
from .statistics import colonies_for_target_coverage, poisson_coverage


def calculate_pool_screen(request, cache_dir=None, progress=None):
    repo=Path(request.get("repo_root",Path(__file__).resolve().parents[3]))
    audit_dir=request.get("audit_dir",str(repo/"docs/codex/parent_union_audit"))
    data=load_pool_catalogue(audit_dir,repo,request.get("selected_sources"),request.get("selected_pools"))
    ids=data["pool_ids"];n=data["summary"]["represented"]
    k=request.get("k",10000);iterations=request.get("iterations",1000);seed=request.get("seed",2026)
    if any(not isinstance(x,(int,float)) or not math.isfinite(x) or int(x)!=x for x in (k,iterations,seed)) or not 0<=k<=10**10 or not 2<=iterations<=10000 or seed<0:
        raise ValueError("Use integer k (0–10^10), iterations (2–10000) and a nonnegative seed")
    k,iterations,seed=int(k),int(iterations),int(seed)
    target=float(request.get("target",.95))
    if not 0<target<1:raise ValueError("Coverage target must be between 0 and 1")
    def constraint(name,default):
        raw=request.get(name)
        if raw is None:return default
        if set(raw)!=set(ids) or any(not np.isfinite(v) or v<0 or v!=int(v) for v in raw.values()):
            raise ValueError(name+" must give a nonnegative integer for each physical pool")
        return np.array([int(raw[p]) for p in ids])
    minimum=constraint("minimum_per_pool",None);maximum=constraint("maximum_per_pool",None)
    mode=request.get("allocation_mode","complexity")
    if mode not in {"equal","complexity","balanced","manual"}:raise ValueError("Unknown primary allocation strategy")
    if mode=="balanced" and k<len(ids):raise ValueError("Balanced allocation requires at least one colony per selected independent pool")
    manual=constraint("manual_allocation",None) if mode=="manual" else None
    if mode=="manual" and manual is None:raise ValueError("Manual mode requires counts for every selected pool")
    if manual is not None:
        k=int(manual.sum())
        if (minimum is not None and np.any(manual<minimum)) or (maximum is not None and np.any(manual>maximum)):
            raise ValueError("Manual allocation violates supplied constraints")
    sizes=sorted({5000,10000,15000,20000,25000,k})
    sizes=[v for v in sizes if (minimum is None or v>=minimum.sum()) and (maximum is None or v<=maximum.sum())]
    if k not in sizes:raise ValueError("Planned budget violates supplied pool constraints")
    cache_request=dict(request,k=k,iterations=iterations,seed=seed,sources=data["sources"])
    engine=b"".join(Path(__file__).with_name(x).read_bytes() for x in ("pool_runner.py","pool_coverage.py","statistics.py"))
    key=hashlib.sha256(engine+json.dumps(cache_request,sort_keys=True).encode()).hexdigest()
    cache=Path(cache_dir)/("pool-screen-"+key+".json") if cache_dir else None
    if cache is not None and cache.exists():return json.loads(cache.read_text())
    out=dict(summary=data["summary"],pools=data["pools"],overlaps=data["overlaps"],observations=data["observations"],
        sources=data["sources"],settings=dict(request,k=k,iterations=iterations,seed=seed,
            allocation_mode="manual" if manual is not None else mode,
            primary_strategy="manual" if manual is not None else mode,
            primary_strategy_label=STRATEGIES["manual" if manual is not None else mode]),
        theory=[],targets=[],practical=[],allocations=[],baselines=[],alleles=[],marginal=[],curves=[],
        gene_access=[],probabilities=[],pool_abundance=[],discovery=[],observation_qc=[],
        limitations=[
          "Pools are modeled as independently constructed insertion libraries on the same genetic background, as specified by the investigator. Pool overlap is overlap in insertion-coordinate/gene content, not evidence about ancestry.",
          "Primary empirical K is total insertion-positive colony/mutant picks. Each pool is sampled from its full observed genomic coordinate catalogue; coding coordinates can recover protein-coding genes, while noncoding/intergenic/unassigned coordinates consume picks without recovering those genes.",
          "Complexity-normalized allocation uses full genomic insertion-coordinate complexity M_j and approximately equalizes K_j/M_j; these are screening weights, not physical abundance estimates.",
          "Each parent observation is normalized by its own read depth; repeated observations are averaged within its physical pool only. Sequencing abundance is not measured viable-colony frequency.",
          "Complexity-normalized and balanced saturation-aware allocation are planning strategies, not physical mixing proportions. The immediate >=1-gene optimizer is an advanced comparison and may overconcentrate screening in one pool.",
          "Parent sequencing is an observed catalogue, not complete physical-library composition. Pool-only means observed only in that pool; nondetection does not prove biological absence from other pools.",
          "Hermes insertion preferences make genes unequally represented. Sequencing abundance reflects growth, extraction, PCR, sequencing, mapping and bottlenecks; very large sensitivity-model requirements are not experimental recommendations.",
          "Observed CDS insertions are an informative coding proxy, not proven loss-of-function mutants. Coordinate identity does not establish clonal independence, and one insertion per clone is assumed.",
          "95% simulation intervals describe sampling variation conditional on assumed probabilities, not sequencing, viability or model uncertainty.",
          "Target estimates are first tested Monte Carlo mean crossings on a local grid, with the previous tested budget reported; they are not guarantees or exact experimental minima.",
          "Exact pool-binomial expectations guide allocation and saturation curves. Gene-level multinomial simulation is a lossless aggregation of per-pool coding-coordinate probabilities; independent-allele outputs draw and globally deduplicate coding coordinates."])
    for t in sorted({.9,.95,.99,.999,target}):
        kk=colonies_for_target_coverage(n,t)
        out["theory"].append(dict(method="Theoretical Poisson / Equal-Gene Benchmark",target=t,k=kk,depth=kk/n,
            represented=n*poisson_coverage(kk,n),pool_allocation="Not defined by theoretical model"))
    curve_sizes=list(range(0,25001,1000))
    for kk in sorted(set(sizes+curve_sizes)):
        c=poisson_coverage(kk,n)
        row=dict(method="Theoretical Poisson / Equal-Gene Benchmark",k=kk,r=1,mean=c,represented=n*c,missed=n*(1-c),lower=None,upper=None)
        out["curves"].append(row)
        if kk in sizes:out["practical"].append(row)
    site_counts=[v["sites"] for v in data["pools"]]
    coding_counts=[v["coding_sites"] for v in data["pools"]]
    # Primary empirical K samples the full observed insertion-positive catalogue.
    complexity_counts=site_counts
    model_items=list(data["gene_probabilities"].items())
    # The interactive planner computes only the selected empirical model; the
    # full two-model analysis remains available to CLI/QC runs. This avoids
    # making one model selection wait for an unrelated sensitivity analysis.
    if request.get("planner_only") and request.get("model_mode") in dict(model_items):
        model_items=[(request["model_mode"],dict(model_items)[request["model_mode"]])]
    for model_index,(model,p) in enumerate(model_items):
        if progress:progress("Pool-stratified strategies: "+model)
        name=LABELS[model];plans={};fits={};allocation_cache={}
        def fit(kk):
            if kk not in fits:fits[kk]=optimize_allocation(p,int(kk),complexity_counts,coding_counts,minimum,maximum)
            return fits[kk]
        def get_allocation(strategy,kk):
            key=(strategy,int(kk))
            if key in allocation_cache:return allocation_cache[key]
            if strategy=="equal":result=integer_allocation(kk,np.ones(len(ids)),minimum,maximum)
            elif strategy=="complexity":result=complexity_allocation(kk,complexity_counts,minimum,maximum)
            elif strategy=="balanced":result=balanced_saturation_allocation(kk,p,data["gene_pool_presence"],minimum,maximum,complexity_counts)
            elif strategy=="optimized":result=fit(kk)["allocation"]
            elif strategy=="site_proportional":result=integer_allocation(kk,site_counts,minimum,maximum)
            elif strategy=="coding_proportional":result=integer_allocation(kk,coding_counts,minimum,maximum)
            elif strategy=="manual":result=manual
            else:raise ValueError("Unknown allocation strategy: "+strategy)
            allocation_cache[key]=result
            return result
        planner_only=bool(request.get("planner_only",False))
        primary_strategies=("complexity","balanced","equal")
        selected_strategy=request.get("allocation_mode","complexity")
        if planner_only and selected_strategy in primary_strategies:
            # The Screen Planner only needs the chosen allocation's curve;
            # changing the strategy reruns this calculation. Full QC/CLI runs
            # continue to calculate all primary allocation comparisons.
            primary_strategies=(selected_strategy,)
        all_strategies=primary_strategies if planner_only else primary_strategies+("site_proportional","coding_proportional","optimized")
        planned_ks=sorted(set(sizes+curve_sizes))
        if manual is not None:planned_ks=sorted(set(planned_ks+[k]))
        for kk in planned_ks:
            for strategy in primary_strategies:
                if strategy=="balanced" and kk<len(ids):continue
                plans[(kk,strategy)]=get_allocation(strategy,kk)
        if not planner_only:
            for kk in sizes:
                for strategy in ("optimized","site_proportional","coding_proportional"):
                    plans[(kk,strategy)]=get_allocation(strategy,kk)
        if manual is not None:plans[(k,"manual")]=manual
        plan_keys=list(plans);allocations=np.array([plans[key] for key in plan_keys])
        if progress:progress("Simulating "+model+" ("+str(len(plans))+" allocations)")
        if request.get("simulate",True):
            simulations=monte_carlo_pools(p,allocations,iterations,seed+100*model_index)
        else:simulations=[None]*len(plans)
        by_key={}
        for (kk,strategy),allocation,mc in zip(plan_keys,allocations,simulations):
            analytic=float(expected_recovery(p,allocation).sum())
            row=dict(method=name,model=model,k=kk,strategy=strategy,r=1,allocation=allocation.tolist(),
                strategy_label=STRATEGIES[strategy],
                exact_expected_genes=analytic,exact_coverage=analytic/n,
                **expected_coordinate_recovery(data["probabilities"][model],allocation),
                optimization_gap_genes=float(fit(kk)["relaxation_gap"]) if strategy=="optimized" else None,
                **(mc or dict(mean=analytic/n,represented=analytic,missed=n-analytic,lower=None,upper=None,
                              sd_genes=None,sd=None,lower_genes=None,upper_genes=None)))
            row["estimate_type"]="Monte Carlo" if mc is not None else "Exact expectation; simulation not run"
            by_key[(kk,strategy)]=row
            if kk in sizes:
                if strategy in all_strategies:out["baselines"].append(row)
                for j,pool in enumerate(ids):
                    out["allocations"].append(dict(method=name,model=model,k=kk,strategy=strategy,pool_id=pool,
                        strategy_label=STRATEGIES[strategy],label=data["pools"][j]["label"],colonies=int(allocation[j]),fraction=float(allocation[j]/kk) if kk else 0,
                        observed_sites=site_counts[j],observed_coding_sites=coding_counts[j],
                        normalized_depth=float(allocation[j]/site_counts[j]) if site_counts[j] else 0.))
            if strategy in primary_strategies and kk in curve_sizes:
                out["curves"].append(row)
            if strategy in primary_strategies and kk in sizes:out["practical"].append(row)
            if strategy=="manual":out["practical"].append(row)
        # Invert exact expected coverage for each primary allocation strategy.
        for t in sorted({.9,.95,.99,target}):
            for strategy in primary_strategies:
                def score(kk):
                    if strategy=="balanced" and kk<len(ids):return 0.
                    return float(expected_recovery(p,get_allocation(strategy,kk)).mean())
                if np.mean(np.any(p>0,axis=1))<t:
                    out["targets"].append(dict(method=name,strategy=strategy,strategy_label=STRATEGIES[strategy],target=t,k=None,status="Target is above the represented-union ceiling"));continue
                low=(int(minimum.sum()) if minimum is not None else (len(ids) if strategy=="balanced" else 0));high=max(low,len(ids),len(p))
                if maximum is not None:high=min(high,int(maximum.sum()))
                while score(high)<t and high<10**10:high=min(10**10,high*2)
                if score(high)<t:
                    out["targets"].append(dict(method=name,strategy=strategy,strategy_label=STRATEGIES[strategy],target=t,k=None,status="Target not reached under constraints"));continue
                while low<high:
                    mid=(low+high)//2
                    if score(mid)>=t:high=mid
                    else:low=mid+1
                allocation=get_allocation(strategy,low);expected=float(expected_recovery(p,allocation).sum())
                out["targets"].append(dict(method=name,strategy=strategy,strategy_label=STRATEGIES[strategy],target=t,k=low,
                    allocation=allocation.tolist(),exact_expected_genes=expected,exact_coverage=expected/n,
                    status="First exact expected-coverage crossing"))
        practical_keys=[key for key in plan_keys if key[0] in sizes and key[1] in set(primary_strategies+("manual",)) and key[0]<=1000000]
        if request.get("simulate",True):
            if progress:progress("Coordinate/allele simulation "+model)
            results=monte_carlo_pool_alleles(data["allele_genes"],data["probabilities"][model],
                np.array([plans[key] for key in practical_keys]),n,iterations,seed+100*model_index+1)
            for (kk,strategy),rows in zip(practical_keys,results):
                out["alleles"].extend(dict(method=name,model=model,k=kk,strategy=strategy,
                    strategy_label=STRATEGIES[strategy],**row) for row in rows)
                # At practical depths, headline gene recovery and allele recovery
                # come from the SAME coordinate draws. Gene-aggregated MC above
                # independently validates these and handles very large targets.
                by_key[(kk,strategy)].update({field:value for field,value in rows[0].items() if field!="r"})
                by_key[(kk,strategy)]["estimate_type"]="Coordinate Monte Carlo"
        for strategy in (list(primary_strategies)+(["manual"] if manual is not None else [])):
            if (k,strategy) not in plans:continue
            allocation=plans[(k,strategy)];full=expected_recovery(p,allocation).sum()
            for j,pool in enumerate(ids):
                without=allocation.copy();without[j]=0
                next_block=allocation.copy();next_block[j]+=100
                pool_expected=float(expected_recovery(p[:,j:j+1],[allocation[j]]).sum())
                nj=int(data["gene_pool_presence"][:,j].sum())
                out["marginal"].append(dict(method=name,pool_id=pool,strategy=strategy,
                    colonies=int(allocation[j]),leave_out_additional_genes=float(full-expected_recovery(p,without).sum()),
                    expected_gain_next_100=float(expected_recovery(p,next_block).sum()-full),
                    pool_represented_genes=nj,pool_expected_genes=pool_expected,
                    pool_saturation_fraction=pool_expected/nj if nj else 0,
                    observed_sites=int(site_counts[j]),
                    normalized_sampling_depth=float(allocation[j]/site_counts[j]) if site_counts[j] else 0.))
    # Pool-specific saturation uses each pool's own represented-gene repertoire N_j.
    # Gene MC supplies intervals; coordinate/allele recovery expectations use exact
    # within-pool categorical probabilities over the full observed catalogue.
    out["pool_saturation"]=[]
    for model_index,(model,gene_p) in enumerate(model_items):
        site_p=data["probabilities"][model]
        for j,pool in enumerate(ids):
            nj=int(data["gene_pool_presence"][:,j].sum());mj=int(site_counts[j])
            depths=[row["colonies"] for row in out["allocations"] if row["model"]==model and row["pool_id"]==pool and row["strategy"] in ("complexity","balanced","equal")]
            grid=saturation_grid(3*mj,depths)
            p_one=gene_p[:,j:j+1]
            if request.get("simulate",True):
                gene_mc=monte_carlo_pools(p_one,grid[:,None],iterations,seed+500+100*model_index+j)
            else:gene_mc=[None]*len(grid)
            allele_values=pool_distinct_allele_expectations(data["allele_genes"],site_p[:,j],grid,n)
            for q,depth in enumerate(grid):
                exact_genes=float(pool_gene_saturation(p_one,int(depth),data["gene_pool_presence"][:,j:j+1])[0]*nj)
                unique_sites=float((-np.expm1(int(depth)*np.log1p(-np.minimum(site_p[:,j],1.)))).sum()) if depth else 0.
                mc=gene_mc[q]
                lower=mc["lower_genes"] if mc else None;upper=mc["upper_genes"] if mc else None
                out["pool_saturation"].append(dict(method=LABELS[model],model=model,pool_id=pool,
                    pool_label=data["pools"][j]["label"],k=int(depth),site_denominator=mj,
                    gene_denominator=nj,expected_unique_sites=unique_sites,
                    insertion_sites_fraction=unique_sites/mj if mj else 0.,
                    expected_genes=exact_genes,gene_fraction=exact_genes/nj if nj else 0.,
                    monte_carlo_genes=mc["represented"]*n if mc else None,
                    lower_genes=lower,upper_genes=upper,
                    lower_gene_fraction=lower/nj if lower is not None and nj else None,
                    upper_gene_fraction=upper/nj if upper is not None and nj else None,
                    expected_genes_ge1_allele=float(allele_values[q,0]),
                    expected_genes_ge2_alleles=float(allele_values[q,1]),
                    expected_genes_ge3_alleles=float(allele_values[q,2]),
                    marginal_genes_per_100=float((pool_gene_saturation(p_one,int(depth)+100,
                        data["gene_pool_presence"][:,j:j+1])[0]-exact_genes/nj)*nj/100) if nj else 0.,
                    normalized_sampling_depth=float(depth/mj) if mj else 0.))
    out["complexity_weights"]=[dict(pool_id=ids[j],label=data["pools"][j]["label"],
        observed_sites=site_counts[j],eligible_coding_sites=coding_counts[j],
        coding_fraction=float(coding_counts[j]/site_counts[j]),
        weight=float(site_counts[j]/sum(site_counts))) for j in range(len(ids))]
    for row in out["marginal"]:
        match=next((x for x in out["pool_saturation"] if x["model"]==row["method"].split("— ")[-1].lower().replace(" ","_")
                    and x["pool_id"]==row["pool_id"] and x["k"]==row["colonies"]),None)
        # Resolve by user-facing model instead of relying on the display string.
        if match is None:
            model_key="equal" if "Equal Insertion" in row["method"] else "abundance"
            match=next((x for x in out["pool_saturation"] if x["model"]==model_key and x["pool_id"]==row["pool_id"] and x["k"]==row["colonies"]),None)
        pool=next(x for x in data["pools"] if x["pool_id"]==row["pool_id"])
        row.update(pool_genomic_sites=pool["sites"],pool_coding_sites=pool["coding_sites"],
            pool_noncoding_sites=pool["noncoding_sites"],pool_coding_fraction=pool["coding_fraction"],
            pool_unique_genes=pool["unique_genes"],pool_unique_coding_alleles=pool["unique_coding_alleles"],
            genes_ge2_pool_alleles=pool["ge2"],genes_ge3_pool_alleles=pool["ge3"])
        if match is not None:
            row.update(expected_unique_sites=match["expected_unique_sites"],
                pool_genes_ge1_distinct=match["expected_genes_ge1_allele"],
                pool_genes_ge2_distinct=match["expected_genes_ge2_alleles"],
                pool_genes_ge3_distinct=match["expected_genes_ge3_alleles"],
                pool_marginal_genes_per_100=match["marginal_genes_per_100"])
    for i,gene in enumerate(data["gene_ids"]):
        access=[ids[j] for j in np.flatnonzero(data["gene_pool_presence"][i])]
        out["gene_access"].append(dict(gene_id=gene,pool_count=len(access),physical_pools=";".join(access)))
        for j,pool in enumerate(ids):
            out["probabilities"].append(dict(gene_id=gene,pool_id=pool,
                coding_length_bp=int(data["gene_lengths"][i]),
                coding_coordinates_in_pool=int(data["gene_coordinate_counts"][i,j]),
                coding_coordinates_global=int(data["gene_coordinate_counts"][i].sum()),
                pools_represented=len(access),
                insertion_coordinates_per_kb=float(data["gene_coordinate_counts"][i,j]/(data["gene_lengths"][i]/1000)),
                p_equal=float(data["gene_probabilities"]["equal"][i,j]),
                p_abundance=float(data["gene_probabilities"]["abundance"][i,j])))
    for j,pool in enumerate(ids):
        positive=data["probabilities"]["abundance"][:,j];positive=positive[positive>0]
        out["pool_abundance"].append(dict(pool_id=pool,probability_sum=float(positive.sum()),
            median=float(np.median(positive)),p99=float(np.quantile(positive,.99)),maximum=float(positive.max()),
            probability_universe="full observed genomic insertion-coordinate catalogue"))
    # Parent-sequencing discovery is a separate diagnostic, not needed to
    # render the default colony-screen planner. Defer its 1,000-ordering
    # calculation during automatic planner bootstrap; an explicit Calculate
    # coverage run includes it for the Parent Sequencing subtab.
    membership=data["observation_reads"]>0;ng=len(data["gene_ids"])
    gene_membership=np.zeros((ng,membership.shape[1]),dtype=bool)
    for j in range(membership.shape[1]):
        g=data["allele_genes"][membership[:,j]];gene_membership[np.unique(g[g>=0]),j]=True
    for j,obs in enumerate(data["observations"]):
        out["observation_qc"].append(dict(**obs,sites=int(membership[:,j].sum()),
            coding_sites=int((membership[:,j]&(data["allele_genes"]>=0)).sum()),genes=int(gene_membership[:,j].sum()),
            uniquely_contributed_sites=int((membership[:,j]&(membership.sum(axis=1)==1)).sum()),
            uniquely_contributed_genes=int((gene_membership[:,j]&(gene_membership.sum(axis=1)==1)).sum()),
            sequencing_reads=float(data["observation_reads"][:,j].sum())))
    if request.get("include_discovery",True):
        rng=np.random.default_rng(seed+999)
        discovery=np.zeros((1000,membership.shape[1],2),dtype=int)
        patterns,weights=np.unique(membership,axis=0,return_counts=True)
        for rep in range(1000):
            order=rng.permutation(membership.shape[1])
            discovery[rep,:,0]=np.cumsum(patterns[:,order],axis=1).astype(bool).T@weights
            discovery[rep,:,1]=np.cumsum(gene_membership[:,order],axis=1).astype(bool).T.sum(axis=1)
        for j in range(membership.shape[1]):
            for metric,c in (("Genomic insertion sites",0),("Protein-coding genes",1)):
                v=discovery[:,j,c]
                out["discovery"].append(dict(observations=j+1,metric=metric,mean=float(v.mean()),
                    lower=float(np.quantile(v,.025)),upper=float(np.quantile(v,.975)),minimum=int(v.min()),maximum=int(v.max())))
    for table in ("theory","targets","practical","baselines","alleles","curves"):
        out[table]=[recovery_denominators(row,n,data["summary"]["annotated"]) for row in out[table]]
    out["seed_map"]=dict(equal_colony=seed,equal_allele=seed+1,abundance_colony=seed+100,abundance_allele=seed+101,discovery=seed+999)
    if cache is not None:
        cache.parent.mkdir(parents=True,exist_ok=True);cache.write_text(json.dumps(out,allow_nan=False))
    return out


def export_result(result,directory,data=None):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    (directory/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False))
    for name,rows in result.items():
        if not isinstance(rows,list) or not rows or not isinstance(rows[0],dict):continue
        columns=list(dict.fromkeys(k for row in rows for k in row))
        with (directory/(name+".csv")).open("w") as handle:
            writer=csv.DictWriter(handle,fieldnames=columns);writer.writeheader()
            writer.writerows({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in row.items()} for row in rows)
    if data is not None:
        with (directory/"coordinate_pool_provenance.csv").open("w") as handle:
            writer=csv.writer(handle);writer.writerow(["chromosome","position","gene_id","physical_pool_ids"])
            for i,coord in enumerate(data["coordinates"]):
                gi=data["allele_genes"][i]
                writer.writerow([*coord,data["gene_ids"][gi] if gi>=0 else "",
                    ";".join(np.array(data["pool_ids"])[data["presence"][i]])])


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--request",required=True);parser.add_argument("--output",required=True);parser.add_argument("--cache-dir")
    args=parser.parse_args();request=json.loads(Path(args.request).read_text())
    result=calculate_pool_screen(request,args.cache_dir,lambda s:print(s,flush=True))
    export_result(result,args.output,load_pool_catalogue(request["audit_dir"],request["repo_root"],request.get("selected_sources")))
