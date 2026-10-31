"""Empirical site catalogues and trimmed-total-reads (TTR) normalization for Hermes raw counts, adapted from TRANSIT's TTRNorm (DeJesus et al. 2015).

Independent implementation of density/trimmed-nonzero-mean scaling. Coordinates
are those already called by YTAB, not inferred TA targets or genome-wide zeros.
"""
import csv
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse
from scipy.special import gammaln
from scipy.stats import trim_mean

from ytab.pipeline.project_config import load_project_config
from ytab.screen_coverage.library import is_parent_sample


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_counts(path):
    """Sum nonnegative raw support at exact 1-based positions, across W/C.

    WIG browser values may be negative to encode strand. We deliberately use the
    canonical hit file's nonnegative Hit count, never sum signed WIG values.
    """
    counts, strands = {}, {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"Chromosome", "Hit position", "Hit count"} <= set(reader.fieldnames or []):
            raise ValueError(f"Not an unnormalized YTAB coordinate hit file: {path}")
        for row in reader:
            position, count = float(row["Hit position"]), float(row["Hit count"])
            if (not math.isfinite(position) or position < 1 or not position.is_integer() or
                    not math.isfinite(count) or count < 0 or not count.is_integer()):
                raise ValueError(f"Expected positive integer coordinates and nonnegative integer raw counts: {path}")
            if not row["Chromosome"].strip():
                raise ValueError(f"Missing chromosome: {path}")
            if count:
                key = (row["Chromosome"], int(position))
                counts[key] = counts.get(key, 0) + int(count)
                strands.setdefault(key, set()).add(row.get("Source", "unspecified"))
    if not counts:
        raise ValueError(f"No positive insertion-support counts: {path}")
    return counts, strands


def ttr_parameters(positive_counts, universe_size, trim=0.05, target=100.0):
    values = np.asarray(positive_counts, dtype=float)
    if (values.ndim != 1 or not len(values) or np.any(~np.isfinite(values)) or
            np.any(values <= 0) or universe_size < len(values) or not 0 <= trim < 0.5 or
            not math.isfinite(target) or target <= 0):
        raise ValueError("TTR requires finite positive counts, a common nonempty universe, 0≤trim<0.5 and target>0")
    theta = len(values) / universe_size
    mu = float(trim_mean(values, proportiontocut=trim))
    factor = target / (theta * mu)
    return dict(universe_sites=universe_size, occupied_sites=len(values), density=theta,
                trimmed_nonzero_mean=mu, trim_fraction_each_tail=trim,
                trimmed_sites_each_tail=math.floor(len(values)*trim), target=target,
                factor=factor, normalized_density_times_trimmed_mean=theta*mu*factor)


def expected_sites_at_depth(counts, depth):
    """Exact rarefaction expectation under count subsampling without replacement.

    Depth is usable insertion-support count, not total sequenced FASTQ reads.
    Groups identical counts before evaluating the hypergeometric zero probability.
    """
    counts = np.asarray(counts, dtype=np.int64)
    total, depth = int(counts.sum()), int(depth)
    if np.any(counts <= 0) or not 0 <= depth <= total:
        raise ValueError("Depth must be between zero and total positive raw support")
    if depth == 0:
        return 0.0
    if depth == 1:
        return 1.0
    if depth == total:
        return float(len(counts))
    values, frequency = np.unique(counts, return_counts=True)
    eligible = total-values >= depth
    v = values[eligible].astype(float)
    log_p0 = (gammaln(total-v+1)-gammaln(total-v-depth+1) +
              gammaln(total-depth+1)-gammaln(total+1))
    present = np.ones(len(values), dtype=float)
    present[eligible] = -np.expm1(np.minimum(log_p0, 0))
    return float(np.dot(frequency, present))


def load_observations(manifest_path):
    manifest_path = Path(manifest_path).resolve()
    manifest = yaml.safe_load(manifest_path.read_text())
    root = (manifest_path.parent / manifest.get("repo_root", ".")).resolve()
    resolve = lambda value, base=root: Path(value) if Path(value).is_absolute() else base / value
    observations, reference_hashes = [], set()
    for source in manifest["sources"]:
        config_path = resolve(source["project_config"]).resolve()
        cfg = load_project_config(config_path)
        base = Path(cfg["_repo_root"])
        project = resolve(cfg["output_project_dir"], base)
        sheet = resolve(cfg["sample_sheet"], base)
        with sheet.open() as handle:
            rows = {r["sample"]: r for r in csv.DictReader(handle)}
        reference_path = project / "config/reference_resolved.json"
        reference = json.loads(reference_path.read_text()) if reference_path.exists() else cfg["reference"]
        fasta = resolve(reference["fasta"], base)
        assembly = sha256(fasta)
        reference_hashes.add(assembly)
        for entry in source["observations"]:
            row = rows[entry["sample"]]
            if str(row.get("include", "true")).strip().lower() not in {"true", "yes", "1"} or not is_parent_sample(row):
                raise ValueError(f"Catalogue accepts included untreated parent/control observations only: {row['sample']}")
            library = entry.get("library_id")
            if not library:
                raise ValueError("An explicit independent library_id is required; repeated observations must share it")
            path = resolve(row.get("hit_file") or ".", base)
            if not path.is_file():
                path = project / "create_hit_file" / row["sample"] / (row["sample"]+"_hits.txt")
            counts, strands = read_counts(path)
            observations.append(dict(observation_id=cfg["project_id"]+"::"+row["sample"],
                sample=row["sample"], display_name=row.get("display_name") or row["sample"],
                project=cfg["project_id"], dataset=source["dataset"], library_id=library,
                pool=row.get("pool_id") or row.get("pool") or row.get("guessed_pool", ""),
                background=row.get("background") or row.get("guessed_background", ""),
                input=str(path.resolve()), input_sha256=sha256(path),
                reference=str(fasta.resolve()), reference_sha256=assembly,
                project_config=str(config_path), sample_sheet_sha256=sha256(sheet),
                counts=counts, strands=strands))
    if not observations or len({o["observation_id"] for o in observations}) != len(observations):
        raise ValueError("Supply nonempty, unique project/sample observations")
    if len(reference_hashes) != 1:
        raise ValueError("Site catalogues require the same reference FASTA checksum")
    for library in {o["library_id"] for o in observations}:
        if len({o["dataset"] for o in observations if o["library_id"] == library}) != 1:
            raise ValueError("An independent library cannot belong to multiple dataset groups")
    return manifest, observations


def write_table(path, fields, rows):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def overlap_rows(groups):
    for (left, a), (right, b) in itertools.combinations(sorted(groups.items()), 2):
        shared = len(a & b)
        yield dict(left=left, right=right, left_sites=len(a), right_sites=len(b), shared_sites=shared,
                   left_only=len(a-b), right_only=len(b-a), union_sites=len(a|b), jaccard=shared/len(a|b))


def run(manifest_path, output, trim=0.05, target=100.0, plots=True, progress=None):
    progress = progress or (lambda message: None)
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory is nonempty; choose a new version to preserve prior results")
    progress("Reading explicit parent observations and raw coordinate counts")
    manifest, observations = load_observations(manifest_path)
    sites = sorted(set().union(*(set(o["counts"]) for o in observations)))
    index = {site: i for i, site in enumerate(sites)}
    progress(f"Catalogue: {len(sites):,} sites; {len(observations)} observations")
    parameters = []
    datasets, libraries, samples = {}, {}, {}
    matrix_rows, matrix_cols, matrix_counts = [], [], []
    for col, observation in enumerate(observations):
        counts = observation["counts"]
        p = ttr_parameters(list(counts.values()), len(sites), trim, target)
        p.update({key: observation[key] for key in ("observation_id","display_name","project","sample","dataset","library_id","pool")})
        p["raw_total_support"] = sum(counts.values())
        p["normalized_total_support"] = p["raw_total_support"]*p["factor"]
        parameters.append(p)
        for collection, name in ((datasets,observation["dataset"]),(libraries,observation["library_id"]),(samples,observation["observation_id"])):
            collection.setdefault(name, set()).update(counts)
        matrix_rows.extend(index[site] for site in counts)
        matrix_cols.extend([col]*len(counts))
        matrix_counts.extend(counts.values())
    matrix = sparse.coo_matrix((np.asarray(matrix_counts,dtype=np.int64),(matrix_rows,matrix_cols)),
                              shape=(len(sites),len(observations))).tocsr()
    normalized = matrix.astype(float).multiply(np.array([p["factor"] for p in parameters])).tocsr()
    output.mkdir(parents=True,exist_ok=True)
    progress("Writing provenance, zero-filled raw matrix and TTR-scaled matrix")
    observation_metadata = [{k:v for k,v in o.items() if k not in {"counts","strands"}} for o in observations]
    write_table(output/"observations.tsv",list(observation_metadata[0]),observation_metadata)
    write_table(output/"ttr_parameters.tsv",list(parameters[0]),parameters)
    sparse.save_npz(output/"raw_counts.npz",matrix)
    sparse.save_npz(output/"ttr_normalized_counts.npz",normalized)
    def site_rows():
        for i,(chrom,position) in enumerate(sites):
            row = matrix.getrow(i)
            members = [observations[c] for c in row.indices]
            yield dict(site_index=i,reference_sha256=observations[0]["reference_sha256"],chromosome=chrom,
                position_1based=position,observations_present=len(members),
                libraries_present=len({m["library_id"] for m in members}),
                datasets_present=";".join(sorted({m["dataset"] for m in members})))
    write_table(output/"sites.tsv.gz",["site_index","reference_sha256","chromosome","position_1based","observations_present","libraries_present","datasets_present"],site_rows())
    def provenance_rows():
        for o in observations:
            for site,count in sorted(o["counts"].items()):
                yield dict(site_index=index[site],observation_id=o["observation_id"],library_id=o["library_id"],
                    dataset=o["dataset"],pool=o["pool"],raw_count=count,source_strands=";".join(sorted(o["strands"][site])))
    write_table(output/"site_provenance.tsv.gz",["site_index","observation_id","library_id","dataset","pool","raw_count","source_strands"],provenance_rows())
    # Dense tabular exports explicitly contain every sample/site zero. NPZ files
    # encode the identical zero-filled matrices sparsely and preserve full precision.
    for name,data in (("raw_counts.tsv.gz",matrix),("ttr_normalized_counts.tsv.gz",normalized)):
        with gzip.open(output/name,"wt",newline="") as handle:
            writer=csv.writer(handle,delimiter="\t")
            writer.writerow(["chromosome","position_1based"]+[o["observation_id"] for o in observations])
            for start in range(0,len(sites),10000):
                dense=data[start:start+10000].toarray()
                writer.writerows([*sites[start+i],*row] for i,row in enumerate(dense))
    overlap_fields=["left","right","left_sites","right_sites","shared_sites","left_only","right_only","union_sites","jaccard"]
    for name,groups in (("dataset",datasets),("library",libraries),("observation",samples)):
        write_table(output/(name+"_overlap.tsv"),overlap_fields,overlap_rows(groups))
    group_stats=[dict(dataset=name,unique_sites=len(group),additional_to_all_other_datasets=len(group-set().union(*(s for n,s in datasets.items() if n!=name)))) for name,group in sorted(datasets.items())]
    write_table(output/"dataset_contributions.tsv",list(group_stats[0]),group_stats)
    # Exact expected catalogue size across all subsets of k independent libraries.
    library_presence=np.zeros(len(sites),dtype=int)
    for group in libraries.values():
        library_presence[[index[s] for s in group]]+=1
    frequency=np.bincount(library_presence)
    j=len(libraries)
    discovery=[dict(independent_libraries=k,expected_union_sites=sum(int(frequency[f])*(1-(math.comb(j-f,k)/math.comb(j,k) if j-f>=k else 0)) for f in range(1,len(frequency)))) for k in range(1,j+1)]
    write_table(output/"library_discovery.tsv",list(discovery[0]),discovery)
    progress("Calculating library overlap and depth-matched expected discovery")
    minimum_depth=min(p["raw_total_support"] for p in parameters)
    depths=sorted({int(minimum_depth*f) for f in (0,.01,.05,.1,.25,.5,.75,1)})
    rarefaction=[]
    for o,p in zip(observations,parameters):
        values=list(o["counts"].values())
        for depth in sorted(set(depths+[p["raw_total_support"]])):
            rarefaction.append(dict(observation_id=o["observation_id"],dataset=o["dataset"],library_id=o["library_id"],
                sampled_support=depth,expected_sites=expected_sites_at_depth(values,depth),
                matched_depth=depth<=minimum_depth))
    write_table(output/"depth_rarefaction.tsv",list(rarefaction[0]),rarefaction)
    summary=dict(status="success",method="TTR",implementation="density * trimmed positive-site mean",
        trim_each_tail=trim,target=target,universe_sites=len(sites),observations=len(observations),
        independent_libraries=len(libraries),datasets=group_stats,minimum_support_depth=minimum_depth,
        reference_sha256=observations[0]["reference_sha256"],manifest=str(Path(manifest_path).resolve()),
        manifest_sha256=sha256(manifest_path),code_sha256=sha256(__file__),
        coordinate_definition="Exact YTAB chromosome + 1-based Hit position; W/C raw counts summed without coordinate shifts",
        matrix_orientation="sites x observations; implicit sparse zeros equal explicit tabular zeros",
        warnings=["Observed catalogue is not all possible Hermes insertion sites or proof of saturation.",
          "Raw support is usable called-insertion count, not total FASTQ sequencing depth.",
          "Read abundance is not measured viable-mutant abundance.",
          "Repeated sequencing/outgrowth observations are not independent induction libraries.",
          "Density correction can scale real differences in library complexity; this is an evaluation, not a validated replacement for fitness CPM or classifier MidLC.",
          "Adding all-zero sites to the common universe changes absolute TTR factors equally, not their ratios."])
    (output/"run_manifest.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    (output/"source_manifest.yaml").write_text(yaml.safe_dump(manifest,sort_keys=False))
    if plots:
        progress("Rendering diagnostic figures")
        draw_plots(output,observations,parameters,rarefaction,discovery)
    progress("Completed; existing project normalization and fitness results were not changed")
    return summary


def draw_plots(output,observations,parameters,rarefaction,discovery):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    labels=[o["display_name"] for o in observations]
    fig,axes=plt.subplots(1,2,figsize=(15,7))
    for ax,key,title in ((axes[0],"density","Observed site density"),(axes[1],"factor","TTR multiplier")):
        ax.barh(range(len(labels)),[p[key] for p in parameters],color="#2f6f9f")
        ax.set_yticks(range(len(labels)),labels)
        ax.invert_yaxis();ax.set_title(title);ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:,.3g}"))
    fig.tight_layout();fig.savefig(output/"ttr_diagnostics.png",dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(15,7))
    for o in observations:
        rows=[r for r in rarefaction if r["observation_id"]==o["observation_id"] and r["matched_depth"]]
        axes[0].plot([r["sampled_support"]/1e6 for r in rows],[r["expected_sites"]/1000 for r in rows],label=o["display_name"])
    axes[0].set(xlabel="Sampled insertion support (millions)",ylabel="Expected observed sites (×1,000)",title="Depth-matched discovery, without replacement")
    axes[0].legend(fontsize=7,loc="upper left",bbox_to_anchor=(0,-.18))
    axes[1].plot([r["independent_libraries"] for r in discovery],[r["expected_union_sites"]/1000 for r in discovery],marker="o")
    axes[1].set(xlabel="Independent libraries included",ylabel="Expected union sites (×1,000)",title="Observed-library accumulation (not a saturation fit)")
    for ax in axes:
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:,.2f}".rstrip("0").rstrip(".") if x else "0"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:,.0f}"))
    fig.tight_layout();fig.savefig(output/"site_discovery.png",dpi=160,bbox_inches="tight");plt.close(fig)
