#!/usr/bin/env python3
"""Run the Hermes window HMM (label-free essentiality) on Hermes libraries of any yeast species.

Inputs: the genome FASTA the reads were aligned to, its GFF annotation, and one or more Hermes hit files
(YTAB *_hits.txt, or any CSV/TSV with chromosome, position and read-count columns). Libraries are fitted
jointly; use independent parent libraries of one strain.

  ytab hermes-hmm OUTPUT --fasta genome.fna --gff genes.gff --hits lib1_hits.txt --hits lib2_hits.txt
  ytab hermes-hmm OUTPUT --project PROJECT --libraries SAMPLE1 SAMPLE2      (genome and GFF from the project config)
With no inputs: Gale 2020 parents Cg1-Cg3 on the C. glabrata reference (the validated configuration).

Checks for other species:
- the FASTA must be the haploid reference used for alignment. A reference with two near-identical
  haplotypes (e.g. full diploid C. albicans) makes almost every 50-mer occur twice, so almost every window
  is masked. The run stops when > 50% of windows are masked (--allow-high-mask to override);
- hit-file chromosome names must match the FASTA (stops if > 10% of insertion sites are on unknown
  chromosomes);
- genes on chromosomes absent from the FASTA (e.g. the C. albicans B haplotype) are excluded and counted.
Ploidy: in a diploid, a heterozygous insertion in an essential gene is viable, so essential genes keep
insertions. Use haploid strains (or haploid derivatives) for essentiality.
Sequence covariate: `--sequence-covariate motif` (default, hermes_hmm_v2) counts the Hermes nTnnnnAn target;
`learned` (hermes_hmm_v3) learns the target-site preference from the libraries' own insertion sites (YTAB hit files
with their W/C strand), so any transposon (Hermes, Ac/Ds, ...) can be modelled; it also writes target_model.tsv.
Writes genes.tsv, windows.npz, fit.json. Reads no essentiality labels.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.essentiality import hermes_hmm as hh
from ytab.essentiality.library_correlation import read_hits
from ytab.normalization.ttr import sha256, write_table

REF = ROOT/"resources/species/glabrata/reference_genome"
FASTA = REF/"GCA_014217725.1_ASM1421772v1_genomic.fna"; GFF = REF/"GCA_014217725.1_ASM1421772v1_genomic.gff"
GALE = ["SRR11561508", "SRR11561507", "SRR11561506"]   # used only by the pre-registered Gale benchmark scripts
CACHE = ROOT/"output/essentiality/hermes_hmm/cache"
MAX_MASKED, MAX_OFF_REFERENCE = .5, .1


def hit_sites(sample, min_reads, project="Gale_2020", path=None):
    """{chromosome: sorted positions with >= min_reads reads} for a project sample or an explicit hit file."""
    path = Path(path) if path else ROOT/"output/projects"/project/"create_hit_file"/sample/f"{sample}_hits.txt"
    return {c:np.asarray(p, dtype=np.int64) for c, (p, _) in read_hits(path, min_reads).items()}


def stranded_sites(sample, min_reads, project="Gale_2020", path=None):
    """{chromosome: (positions, strands)}: each distinct (position, strand) with >= min_reads reads, from a YTAB hit file."""
    import csv
    path = Path(path) if path else ROOT/"output/projects"/project/"create_hit_file"/sample/f"{sample}_hits.txt"
    rows = csv.DictReader(path.read_text().splitlines(), delimiter="\t")
    if not {"Chromosome", "Source", "Hit position", "Hit count"} <= set(rows.fieldnames or []):
        raise ValueError(f"{path}: the learned target model needs a YTAB hit file (Chromosome, Source, Hit position, Hit count)")
    counts = {}
    for r in rows:
        key = (r["Chromosome"], int(r["Hit position"]), r["Source"]); counts[key] = counts.get(key, 0)+int(r["Hit count"])
    return {k for k, n in counts.items() if n >= min_reads}


def learn_target(samples, paths, settings, project, fasta):
    """Target-site model from the distinct sites of all libraries together (each site once)."""
    sites = set().union(*(stranded_sites(s, settings.min_reads, project, p) for s, p in zip(samples, paths)))
    by_chrom = {}
    for c, pos, strand in sorted(sites):by_chrom.setdefault(c, ([], []))[0].append(pos); by_chrom[c][1].append(strand)
    seqs = {n:hh.encode(s) for n, s in hh.read_fasta(fasta).items()}
    return hh.learn_target_model(seqs, {c:(np.array(p), np.array(s)) for c, (p, s) in by_chrom.items()})


def write_target_model(path, model):
    rows = [dict(offset=int(j), **{b:float(model["P"][k, i]) for i, b in enumerate("ACGT")},
                 preferred=("ACGT"[int(np.argmax(model["logodds"][k, :4]))]), information_bits=float(model["information"][k]))
            for k, j in enumerate(model["offsets"])]
    write_table(path, list(rows[0]), rows)


def features(settings, genes, fasta=FASTA, target=None):
    """Window covariates and mappability for this genome; cached per genome (FASTA hash), settings and target model."""
    CACHE.mkdir(parents=True, exist_ok=True)
    fasta = Path(fasta); tag = f"w{settings.window}_k{settings.kmer}_u{settings.unmappable_fraction}"
    if target is not None:tag += "_t"+hashlib.sha256(target["logodds"].tobytes()).hexdigest()[:12]
    legacy = CACHE/f"features_{tag}.npz"                      # C. glabrata cache written before genome-keyed caching
    key = legacy if fasta.resolve() == FASTA.resolve() and legacy.exists() else CACHE/f"features_{hashlib.sha256(fasta.read_bytes()).hexdigest()[:16]}_{tag}.npz"
    if key.exists():
        z = np.load(key, allow_pickle=True); return z["features"].item()
    seqs = {n:hh.encode(s) for n, s in hh.read_fasta(fasta).items()}
    cds = {}
    for g in genes:cds.setdefault(g["chrom"], []).extend(g["blocks"])
    f = hh.window_features(seqs, cds, settings, target); np.savez(key, features=np.array(f, dtype=object)); return f


def design(feats, libraries, settings):
    chroms, names = [], sorted(feats)
    learned = all("target" in feats[c] for c in names)
    seqcov = (lambda f:f["target"]) if learned else (lambda f:np.log1p(f["motif"]))
    allX = np.vstack([np.column_stack([seqcov(feats[c]), feats[c]["at"], feats[c]["outside_cds"]]) for c in names])
    mean, sd = allX.mean(0), allX.std(0)
    for c in names:
        f = feats[c]; n = len(f["start"])
        X = (np.column_stack([seqcov(f), f["at"], f["outside_cds"]])-mean)/sd
        Y = np.column_stack([hh.window_counts(lib.get(c, []), f["start"], settings.window, n) for lib in libraries])
        keep = ~f["masked"] & (f["length"] == settings.window)
        chroms.append(dict(name=c, X=X, Y=Y, keep=keep))
    first = "log mean learned target propensity" if learned else "log1p(nTnnnnAn count)"
    return chroms, dict(mean=mean, sd=sd, covariates=[first, "AT fraction", "fraction outside CDS"])


def check_inputs(feats, genes, libraries, names, settings, allow_high_mask=False):
    """Species-agnostic input checks; returns (genes on reference chromosomes, notes). Raises ValueError on mismatch."""
    reference = set(feats); notes = []
    off = sorted({g["chrom"] for g in genes} - reference)
    kept = [g for g in genes if g["chrom"] in reference]
    if not kept:raise ValueError(f"No annotated gene lies on a FASTA chromosome; GFF chromosomes {off[:5]} vs FASTA {sorted(reference)[:5]}")
    if off:notes.append(f"{len(genes)-len(kept)} genes on {len(off)} chromosomes absent from the FASTA were excluded (e.g. {off[:3]}); "
                        "expected for a haploid reference of a diploid annotation")
    for name, lib in zip(names, libraries):
        total = sum(len(v) for v in lib.values()); outside = sum(len(v) for c, v in lib.items() if c not in reference)
        if not total:raise ValueError(f"{name}: no insertion sites")
        if outside/total > MAX_OFF_REFERENCE:
            raise ValueError(f"{name}: {outside/total:.0%} of insertion sites are on chromosomes not in the FASTA "
                             f"(hit-file chromosomes {sorted(set(lib)-reference)[:5]}; FASTA {sorted(reference)[:5]}). Use the reference the reads were aligned to")
        if outside:notes.append(f"{name}: {outside} sites ({outside/total:.1%}) on chromosomes absent from the FASTA were ignored")
    full = [c for c in reference if len(feats[c]["start"])]
    windows = sum(int((feats[c]["length"] == settings.window).sum()) for c in full)
    masked = sum(int((feats[c]["masked"] & (feats[c]["length"] == settings.window)).sum()) for c in full)
    if windows and masked/windows > MAX_MASKED and not allow_high_mask:
        raise ValueError(f"{masked/windows:.0%} of windows are unmappable (non-unique {settings.kmer}-mers). The FASTA probably contains two "
                         "near-identical haplotypes (e.g. a diploid assembly); use the haploid reference the reads were aligned to "
                         "(for C. albicans: the HapA FASTA). --allow-high-mask overrides")
    return kept, notes


def run(output, samples, settings, label, *, fasta=FASTA, gff=GFF, hits=None, project="Gale_2020", allow_high_mask=False,
        sequence_covariate="motif", target=None):
    """samples: names; hits: matching hit-file paths (default: the project's create_hit_file outputs).
    sequence_covariate: "motif" (v2) or "learned" (v3; `target` reuses a model learned elsewhere, else it is learned
    from these libraries)."""
    start = time.time(); genes = ccg.annotations(Path(gff).resolve())
    paths = list(hits) if hits else [None]*len(samples)
    if sequence_covariate not in ("motif", "learned"):raise ValueError("sequence_covariate must be motif or learned")
    if sequence_covariate == "learned" and target is None:target = learn_target(samples, paths, settings, project, fasta)
    feats = features(settings, genes, fasta, target if sequence_covariate == "learned" else None)
    libraries = [hit_sites(s, settings.min_reads, project, p) for s, p in zip(samples, paths)]
    genes, notes = check_inputs(feats, genes, libraries, samples, settings, allow_high_mask)
    chroms, scaling = design(feats, libraries, settings)
    fitres = hh.fit(chroms, settings)
    if len(fitres["trace"]) >= settings.max_iterations:
        notes.append(f"EM did not converge within {settings.max_iterations} iterations (relative tolerance {settings.tolerance:g}). "
                     "Check fit.json: rates of neighbouring states that are nearly equal mean the data support fewer density classes; "
                     "P_ES_run is usually still usable, but rerun with more iterations before relying on GD/GA classes")
    G, offsets = fitres["gamma"], np.cumsum([0]+[len(c["Y"]) for c in chroms])
    gamma = {c["name"]:G[offsets[i]:offsets[i+1]] for i, c in enumerate(chroms)}; keep = {c["name"]:c["keep"] for c in chroms}
    chain = {c["name"]:fitres["chains"][i] for i, c in enumerate(chroms)}
    rows = hh.gene_table(genes, {c["name"] for c in chroms}, gamma, keep, settings, chain, fitres["A"])
    output.mkdir(parents=True, exist_ok=True)
    if sequence_covariate == "learned":write_target_model(output/"target_model.tsv", target)
    write_table(output/"genes.tsv", list(dict.fromkeys(k for r in rows for k in r)), rows)
    np.savez_compressed(output/"windows.npz", **{f"{c['name']}_gamma":gamma[c["name"]] for c in chroms}, **{f"{c['name']}_Y":c["Y"] for c in chroms},
                        **{f"{c['name']}_keep":c["keep"] for c in chroms})
    calls = {s:sum(r["call"] == s for r in rows) for s in hh.STATES+("uncertain", "N/A")}
    T = sum(len(c["Y"]) for c in chroms); kept = sum(int(c["keep"].sum()) for c in chroms)
    summary = dict(model=hh.MODEL, revision=hh.REVISION_LEARNED if sequence_covariate == "learned" else hh.REVISION, label=label, samples=samples,
                   sequence_covariate=sequence_covariate,
                   target_model=dict(sites=target["sites"], flank=target["flank"], offsets=[int(target["offsets"][0]), int(target["offsets"][-1])],
                                     information_bits=target["information"]) if sequence_covariate == "learned" else None,
                   inputs=dict(fasta=str(fasta), fasta_sha256=sha256(Path(fasta)), gff=str(gff), gff_sha256=sha256(Path(gff)),
                               hits=[str(p) for p in paths] if hits else f"output/projects/{project}/create_hit_file"),
                   notes=notes, settings=asdict(settings),
                   sites_per_library=[int(sum(len(v) for v in lib.values())) for lib in libraries],
                   windows=T, kept_windows=kept, masked_fraction=1-kept/T, covariate_scaling=scaling,
                   r=fitres["r"], a=fitres["a"], b=fitres["b"], psi=fitres["psi"], A=fitres["A"], pi=fitres["pi"],
                   window_state_fraction=G[np.concatenate([c["keep"] for c in chroms])].mean(0), loglik=fitres["loglik"],
                   iterations=fitres["iterations"], converged=fitres["converged"], trace_tail=fitres["trace"][-5:],
                   monotone=bool(np.all(np.diff(fitres["trace"]) > -1e-7*abs(fitres["trace"][0]))),
                   gene_calls=calls, seconds=time.time()-start, implementation_sha256=sha256(Path(hh.__file__)))
    (output/"fit.json").write_text(json.dumps(hh._json_ready(summary), indent=2)+"\n")
    return summary, rows


def project_reference(project):
    config = yaml.safe_load((ROOT/"output/projects"/project/"config/project.yaml").read_text())
    return ROOT/config["reference"]["fasta"], ROOT/config["reference"]["gff"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("output", type=Path)
    p.add_argument("--fasta", type=Path, help="Haploid genome FASTA the reads were aligned to")
    p.add_argument("--gff", type=Path, help="GFF annotation of that genome")
    p.add_argument("--hits", type=Path, action="append", help="Hermes hit file; repeat for each library")
    p.add_argument("--project", help="YTAB project whose create_hit_file outputs and reference to use")
    p.add_argument("--libraries", nargs="+", help="Independent parent libraries of one strain in --project (sample IDs)")
    p.add_argument("--window", type=int, default=100); p.add_argument("--min-reads", type=int, default=1)
    p.add_argument("--each", action="store_true", help="also fit each library alone")
    p.add_argument("--allow-high-mask", action="store_true", help="run even if > 50%% of windows are unmappable")
    p.add_argument("--sequence-covariate", choices=("motif", "learned"), default="motif",
                   help="motif: Hermes nTnnnnAn count (v2); learned: target preference learned from the libraries' sites (v3, any transposon)")
    a = p.parse_args(); settings = hh.Settings(window=a.window, min_reads=a.min_reads)
    if a.hits:
        if not (a.fasta and a.gff):p.error("--hits needs --fasta and --gff")
        samples, fasta, gff, hits = [h.name.removesuffix("_hits.txt") for h in a.hits], a.fasta, a.gff, a.hits
    else:
        if not (a.project and a.libraries):p.error("give --project and --libraries, or --hits files with --fasta and --gff")
        samples = a.libraries
        pf, pg = project_reference(a.project)  # every project's own reference, from its project.yaml
        fasta, gff, hits = a.fasta or pf, a.gff or pg, None
    try:
        s, _ = run(a.output/"joint", samples, settings, "joint", fasta=fasta, gff=gff, hits=hits, project=a.project, allow_high_mask=a.allow_high_mask,
                   sequence_covariate=a.sequence_covariate)
    except ValueError as error:p.exit(1, f"ERROR: {error}\n")
    print(json.dumps(hh._json_ready({k:v for k, v in s.items() if k not in ("settings", "covariate_scaling")}), indent=1))
    for note in s["notes"]:print(f"NOTE: {note}", file=sys.stderr)
    if a.each:
        for i, lib in enumerate(samples):
            s, _ = run(a.output/lib, [lib], settings, lib, fasta=fasta, gff=gff, hits=[hits[i]] if hits else None, project=a.project,
                       allow_high_mask=a.allow_high_mask, sequence_covariate=a.sequence_covariate)
            print(lib, s["gene_calls"], s["r"], s["iterations"])


if __name__ == "__main__":main()
