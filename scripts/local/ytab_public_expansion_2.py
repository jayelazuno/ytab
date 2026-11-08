#!/usr/bin/env python3
"""Transposon chemistry and the density rule beyond the four species: eight new public datasets, 84 runs.

Implements output/public_expansion_2/PREREGISTRATION.md (frozen; every phase refuses to run if it changed).
  libraries: the library table, merged hit files for multi-run libraries, insertion density per library and the
             per-species insertable-site ceiling (TTAA for piggyBac, every position for the hAT transposons) -> H1
  fit:       HMM v3 fits on the pre-registered density sets per species (reads no labels)
  conditions: treated-vs-control resampling under total reads, TTR and BGC on the two new condition designs -> H3
  evaluate:  reads the truth sets; writes evaluation.tsv, per_library.tsv and result.json -> H2

The library rules follow R24: one library per (project, pool, condition, culture replicate), its runs merged, and
density as the median insertion sites per coding kb over the genes on the mapped chromosomes.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src")); sys.path.append(str(ROOT/"scripts/local"))
from ytab.normalization.ttr import sha256, write_table   # noqa: E402
from ytab.mapping import transposons   # noqa: E402

BENCH = ROOT/"output/public_expansion_2"
THRESHOLD_KB = 36          # R16's essentiality density recommendation, the number H2 tests out of lineage
PROJECTS = {               # project: (species, transposon)
    "Patterson_2018_lipolytica":("lipolytica", "hermes"),
    "Billmyre_2025_neoformans":("neoformans", "ac_ds"),
    "Gao_2021_auris":("auris", "piggybac"),
    "Song_2026_auris":("auris", "piggybac"),
    "Chow_2024_glabrata":("glabrata_CBS138", "piggybac"),
    "Yao_2026_phaffii":("phaffii", "piggybac"),
    "Nickels_2026_glabrata_tolerance":("glabrata", "hermes"),
    "Nickels_2026_glabrata_FKS2":("glabrata", "hermes"),
}
TARGET_SITE = {"piggybac":"TTAA", "hermes":None, "ac_ds":None}   # None: the transposon inserts at (nearly) any position
# H1's comparison set: the C. glabrata Hermes libraries already measured, from R24's frozen library table.
HERMES_GLABRATA = ROOT/"output/public_expansion/libraries.tsv"


def check_preregistration():
    frozen = (BENCH/"PREREGISTRATION.sha256").read_text().split()
    if sha256(BENCH/"PREREGISTRATION.md") != frozen[0]:
        raise SystemExit("PREREGISTRATION.md changed after it was frozen; refusing to run")
    return dict(sha256=frozen[0], frozen_utc=frozen[1])


def reference(project):
    cfg = yaml.safe_load((ROOT/"output/projects"/project/"config/project.yaml").read_text())["reference"]
    return ROOT/cfg["fasta"], ROOT/cfg["gff"]


def run_hits(project, run):
    return ROOT/"output/projects"/project/"create_hit_file"/run/f"{run}_hits.txt"


def library_table(projects=None):
    """One library per (project, pool, condition, culture replicate); its runs are the sequencing runs of it."""
    libs = {}
    for project, (species, transposon) in PROJECTS.items():
        if projects and project not in projects:continue
        sheet = ROOT/"output/projects"/project/"config/sample_sheet.csv"
        if not sheet.exists():continue
        for row in csv.DictReader(open(sheet)):
            if row.get("include", "True") != "True":continue
            if not run_hits(project, row["sample"]).exists():continue      # not mapped yet; reported by `status`
            key = (project, row["pool"], row["condition"], row["culture_replicate"])
            # A density is a property of the library as built, so only the arms that were not selected can carry
            # one: `parent`, and the `control` arm of a treated-vs-control design. The role column already says
            # which, so no condition name is matched (R24 rule 4 matched "no_drug"; this does not).
            lib = libs.setdefault(key, dict(library="", project=project, species=species, transposon=transposon,
                                            background=row["background"], pool=row["pool"], condition=row["condition"],
                                            role=row["library_role"],
                                            unselected=row["library_role"] in ("parent", "control"), runs=[]))
            lib["runs"].append(row["sample"])
    for lib in libs.values():
        lib["runs"].sort()
        lib["library"] = lib["runs"][0] if len(lib["runs"]) == 1 else f"{lib['pool']}_{lib['condition']}_lib{lib['runs'][0]}"
    return sorted(libs.values(), key=lambda l:(l["project"], l["library"]))


def merged_hits(lib):
    """One hit file per library; runs summed per (chromosome, strand, position), as R24 rule 2."""
    if len(lib["runs"]) == 1:return run_hits(lib["project"], lib["runs"][0])
    out = BENCH/"hits"/lib["project"]/f"{lib['library']}_hits.txt"
    if out.exists():return out
    out.parent.mkdir(parents=True, exist_ok=True); rows, header = {}, None
    for run in lib["runs"]:
        with open(run_hits(lib["project"], run)) as f:
            r = csv.reader(f, delimiter="\t"); header = next(r)
            pos, cnt = header.index("Hit position"), header.index("Hit count")
            for row in r:
                key = (row[0], row[1], row[pos]); base = rows.setdefault(key, row[:])
                if base is not row:base[cnt] = str(int(base[cnt])+int(row[cnt]))
    with open(out, "w") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(header); w.writerows(rows.values())
    return out


def coding_genes(project):
    """The annotated genes on the chromosomes the reads were mapped to, with their CDS blocks."""
    from ytab.essentiality import control_conditioned_gumbel as ccg
    fasta, gff = reference(project)
    chroms = {l[1:].split()[0] for l in open(fasta) if l.startswith(">")}
    return [g for g in ccg.annotations(gff) if g["blocks"] and g["chrom"] in chroms]


def density(hit_file, genes, genome=None, site=None):
    """Median insertion sites per coding kb over genes: distinct hit positions inside the CDS blocks, either strand
    (as R24 rule 5 and ytab_condition_density.density).

    A transposon that inserts only at a fixed target site can have no insertion anywhere else, so for those a
    position off the motif is an artefact of read processing and is not counted. The motif comes from the
    transposon's own profile, so the rule holds for any transposon that declares one and changes nothing for the
    hAT elements, which declare none.
    """
    sites = {}
    with open(hit_file) as handle:
        r = csv.reader(handle, delimiter="\t"); h = next(r); pos = h.index("Hit position")
        for row in r:
            if site and not transposons.on_target(genome.get(row[0], ""), int(row[pos]), site):continue
            sites.setdefault(row[0], set()).add(int(row[pos]))
    arrays = {c:np.array(sorted(s)) for c, s in sites.items()}; per_kb = []
    for g in genes:
        a = arrays.get(g["chrom"]); length = sum(b-s+1 for s, b in g["blocks"])
        n = sum(int(np.searchsorted(a, b, "right")-np.searchsorted(a, s, "left")) for s, b in g["blocks"]) if a is not None else 0
        per_kb.append(n/length*1000)
    return float(np.median(per_kb))


def sequences(fasta):
    name, chunks, out = None, [], {}
    for line in open(fasta):
        if line.startswith(">"):
            if name:out[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        else:chunks.append(line.strip())
    if name:out[name] = "".join(chunks).upper()
    return out


def ceiling(project, genes, motif):
    """The insertable-site ceiling in the unit the densities are counted in — distinct *recordable* hit positions per
    coding kb if every insertable position in the CDS were occupied.

    A hit file records one position per read end, so an occupied target site appears as two positions: measured in
    these libraries (`DEVIATIONS.md`, 2026-10-31), a piggyBac insertion in a TTAA starting at q is recorded as q on
    the W strand and q+3 on the C strand, and the density counts both. The piggyBac ceiling is therefore the size of
    {q, q+3 : TTAA at q} inside the CDS, not the number of TTAA sites. For a transposon with no target restriction
    every position is recordable, so the ceiling is 1000 per kb."""
    if motif is None:return 1000.
    fasta, _ = reference(project); genome = sequences(fasta); per_kb = []
    duplication = len(motif)-1                      # TTAA: the two records are 3 apart
    for g in genes:
        seq = genome.get(g["chrom"], "")
        length = sum(b-s+1 for s, b in g["blocks"]); recordable = set()
        for s, b in g["blocks"]:
            block = seq[s-1:b]                      # GFF coordinates are 1-based inclusive
            start = block.find(motif)
            while start != -1:
                recordable |= {s+start, s+start+duplication}
                start = block.find(motif, start+1)
        inside = lambda p:any(lo <= p <= hi for lo, hi in g["blocks"])
        per_kb.append(sum(1 for p in recordable if inside(p))/length*1000)
    return float(np.median(per_kb))


def hermes_glabrata_median():
    """H1's comparison number: the median density of the C. glabrata Hermes libraries measured in R24."""
    if not HERMES_GLABRATA.exists():return None, 0
    # Restricted to R24's own unselected libraries (its rule 4), so both arms of H1 measure the same thing.
    rows = [r for r in csv.DictReader(open(HERMES_GLABRATA), delimiter="\t")
            if r["species"].startswith("glabrata") and r["unselected"] == "True"]
    values = [float(r["sites_per_coding_kb"]) for r in rows]
    return (float(np.median(values)), len(values)) if values else (None, 0)


def status():
    """What is mapped and what is still missing, per project: the pull from the HPC run is incremental."""
    for project in PROJECTS:
        sheet = ROOT/"output/projects"/project/"config/sample_sheet.csv"
        if not sheet.exists():print(f"{project:34s} no sample sheet"); continue
        rows = [r for r in csv.DictReader(open(sheet)) if r.get("include", "True") == "True"]
        have = [r for r in rows if run_hits(project, r["sample"]).exists()]
        print(f"{project:34s} {len(have):3d}/{len(rows):3d} runs with hit files" + ("" if len(have) == len(rows) else "   INCOMPLETE"))


def libraries(projects=None):
    prereg = check_preregistration()
    libs = library_table(projects)
    if not libs:raise SystemExit("No library has a hit file yet; run `status` to see what is still missing")
    genes, ceilings, genomes, target = {}, {}, {}, {}
    for lib in libs:
        if lib["project"] not in genes:
            genes[lib["project"]] = coding_genes(lib["project"])
            ceilings[lib["project"]] = ceiling(lib["project"], genes[lib["project"]], TARGET_SITE[lib["transposon"]])
            target[lib["project"]] = transposons.target_site(lib["transposon"], ROOT)
            genomes[lib["project"]] = sequences(reference(lib["project"])[0]) if target[lib["project"]] else None
        lib["hits"] = str(merged_hits(lib).relative_to(ROOT))
        lib["coding_genes"] = len(genes[lib["project"]])
        lib["sites_per_coding_kb"] = density(ROOT/lib["hits"], genes[lib["project"]],
                                             genomes[lib["project"]], target[lib["project"]])
        lib["insertable_sites_per_coding_kb"] = ceilings[lib["project"]]
        lib["share_of_ceiling"] = lib["sites_per_coding_kb"]/lib["insertable_sites_per_coding_kb"]
        lib["meets_threshold"] = lib["sites_per_coding_kb"] >= THRESHOLD_KB
        print(f"{lib['project']:34s} {lib['library']:28s} runs {len(lib['runs'])}  "
              f"{lib['sites_per_coding_kb']:7.2f} / {lib['insertable_sites_per_coding_kb']:7.2f} per coding kb  "
              f"({lib['share_of_ceiling']*100:5.1f}% of ceiling)", flush=True)
    hermes_median, hermes_n = hermes_glabrata_median()
    rows = [dict(lib, runs=" ".join(lib["runs"])) for lib in libs]
    write_table(BENCH/"libraries.tsv", list(rows[0]), rows)
    (BENCH/"libraries.json").write_text(json.dumps(dict(
        preregistration=prereg, script_sha256=sha256(Path(__file__)), threshold_sites_per_coding_kb=THRESHOLD_KB,
        hermes_glabrata_reference=dict(source=str(HERMES_GLABRATA.relative_to(ROOT)), libraries=hermes_n, median=hermes_median),
        complete=projects is None, libraries=libs), indent=1)+"\n")
    h1_summary(libs, hermes_median, hermes_n)


def h1_summary(libs, hermes_median, hermes_n):
    """H1 as registered: are the C. glabrata piggyBac libraries sparser than that species' Hermes libraries, and
    below the 36 sites per coding kb the essentiality recommendation asks for?"""
    pb = [l for l in libs if l["project"] == "Chow_2024_glabrata" and l["unselected"]]
    if not pb or hermes_median is None:
        print("\nH1 not evaluable yet (needs Chow_2024_glabrata and R24's library table)"); return
    median = float(np.median([l["sites_per_coding_kb"] for l in pb]))
    below_hermes, below_threshold = median < hermes_median, median < THRESHOLD_KB
    print(f"\nH1  C. glabrata piggyBac libraries (n={len(pb)}): median {median:.2f} sites per coding kb")
    print(f"    C. glabrata Hermes libraries (n={hermes_n}, R24):  median {hermes_median:.2f}")
    print(f"    below the Hermes median: {below_hermes}; below {THRESHOLD_KB}: {below_threshold}")
    print(f"    H1 {'supported' if below_hermes and below_threshold else 'not supported'}")
    print(f"    ceiling: piggyBac {pb[0]['insertable_sites_per_coding_kb']:.2f} per coding kb (TTAA), "
          f"Hermes 1000.00 (every position)")
    per_project = {}
    for l in libs:
        if not l["unselected"]:continue
        per_project.setdefault((l["project"], l["species"], l["transposon"]), []).append(l["sites_per_coding_kb"])
    print("\nDensity per project (median over its unselected libraries):")
    for (project, species, transposon), values in sorted(per_project.items()):
        print(f"    {project:34s} {species:18s} {transposon:9s} n={len(values):3d}  {np.median(values):7.2f}")
    (BENCH/"h1.json").write_text(json.dumps(dict(
        piggybac_glabrata=dict(libraries=len(pb), median=median, ceiling=pb[0]["insertable_sites_per_coding_kb"]),
        hermes_glabrata=dict(libraries=hermes_n, median=hermes_median, ceiling=1000.),
        threshold=THRESHOLD_KB, below_hermes_median=below_hermes, below_threshold=below_threshold,
        supported=bool(below_hermes and below_threshold)), indent=1)+"\n")


# ---------------------------------------------------------------- H3: the two new condition designs

CONDITION_PROJECTS = ("Song_2026_auris", "Billmyre_2025_neoformans")
CONDITION_NORMALIZATIONS = ("totreads", "ttr", "betageom")
H3_DIFFERENCE, H3_JACCARD = 0.25, 0.90       # R25's criteria, as registered


def load_libraries():
    return json.loads((BENCH/"libraries.json").read_text())["libraries"]


def condition_design(libs):
    """Rows (treated library, control library): each treated library against every control library of its own
    background and pool, so a dose series is compared within the pool it came from."""
    rows = []
    for t in libs:
        if t["role"] != "treated":continue
        same = [c for c in libs if c["background"] == t["background"] and c["pool"] == t["pool"]]
        controls = [c for c in same if c["role"] == "control"] or [c for c in same if c["role"] == "parent"]
        rows += [(t, c) for c in controls]
    return rows


def write_csv(path, rows):
    with open(path, "w") as f:
        w = csv.DictWriter(f, list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)


def derived_project(project, libs):
    """A copy of the project's config whose samples are the libraries (merged hit files) and whose design is
    condition_design(); the project itself is never changed. As R24."""
    import yaml as _yaml
    out = BENCH/"projects"/project; config = out/"config"; config.mkdir(parents=True, exist_ok=True)
    cfg = _yaml.safe_load((ROOT/"output/projects"/project/"config/project.yaml").read_text())
    cfg.update(repo_root="../../../../..", output_project_dir=str(out.relative_to(ROOT)),
               output_export_dir=str((out/"exports").relative_to(ROOT)),
               sample_sheet=str((config/"sample_sheet.csv").relative_to(ROOT)), samples=[])
    (config/"project.yaml").write_text(_yaml.safe_dump(cfg, sort_keys=False))
    (config/"reference_resolved.json").write_text((ROOT/"output/projects"/project/"config/reference_resolved.json").read_text())
    write_csv(config/"sample_sheet.csv",
              [dict(sample=l["library"], library_role=l["role"], condition=l["condition"], background=l["background"],
                    pool=l["pool"], display_name=l["library"], include="True", hit_file=l["hits"],
                    runs=" ".join(l["runs"])) for l in libs])
    pairs = condition_design(libs)
    write_csv(config/"comparison_design.csv",
              [dict(comparison_id=f"{t['library']}_vs_{c['library']}", pair_id=f"{t['background']}_{t['pool']}",
                    parent_sample=c["library"], treated_sample=t["library"], background=t["background"], pool=t["pool"],
                    parent_condition=c["condition"], treated_condition=t["condition"], include="true", notes="")
               for t, c in pairs])
    return config/"project.yaml", sorted({l["library"] for t, c in pairs for l in (t, c)})


def conditions():
    """H3: treated-vs-control resampling under each of the three normalizations, on each new condition design."""
    import subprocess
    prereg = check_preregistration(); libs = load_libraries(); runs = {}
    for project in CONDITION_PROJECTS:
        mine = [l for l in libs if l["project"] == project]
        if not mine:print(f"{project}: no libraries yet, skipped"); continue
        config, samples = derived_project(project, mine); base = BENCH/"projects"/project
        runs[project] = {}
        for method in CONDITION_NORMALIZATIONS:
            norm, tvc = base/f"normalization/{method}_together", base/f"treated_vs_control/resampling_{method}"
            steps = [("normalize", norm, [ROOT/"scripts/local/ytab_normalize_project.py", "--project-config", config,
                                          "--method", method, "--output", norm, "--normalize-together",
                                          *[x for s in samples for x in ("--sample", s)]]),
                     ("resampling", tvc, [ROOT/"scripts/local/ytab_treated_vs_control.py", "--project-config", config,
                                          "--method", "resampling", "--normalization-run", norm,
                                          "--normalization", method, "--output", tvc])]
            for name, out, cmd in steps:
                if (out/"manifest.json").exists() or (out/"genes.tsv").exists():continue
                log = base/f"{method}_{name}.log"; log.parent.mkdir(parents=True, exist_ok=True)
                print(project, method, name, flush=True)
                with open(log, "w") as handle:
                    r = subprocess.run([sys.executable]+[str(c) for c in cmd], stdout=handle,
                                       stderr=subprocess.STDOUT, cwd=ROOT)
                if r.returncode:raise SystemExit(f"{project} {method} {name} failed; see {log}")
            runs[project][method] = dict(normalization=str(norm.relative_to(ROOT)), tests=str(tvc.relative_to(ROOT)))
    (BENCH/"conditions.json").write_text(json.dumps(dict(preregistration=prereg, script_sha256=sha256(Path(__file__)),
                                                         normalizations=list(CONDITION_NORMALIZATIONS), runs=runs), indent=1)+"\n")
    h3_summary(runs)


def depleted(genes_tsv):
    """The genes a treated-vs-control run calls depleted. `call` is the run's own verdict per gene and reads
    "depleted", "enriched", "not significant" or "not tested"; a comparison writes one row per gene per comparison,
    so a gene depleted in any comparison of the design counts once."""
    with open(genes_tsv) as handle:
        return {row["ORF"] for row in csv.DictReader(handle, delimiter="\t") if row["call"].strip() == "depleted"}


def h3_summary(runs):
    """H3 as registered: does the depleted-gene count under total reads differ from BGC by at least 25% in at least
    one design, and do the BGC calls agree with TTR at Jaccard >= 0.90?"""
    rows = []
    for project, methods in runs.items():
        calls = {}
        for method, where in methods.items():
            path = ROOT/where["tests"]/"genes.tsv"
            if path.exists():calls[method] = depleted(path)
        if {"totreads", "betageom", "ttr"} - set(calls):print(f"{project}: incomplete, H3 not evaluable"); continue
        n_tot, n_bgc = len(calls["totreads"]), len(calls["betageom"])
        difference = abs(n_tot-n_bgc)/max(n_bgc, 1)
        union = calls["betageom"] | calls["ttr"]
        jaccard = len(calls["betageom"] & calls["ttr"])/len(union) if union else 1.
        rows.append(dict(project=project, depleted_totreads=n_tot, depleted_ttr=len(calls["ttr"]), depleted_bgc=n_bgc,
                         relative_difference=difference, meets_difference=difference >= H3_DIFFERENCE,
                         jaccard_bgc_ttr=jaccard, meets_jaccard=jaccard >= H3_JACCARD))
        print(f"{project}: depleted totreads {n_tot}, TTR {len(calls['ttr'])}, BGC {n_bgc}; "
              f"|totreads-BGC|/BGC = {difference:.3f}; Jaccard(BGC,TTR) = {jaccard:.3f}")
    if not rows:return
    supported = any(r["meets_difference"] for r in rows) and all(r["meets_jaccard"] for r in rows)
    print(f"H3 {'supported' if supported else 'not supported'}")
    (BENCH/"h3.json").write_text(json.dumps(dict(criteria=dict(relative_difference=H3_DIFFERENCE, jaccard=H3_JACCARD),
                                                 designs=rows, supported=supported), indent=1)+"\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("phase", choices=["status", "libraries", "conditions"])
    p.add_argument("--projects", nargs="*", help="limit to these projects (a partial pull); recorded in libraries.json")
    a = p.parse_args()
    if a.phase == "status":status(); return
    if a.phase == "conditions":conditions(); return
    libraries(a.projects)


if __name__ == "__main__":
    main()
