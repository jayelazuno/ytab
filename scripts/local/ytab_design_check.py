#!/usr/bin/env python3
"""Which analyses a project's experimental design supports, which libraries each uses, and why the others are not possible.

Reads the project's sample sheet (included samples: role parent / control / treated, condition, strain = background, pool)
and applies each analysis's own minimum design (docs/design_guide.md):

  essentiality (HMM)        >= 1 unselected parent library of a strain; denser = better (pool independent inductions)
  treated vs control:
    resampling              >= 1 control and >= 1 treated library of a strain and condition, from the same pool
    z-score (Gale 2020)     as resampling, and >= 2 controls in the strain (its noise model compares two controls)
    control-cond. Gumbel    >= 2 controls (depletion only)
  multi-condition           >= 2 conditions (one the reference) in a pool, or strains made as separate libraries
  genetic interactions      2 strains x (control, treated) in the same condition, >= 2 control/treated pairs per strain
  library correlation       >= 2 libraries

Sequencing runs of one library (same strain, pool, condition and culture replicate) are flagged: merge them first
(`ytab merge-runs`), or they would be counted as replicates.

  ytab design-check --project-config output/projects/P/config/project.yaml          # report
  ytab design-check --project-config ... --json                                     # machine-readable
"""
import argparse
import csv
from collections import defaultdict
from itertools import combinations
import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]


def resolve(path, base):
    path = Path(path)
    return path if path.is_absolute() else (base/path).resolve()


def role(r):
    value = str(r.get("library_role") or "").strip().lower()
    return value if value in ("parent", "control", "treated") else ("parent" if str(r.get("condition", "")).lower() in ("parent", "") else "treated")


def check(rows, project, config):
    included = [r for r in rows if str(r.get("include", "True")).strip().lower() not in ("false", "0", "no")]
    cfg = f"--project-config {config}"
    out = []
    add = lambda analysis, status, samples, why, command="":out.append(dict(analysis=analysis, status=status, samples=samples, why=why, command=command))

    runs = defaultdict(list)
    for r in included:runs[(r.get("background", ""), r.get("pool", ""), r.get("condition", ""), r.get("culture_replicate", ""))].append(r["sample"])
    for (bg, pool, cond, rep), ss in runs.items():
        if len(ss) > 1:
            add("Merge sequencing runs", "needed", ss, f"{len(ss)} runs share strain {bg}, pool {pool}, condition {cond}, replicate {rep or '-'}: if they are lanes or "
                "digests of one library, merge them; if they are separate cultures, give each its own culture_replicate",
                f"ytab merge-runs {cfg} --library {pool}_{cond}_lib{ss[0]} --runs {' '.join(ss)}")

    strains = defaultdict(list)
    for r in included:strains[r.get("background", "")].append(r)
    for bg, rs in strains.items():
        parents = [r["sample"] for r in rs if role(r) == "parent"]
        unselected = parents or [r["sample"] for r in rs if role(r) == "control"]
        if unselected:
            pools = len({r["pool"] for r in rs if r["sample"] in unselected})
            add(f"Essentiality (HMM), strain {bg}", "ready", unselected,
                f"{len(unselected)} unselected {'parent' if parents else 'control'} libraries from {pools} pool(s); aim for >= 36 insertion sites per coding kb "
                "(Hermes; Ac/Ds needs more), reached by pooling independent inductions", f"ytab essentiality {cfg} --libraries {' '.join(unselected)}")
        else:
            add(f"Essentiality (HMM), strain {bg}", "not possible", [], "no unselected (parent or control) library of this strain")

    pairs = defaultdict(lambda:defaultdict(lambda:{"control":[], "treated":[]}))   # (strain, condition) -> pool -> samples
    for bg, rs in strains.items():
        for t in (r for r in rs if role(r) == "treated"):
            same = [r for r in rs if r["pool"] == t["pool"] and role(r) != "treated"]
            controls = [r for r in same if role(r) == "control"] or same
            g = pairs[(bg, t["condition"])][t["pool"]]
            g["treated"].append(t["sample"]); g["control"] += [c["sample"] for c in controls if c["sample"] not in g["control"]]
    for (bg, cond), pools in sorted(pairs.items()):
        controls = sorted({s for p in pools.values() for s in p["control"]}); treated = sorted({s for p in pools.values() for s in p["treated"]})
        name = f"Treated vs control, strain {bg}, {cond}"
        if not controls:
            others = [r["sample"] for r in strains[bg] if role(r) != "treated"]
            add(f"{name} (resampling)", "possible with caution" if others else "not possible", others + treated,
                "no control or parent in the treated library's pool; controls from other pools share no mutants, so only much denser libraries "
                "detect effects (resampling with independent libraries)" if others else "no control or parent of this strain")
            continue
        add(f"{name} (resampling)", "ready", controls + treated, f"{len(controls)} control(s) vs {len(treated)} treated from {len(pools)} pool(s); "
            "replicates per side add power", "ytab comparison-design --controls all ...; ytab treated-vs-control --method resampling ...")
        all_controls = {s for (b, _), ps in pairs.items() if b == bg for p in ps.values() for s in p["control"]}
        if len(all_controls) >= 2:
            add(f"{name} (z-score)", "ready", controls + treated, f"{len(all_controls)} controls in strain {bg} for the noise model", f"ytab fitness {cfg} --comparisons ...")
        else:
            add(f"{name} (z-score)", "not possible", controls + treated, f"needs >= 2 controls in strain {bg} for its noise model (has {len(all_controls)}); use resampling")
        add(f"{name} (control-conditioned Gumbel)", "ready" if len(controls) >= 2 else "not possible", controls + treated,
            "depletion only" if len(controls) >= 2 else f"needs >= 2 controls (has {len(controls)})")

    for bg, rs in strains.items():
        for pool in sorted({r["pool"] for r in rs}):
            conds = sorted({r["condition"] for r in rs if r["pool"] == pool})
            if len(conds) >= 2:
                ref = next((r["condition"] for r in rs if r["pool"] == pool and role(r) == "control"), None) or next(
                    (r["condition"] for r in rs if r["pool"] == pool and role(r) == "parent"), conds[0])
                add(f"Multi-condition, strain {bg}, pool {pool}", "ready", [r["sample"] for r in rs if r["pool"] == pool],
                    f"{len(conds)} conditions ({', '.join(conds)}); reference {ref}", f"ytab multicondition-project {cfg} --reference-condition {ref} --sample ... ")
    by_condition = defaultdict(set)
    for r in included:by_condition[r["condition"]].add(r.get("background", ""))
    for cond, bgs in sorted(by_condition.items()):
        if len(bgs) >= 2:
            add(f"Strain comparison, {cond}", "ready", [r["sample"] for r in included if r["condition"] == cond],
                f"strains {', '.join(sorted(bgs))} in the same condition: multi-condition with strains as conditions; --independent-libraries when each "
                "strain is its own library", f"ytab multicondition-project {cfg} --independent-libraries --sample P::LIB=STRAIN ...")

    for (a, b) in combinations(sorted(strains), 2):
        for cond in sorted({c for (s, c) in pairs if s == a} & {c for (s, c) in pairs if s == b}):
            n = {s:sum(min(len(p["control"]), len(p["treated"])) for p in pairs[(s, cond)].values()) for s in (a, b)}
            ok = all(v >= 2 for v in n.values())
            add(f"Genetic interactions, {b} vs {a}, {cond}", "ready" if ok else "not possible",
                [x for s in (a, b) for p in pairs[(s, cond)].values() for x in p["control"] + p["treated"]],
                f"control/treated pairs: {a} {n[a]}, {b} {n[b]}" + ("" if ok else "; each strain needs >= 2 (two pools each split into control and treated, "
                "or replicate control/treated cultures of one pool)"), "ytab genetic-interactions --reference-strain A --query-strain B --sample P::LIB=STRAIN:control ..." if ok else "")
    if not any(r["analysis"].startswith("Genetic interactions") for r in out) and len(strains) >= 2:
        add("Genetic interactions", "not possible", [], "no treated condition shared by two strains")
    if len(strains) < 2:add("Genetic interactions", "not possible", [], "one strain only; needs a reference and a query strain (e.g. wild type and a deletion)")
    add("Library correlation", "ready" if len(included) >= 2 else "not possible", [r["sample"] for r in included], f"{len(included)} libraries",
        f"ytab library-correlation --project {project} --output output/projects/{project}/library_correlation/RUN")
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-config", required=True, type=Path)
    p.add_argument("--json", action="store_true")
    a = p.parse_args()
    cfg = yaml.safe_load(a.project_config.read_text()); base = resolve(cfg.get("repo_root", "."), a.project_config.parent)
    sheet = resolve(cfg.get("sample_sheet") or a.project_config.parent/"sample_sheet.csv", base)
    rows = list(csv.DictReader(open(sheet, newline="", encoding="utf-8")))
    config = str(a.project_config.resolve().relative_to(ROOT)) if a.project_config.resolve().is_relative_to(ROOT) else str(a.project_config)
    result = check(rows, cfg.get("project_id", ""), config)
    if a.json:print(json.dumps(result, indent=1)); return
    print(f"Experimental design of {cfg.get('project_id', '')}: what each analysis can use (docs/design_guide.md)\n")
    for r in result:
        mark = {"ready":"OK ", "needed":"!! ", "not possible":"-- "}.get(r["status"], "?  ")
        samples = " ".join(r["samples"][:12]) + (f" (+{len(r['samples'])-12})" if len(r["samples"]) > 12 else "")
        print(f"{mark}{r['analysis']}: {r['status']}. {r['why']}")
        if samples:print(f"     libraries: {samples}")
        if r["command"]:print(f"     {r['command']}")


if __name__ == "__main__":main()
