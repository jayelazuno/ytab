#!/usr/bin/env python3
"""Hermes multi-condition tests: site-restricted permutation ANOVA and all-pairs resampling.

For one library (pool) split into several conditions (e.g. control, FLU, Zn, other drugs), or several
libraries each split into the same conditions (--stratify-libraries). Input is a YTAB `normalize` output;
BGC (betageom) is recommended. Assign every sample to a condition with --sample OBSERVATION_ID=CONDITION.

  anova      does a gene respond to any condition?  (q <= FDR -> "varies"; top_condition = largest |log2FC|)
  all-pairs  which conditions differ, and in which direction? (resampling v3 for every pair, one FDR family)

If the ANOVA replicate-noise warning is printed, use the all-pairs results: ANOVA p-values are
anti-conservative when one condition's replicates are much noisier (pre-registered evaluation:
output/essentiality/multicondition/benchmark_v2/RESULT.md).

Example:
  ytab multicondition --normalized output/normalization/RUN --normalization betageom \\
      --gff genes.gff --reference genome.fna --output output/multicondition/RUN \\
      --sample P::ctrl1=control --sample P::ctrl2=control --sample P::flu1=FLU --sample P::flu2=FLU \\
      --sample P::zn1=Zn --sample P::zn2=Zn --reference-condition control
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"src"))
from ytab.essentiality.control_conditioned_gumbel import NONLINEAR_NORMALIZATIONS, NORMALIZATIONS
from ytab.essentiality.multicondition import Settings, run


def sample(text):
    observation, sep, condition = text.rpartition("=")
    if not sep or not observation or not condition:raise argparse.ArgumentTypeError(f"expected OBSERVATION_ID=CONDITION, got {text!r}")
    return observation, condition


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--normalized", required=True, help="YTAB normalize output directory, or one method folder inside it")
    parser.add_argument("--normalization", choices=NORMALIZATIONS+NONLINEAR_NORMALIZATIONS, help="Method folder when --normalized is a run directory (betageom recommended)")
    for name in ("gff", "reference", "output"):
        parser.add_argument("--"+name, required=True)
    parser.add_argument("--sample", type=sample, action="append", required=True, help="OBSERVATION_ID=CONDITION; repeat for every sample")
    parser.add_argument("--reference-condition", default="", help="Condition used for recentering and log2FC (default: the first condition listed)")
    parser.add_argument("--test", choices=["both", "anova", "all-pairs"], default="both")
    parser.add_argument("--stratify-libraries", action="store_true", help="Several libraries (pools), each split into the conditions: permute only within library_id")
    parser.add_argument("--independent-libraries", action="store_true",
                        help="Conditions are separate, independently made libraries (e.g. one per strain): permute among all samples at each site, as TRANSIT")
    defaults = Settings()
    parser.add_argument("--recenter", choices=["iterative-median", "trimmed-sum", "none"], default=defaults.recenter,
                        help="Per-condition rescaling to the reference (iterative-median is the validated default)")
    for name in ("min_sites", "seed"):
        parser.add_argument("--"+name.replace("_", "-"), type=int, default=getattr(defaults, name))
    for name in ("fdr", "pseudocount", "nterm", "cterm", "asymmetry_warning"):
        parser.add_argument("--"+name.replace("_", "-"), type=float, default=getattr(defaults, name))
    parser.add_argument("--include-stop", action="store_true")
    args = parser.parse_args()
    settings = Settings(reference=args.reference_condition, recenter=args.recenter, min_sites=args.min_sites, seed=args.seed, fdr=args.fdr,
                        pseudocount=args.pseudocount, nterm=args.nterm, cterm=args.cterm, asymmetry_warning=args.asymmetry_warning,
                        ignore_stop=not args.include_stop)
    tests = ("anova", "all-pairs") if args.test == "both" else (args.test,)
    try:
        result = run(args.normalized, args.gff, args.reference, args.output, args.sample, settings, args.normalization, args.stratify_libraries, tests, independent=args.independent_libraries)
    except (ValueError, KeyError, OSError) as error:parser.exit(1, f"ERROR: {error}\n")
    d = result["diagnostics"]
    print(json.dumps(dict(status=result["status"], analysis=result["analysis"], revision=result["algorithm_revision"],
                          normalization=result["normalization"]["method"], conditions=result["conditions"], reference=result["reference_condition"],
                          anova_calls=(d.get("anova") or {}).get("calls"),
                          all_pairs_calls={k:v["calls"] for k, v in (d.get("all_pairs") or {}).items() if isinstance(v, dict) and "calls" in v},
                          files=result["files"], recommendation=result["recommendation"]), indent=2))
    warning = (d.get("anova") or {}).get("warning")
    if warning:print(f"WARNING: {warning}. Use the all-pairs results.", file=sys.stderr)


if __name__ == "__main__":main()
