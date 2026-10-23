#!/usr/bin/env python3
"""Hermes genetic interactions, adapted from TRANSIT's GI method (DeJesus et al. 2017): does a gene's treatment
response differ between two strains?

delta = log2FC(query) - log2FC(reference); interaction classes Aggravating / Alleviating / Suppressive.
Each strain needs >= 2 control/treated pairs: replicate cultures split from one pool (paired in order),
or several pools each split into control and treated. Input is a YTAB `normalize` output (BGC recommended)
whose catalogue gives each sample's library (pool). Validated in a pre-registered simulation
(output/essentiality/gi/benchmark/RESULT.md).

Example:
  ytab genetic-interactions --normalized output/normalization/RUN --normalization betageom \\
      --gff genes.gff --reference genome.fna --output output/gi/RUN \\
      --reference-strain WT --query-strain mutA \\
      --sample P::wt_c1=WT:control --sample P::wt_t1=WT:treated --sample P::wt_c2=WT:control --sample P::wt_t2=WT:treated \\
      --sample P::mu_c1=mutA:control --sample P::mu_t1=mutA:treated --sample P::mu_c2=mutA:control --sample P::mu_t2=mutA:treated
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"src"))
from ytab.essentiality.control_conditioned_gumbel import NONLINEAR_NORMALIZATIONS, NORMALIZATIONS
from ytab.essentiality.genetic_interactions import Settings, run


def sample(text):
    observation, sep, rest = text.rpartition("=")
    strain, sep2, condition = rest.rpartition(":")
    if not (sep and sep2 and observation and strain and condition in ("control", "treated")):
        raise argparse.ArgumentTypeError(f"expected OBSERVATION_ID=STRAIN:control|treated, got {text!r}")
    return observation, strain, condition


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--normalized", required=True); parser.add_argument("--normalization", choices=NORMALIZATIONS+NONLINEAR_NORMALIZATIONS)
    for name in ("gff", "reference", "output"):parser.add_argument("--"+name, required=True)
    parser.add_argument("--sample", type=sample, action="append", required=True, help="OBSERVATION_ID=STRAIN:control|treated; repeat")
    parser.add_argument("--reference-strain", required=True, help="Strain A (e.g. wild type)")
    parser.add_argument("--query-strain", required=True, help="Strain B (e.g. the mutant)")
    d = Settings()
    for name in ("pseudocount", "fdr", "min_delta", "clone_warning", "nterm", "cterm"):
        parser.add_argument("--"+name.replace("_", "-"), type=float, default=getattr(d, name))
    parser.add_argument("--min-sites", type=int, default=d.min_sites)
    parser.add_argument("--empirical-null", action="store_true", help="Efron empirical null (inflate-only); off by default (not used in validation)")
    parser.add_argument("--rope", type=float, default=d.rope,
                        help="Call only interactions with |delta| demonstrably > ROPE (TREAT test, McCarthy & Smyth 2009; 0.5 is a common choice). 0 = test delta != 0")
    parser.add_argument("--correction", choices=["bh", "holm"], default=d.correction,
                        help="bh: false discovery rate (default, validated); holm: family-wise error rate (NOT validated: 8.3%% family-wise error in simulation)")
    parser.add_argument("--confidence", type=float, default=d.confidence, help="Level of the delta interval (delta_lower/upper)")
    parser.add_argument("--include-stop", action="store_true")
    a = parser.parse_args()
    settings = Settings(pseudocount=a.pseudocount, fdr=a.fdr, min_delta=a.min_delta, min_sites=a.min_sites, clone_warning=a.clone_warning,
                        empirical_null=a.empirical_null, rope=a.rope, correction=a.correction, confidence=a.confidence, nterm=a.nterm, cterm=a.cterm, ignore_stop=not a.include_stop)
    try:
        m = run(a.normalized, a.gff, a.reference, a.output, a.sample, a.reference_strain, a.query_strain, settings, a.normalization)
    except (ValueError, KeyError, OSError) as error:parser.exit(1, f"ERROR: {error}\n")
    d = m["diagnostics"]
    print(json.dumps(dict(status=m["status"], analysis=m["analysis"], revision=m["algorithm_revision"], reference=d["reference_strain"],
                          query=d["query_strain"], tested=d["tested"], calls=d["calls"], site_variance_scale=d["site_variance_scale"]), indent=2))
    if d["warning"]:print(f"WARNING: {d['warning']}", file=sys.stderr)
    if a.correction == "holm":print("WARNING: --correction holm is not validated (simulated family-wise error 8.3% vs 5%); prefer bh", file=sys.stderr)


if __name__ == "__main__":main()
