#!/usr/bin/env python3
"""Experimental control-conditioned Hermes Gumbel on combined replicate tracks.

Input is a YTAB `normalize` output (ttr default, totreads or nzmean); no normalization is run here.
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
from ytab.essentiality.control_conditioned_gumbel import NORMALIZATIONS,Settings,run

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--normalized",required=True,help="YTAB normalize output directory, or one method folder inside it (e.g. .../ttr)")
    parser.add_argument("--normalization",choices=NORMALIZATIONS,help="Method folder to use when --normalized is a run directory (default ttr); must match when it is a method folder")
    for name in ("gff","reference","output"):
        parser.add_argument("--"+name,required=True)
    parser.add_argument("--control",action="append",required=True,help="Matched control observation ID; repeat for each control (>=2); combined into the null universe")
    parser.add_argument("--treatment",action="append",required=True,help="Treatment observation ID; repeat for each replicate; combined for the calls")
    defaults=Settings()
    parser.add_argument("--inference",choices=["exact","gumbel-mc"],default=defaults.inference,
                        help="exact (default): deterministic exact conditional tail p, Gumbel reported as a diagnostic; gumbel-mc: fitted Gumbel p with Monte Carlo-aware adequacy gates (20000 simulations by default)")
    for name in ("calibration_sims","validation_sims"):
        parser.add_argument("--"+name.replace("_","-"),type=int,default=None,help="default 1000 for exact (diagnostic only), 20000 for gumbel-mc")
    for name in ("min_control_reads","min_sites","seed","nodes"):
        parser.add_argument("--"+name.replace("_","-"),type=int,default=getattr(defaults,name))
    for name in ("alpha_site","fdr","tail_tolerance","ks_confidence","nterm","cterm"):
        parser.add_argument("--"+name.replace("_","-"),type=float,default=getattr(defaults,name))
    parser.add_argument("--eligibility",choices=["control-support","detectable"],default=defaults.eligibility)
    parser.add_argument("--treatment-dispersion",choices=["replicates","control"],default=defaults.treatment_dispersion,
                        help="replicates: fit from treatment replicate disagreement (needs >=2); control: assume control dispersion explicitly")
    parser.add_argument("--dispersion-model",choices=["linear-quadratic","quadratic"],default=defaults.dispersion_model,
                        help="linear-quadratic: Var=mu*(1+linear)+phi*mu^2 per condition; quadratic: linear term fixed at 0")
    parser.add_argument("--include-stop",action="store_true")
    args=parser.parse_args()
    values={k:getattr(args,k) for k in Settings.__dataclass_fields__ if k!="ignore_stop"}
    values["ignore_stop"]=not args.include_stop
    try:
        result=run(args.normalized,args.gff,args.reference,args.output,args.control,args.treatment,Settings(**values),args.normalization)
    except (ValueError,KeyError,OSError) as error:parser.exit(1,f"ERROR: {error}\n")
    diagnostics=result["diagnostics"]
    print(json.dumps(dict(status=result["status"],analysis=result["analysis"],normalization=result["normalization"]["method"],inference=diagnostics["inference"],
        treatments=[t["observation_id"] for t in result["treatments"]],exposures=result["exposures"]["values"],
        background={k:diagnostics["background"][k] for k in ("m","tau","phi","boundary")},dispersion=diagnostics["dispersion"],
        **{k:diagnostics[k] for k in ("observed_depleted_sites","expected_depleted_sites_null","status_counts","family_size","calls")}),indent=2))

if __name__=="__main__":main()
