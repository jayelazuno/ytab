#!/usr/bin/env python3
import argparse, copy, sys, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from ytab.pipeline.treated_vs_parent_runner import (load_project_for_treated_vs_parent,load_comparison_design,get_included_comparisons,resolve_treated_vs_parent_inputs,inspect_treated_vs_parent_script,verify_cpm_normalization,build_fitness_cache_context,fitness_cache_reusable)
p=argparse.ArgumentParser();p.add_argument("--project-config",required=True);a=p.parse_args()
config=load_project_for_treated_vs_parent(Path(a.project_config));design=Path(config["output_project_dir"])/"config/comparison_design.csv"
rows=get_included_comparisons(load_comparison_design(design));inspection=verify_cpm_normalization(inspect_treated_vs_parent_script(ROOT/"scripts/ytab_treated_vs_parent_screen.R"))
def context(selected, annotation=None, params=None):
    resolved=resolve_treated_vs_parent_inputs(config,selected,"raw-summary")
    return build_fitness_cache_context(config,"H2O2_vs_parent",selected,resolved,inspection,annotation,params)
all4=context(rows);backgrounds={x["background"] for x in rows};groups=[[x for x in rows if x["background"]==b] for b in sorted(backgrounds)]
variants=[all4]+[context(x) for x in groups]+[context([rows[0],rows[-1]])]
assert len({x["cache_signature"] for x in variants})==len(variants)
assert context(list(reversed(rows)))["cache_signature"]==all4["cache_signature"]
changed=copy.deepcopy(all4);key=next(iter(changed["input_summary_hashes"]));changed["input_summary_hashes"][key]="changed";assert changed!=all4
assert context(rows,params={"pseudocount":"0.6"})["cache_signature"]!=all4["cache_signature"]
assert context(rows,params={"z_threshold_quantile":"0.95"})["cache_signature"]!=all4["cache_signature"]
with tempfile.TemporaryDirectory() as td:
    one=Path(td)/"one.csv";two=Path(td)/"two.csv";one.write_text("feature_id,label\nx,a\n");two.write_text("feature_id,label\nx,b\n")
    c1=context(rows,{"prediction_table":one,"target_tag":"one"});c2=context(rows,{"prediction_table":two,"target_tag":"two"})
    assert c1["cache_signature"]!=all4["cache_signature"] and c1["cache_signature"]!=c2["cache_signature"]
manifest={"status":"success","cache_signature":all4["cache_signature"],"stable_result_table":next(iter(config.get("_config_path","") and []),"")}
manifest["stable_result_table"]=str(next(iter(resolve_treated_vs_parent_inputs(config,rows)["sample_files"].values())))
assert fitness_cache_reusable(manifest,all4) and not fitness_cache_reusable(manifest,all4,True)
assert "MidLC" not in str(all4) and all4["scientific_parameters"]["input_mode"]=="raw-summary"
print("PASS")
