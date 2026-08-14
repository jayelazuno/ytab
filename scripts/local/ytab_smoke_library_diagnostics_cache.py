#!/usr/bin/env python3
"""Verify LibraryDiagnostics cache signatures without executing diagnostics."""
import argparse, copy, json, sys, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from ytab.pipeline.library_diagnostics_runner import build_cache_context,cache_reusable,get_included_samples,load_project_for_library_diagnostics
def main():
 p=argparse.ArgumentParser();p.add_argument("--project-config",required=True,type=Path);a=p.parse_args();config=load_project_for_library_diagnostics(a.project_config);rows=get_included_samples(config);names=[str(x["sample"]) for x in rows];condition={str(x["sample"]):str(x.get("guessed_condition","")).lower() for x in rows};parents=[x for x in names if condition[x]=="parent"];treated=[x for x in names if condition[x]=="treated"]
 custom=[parents[0],treated[0]];contexts=[build_cache_context(config,x) for x in (names,parents,treated,custom)];assert len({x["cache_signature"] for x in contexts})==4;assert build_cache_context(config,list(reversed(names)))["cache_signature"]==contexts[0]["cache_signature"]
 with tempfile.TemporaryDirectory() as directory:
  manifest={**contexts[0],"status":"success","output_dir":directory};assert cache_reusable(manifest,contexts[0]);assert not cache_reusable(manifest,contexts[1]);changed=copy.deepcopy(contexts[0]);first=next(iter(changed["input_hit_file_hashes"]));changed["input_hit_file_hashes"][first]="changed";changed["cache_signature"]="changed";assert not cache_reusable(manifest,changed);parameter=build_cache_context(config,names,scientific_parameters={"option":"changed"});assert not cache_reusable(manifest,parameter);assert not cache_reusable(manifest,contexts[0],force=True)
 assert contexts[0]["run_id"]=="all_samples" and contexts[1]["run_id"]=="parents" and contexts[2]["run_id"]=="treated" and contexts[3]["run_id"].startswith("custom_");print("PASS")
if __name__=="__main__":main()
