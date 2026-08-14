#!/usr/bin/env bash
set -euo pipefail
script_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
repo_root="$(CDPATH= cd -- "$script_dir/../.." && pwd)"
host="127.0.0.1"; port="3838"; project_config=""; no_browser=0; skip_check=0
usage(){ echo "Usage: $0 [--project-config FILE] [--host HOST] [--port PORT] [--no-browser] [--skip-env-check]"; }
while (($#)); do case "$1" in --project-config) project_config="$2";shift 2;;--host) host="$2";shift 2;;--port) port="$2";shift 2;;--no-browser) no_browser=1;shift;;--skip-env-check) skip_check=1;shift;;-h|--help) usage;exit 0;;*) echo "Unknown option: $1" >&2;usage;exit 2;;esac;done
command -v Rscript >/dev/null || { echo "Rscript not found. Activate the ytab-local environment." >&2;exit 1; }
[[ -f "$repo_root/app/shiny/app.R" ]] || { echo "Missing app/shiny/app.R" >&2;exit 1; }
if (( ! skip_check )); then python_exec="$(command -v python || command -v python3 || true)"; [[ -n "$python_exec" ]] || { echo "Python not found. Activate ytab-local." >&2;exit 1; }; "$python_exec" "$script_dir/ytab_check_local_env.py"; fi
args=("$repo_root/app/shiny/run_app.R" --host "$host" --port "$port")
if [[ -n "$project_config" ]]; then [[ "$project_config" = /* ]] || project_config="$repo_root/$project_config"; [[ -f "$project_config" ]] || { echo "Project config not found: $project_config" >&2;exit 1; };args+=(--project-config "$project_config");fi
if command -v lsof >/dev/null && lsof -iTCP:"$port" -sTCP:LISTEN -t >/dev/null 2>&1; then
  original="$port"; while lsof -iTCP:"$port" -sTCP:LISTEN -t >/dev/null 2>&1; do port=$((port+1));done
  args=("$repo_root/app/shiny/run_app.R" --host "$host" --port "$port"); [[ -n "$project_config" ]] && args+=(--project-config "$project_config"); echo "Port $original occupied; using $port."
fi
(( no_browser )) && args+=(--no-browser)
echo "YTAB app URL: http://$host:$port"; exec Rscript "${args[@]}"
