#!/usr/bin/env bash
set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
repo_root="$(CDPATH= cd -- "$script_dir/../.." && pwd)"
env_file="$repo_root/environment.local.yml"
env_name="$(awk -F': *' '/^name:/ {print $2; exit}' "$env_file")"
env_name="${env_name:-ytab-local}"
project_config="output/projects/H2O2_screen_v1/config/project.yaml"
launch_app=0
skip_validation=0
force_update=1

usage() {
  cat <<EOF
Usage: $0 [--project-config FILE] [--launch-app] [--skip-validation] [--no-update]

Creates or updates the local YTAB conda environment from environment.local.yml.

Behavior:
  1. Uses mamba if available.
  2. If conda exists but mamba is missing, installs mamba into conda base.
  3. If neither conda nor mamba exists, installs a local micromamba bootstrap
     under .tools/micromamba and creates the environment with it.
  4. Validates Python and R packages needed by the Shiny app.

Options:
  --project-config FILE   Project config used with --launch-app.
  --launch-app            Launch the Shiny app after setup.
  --skip-validation       Skip package validation.
  --no-update             Do not update an existing environment.
  -h, --help              Show this help.
EOF
}

while (($#)); do
  case "$1" in
    --project-config) project_config="$2"; shift 2 ;;
    --launch-app) launch_app=1; shift ;;
    --skip-validation) skip_validation=1; shift ;;
    --no-update) force_update=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ ! -f "$env_file" ]]; then
  echo "Missing environment file: $env_file" >&2
  exit 1
fi

install_local_micromamba() {
  local platform arch target bootstrap_dir tmp_dir
  platform="$(uname -s)"
  arch="$(uname -m)"
  case "$platform:$arch" in
    Darwin:arm64) target="osx-arm64" ;;
    Darwin:x86_64) target="osx-64" ;;
    Linux:x86_64) target="linux-64" ;;
    Linux:aarch64|Linux:arm64) target="linux-aarch64" ;;
    *) echo "Unsupported platform for automatic micromamba install: $platform $arch" >&2; exit 1 ;;
  esac
  command -v curl >/dev/null || { echo "curl is required to bootstrap micromamba." >&2; exit 1; }
  command -v tar >/dev/null || { echo "tar is required to bootstrap micromamba." >&2; exit 1; }
  bootstrap_dir="$repo_root/.tools/micromamba"
  mkdir -p "$bootstrap_dir/bin"
  tmp_dir="$(mktemp -d)"
  trap 'rm -rf "$tmp_dir"' RETURN
  echo "mamba/conda not found. Installing local micromamba bootstrap for $target..."
  curl -Ls "https://micro.mamba.pm/api/micromamba/${target}/latest" -o "$tmp_dir/micromamba.tar.bz2"
  tar -xjf "$tmp_dir/micromamba.tar.bz2" -C "$tmp_dir" bin/micromamba
  mv "$tmp_dir/bin/micromamba" "$bootstrap_dir/bin/micromamba"
  chmod +x "$bootstrap_dir/bin/micromamba"
  echo "$bootstrap_dir/bin/micromamba"
}

mamba_cmd=""
root_prefix=""

if command -v mamba >/dev/null 2>&1; then
  mamba_cmd="$(command -v mamba)"
elif command -v conda >/dev/null 2>&1; then
  echo "mamba not found. Installing mamba into the existing conda base environment..."
  conda install -n base -c conda-forge mamba -y
  mamba_cmd="$(command -v mamba || true)"
  if [[ -z "$mamba_cmd" ]]; then
    mamba_cmd="$(conda run -n base which mamba)"
  fi
else
  mamba_cmd="$(install_local_micromamba)"
  root_prefix="$repo_root/.micromamba"
  mkdir -p "$root_prefix"
fi

run_mamba() {
  if [[ -n "$root_prefix" ]]; then
    "$mamba_cmd" --root-prefix "$root_prefix" "$@"
  else
    "$mamba_cmd" "$@"
  fi
}

detect_env_prefix() {
  local prefix
  if [[ -n "$root_prefix" ]]; then
    prefix="$root_prefix/envs/$env_name"
    [[ -d "$prefix" ]] && { echo "$prefix"; return 0; }
  fi
  if command -v conda >/dev/null 2>&1; then
    prefix="$(conda env list | awk -v env="$env_name" '$1 == env {print $NF; exit}')"
    [[ -n "$prefix" && -d "$prefix" ]] && { echo "$prefix"; return 0; }
  fi
  prefix="$(run_mamba env list | awk -v env="$env_name" '$1 == env {print $NF; exit}')"
  [[ -n "$prefix" && -d "$prefix" ]] && { echo "$prefix"; return 0; }
  return 1
}

if run_mamba env list | awk '{print $1}' | grep -Fxq "$env_name"; then
  if (( force_update )); then
    echo "Updating existing environment: $env_name"
    run_mamba env update -n "$env_name" -f "$env_file"
  else
    echo "Environment already exists; skipping update: $env_name"
  fi
else
  echo "Creating environment: $env_name"
  run_mamba env create -f "$env_file"
fi

env_prefix="$(detect_env_prefix)" || {
  echo "Could not locate environment prefix for $env_name after setup." >&2
  exit 1
}
env_python="$env_prefix/bin/python"
env_rscript="$env_prefix/bin/Rscript"
[[ -x "$env_python" ]] || { echo "Missing environment Python: $env_python" >&2; exit 1; }
[[ -x "$env_rscript" ]] || { echo "Missing environment Rscript: $env_rscript" >&2; exit 1; }

if (( ! skip_validation )); then
  echo "Validating $env_name..."
  "$env_python" - <<'PY'
import importlib
for name in ["pandas", "yaml", "Bio", "pysam", "gffutils"]:
    importlib.import_module(name)
print("Python packages OK")
PY
  "$env_rscript" -e 'pkgs <- c("shiny","bslib","DT","jsonlite","yaml","processx","ggplot2","dplyr","readr","ggrepel","slider","minpack.lm"); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if(length(missing)) stop(paste("Missing R packages:", paste(missing, collapse=", "))); cat("R packages OK\n")'
fi

cat <<EOF

YTAB environment is ready: $env_name

Activate it with:
  mamba activate $env_name

Then launch the app with:
  ./scripts/local/ytab_launch_app.sh --project-config $project_config
EOF

if [[ -n "$root_prefix" ]]; then
  cat <<EOF

This machine used the local micromamba bootstrap. If your shell does not know
micromamba yet, initialize it with:
  eval "\$($mamba_cmd shell hook -s bash -r $root_prefix)"
  micromamba activate $env_name
EOF
fi

if (( launch_app )); then
  echo
  echo "Launching YTAB..."
  env PATH="$env_prefix/bin:$PATH" "$repo_root/scripts/local/ytab_launch_app.sh" \
    --project-config "$project_config"
fi
