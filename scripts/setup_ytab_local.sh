#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ENV_FILE="${REPO_ROOT}/environment.local.yml"
ENV_NAME="ytab-local"

if command -v mamba >/dev/null 2>&1; then
  CONDA_TOOL="mamba"
elif command -v conda >/dev/null 2>&1; then
  CONDA_TOOL="conda"
else
  echo "ERROR: Neither mamba nor conda was found on PATH." >&2
  echo "Install Miniforge, restart your shell, and run this script again." >&2
  echo "Miniforge: https://github.com/conda-forge/miniforge" >&2
  exit 1
fi

cd "${REPO_ROOT}"
if "${CONDA_TOOL}" env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
  echo "Updating existing ${ENV_NAME} environment with ${CONDA_TOOL}..."
  "${CONDA_TOOL}" env update --name "${ENV_NAME}" --file "${ENV_FILE}" --prune
else
  echo "Creating ${ENV_NAME} environment with ${CONDA_TOOL}..."
  "${CONDA_TOOL}" env create --file "${ENV_FILE}"
fi

echo
echo "Setup complete. Run:"
echo "  mamba activate ${ENV_NAME}"
if [[ "${CONDA_TOOL}" == "conda" ]]; then
  echo "  (or: conda activate ${ENV_NAME})"
fi
echo "  python scripts/local/ytab_check_local_env.py"
echo "  python scripts/local/ytab_smoke_reference_discovery.py"
echo "  Rscript app/shiny/run_app.R"
