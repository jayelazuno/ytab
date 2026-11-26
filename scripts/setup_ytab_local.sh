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

# A lock for this platform, if one exists, is used in preference to the unpinned specification: it installs an
# environment that is known to work, by download rather than by solving, and gives the same one to everybody.
# environment.local.yml names 37 packages across two channels with almost no versions, so solving it is slow
# and the result drifts between machines and over time.
#
# Locks are per platform; envs/README.md says how to make the others. Set YTAB_SOLVE=1 to ignore them and
# solve from the specification, which is what to do when adding a package.
case "$(uname -s)-$(uname -m)" in
  Darwin-arm64) LOCK_PLATFORM="osx-arm64" ;;
  Darwin-*)     LOCK_PLATFORM="osx-64" ;;
  Linux-aarch64) LOCK_PLATFORM="linux-aarch64" ;;
  Linux-*)      LOCK_PLATFORM="linux-64" ;;
  *)            LOCK_PLATFORM="" ;;
esac
LOCK_FILE="${REPO_ROOT}/envs/${ENV_NAME}.${LOCK_PLATFORM}.lock"

if [[ -z "${YTAB_SOLVE:-}" && -n "${LOCK_PLATFORM}" && -f "${LOCK_FILE}" ]]; then
  if "${CONDA_TOOL}" env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
    echo "The ${ENV_NAME} environment already exists; leaving it alone."
    echo "To rebuild it from the lock:  conda env remove -n ${ENV_NAME} && bash scripts/setup_ytab_local.sh"
  else
    echo "Creating ${ENV_NAME} from envs/${ENV_NAME}.${LOCK_PLATFORM}.lock (no solve)..."
    conda create --yes --name "${ENV_NAME}" --file "${LOCK_FILE}"
  fi
  echo
  echo "Setup complete. Run:"
  echo "  mamba activate ${ENV_NAME}"
  exit 0
fi

if [[ -n "${LOCK_PLATFORM}" && ! -f "${LOCK_FILE}" ]]; then
  echo "No lock for ${LOCK_PLATFORM}; solving from environment.local.yml. See envs/README.md to make one."
fi

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
