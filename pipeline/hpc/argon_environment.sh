#!/bin/bash
# Shared personal Argon environment setup; source from a launcher/batch job.
YTAB="${YTAB:-/Users/jayelazuno/workspace/ytab}"
YTAB_ENV_MODE="${YTAB_ENV_MODE:-conda}"
YTAB_CONDA_SH="${YTAB_CONDA_SH:-/old_Users/jayelazuno/miniforge3/etc/profile.d/conda.sh}"
YTAB_CONDA_ENV="${YTAB_CONDA_ENV:-tnseq}"
PYTHON_BIN="${PYTHON_BIN:-python}"
[[ -d "$YTAB" ]] || { echo "ERROR: YTAB checkout missing: $YTAB" >&2; return 2; }
YTAB="$(cd "$YTAB" && pwd -P)"
case "$YTAB_ENV_MODE" in
    conda)
        [[ -s "$YTAB_CONDA_SH" ]] || { echo "ERROR: Conda initialization missing." >&2; return 2; }
        source "$YTAB_CONDA_SH"
        conda activate "$YTAB_CONDA_ENV"
        ;;
    modules)
        if ! type module >/dev/null 2>&1; then
            if [[ -f /etc/profile.d/modules.sh ]]; then source /etc/profile.d/modules.sh
            elif [[ -f /usr/share/lmod/lmod/init/bash ]]; then source /usr/share/lmod/lmod/init/bash
            else echo "ERROR: Environment modules unavailable." >&2; return 2
            fi
        fi
        [[ -n "${YTAB_MODULES:-}" ]] || { echo "ERROR: Set YTAB_MODULES." >&2; return 2; }
        read -r -a requested_modules <<< "$YTAB_MODULES"
        module load "${requested_modules[@]}"
        ;;
    *) echo "ERROR: YTAB_ENV_MODE must be conda or modules." >&2; return 2 ;;
esac
command -v "$PYTHON_BIN" >/dev/null
export YTAB YTAB_ENV_MODE YTAB_CONDA_SH YTAB_CONDA_ENV PYTHON_BIN
