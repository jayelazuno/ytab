#!/bin/bash
# Run on the login node with bash, not qsub; submits the stage jobs itself.
set -euo pipefail
YTAB="${YTAB:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)}"
source "$YTAB/pipeline/hpc/argon_environment.sh"
exec "$PYTHON_BIN" "$YTAB/pipeline/hpc/argon_pipeline.py" submit "$@"
