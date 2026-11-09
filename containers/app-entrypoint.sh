#!/bin/sh
# The image carries only the environment; the repository is bind-mounted at /opt/ytab, so a missing mount is the
# most likely mistake and should say so rather than fail somewhere inside R.
set -e
if [ ! -f /opt/ytab/app/shiny/run_app.R ]; then
  echo "ERROR: the YTAB repository is not mounted at /opt/ytab." >&2
  echo "Run the image from your clone, for example:" >&2
  echo '  docker run --rm -p 3838:3838 -v "$PWD:/opt/ytab" ghcr.io/jayelazuno/ytab-app:v0.1' >&2
  exit 1
fi
exec micromamba run -n ytab-local Rscript app/shiny/run_app.R --host 0.0.0.0 --port "${YTAB_PORT:-3838}" --no-browser "$@"
