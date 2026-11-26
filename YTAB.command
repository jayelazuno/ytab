#!/bin/bash
# Double-click this file to start YTAB. It needs Docker and nothing else: R, Shiny, Python and the aligners all live
# inside the image, and your data never leaves this computer.
#
# macOS marks a downloaded file as unopenable until it is allowed. If double-clicking does nothing, open Terminal in
# this folder and run:  chmod +x YTAB.command
#
# Windows: use YTAB.bat instead.

set -u
IMAGE="ghcr.io/jayelazuno/ytab-app:v0.1"
cd "$(dirname "$0")" || exit 1
REPO="$(pwd -P)"

say() { printf '\n%s\n' "$*"; }
fail() { printf '\nERROR: %s\n\n' "$*"; printf 'Press return to close this window. '; read -r _; exit 1; }

printf 'YTAB — Yeast Transposon Analysis Browser\n%s\n' "----------------------------------------"

[ -f "$REPO/app/shiny/run_app.R" ] || fail "This file must stay inside the YTAB folder; it is looking in $REPO."

# Any container runtime will do. Docker Desktop is the usual one, but it needs a paid subscription in larger
# organizations, so Podman and Colima are tried as well: they are open source and run this same image.
RUNTIME=""
for candidate in ${YTAB_RUNTIME:-} docker podman; do
    [ -n "$candidate" ] || continue
    command -v "$candidate" >/dev/null 2>&1 || continue
    if "$candidate" info >/dev/null 2>&1; then RUNTIME="$candidate"; break; fi
    FOUND_BUT_STOPPED="$candidate"
done

if [ -z "$RUNTIME" ] && [ -n "${FOUND_BUT_STOPPED:-}" ]; then
  fail "$FOUND_BUT_STOPPED is installed but not running.
Start it and double-click this again:
  Docker Desktop  — open the app and wait for it to say it is running
  Colima          — run 'colima start'
  Podman          — run 'podman machine start'"
fi

# Offering to install the runtime, rather than only naming it.
#
# Colima is the one offered: it is open source and carries no subscription, where Docker Desktop is free only for
# personal, educational and small-business use and would otherwise be installed by default on a laboratory machine
# that never agreed to its terms. Docker Desktop stays a documented choice, never an automatic one.
#
# The install is offered, never silent. It puts software on someone's computer, and that is their decision.
ask() {
  printf '\n%s ' "$1"
  if [ ! -t 0 ]; then printf '(no terminal to answer on: assuming no)\n'; return 1; fi
  read -r reply
  case "$reply" in [yY]|[yY][eE][sS]) return 0 ;; *) return 1 ;; esac
}

install_colima() {
  if ! command -v brew >/dev/null 2>&1; then
    say "Colima is installed with Homebrew, the package manager for macOS, which is not on this computer."
    printf '  Homebrew installs into /opt/homebrew and will ask for your password.\n'
    printf '  It is widely used and open source: https://brew.sh\n'
    ask "Install Homebrew now? [y/N]" || return 1
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" || return 1
    # A fresh install is not on PATH until a new shell starts, so this one is told where to look.
    for prefix in /opt/homebrew /usr/local; do
      [ -x "$prefix/bin/brew" ] && eval "$("$prefix/bin/brew" shellenv)" && break
    done
    command -v brew >/dev/null 2>&1 || { say "Homebrew did not end up on PATH. Open a new Terminal and run this again."; return 1; }
  fi
  say "Installing Colima and the docker command (a few minutes)…"
  brew install colima docker || return 1
  say "Starting Colima…"
  colima start || return 1
}

if [ -z "$RUNTIME" ]; then
  say "YTAB runs inside a container, and this computer has no container runtime yet."
  printf '  Colima        open source, no subscription  (offered below)\n'
  printf '  Docker Desktop  https://docs.docker.com/get-started/get-docker/\n'
  printf '                  free for personal, educational and small-business use;\n'
  printf '                  larger organizations need a subscription\n'
  printf '  Podman Desktop  https://podman-desktop.io\n'
  if ask "Install Colima now? [y/N]"; then
    install_colima || fail "The Colima install did not finish. Install a runtime from the list above and double-click this again."
    command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1 && RUNTIME="docker"
  fi
fi

[ -n "$RUNTIME" ] || fail "No container runtime is available. Install one of the above, then double-click this again."

say "Using $RUNTIME"

if ! "$RUNTIME" image inspect "$IMAGE" >/dev/null 2>&1; then
  say "First run: downloading YTAB (about 1.4 GB). This happens once and takes a few minutes…"
  "$RUNTIME" pull "$IMAGE" || fail "The download failed. Check your internet connection and try again."
fi

# 3838 is the usual port; step along if something else is already using it.
PORT=3838
while [ "$PORT" -lt 3858 ] && (command -v lsof >/dev/null 2>&1 && lsof -nP -iTCP:"$PORT" -sTCP:LISTEN >/dev/null 2>&1); do
  PORT=$((PORT + 1))
done
NAME="ytab-$PORT"
"$RUNTIME" rm -f "$NAME" >/dev/null 2>&1

say "Starting YTAB on http://localhost:$PORT"
# Your home folder is handed to YTAB read-only at /host, so "Start from FASTQ" can browse to your own sequencing
# files. YTAB only ever reads them; everything it writes goes to this YTAB folder. Set YTAB_FILES to share a
# narrower folder instead, for example: YTAB_FILES=~/sequencing open YTAB.command
"$RUNTIME" run -d --name "$NAME" -p "$PORT":3838 -v "$REPO:/opt/ytab" -v "${YTAB_FILES:-$HOME}:/host:ro" "$IMAGE" >/dev/null \
  || fail "YTAB could not start. Run this for the reason:  $RUNTIME logs $NAME"

printf 'Waiting for it to be ready'
for _ in $(seq 1 60); do
  if curl -fsS --max-time 3 "http://localhost:$PORT/" >/dev/null 2>&1; then
    printf '\n'
    (command -v open >/dev/null 2>&1 && open "http://localhost:$PORT") || \
    (command -v xdg-open >/dev/null 2>&1 && xdg-open "http://localhost:$PORT") || \
    say "Open this in your browser: http://localhost:$PORT"
    break
  fi
  printf '.'; sleep 2
done

cat <<EOF

YTAB is running at http://localhost:$PORT
Your results are saved in this folder, under output/.

To analyse your own reads, use "Start from FASTQ" and click "Browse local folders": your files are there under
"My files", read-only. YTAB only reads them; everything it writes stays in this YTAB folder.

Leave this window open while you work.
To stop YTAB, close this window or press Control-C.
EOF

stop() { printf '\nStopping YTAB…\n'; "$RUNTIME" rm -f "$NAME" >/dev/null 2>&1; exit 0; }
trap stop INT TERM HUP
# Follow the app's log, so the window shows what it is doing and stays open while it runs.
"$RUNTIME" logs -f "$NAME" 2>&1 | sed 's/^/  /'
stop
