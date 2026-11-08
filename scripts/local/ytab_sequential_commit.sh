#!/usr/bin/env bash
# Commit a change set as a sequence of coherent commits, one subject at a time, then push once.
#
#   bash scripts/local/ytab_sequential_commit.sh            # show what each commit would contain, change nothing
#   bash scripts/local/ytab_sequential_commit.sh --commit   # make the commits, do not push
#   bash scripts/local/ytab_sequential_commit.sh --commit --push
#
# Every path belongs to exactly one group, so no commit carries part of another's subject. Each group is checked
# before it is committed: the app must parse and the smoke tests named for that group must pass. A failing check
# stops the run with the earlier commits intact, so nothing half-checked is pushed.
set -euo pipefail

cd "$(dirname "$0")/../.."
ROOT="$PWD"
R_BIN="${YTAB_R:-Rscript}"
PUSH_URL="https://github.com/jayelazuno/ytab.git"

MODE="preview"; PUSH="no"
for arg in "$@"; do
  case "$arg" in
    --commit) MODE="commit" ;;
    --push) PUSH="yes" ;;
    *) echo "Unknown option: $arg" >&2; exit 2 ;;
  esac
done

# group | subject | paths (space separated) | smoke scripts to run first (repo-root argument)
# NB: not `GROUPS` — bash reserves that name for the caller's Unix group ids and ignores the assignment.
COMMIT_GROUPS=(
"Cite the published method behind each adapted statistic|docs/references.md src/ytab/benchmarks src/ytab/core src/ytab/essentiality src/ytab/fitness src/ytab/normalization src/ytab/qc|"
"Gene Explorer: draw another project's tracks for the same gene|src/ytab/domain_explorer scripts/local/ytab_draw_gene_domain_insertions.py|"
"Compare a project with other projects, tab by tab|app/shiny|ytab_smoke_qc_module.R ytab_smoke_app_ui_release.R ytab_smoke_normalization_ui.R"
"Smoke tests for the comparison and the renamed Normalization tab|scripts/local/ytab_smoke_*.R|ytab_smoke_app_visualizations.R"
)

parse_check() {
  "$R_BIN" -e 'ok <- TRUE
    for (f in c(list.files("app/shiny/R", full.names = TRUE), "app/shiny/app.R"))
      tryCatch(parse(f), error = function(e) { cat("PARSE FAILED:", f, conditionMessage(e), "\n"); ok <<- FALSE })
    if (!ok) quit(status = 1)' >/dev/null
}

run_smokes() {
  local list="$1"
  [ -z "$list" ] && return 0
  for smoke in $list; do
    printf '    %-40s' "$smoke"
    if "$R_BIN" "scripts/local/$smoke" "$ROOT" >/dev/null 2>&1; then echo "pass"; else echo "FAILED"; return 1; fi
  done
}

echo "Repository: $ROOT"
echo "Mode: $MODE${PUSH:+, push: $PUSH}"
echo

for entry in "${COMMIT_GROUPS[@]}"; do
  IFS='|' read -r subject paths smokes <<<"$entry"
  # shellcheck disable=SC2086
  files=$(git status --porcelain -- $paths | awk '{print $NF}')
  if [ -z "$files" ]; then echo "— $subject: nothing to commit"; continue; fi
  echo "— $subject"
  echo "$files" | sed 's/^/    /'
  if [ "$MODE" != "commit" ]; then echo; continue; fi
  # shellcheck disable=SC2086
  git add -A -- $paths
  parse_check
  run_smokes "$smokes"
  git commit -q -m "$subject" -m "Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
  echo "    committed $(git rev-parse --short HEAD)"
  echo
done

if [ "$MODE" = "commit" ]; then
  left=$(git status --porcelain | grep -v '^??' || true)
  if [ -n "$left" ]; then
    echo "Tracked changes outside every group, left uncommitted:"
    echo "$left" | sed 's/^/    /'
  fi
fi

if [ "$PUSH" = "yes" ]; then
  if [ "$MODE" != "commit" ]; then echo "Nothing to push in preview mode." >&2; exit 2; fi
  echo "Pushing to $PUSH_URL"
  git -c credential.helper= -c credential.helper='!gh auth git-credential' push "$PUSH_URL" main
fi
