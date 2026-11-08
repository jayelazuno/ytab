#!/usr/bin/env bash
# Commit a mapped project's outputs one folder at a time and push in batches. For use on Argon, where a finished
# dataset is thousands of files: a single commit (and a single push) of all of it is slow and fails as a whole, while
# one commit per sample folder pushes in minutes and resumes where it stopped.
#
#   bash scripts/local/ytab_commit_projects.sh Song2026                 # list what would be committed, change nothing
#   bash scripts/local/ytab_commit_projects.sh Song2026 --commit        # commit, push every 25 commits
#   bash scripts/local/ytab_commit_projects.sh Song2026 Chow2024 --commit --push-every 10
#   bash scripts/local/ytab_commit_projects.sh --all --commit           # every project with uncommitted outputs
#
# Re-run it after an interrupted run: it only ever looks at what is still uncommitted, so finished folders are
# skipped. Nothing is deleted and nothing is force-pushed.
set -euo pipefail

cd "$(dirname "$0")/../.."
BRANCH="${YTAB_BRANCH:-main}"
REMOTE="${YTAB_REMOTE:-origin}"
PUSH_EVERY=25
MODE="preview"
PROJECTS=()

while [ $# -gt 0 ]; do
  case "$1" in
    --commit) MODE="commit" ;;
    --all) PROJECTS=(ALL) ;;
    --push-every) PUSH_EVERY="$2"; shift ;;
    --branch) BRANCH="$2"; shift ;;
    --remote) REMOTE="$2"; shift ;;
    -*) echo "Unknown option: $1" >&2; exit 2 ;;
    *) PROJECTS+=("$1") ;;
  esac
  shift
done

if [ "${#PROJECTS[@]}" -eq 0 ]; then
  echo "Usage: bash scripts/local/ytab_commit_projects.sh <PROJECT_ID>... [--commit] [--push-every N]" >&2
  echo "       bash scripts/local/ytab_commit_projects.sh --all --commit" >&2
  exit 2
fi

if [ "${PROJECTS[0]}" = "ALL" ]; then
  mapfile -t PROJECTS < <(git status --porcelain --untracked-files=all -- output/projects output/exports |
    awk '{print $NF}' | awk -F/ '{print $3}' | sort -u)
  [ "${#PROJECTS[@]}" -eq 0 ] && { echo "No project has uncommitted outputs."; exit 0; }
fi

commits=0
pushed=0

push_now() {
  echo "  pushing $commits commit(s) to $REMOTE/$BRANCH"
  git push "$REMOTE" "$BRANCH"
  pushed=$((pushed + commits)); commits=0
}

# Every folder of a project that still holds uncommitted work, deepest first so a sample's own folder is its own
# commit rather than being swallowed by its parent.
folders_with_changes() {
  local project="$1"
  git status --porcelain --untracked-files=all -- \
      "output/projects/$project" "output/exports/$project" 2>/dev/null |
    awk '{print $NF}' | sed 's/^"//; s/"$//' |
    while read -r path; do [ -n "$path" ] && dirname "$path"; done |
    sort -u
}

for project in "${PROJECTS[@]}"; do
  mapfile -t folders < <(folders_with_changes "$project")
  if [ "${#folders[@]}" -eq 0 ]; then echo "— $project: nothing to commit"; continue; fi
  echo "— $project: ${#folders[@]} folder(s) with uncommitted outputs"
  if [ "$MODE" != "commit" ]; then printf '    %s\n' "${folders[@]}" | head -40
    [ "${#folders[@]}" -gt 40 ] && echo "    … and $(( ${#folders[@]} - 40 )) more"
    continue
  fi
  for folder in "${folders[@]}"; do
    git add -A -- "$folder"
    git diff --cached --quiet && continue          # everything in it was ignored
    git commit -q -m "Public Tn-seq datasets: $folder" \
                  -m "Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
    commits=$((commits + 1))
    printf '    %s  %s\n' "$(git rev-parse --short HEAD)" "$folder"
    [ "$commits" -ge "$PUSH_EVERY" ] && push_now
  done
done

if [ "$MODE" = "commit" ]; then
  [ "$commits" -gt 0 ] && push_now
  echo "Pushed $pushed commit(s) to $REMOTE/$BRANCH."
  left=$(git status --porcelain --untracked-files=all -- output/projects output/exports | head -5)
  [ -n "$left" ] && { echo "Still uncommitted under output/projects or output/exports:"; echo "$left" | sed 's/^/    /'; }
fi
