#!/usr/bin/env bash
set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
repo_root="$(CDPATH= cd -- "$script_dir/../.." && pwd)"
remote="origin"
branch=""

usage() {
  cat <<EOF
Usage: $0 [--remote REMOTE] [--branch BRANCH]

Synchronizes the local checkout before starting YTAB app development.

This script is intentionally conservative:
  - refuses to run if local staged, unstaged, or untracked files exist;
  - fetches from the selected remote;
  - pulls with rebase from the selected branch;
  - does not stage, commit, push, clean, reset, or delete files.

Options:
  --remote REMOTE   Git remote to sync from. Default: origin.
  --branch BRANCH   Branch to pull. Default: current checked-out branch.
  -h, --help        Show this help.
EOF
}

while (($#)); do
  case "$1" in
    --remote) remote="$2"; shift 2 ;;
    --branch) branch="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

cd "$repo_root"

git rev-parse --is-inside-work-tree >/dev/null 2>&1 || {
  echo "Not inside a Git work tree: $repo_root" >&2
  exit 1
}

git remote get-url "$remote" >/dev/null 2>&1 || {
  echo "Git remote not found: $remote" >&2
  echo "Available remotes:" >&2
  git remote -v >&2 || true
  exit 1
}

if [[ -z "$branch" ]]; then
  branch="$(git symbolic-ref --quiet --short HEAD || true)"
  if [[ -z "$branch" ]]; then
    echo "Detached HEAD. Re-run with --branch <branch> after checking out a branch." >&2
    exit 1
  fi
fi

status="$(git status --porcelain=v1)"
if [[ -n "$status" ]]; then
  cat >&2 <<EOF
Local changes exist. Commit, stash, or discard them before syncing.

Current status:
$status

No fetch or pull was run.
EOF
  exit 1
fi

echo "Fetching $remote..."
git fetch --prune "$remote"

echo "Pulling $remote/$branch with rebase..."
git pull --rebase "$remote" "$branch"

echo
echo "Sync complete. Current branch:"
git branch --show-current
echo
echo "Current status:"
git status --short
