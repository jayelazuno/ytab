#!/usr/bin/env bash
set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
repo_root="$(CDPATH= cd -- "$script_dir/../.." && pwd)"
image="ytab-app:local"

usage() {
  echo "Usage: $0 [--image IMAGE]"
}

while (($#)); do
  case "$1" in
    --image) image="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
  esac
done

command -v docker >/dev/null || {
  echo "Docker is required to build the YTAB app image." >&2
  exit 1
}

docker build \
  -f "$repo_root/containers/app.Dockerfile" \
  -t "$image" \
  "$repo_root"

echo "Built $image"
