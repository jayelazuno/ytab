#!/usr/bin/env bash
set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
repo_root="$(CDPATH= cd -- "$script_dir/../.." && pwd)"
image="ytab-app:local"
port="3838"
project_config=""
container_name="ytab-app"
detach=0
restart_policy="unless-stopped"

usage() {
  echo "Usage: $0 [--image IMAGE] [--port PORT] [--project-config FILE] [--name CONTAINER_NAME] [--detach] [--restart-policy POLICY]"
}

while (($#)); do
  case "$1" in
    --image) image="$2"; shift 2 ;;
    --port) port="$2"; shift 2 ;;
    --project-config) project_config="$2"; shift 2 ;;
    --name) container_name="$2"; shift 2 ;;
    --detach) detach=1; shift ;;
    --restart-policy) restart_policy="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
  esac
done

command -v docker >/dev/null || {
  echo "Docker is required to run the YTAB app container." >&2
  exit 1
}

docker_args=(
  run
  --name "$container_name"
  -p "$port:3838"
  -v "$repo_root/output:/opt/ytab/output:ro"
)

if (( detach )); then
  docker_args+=(-d --restart "$restart_policy")
else
  docker_args+=(--rm)
fi

if [[ -n "$project_config" ]]; then
  [[ "$project_config" = /* ]] || project_config="$repo_root/$project_config"
  [[ -f "$project_config" ]] || {
    echo "Project config not found: $project_config" >&2
    exit 1
  }
  case "$project_config" in
    "$repo_root"/*) container_project="/opt/ytab/${project_config#"$repo_root"/}" ;;
    *) echo "Project config must be inside the repository when using the Docker helper." >&2; exit 1 ;;
  esac
  docker_args+=(-e "YTAB_PROJECT_CONFIG=$container_project")
fi

docker_args+=("$image")

echo "YTAB app URL: http://127.0.0.1:$port"
exec docker "${docker_args[@]}"
