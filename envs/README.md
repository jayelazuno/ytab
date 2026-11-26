# Locked environments

`environment.local.yml` lists what YTAB needs; it is unpinned, so every `setup.sh` run solves it afresh. That
solve is slow, and two people running it a month apart do not get the same environment.

A lock file records the exact build of every package, with its hash, so an install is a download rather than a
solve. It is produced from an environment that is known to work:

```sh
conda list --explicit -p "$(conda info --base)/envs/ytab-local" > envs/ytab-local.<platform>.lock
```

and consumed with:

```sh
conda create --name ytab-local --file envs/ytab-local.osx-64.lock
```

**A lock is per platform.** Present here:

| file | packages | how it was made |
|---|---|---|
| `ytab-local.osx-64.lock` | 348 | `conda list --explicit --md5` of a working environment |
| `ytab-local.osx-arm64.lock` | 345 | solved for that platform with `--subdir` |
| `ytab-local.linux-64.lock` | 359 | solved for that platform with `--subdir` |

`win-64` is not here: the Windows path is the container, which carries its own environment.

conda solves for another platform without being on it, so the other locks were made from this machine:

```sh
conda create --dry-run --json --subdir linux-64 -n probe \
  -c conda-forge -c bioconda --override-channels --file specs.txt
```

where `specs.txt` is the `dependencies:` of `environment.local.yml`, one per line. The JSON's `FETCH`
entries carry each package's URL and md5; `LINK` entries without a `FETCH` are noarch packages already in
the local cache, whose URL and hash are the same on every platform and can be taken from another lock.

**A cross-solved lock has not been installed on its own platform.** It resolves and every package is pinned
and hashed, but `osx-arm64` and `linux-64` here were produced on an Intel Mac. Install one on the platform it
names before relying on it; regenerating it there with `conda list --explicit --md5` afterwards is better
still, because that records what actually worked.

The container image carries its own solved environment and does not use these files. They are for the local
install described in `docs/setup_local.md`, and for cluster installs where a container is not available.
