# YTAB Shiny app runtime.
#
# The image carries the environment only - R, Shiny, Bioconductor, Python and the aligner tools from
# environment.local.yml - and no YTAB code or data. The repository is bind-mounted at /opt/ytab when the container
# runs, so testers see their own clone (including the demo projects under output/) and pick up code fixes with
# `git pull` instead of a rebuild. Copying the code in as well would only shadow the mount and add ~800 MB.
#
#   docker run --rm -p 3838:3838 -v "$PWD:/opt/ytab" ghcr.io/jayelazuno/ytab-app:v0.1
#
# Built for linux/amd64: several Bioconductor and bioconda packages this environment needs have no linux-aarch64
# build, so Apple Silicon runs it under emulation (slower to start, fine for the interface).
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="YTAB Shiny app"
LABEL org.opencontainers.image.description="Yeast Transposon Analysis Browser Shiny runtime (environment only; mount the repository at /opt/ytab)"
LABEL org.opencontainers.image.source="https://github.com/jayelazuno/ytab"

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      ca-certificates \
      fonts-dejavu-core \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      procps \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Only the environment file is copied, so editing YTAB code never invalidates this layer.
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.local.yml /tmp/ytab-environment.local.yml
RUN micromamba create -y -f /tmp/ytab-environment.local.yml \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/app-entrypoint.sh /usr/local/bin/ytab-app
RUN chmod +x /usr/local/bin/ytab-app

WORKDIR /opt/ytab
ENV YTAB_DOCKER=1
ENV PYTHONUNBUFFERED=1
EXPOSE 3838

ENTRYPOINT ["/usr/local/bin/ytab-app"]
