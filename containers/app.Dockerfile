FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="YTAB Shiny app"
LABEL org.opencontainers.image.description="Yeast Transposon Analysis Browser Shiny runtime"

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
WORKDIR /opt/ytab

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.local.yml /tmp/ytab-environment.local.yml
RUN micromamba create -y -f /tmp/ytab-environment.local.yml \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER app /opt/ytab/app
COPY --chown=$MAMBA_USER:$MAMBA_USER resources /opt/ytab/resources
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts /opt/ytab/scripts
COPY --chown=$MAMBA_USER:$MAMBA_USER src /opt/ytab/src
COPY --chown=$MAMBA_USER:$MAMBA_USER pipeline /opt/ytab/pipeline
COPY --chown=$MAMBA_USER:$MAMBA_USER README.md environment.local.yml /opt/ytab/

ENV YTAB_DOCKER=1
ENV PYTHONUNBUFFERED=1
EXPOSE 3838

CMD ["micromamba", "run", "-n", "ytab-local", "Rscript", "app/shiny/run_app.R", "--host", "0.0.0.0", "--port", "3838", "--no-browser"]
