# Minimal image for the MANGO batch WDL workflow (Cromwell / AnVIL).
#
# It contains only what the headless pipeline needs: modkit, samtools/htslib
# (tabix), Python, and the MANGO package (which provides the `mango-run` CLI and
# its matplotlib-based figures). The interactive app's optional visualization
# extras (methylartist, modbamtools) are intentionally excluded so the image is
# small and the conda solve is robust.
#
# Build from the repository root:
#   docker build -f wdl/Dockerfile -t <image> .
FROM mambaorg/micromamba:2.3.3

USER root
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl tar gzip \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# External binaries the pipeline shells out to (modkit, samtools, tabix).
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      python=3.11 ont-modkit samtools htslib pip \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
WORKDIR /app

# Install the package (non-editable) so the `mango-run` console script and the
# core/app modules are importable. MAMBA_DOCKERFILE_ACTIVATE makes pip the
# base-env pip during the build.
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-cache-dir .

ENV MPLBACKEND=Agg
CMD ["mango-run", "--help"]
