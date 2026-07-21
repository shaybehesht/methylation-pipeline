# MANGO as an interactive application inside an AnVIL / Terra Cloud Environment.
#
# This image layers Jupyter + jupyter-server-proxy and the Google Cloud SDK on
# top of the same conda environment the app already uses (modkit, samtools,
# htslib, methylartist, modbamtools, and the Python deps from environment.yml).
# Terra's Leonardo launches it as a Jupyter server; MANGO then appears in the
# Jupyter Launcher and is served at ``.../proxy/mango/`` via the
# ``jupyter_serverproxy_servers`` entry point declared in pyproject.toml.
#
# Build from the repository root so the whole project is in the build context:
#   docker build -f anvil/Dockerfile -t <registry>/mango-anvil:latest .
FROM mambaorg/micromamba:2.3.3

USER root
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml pyproject.toml /app/
WORKDIR /app
# Analysis tools + app (environment.yml runs `pip install -e .`, which also
# registers the jupyter-server-proxy entry point), plus the Jupyter/proxy layer
# and gsutil/gcloud for reading requester-pays workspace buckets.
RUN micromamba install -y -n base -f environment.yml \
    && micromamba install -y -n base -c conda-forge \
       jupyterlab jupyter-server-proxy google-cloud-sdk \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

ENV PATH="/opt/conda/bin:${PATH}" \
    PYTHONPATH=/app \
    METHYL_TRIO_DATA_ROOT=/home/mambauser \
    METHYL_TRIO_REFERENCE_CACHE=/home/mambauser/.cache/methyl-trio/references

EXPOSE 8888
# Leonardo overrides the command in most configurations; this default lets the
# image also be run directly (`docker run -p 8888:8888 ...`) for local testing.
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", \
     "--ServerApp.allow_origin=*", "--ServerApp.token=", "--ServerApp.base_url=/"]
