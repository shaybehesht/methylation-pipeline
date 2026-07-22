# MANGO as an interactive application inside an AnVIL / Terra Cloud Environment.
#
# Layers Jupyter + jupyter-server-proxy and the Google Cloud SDK on top of the
# analysis tools (modkit, samtools/htslib, methylartist) and the MANGO package.
# Terra's Leonardo launches it as a Jupyter server; MANGO then appears in the
# Jupyter Launcher and is served at ``.../proxy/mango/`` via the
# ``jupyter_serverproxy_servers`` entry point declared in pyproject.toml.
#
# Note: the source is copied before the editable install so the Streamlit pages
# and console/proxy entry points resolve. modbamtools is not included (it is not
# a conda package here and its PyPI build caps at Python 3.10; it is unused).
#
# Build from the repository root so the whole project is in the build context:
#   docker build -f anvil/Dockerfile -t <registry>/mango-anvil:latest .
FROM mambaorg/micromamba:2.3.3

USER root
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl tar gzip \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Analysis binaries + Jupyter/proxy layer + gsutil/gcloud (requester-pays reads).
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      python=3.11 ont-modkit samtools htslib methylartist pip \
      jupyterlab jupyter-server-proxy google-cloud-sdk \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
WORKDIR /app

# Editable install (source present) registers the mango console/proxy entry
# points and keeps app/pages/ next to streamlit_app.py.
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-cache-dir -e .

ENV PATH="/opt/conda/bin:${PATH}" \
    PYTHONPATH=/app \
    MPLBACKEND=Agg \
    METHYL_TRIO_DATA_ROOT=/home/mambauser \
    METHYL_TRIO_REFERENCE_CACHE=/home/mambauser/.cache/methyl-trio/references

EXPOSE 8888
# Leonardo overrides the command in most configurations; this default lets the
# image also be run directly (`docker run -p 8888:8888 ...`) for local testing.
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", \
     "--ServerApp.allow_origin=*", "--ServerApp.token=", "--ServerApp.base_url=/"]
