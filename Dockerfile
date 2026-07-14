FROM mambaorg/micromamba:2.3.3

USER root
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml pyproject.toml /app/
WORKDIR /app
RUN micromamba install -y -n base -f environment.yml \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
USER $MAMBA_USER
# Reference FASTA and matching GENCODE/CpG annotations are downloaded once at
# runtime into the host-mounted cache rather than baked into the image.
ENV PATH="/opt/conda/bin:${PATH}" PYTHONPATH=/app \
    METHYL_TRIO_DATA_ROOT=/data \
    METHYL_TRIO_REFERENCE_CACHE=/reference-cache
EXPOSE 8501
HEALTHCHECK CMD curl -f http://localhost:8501/_stcore/health || exit 1
CMD ["streamlit", "run", "app/streamlit_app.py", "--server.address=0.0.0.0", "--server.port=8501"]
