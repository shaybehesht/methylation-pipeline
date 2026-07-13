FROM mambaorg/micromamba:2.3.3

USER root
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml pyproject.toml /app/
WORKDIR /app
RUN micromamba install -y -n base -f environment.yml \
    && micromamba clean --all --yes

ARG GENCODE_RELEASE=49
ARG GENOME_BUILD=GRCh38
USER root
RUN mkdir -p /app/annotations \
    && curl -fsSL "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz" \
       -o /app/annotations/gencode.annotation.gtf.gz \
    && curl -fsSL "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz" \
       -o /app/annotations/cpgIslandExt.txt.gz \
    && gzip -dc /app/annotations/cpgIslandExt.txt.gz \
       | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5}' \
       > /app/annotations/cpg_islands.bed \
    && printf '%s\n' "${GENOME_BUILD}" > /app/annotations/genome-build.txt

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
USER $MAMBA_USER
ENV PATH="/opt/conda/bin:${PATH}" PYTHONPATH=/app
EXPOSE 8501
HEALTHCHECK CMD curl -f http://localhost:8501/_stcore/health || exit 1
CMD ["streamlit", "run", "app/streamlit_app.py", "--server.address=0.0.0.0", "--server.port=8501"]
