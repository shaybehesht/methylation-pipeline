#!/usr/bin/env bash
# One-shot setup for MANGO (Methylation Analysis for Novel Genomic Outcomes).
#
#   bash setup.sh            # create/update the conda env and install the app
#
# Creates a conda env named "methyl-trio-ui" with the external tools
# (modkit/samtools/htslib) and all Python dependencies, then registers the
# `methyl-trio` launch command. Requires conda or mamba/micromamba.
set -euo pipefail

ENV_NAME="${METHYL_TRIO_ENV:-methyl-trio-ui}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if command -v mamba >/dev/null 2>&1; then
  CONDA=mamba
elif command -v conda >/dev/null 2>&1; then
  CONDA=conda
else
  echo "ERROR: conda (or mamba) is required. Install Miniforge: https://github.com/conda-forge/miniforge"
  exit 1
fi

echo "== Creating/updating env '$ENV_NAME' from environment.yml (using $CONDA) =="
if $CONDA env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  $CONDA env update -n "$ENV_NAME" -f "$here/environment.yml"
else
  $CONDA env create -n "$ENV_NAME" -f "$here/environment.yml"
fi

echo
echo "== Verifying tools in '$ENV_NAME' =="
$CONDA run -n "$ENV_NAME" bash -lc '
  python -c "import streamlit, pysam, pandas, plotly, paramiko; print(\"python deps OK\")"
  for t in modkit samtools tabix; do printf "%s: " "$t"; command -v "$t" || echo "(missing — install via conda)"; done
' || true

echo
echo "=========================================================="
echo "Done. To launch:"
echo "  conda activate $ENV_NAME"
echo "  cd \"$here\""
echo "  methyl-trio        # or: ./run.sh"
echo "=========================================================="
