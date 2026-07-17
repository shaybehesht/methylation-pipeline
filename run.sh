#!/usr/bin/env bash
# Launch the MANGO (Methylation Analysis for Novel Genomic Outcomes) app from anywhere.
#   ./run.sh                 # defaults
#   ./run.sh --server.port 8502
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$here${PYTHONPATH:+:$PYTHONPATH}"
exec streamlit run "$here/app/streamlit_app.py" "$@"
