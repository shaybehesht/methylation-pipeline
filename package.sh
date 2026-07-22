#!/usr/bin/env bash
# Build a distributable archive of MANGO to hand to a colleague (e.g. by email
# or drive) when they do not have access to the Git repository.
#
#   bash package.sh                 # writes ./mango.zip from the current branch
#   bash package.sh /tmp/out.zip    # custom output path
#
# Only version-controlled files are included, so runs/, data/, caches, and .git
# are automatically excluded. The colleague unzips it and runs setup.sh.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="${1:-$here/mango.zip}"
ref="$(git -C "$here" rev-parse --abbrev-ref HEAD)"

git -C "$here" archive --format=zip --prefix=mango/ -o "$out" "$ref"

echo "Wrote $out  (from branch: $ref)"
echo
echo "Tell your colleague to:"
echo "  unzip $(basename "$out") && cd mango"
echo "  bash setup.sh"
echo "  conda activate methyl-trio-ui"
echo "  python -m streamlit run app/streamlit_app.py     # or: ./run.sh"
