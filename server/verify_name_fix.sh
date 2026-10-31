#!/usr/bin/env bash
#
# verify_name_fix.sh - validate the product-name specificity filters on the
# existing end-to-end run without touching the original annotation files.
#
# The script syncs the latest code from GitHub main, reruns rank so the
# predictions carry GO depth, reruns annotate, and then scans changes.tsv for
# shallow generic terms. The annotate stage applies a depth gate and a
# propagated-ancestor filter, so generic terms such as membrane are excluded.
# It is idempotent: previous rank and annotate outputs are moved aside with a
# timestamp rather than deleted.
#
# Override the defaults with environment variables:
#   PUFPROBE_ROOT, E2E_WORKDIR, DBDIR, E2E_CONFIG, PUFPROBE_ENV
#
set -euo pipefail

REPO_ROOT="${PUFPROBE_ROOT:-/data/liushiwei/PUFprobe}"
WD="${E2E_WORKDIR:-$REPO_ROOT/test_e2e/run}"
DBDIR="${DBDIR:-$REPO_ROOT/dbs/dbdir}"
CONFIG="${E2E_CONFIG:-$REPO_ROOT/config_full.json}"
ENV_NAME="${PUFPROBE_ENV:-pufprobe}"

cd "$REPO_ROOT"
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# 1. Sync the latest code from the GitHub main branch. --------------------
rm -rf /tmp/puf_extract && mkdir -p /tmp/puf_extract
curl -fL -o /tmp/puf.tar.gz \
  https://codeload.github.com/LLQ95/PUFprobe/tar.gz/refs/heads/main
tar -xzf /tmp/puf.tar.gz -C /tmp/puf_extract --strip-components=1
cp -a /tmp/puf_extract/src/pufprobe/. "$REPO_ROOT/src/pufprobe/"
cp -a /tmp/puf_extract/server/. "$REPO_ROOT/server/"
find "$REPO_ROOT/src" "$REPO_ROOT/server" -name __pycache__ -type d \
  -prune -exec rm -rf {} +
python -c "import pufprobe; print('synced pufprobe', pufprobe.__version__)"

# 2. Back up the current outputs and rerun rank and annotate. -------------
stamp="$(date +%Y%m%d_%H%M%S)"
[ -d "$WD/07_rank" ] && mv "$WD/07_rank" "$WD/07_rank_pre091_$stamp"
[ -d "$WD/09_annotate" ] && mv "$WD/09_annotate" \
  "$WD/09_annotate_pre091_$stamp"

pufprobe --workdir "$WD" --config "$CONFIG" rank --dbdir "$DBDIR" --force
pufprobe --workdir "$WD" --config "$CONFIG" annotate --force

# 3. Inspect the edits and confirm no generic term leaked through. --------
echo "=== changes.tsv ==="
column -t -s $'\t' "$WD/09_annotate/changes.tsv" || true

echo "=== shallow generic-term scan (expect no match) ==="
if grep -E 'catalytic activity|membrane|molecular_function|biological_process|cellular_component' \
  "$WD/09_annotate/changes.tsv"; then
  echo "WARNING: a shallow generic GO term is still present in a product name."
else
  echo "OK: no shallow generic GO term found in product names."
fi
