#!/usr/bin/env bash
#
# 07_fetch_reference_proteomes.sh — build a cross-species hypothetical protein
# dataset from UniProt reference proteomes.
#
# Proteomes are sampled by phylum so the major bacterial lineages are covered,
# then each sampled proteome is downloaded and filtered to its uncharacterized
# proteins. The download is resumable; completed proteomes are skipped on re-run.
#
# The dataset lands under $CROSS_SPECIES_ROOT (default: $DB_ROOT/cross_species).
# Override the size or sampling with environment variables, for example:
#   MAX_PROTEOMES=300 MIN_PER_PHYLUM=3 bash server/07_fetch_reference_proteomes.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="${PUFPROBE_SERVER_CONFIG:-$SCRIPT_DIR/configs/server_config.yaml}"
set -a
# shellcheck source=configs/server_config.yaml
source "$CONFIG_FILE"
set +a

mkdir -p "$LOG_DIR"
CROSS_SPECIES_ROOT="${CROSS_SPECIES_ROOT:-$DB_ROOT/cross_species}"
MAX_PROTEOMES="${MAX_PROTEOMES:-200}"
MIN_PER_PHYLUM="${MIN_PER_PHYLUM:-2}"
SAMPLING_SEED="${SAMPLING_SEED:-42}"

LOG="$LOG_DIR/07_fetch_reference_proteomes_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log() { echo "[07_fetch $(date +%H:%M:%S)] $*"; }
die() { echo "[07_fetch ERROR] $*" >&2; exit 1; }

ENV_PREFIX="/data/liushiwei/anaconda3/envs/${PUFPROBE_ENV}"
PY="$ENV_PREFIX/bin/python"
[ -x "$PY" ] || die "python not found at $PY; run 00_setup_env.sh first."

log "Cross-species root: $CROSS_SPECIES_ROOT"
log "Target proteomes: $MAX_PROTEOMES (minimum $MIN_PER_PHYLUM per phylum)"

"$PY" "$SCRIPT_DIR/fetch_reference_proteomes.py" \
    --output "$CROSS_SPECIES_ROOT" \
    --max-proteomes "$MAX_PROTEOMES" \
    --min-per-phylum "$MIN_PER_PHYLUM" \
    --seed "$SAMPLING_SEED"

log "Cross-species dataset written to $CROSS_SPECIES_ROOT"
log "Next: bash server/08_embed_cluster.sh"
