#!/usr/bin/env bash
#
# 10_build_network.sh — build and plot the cross-species hypothetical-protein
# family network.
#
# Family representatives are linked into a k-nearest-neighbour graph from their
# ESM-2 embeddings (optionally combined with ProstT5), Leiden communities define
# the macro-clusters, and R renders the community network, a core-family
# detail network and a ranked community-size bar chart.
#
# Inputs:  $CROSS_SPECIES_ROOT/clusters/ (from 08) and structure_search/ (09)
# Outputs: $CROSS_SPECIES_ROOT/network/ (tables and PDF/SVG figures)
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
CLUSTER_DIR="$CROSS_SPECIES_ROOT/clusters"
STRUCTURE_DIR="$CROSS_SPECIES_ROOT/structure_search"
NETWORK_DIR="$CROSS_SPECIES_ROOT/network"
CS_CORE_MIN_SIZE="${CS_CORE_MIN_SIZE:-2}"
RSCRIPT="${RSCRIPT:-$(command -v Rscript || true)}"

LOG="$LOG_DIR/10_build_network_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log() { echo "[10_network $(date +%H:%M:%S)] $*"; }
die() { echo "[10_network ERROR] $*" >&2; exit 1; }

ENV_PREFIX="/data/liushiwei/anaconda3/envs/${PUFPROBE_ENV}"
PY="$ENV_PREFIX/bin/python"
[ -x "$PY" ] || die "python not found at $PY; run 00_setup_env.sh first."
[ -f "$CLUSTER_DIR/representative_embeddings.npz" ] || \
    die "missing $CLUSTER_DIR/representative_embeddings.npz; run 08_embed_cluster.sh."

mkdir -p "$NETWORK_DIR"

PROSTT5_FLAG=()
if [ "${CS_PROSTT5:-0}" = "1" ] && [ -f "$CLUSTER_DIR/prostt5_embeddings.npz" ]; then
    PROSTT5_FLAG=(--prostt5)
    log "Combining ESM-2 with ProstT5 embeddings."
fi

STRUCTURE_FLAG=()
if [ -f "$STRUCTURE_DIR/structure_hits.tsv" ]; then
    STRUCTURE_FLAG=(--structure-dir "$STRUCTURE_DIR")
fi

log "Building family network (k=$CS_NETWORK_K, resolution=$CS_LEIDEN_RESOLUTION)"
"$PY" "$SCRIPT_DIR/cross_species_network.py" \
    --cluster-dir "$CLUSTER_DIR" \
    "${STRUCTURE_FLAG[@]}" \
    --output "$NETWORK_DIR" \
    --k "$CS_NETWORK_K" \
    --resolution "$CS_LEIDEN_RESOLUTION" \
    --seed "$CS_LEIDEN_SEED" \
    --core-min-size "$CS_CORE_MIN_SIZE" \
    "${PROSTT5_FLAG[@]}"

if [ -z "$RSCRIPT" ]; then
    log "Rscript not found; tables written but figures skipped."
    log "Install R (>=4.6) with ggplot2/dplyr/readr/ggrepel/patchwork, then run:"
    log "  Rscript server/plot_network.R $NETWORK_DIR"
else
    log "Rendering network figures with R"
    "$RSCRIPT" "$SCRIPT_DIR/plot_network.R" "$NETWORK_DIR" "$NETWORK_DIR"
fi

log "Network tables and figures written to $NETWORK_DIR"
log "Next: bash server/11_baseline_compare.sh"
