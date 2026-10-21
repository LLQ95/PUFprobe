#!/usr/bin/env bash
#
# 08_embed_cluster.sh — cluster cross-species hypothetical proteins into
# families and (optionally) embed family representatives with ESM-2.
#
# The sequence families are built with MMseqs2 linclust, which scales linearly
# and handles hundreds of thousands of proteins. ESM-2 embeddings are computed
# only for the family representatives, so the GPU cost tracks the number of
# families. Both the clustering and each embedding shard are resumable.
#
# Inputs:  $CROSS_SPECIES_ROOT/all_hypothetical_proteins.tsv  (from 07)
# Outputs: $CROSS_SPECIES_ROOT/clusters/
#
# Re-cluster after changing the cutoffs by deleting the linclust output, or
# skip embeddings with CS_EMBED_REPRESENTATIVES=0:
#   CS_EMBED_REPRESENTATIVES=0 bash server/08_embed_cluster.sh
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
PROTEINS_TSV="$CROSS_SPECIES_ROOT/all_hypothetical_proteins.tsv"

LOG="$LOG_DIR/08_embed_cluster_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log() { echo "[08_cluster $(date +%H:%M:%S)] $*"; }
die() { echo "[08_cluster ERROR] $*" >&2; exit 1; }

ENV_PREFIX="/data/liushiwei/anaconda3/envs/${PUFPROBE_ENV}"
PY="$ENV_PREFIX/bin/python"
[ -x "$PY" ] || die "python not found at $PY; run 00_setup_env.sh first."
command -v mmseqs >/dev/null 2>&1 || die "mmseqs not found on PATH; run 06_update_software.sh."

[ -f "$PROTEINS_TSV" ] || die "missing $PROTEINS_TSV; run 07_fetch_reference_proteomes.sh first."

mkdir -p "$CLUSTER_DIR"
log "Clustering hypothetical proteins from $PROTEINS_TSV"
log "linclust: min-seq-id=$CS_MIN_SEQ_ID coverage=$CS_COVERAGE cov-mode=$CS_COV_MODE cluster-mode=$CS_CLUSTER_MODE threads=$THREADS"

EMBED_FLAG=()
if [ "${CS_EMBED_REPRESENTATIVES:-1}" = "0" ]; then
    EMBED_FLAG=(--skip-embedding)
    log "Representative ESM-2 embedding disabled (CS_EMBED_REPRESENTATIVES=0)."
fi

"$PY" "$SCRIPT_DIR/cross_species_cluster.py" \
    --input "$PROTEINS_TSV" \
    --output "$CLUSTER_DIR" \
    --min-seq-id "$CS_MIN_SEQ_ID" \
    --coverage "$CS_COVERAGE" \
    --cov-mode "$CS_COV_MODE" \
    --cluster-mode "$CS_CLUSTER_MODE" \
    --threads "$THREADS" \
    --embed-model "$CS_EMBED_MODEL" \
    --embed-batch "$CS_EMBED_BATCH" \
    ${EMBED_FLAG[@]+"${EMBED_FLAG[@]}"}

log "Family clustering written to $CLUSTER_DIR"
log "  families.tsv              one row per family"
log "  family_members.tsv        accession -> family with taxonomy"
log "  representatives.faa       family representative sequences"
log "  representative_embeddings.npz  (when embedding enabled)"
log "Next: bash server/09_fetch_alphafold.sh"
