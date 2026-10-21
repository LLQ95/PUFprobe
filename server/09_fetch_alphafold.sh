#!/usr/bin/env bash
#
# 09_fetch_alphafold.sh — obtain structures for family representatives and run
# a structural search.
#
# Representative sequences are first matched to AlphaFold DB precomputed models
# (no GPU cost). Accessions absent from AlphaFold DB are optionally folded with
# ESMFold on the local GPU (CS_ESMFOLD_FALLBACK=1). Every available model is
# then searched against the foldseek structural database, with hits mapped to
# Pfam through SIFTS.
#
# Inputs:  $CROSS_SPECIES_ROOT/clusters/representatives.faa   (from 08)
# Outputs: $CS_ALPHAFOLD_DIR/*.pdb, $CROSS_SPECIES_ROOT/structure_search/
#
# Smoke test with CS_MAX_REPS=50. Downloads and the foldseek table are
# resumable.
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
CS_ALPHAFOLD_DIR="${CS_ALPHAFOLD_DIR:-$CROSS_SPECIES_ROOT/structures}"
STRUCTURE_OUT="$CROSS_SPECIES_ROOT/structure_search"
REPS="$CROSS_SPECIES_ROOT/clusters/representatives.faa"
SIFTS_FILE="${SIFTS_FILE:-$DBDIR/pdb_chain_pfam.tsv.gz}"
CS_FOLDSEEK_DB="${CS_FOLDSEEK_DB:-foldseek_db}"
CS_WORKERS="${CS_WORKERS:-8}"
CS_MAX_REPS="${CS_MAX_REPS:-0}"
CS_ESMFOLD_FALLBACK="${CS_ESMFOLD_FALLBACK:-1}"

LOG="$LOG_DIR/09_fetch_alphafold_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log() { echo "[09_struct $(date +%H:%M:%S)] $*"; }
die() { echo "[09_struct ERROR] $*" >&2; exit 1; }

ENV_PREFIX="/data/liushiwei/anaconda3/envs/${PUFPROBE_ENV}"
PY="$ENV_PREFIX/bin/python"
[ -x "$PY" ] || die "python not found at $PY; run 00_setup_env.sh first."
command -v foldseek >/dev/null 2>&1 || die "foldseek not found on PATH; run 06_update_software.sh."

[ -f "$REPS" ] || die "missing $REPS; run 08_embed_cluster.sh first."

mkdir -p "$CS_ALPHAFOLD_DIR" "$STRUCTURE_OUT"

SIFTS_FLAG=()
if [ -f "$SIFTS_FILE" ]; then
    SIFTS_FLAG=(--sifts "$SIFTS_FILE")
    log "SIFTS mapping: $SIFTS_FILE"
else
    log "SIFTS mapping not found at $SIFTS_FILE; hits will be kept as raw PDB targets."
fi

FALLBACK_FLAG=()
if [ "$CS_ESMFOLD_FALLBACK" = "1" ]; then
    FALLBACK_FLAG=(--esmfold-fallback)
    log "ESMFold fallback enabled for representatives absent from AlphaFold DB."
fi

MAX_FLAG=()
if [ "$CS_MAX_REPS" -gt 0 ] 2>/dev/null; then
    MAX_FLAG=(--max-reps "$CS_MAX_REPS")
    log "SMOKE: limited to $CS_MAX_REPS representatives."
fi

log "Fetching structures for representatives from $REPS"
"$PY" "$SCRIPT_DIR/cross_species_structures.py" \
    --reps "$REPS" \
    --struct-dir "$CS_ALPHAFOLD_DIR" \
    --out-dir "$STRUCTURE_OUT" \
    --dbdir "$DBDIR" \
    --foldseek-db "$CS_FOLDSEEK_DB" \
    --workers "$CS_WORKERS" \
    --threads "$THREADS" \
    --min-tm "$CS_FOLDSEEK_TM" \
    --max-evalue "$CS_FOLDSEEK_STRUCT_EVALUE" \
    "${SIFTS_FLAG[@]}" \
    "${FALLBACK_FLAG[@]}" \
    "${MAX_FLAG[@]}"

log "Structures written to $CS_ALPHAFOLD_DIR"
log "Structure search outputs written to $STRUCTURE_OUT"

# Optional ProstT5 structure-aware embeddings (separate from the structure search).
if [ "${CS_PROSTT5:-0}" = "1" ]; then
    log "Computing ProstT5 structure-aware embeddings (model: $CS_PROSTT5_MODEL)"
    "$PY" "$SCRIPT_DIR/cross_species_prostt5.py" \
        --reps "$REPS" \
        --output "$CROSS_SPECIES_ROOT/clusters" \
        --model "$CS_PROSTT5_MODEL" \
        --batch-size "$CS_PROSTT5_BATCH" \
        || log "ProstT5 step skipped (model/runtime unavailable); continuing."
fi

log "Next: bash server/10_build_network.sh"
