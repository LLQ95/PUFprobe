#!/usr/bin/env bash
#
# 15_run_hpdb_gpu.sh - end-to-end HPDB build on a GPU server.
#
# This script runs the complete HPDB workflow for the stratified panel of
# reference pathogens:
#
#   1. Prepare the inference-only ESMFold environment (idempotent).
#   2. Download the reference pathogen genomes from NCBI (cached manifest,
#      resumable).
#   3. Run the eight PUFprobe stages on GPU with real ESM-2 embeddings and
#      ESMFold structures:
#        prepare -> embed -> fold -> search -> context -> rank -> report
#        -> annotate
#   4. Build the queryable SQLite database.
#   5. Optionally build the MMseqs2 sequence index (BUILD_INDEX=1).
#
# Every stage is guarded by its _DONE marker and the NCBI download is cached,
# so an interrupted run can be restarted with the same command. Completed
# stages are skipped.
#
# Usage:
#   conda run -n pufprobe bash server/15_run_hpdb_gpu.sh
#
# Useful overrides:
#   THREADS=64                 external tool threads
#   HPDB_ROOT=/path/to/hpdb    root for data, workdir and the database
#   SEARCH_CHANNELS=mmseqs,foldseek  override channel autodetection
#   EXTRA_CHANNELS=eggnog      comma-separated channels added to autodetection
#   BUILD_INDEX=1              also build the MMseqs2 sequence index
#   NCBI_API_KEY=...           higher NCBI rate limit (optional)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="${PUFPROBE_SERVER_CONFIG:-$SCRIPT_DIR/configs/server_config.yaml}"
set -a
# shellcheck source=configs/server_config.yaml
source "$CONFIG_FILE"
set +a

# ---------------------------------------------------------------------------
# Directory layout
# ---------------------------------------------------------------------------
HPDB_ROOT="${HPDB_ROOT:-$DB_ROOT/hpdb}"
HPDB_DATA="$HPDB_ROOT/data"
HPDB_WORKDIR="$HPDB_ROOT/workdir"
HPDB_DB="$HPDB_ROOT/hpdb.sqlite"
MANIFEST="$HPDB_DATA/pathogen_manifest.json"
RUN_CFG_DIR="$HPDB_ROOT/.config"
RUN_CONFIG="$RUN_CFG_DIR/hpdb_runtime.json"
HPDB_LOG_DIR="$HPDB_ROOT/logs"

mkdir -p "$HPDB_ROOT" "$HPDB_DATA" "$HPDB_WORKDIR" "$RUN_CFG_DIR" "$HPDB_LOG_DIR"
MASTER_LOG="$HPDB_LOG_DIR/15_run_hpdb_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$MASTER_LOG") 2>&1

log() { echo "[15_hpdb $(date +%H:%M:%S)] $*"; }
die() { echo "[15_hpdb ERROR] $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 0. Preconditions
# ---------------------------------------------------------------------------
conda env list | awk '{print $1}' | grep -qx "$PUFPROBE_ENV" \
    || die "conda environment '$PUFPROBE_ENV' not found; run server/00_setup_env.sh first."
CONDA_RUN=(conda run --no-capture-output -n "$PUFPROBE_ENV")

command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1 \
    || die "No usable GPU detected. This script runs ESMFold on a GPU; check nvidia-smi or the driver."
log "GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -1)"

AVAIL_GB="$(df -Pk "$HPDB_ROOT" | awk 'NR==2 {printf "%.0f", $4/1024/1024}')"
[ "$AVAIL_GB" -ge 20 ] \
    || die "Only ${AVAIL_GB} GB free at $HPDB_ROOT; need at least 20 GB for weights and structures."
log "Free space: ${AVAIL_GB} GB; threads: $THREADS"

export TORCH_HOME

# ---------------------------------------------------------------------------
# 1. ESMFold inference environment (idempotent)
# ---------------------------------------------------------------------------
log "=== Preparing ESMFold inference environment ==="
"${CONDA_RUN[@]}" bash "$SCRIPT_DIR/install_esmfold_inference.sh"

# ---------------------------------------------------------------------------
# 2. Download reference pathogen genomes (resumable)
# ---------------------------------------------------------------------------
log "=== Downloading reference pathogen genomes from NCBI ==="
"${CONDA_RUN[@]}" python "$SCRIPT_DIR/ncbi_pathogen_sampling.py" --output "$HPDB_DATA"
[ -f "$MANIFEST" ] || die "Sampling manifest missing after download step: $MANIFEST"

# ---------------------------------------------------------------------------
# 3. Runtime JSON config
# ---------------------------------------------------------------------------
"${CONDA_RUN[@]}" python - "$RUN_CONFIG" "$ESM_BATCH_SIZE" "$MIN_PLDDT" "$HIGH_CONFIDENCE" "$MIN_CONFIDENCE_ANNOTATE" <<'PY'
import json
import sys
path, batch, min_plddt, high_conf, min_conf = sys.argv[1:6]
cfg = {
    "embed": {"batch_size": int(batch)},
    "fold": {"min_plddt": float(min_plddt), "allow_cpu_fold": False},
    "rank": {
        "high_confidence": float(high_conf),
        "min_confidence_annotate": float(min_conf),
    },
}
with open(path, "w", encoding="utf-8") as handle:
    json.dump(cfg, handle, indent=2)
print(f"runtime config written: {path}")
PY

# ---------------------------------------------------------------------------
# 4. Autodetect search channels from the databases that are present
# ---------------------------------------------------------------------------
# Channel selection. Set SEARCH_CHANNELS to override autodetection, e.g.
# SEARCH_CHANNELS=mmseqs,foldseek to skip the slow serial hhblits channel.
CHANNELS=()
if [ -n "${SEARCH_CHANNELS:-}" ]; then
    for ch in $(echo "$SEARCH_CHANNELS" | tr ',' ' '); do
        CHANNELS+=("$ch")
    done
else
    if [ -f "$DBDIR/uniref90.index" ] || [ -f "$DBDIR/swissprot.index" ]; then
        CHANNELS+=(mmseqs)
    fi
    if [ -f "$DBDIR/foldseek_db.index" ]; then
        CHANNELS+=(foldseek)
    fi
    if [ -f "$DBDIR/hhblits_db_hhm.ffdata" ]; then
        CHANNELS+=(hhblits)
    fi
fi
if [ -n "${EXTRA_CHANNELS:-}" ]; then
    for extra in $(echo "$EXTRA_CHANNELS" | tr ',' ' '); do
        CHANNELS+=("$extra")
    done
fi
CHANNEL_LIST="$(IFS=,; echo "${CHANNELS[*]}")"
log "Search channels: ${CHANNEL_LIST:-none (databases absent -> search skipped)}"

# Context options. Optional CARD/VFDB references are added only when present.
CONTEXT_OPTS=(--reference-fna "$HPDB_DATA/genomes")
[ -f "$CARD_FA" ] && CONTEXT_OPTS+=(--card "$CARD_FA")
[ -f "$VFDB_FA" ] && CONTEXT_OPTS+=(--vfdb "$VFDB_FA")

run_stage() {
    local name="$1"; shift
    local stagelog="$HPDB_LOG_DIR/stage_${name}_$(date +%Y%m%d_%H%M%S).log"
    log "=== stage $name START (log: $stagelog) ==="
    local t0=$SECONDS
    set +e
    "${CONDA_RUN[@]}" pufprobe --workdir "$HPDB_WORKDIR" --threads "$THREADS" \
        --config "$RUN_CONFIG" "$name" "$@" 2>&1 | tee "$stagelog"
    local rc=${PIPESTATUS[0]}
    set -e
    local dt=$((SECONDS - t0))
    [ "$rc" -eq 0 ] || die "stage $name FAILED (exit $rc) after ${dt}s. Inspect: $stagelog"
    log "=== stage $name OK (${dt}s) ==="
}

# ---------------------------------------------------------------------------
# 5. Eight pipeline stages
# ---------------------------------------------------------------------------
run_stage prepare --input "$HPDB_DATA/genomes"
run_stage embed
run_stage fold
if [ "${#CHANNELS[@]}" -gt 0 ]; then
    run_stage search --channels "$CHANNEL_LIST" --dbdir "$DBDIR"
else
    log "No search databases present; skipping search (run server/01_download_dbs.sh)."
fi
run_stage context "${CONTEXT_OPTS[@]}"
run_stage rank --dbdir "$DBDIR"
run_stage report --formats json,html,md
run_stage annotate --min-confidence "$MIN_CONFIDENCE_ANNOTATE"

# ---------------------------------------------------------------------------
# 6. Build the SQLite database
# ---------------------------------------------------------------------------
log "=== Building HPDB SQLite database ==="
"${CONDA_RUN[@]}" pufprobe --workdir "$HPDB_WORKDIR" \
    db build --output "$HPDB_DB" --manifest "$MANIFEST" --force

# ---------------------------------------------------------------------------
# 7. Optional MMseqs2 sequence index
# ---------------------------------------------------------------------------
if [ "${BUILD_INDEX:-0}" = "1" ]; then
    log "=== Building MMseqs2 sequence index ==="
    MMSEQS_THREADS="$THREADS" bash "$SCRIPT_DIR/14_build_sequence_index.sh" "$HPDB_DB"
else
    log "Skipping sequence index (set BUILD_INDEX=1 to build it)."
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log "Workflow finished. Markers:"
find "$HPDB_WORKDIR" -maxdepth 2 -name "_DONE" | sort | sed 's/^/  marker: /'
log "Database:   $HPDB_DB"
log "Reports:    $HPDB_WORKDIR/08_report/"
log "Master log: $MASTER_LOG"
cat <<EOF
Verify the build with:
  conda run -n $PUFPROBE_ENV pufprobe --workdir "$HPDB_WORKDIR" db query \\
    --db "$HPDB_DB" --summary
Browse locally with:
  datasette serve --host 127.0.0.1 --port 8001 --immutable "$HPDB_DB"
EOF
