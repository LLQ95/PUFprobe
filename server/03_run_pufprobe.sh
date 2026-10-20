#!/usr/bin/env bash
#
# 03_run_pufprobe.sh — run the eight PUFprobe stages on annotated genomes.
#
#   prepare -> embed -> fold -> search -> context -> rank -> report -> annotate
#
# The annotation source (bakta or prokka) is selected by ANNOTATION_SOURCE.
# GPU/CPU mode is auto-detected; for the first smoke run keep FORCE_FALLBACK=1
# so the whole pipeline finishes in minutes without any model download.
#
# Every stage is idempotent (per-stage _DONE markers); rerun after a failure
# and completed stages are skipped.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="${PUFPROBE_SERVER_CONFIG:-$SCRIPT_DIR/configs/server_config.yaml}"
set -a
# shellcheck source=configs/server_config.yaml
source "$CONFIG_FILE"
set +a

mkdir -p "$LOG_DIR" "$RUNS_ROOT" "$RESULTS_DIR"
RUN_CFG_DIR="$RUNS_ROOT/.config"
mkdir -p "$RUN_CFG_DIR"
LOG="$LOG_DIR/03_run_${BATCH_NAME}_${ANNOTATION_SOURCE}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log() { echo "[03_run $(date +%H:%M:%S)] $*"; }
die() { echo "[03_run ERROR] $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 0. Preconditions
# ---------------------------------------------------------------------------
conda env list | awk '{print $1}' | grep -qx "$PUFPROBE_ENV" \
    || die "conda environment '$PUFPROBE_ENV' not found; run 00_setup_env.sh first."
CONDA_RUN=(conda run --no-capture-output -n "$PUFPROBE_ENV")

INPUT_DIR="$ANNOT_ROOT/$ANNOTATION_SOURCE"
[ -d "$INPUT_DIR" ] || die "No annotations at $INPUT_DIR. Run 02_annotate.sh first (ANNOTATION_SOURCE=$ANNOTATION_SOURCE)."
find "$INPUT_DIR" -maxdepth 2 \( -name "*.gff" -o -name "*.gff3" -o -name "*.gbk" -o -name "*.gbff" \) | grep -q . \
    || die "No GFF/GBK files found under $INPUT_DIR."

WD="$RUNS_ROOT/${BATCH_NAME}_${ANNOTATION_SOURCE}"
mkdir -p "$WD"
export TORCH_HOME

# ---------------------------------------------------------------------------
# 1. Mode selection (GPU vs CPU vs fallback)
# ---------------------------------------------------------------------------
MODE="fallback"
EMBED_OPTS=()
if [ "$FORCE_FALLBACK" = "1" ]; then
    MODE="fallback (forced; no torch models)"
    GLOBAL=(--cpu)
    EMBED_OPTS=(--fallback-only)
elif command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    MODE="GPU ($(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1))"
    GLOBAL=()
else
    GLOBAL=(--cpu)
    if [ "$USE_CPU_ESM" = "1" ]; then
        MODE="CPU with small ESM-2 model (slow)"
        EMBED_OPTS=(--model esm2_t12_35M_UR50D)
    else
        MODE="CPU fallback features"
        EMBED_OPTS=(--fallback-only)
    fi
fi

# ---------------------------------------------------------------------------
# 2. Runtime JSON config (threads are a CLI flag; thresholds live here)
# ---------------------------------------------------------------------------
RUN_CONFIG="$RUN_CFG_DIR/${BATCH_NAME}_${ANNOTATION_SOURCE}.json"
MIN_CONF="${MIN_CONFIDENCE_ANNOTATE:-$HIGH_CONFIDENCE}"
"${CONDA_RUN[@]}" python - "$RUN_CONFIG" "$ESM_BATCH_SIZE" "$MIN_PLDDT" "$HIGH_CONFIDENCE" "$MIN_CONF" "$FOLD_ON_CPU" <<'PY'
import json
import sys
path, batch, min_plddt, high_conf, min_conf, fold_on_cpu = sys.argv[1:7]
cfg = {
    "embed": {"batch_size": int(batch)},
    "fold": {"min_plddt": float(min_plddt), "allow_cpu_fold": fold_on_cpu == "1"},
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
# 3. Channels that actually have databases present
# ---------------------------------------------------------------------------
CHANNELS=()
if [ -f "$DBDIR/uniref90.index" ] || [ -f "$DBDIR/swissprot.index" ]; then
    CHANNELS+=(mmseqs)
fi
if [ -f "$DBDIR/foldseek_db.index" ]; then
    CHANNELS+=(foldseek)
fi
if [ -f "$DBDIR/hhblits_db_hhm.ffdata" ]; then
    CHANNELS+=(hhblits)
fi

CONTEXT_OPTS=()
[ -f "$CARD_FA" ] && CONTEXT_OPTS+=(--card "$CARD_FA")
[ -f "$VFDB_FA" ] && CONTEXT_OPTS+=(--vfdb "$VFDB_FA")
if [ -n "$CONTIGS_DIR" ] && [ -d "$CONTIGS_DIR" ]; then
    CONTEXT_OPTS+=(--reference-fna "$CONTIGS_DIR")
fi

log "Working directory: $WD"
log "Compute mode: $MODE"
log "Search channels: ${CHANNELS[*]:-none (databases absent -> skipped)}"
log "Context options: ${CONTEXT_OPTS[*]:-none}"

run_stage() {
    local name="$1"; shift
    local stagelog="$LOG_DIR/03_${name}_$(date +%Y%m%d_%H%M%S).log"
    log "=== stage $name START (log: $stagelog) ==="
    local t0=$SECONDS
    set +e
    "${CONDA_RUN[@]}" pufprobe --workdir "$WD" --threads "$THREADS" \
        --config "$RUN_CONFIG" "${GLOBAL[@]}" "$name" "$@" 2>&1 | tee "$stagelog"
    local rc=${PIPESTATUS[0]}
    set -e
    local dt=$((SECONDS - t0))
    if [ "$rc" -ne 0 ]; then
        die "stage $name FAILED (exit $rc) after ${dt}s. Inspect: $stagelog"
    fi
    log "=== stage $name OK (${dt}s) ==="
}

# ---------------------------------------------------------------------------
# 4. Eight stages
# ---------------------------------------------------------------------------
run_stage prepare --input "$INPUT_DIR"
run_stage embed "${EMBED_OPTS[@]}"
# fold skips gracefully when torch/ESMFold is unavailable (all families -> skipped).
run_stage fold
if [ "${#CHANNELS[@]}" -gt 0 ]; then
    run_stage search --channels "$(IFS=,; echo "${CHANNELS[*]}")" --dbdir "$DBDIR"
else
    log "No search databases present; skipping search stage (run 01_download_dbs.sh)."
fi
run_stage context ${CONTEXT_OPTS[@]+"${CONTEXT_OPTS[@]}"}
run_stage rank
run_stage report --formats json,html,md
run_stage annotate --min-confidence "$MIN_CONF"

log "Pipeline finished. Artifacts:"
find "$WD" -maxdepth 2 -name "_DONE" | sort | sed 's/^/  marker: /'
log "Reports:        $WD/08_report/"
log "Backfilled:     $WD/09_annotate/"
log "Next: bash server/04_summary.sh"
