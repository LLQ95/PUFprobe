#!/usr/bin/env bash
#
# 04_summary.sh — collect per-strain/family counts and wet-lab priority tables.
#
# Outputs (results/<batch>_<source>/):
#   per_sample_summary.csv         Per strain: predicted proteins, functional terms, high-confidence terms
#   channel_hits.csv               Number of hits per evidence channel
#   family_overview.csv            Family size, fold status, pLDDT, priority
#   high_confidence_predictions.csv  All high-confidence functional predictions
#   wetlab_priorities_topN.csv     Wet-lab priority for anchorless families
#   summary.md                     Human-readable summary
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="${PUFPROBE_SERVER_CONFIG:-$SCRIPT_DIR/configs/server_config.yaml}"
set -a
# shellcheck source=configs/server_config.yaml
source "$CONFIG_FILE"
set +a

mkdir -p "$LOG_DIR" "$RESULTS_DIR"
LOG="$LOG_DIR/04_summary_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log() { echo "[04_summary $(date +%H:%M:%S)] $*"; }
die() { echo "[04_summary ERROR] $*" >&2; exit 1; }

WD="$RUNS_ROOT/${BATCH_NAME}_${ANNOTATION_SOURCE}"
[ -f "$WD/07_rank/predictions.json" ] \
    || die "No completed run at $WD (run 03_run_pufprobe.sh first)."

OUT="$RESULTS_DIR/${BATCH_NAME}_${ANNOTATION_SOURCE}"
conda run --no-capture-output -n "$PUFPROBE_ENV" \
    python "$SCRIPT_DIR/summarize.py" \
        --workdir "$WD" \
        --out "$OUT" \
        --top "$WETLAB_TOP_N" \
        --high-confidence "$HIGH_CONFIDENCE"

log "Summary tables written to $OUT"
log "Open $OUT/summary.md for the human-readable overview."
