#!/usr/bin/env bash
#
# 11_baseline_compare.sh — Swiss-Prot hold-out benchmark.
#
# Reviewed bacterial proteins with curated Pfam annotations are treated as
# unknown: the self-hit is removed, then the sequence channel (MMseqs2), the
# structure channel (AlphaFold + Foldseek + SIFTS) and their fusion are scored
# against the hidden Pfam labels. Results are broken down by phylum and by
# sequence identity so the structural gain in the remote-homology regime is
# visible. R renders the gain figures.
#
# Outputs: $CROSS_SPECIES_ROOT/benchmark/
#
# Environment overrides:
#   CS_BENCHMARK_SAMPLE  number of hold-out proteins (default 20000; use ~50 for smoke)
#   CS_NO_STRUCTURES=1   sequence channel only (no AlphaFold/Foldseek)
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
BENCHMARK_DIR="$CROSS_SPECIES_ROOT/benchmark"
SIFTS_FILE="${SIFTS_FILE:-$DBDIR/pdb_chain_pfam.tsv.gz}"
CS_BENCHMARK_SAMPLE="${CS_BENCHMARK_SAMPLE:-20000}"
CS_BENCHMARK_PAGES="${CS_BENCHMARK_PAGES:-120}"
CS_WORKERS="${CS_WORKERS:-8}"
RSCRIPT="${RSCRIPT:-$(command -v Rscript || true)}"

LOG="$LOG_DIR/11_baseline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log() { echo "[11_baseline $(date +%H:%M:%S)] $*"; }
die() { echo "[11_baseline ERROR] $*" >&2; exit 1; }

ENV_PREFIX="/data/liushiwei/anaconda3/envs/${PUFPROBE_ENV}"
PY="$ENV_PREFIX/bin/python"
[ -x "$PY" ] || die "python not found at $PY; run 00_setup_env.sh first."
[ -f "$DBDIR/swissprot.index" ] || \
    die "missing MMseqs2 Swiss-Prot database ($DBDIR/swissprot); run 01_download_dbs.sh."

mkdir -p "$BENCHMARK_DIR"

STRUCTURE_FLAG=()
if [ "${CS_NO_STRUCTURES:-0}" = "1" ]; then
    STRUCTURE_FLAG=(--no-structures)
    log "Running sequence-only benchmark (structure channel disabled)."
else
    [ -f "$DBDIR/foldseek_db.index" ] || \
        die "missing Foldseek database ($DBDIR/foldseek_db); run 01_download_dbs.sh or set CS_NO_STRUCTURES=1."
    [ -s "$SIFTS_FILE" ] || die "missing SIFTS mapping at $SIFTS_FILE."
fi

log "Hold-out benchmark: sample=$CS_BENCHMARK_SAMPLE threads=$THREADS workers=$CS_WORKERS"
"$PY" "$SCRIPT_DIR/baseline_holdout.py" \
    --output "$BENCHMARK_DIR" \
    --dbdir "$DBDIR" \
    --swissprot-db "$DBDIR/swissprot" \
    --foldseek-db "foldseek_db" \
    --sifts "$SIFTS_FILE" \
    --sample "$CS_BENCHMARK_SAMPLE" \
    --max-pages "$CS_BENCHMARK_PAGES" \
    --threads "$THREADS" \
    --workers "$CS_WORKERS" \
    "${STRUCTURE_FLAG[@]}"

if [ -z "$RSCRIPT" ]; then
    log "Rscript not found; tables written but figures skipped."
    log "Install R (>=4.6) with ggplot2/dplyr/readr, then run:"
    log "  Rscript server/plot_benchmark.R $BENCHMARK_DIR"
else
    log "Rendering benchmark figures with R"
    "$RSCRIPT" "$SCRIPT_DIR/plot_benchmark.R" "$BENCHMARK_DIR" "$BENCHMARK_DIR"
fi

log "Benchmark tables and figures written to $BENCHMARK_DIR"
log "See benchmark_summary.tsv and metrics_by_identity.tsv for the structure-channel gain."
