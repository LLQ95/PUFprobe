#!/usr/bin/env bash
#
# run_cross_species.sh — one-stop driver for the cross-species expansion.
#
# Runs the stages in order. Every stage is idempotent and resumable, so the
# driver can be re-launched after an interruption; completed stages are skipped.
#
#   07  fetch reference proteomes and their hypothetical proteins
#   08  cluster families (MMseqs2 linclust) and embed representatives (ESM-2)
#   09  fetch AlphaFold models, Foldseek search + SIFTS, optional ProstT5
#   10  build the family kNN network, Leiden communities, render figures
#   11  Swiss-Prot hold-out benchmark of sequence/structure/fusion gain
#
# Run 06_update_software.sh once beforehand to install igraph/leidenalg/pandas
# and refresh Foldseek, MMseqs2 and the structure databases.
#
# Environment:
#   SMOKE=1  small end-to-end smoke run (few proteomes, few reps, small sample)
#   CS_PROSTT5=1  additionally compute ProstT5 embeddings (needs transformers)
#   CS_NO_STRUCTURES=1  skip the structure channel in stage 11
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

log() { echo "[run_cross_species $(date +%H:%M:%S)] $*"; }
die() { echo "[run_cross_species ERROR] $*" >&2; exit 1; }

if [ "${SMOKE:-0}" = "1" ]; then
    export MAX_PROTEOMES="${MAX_PROTEOMES:-6}"
    export MIN_PER_PHYLUM="${MIN_PER_PHYLUM:-1}"
    export CS_MAX_REPS="${CS_MAX_REPS:-200}"
    export CS_BENCHMARK_SAMPLE="${CS_BENCHMARK_SAMPLE:-200}"
    export CS_BENCHMARK_PAGES="${CS_BENCHMARK_PAGES:-8}"
    log "SMOKE mode: reduced sampling and representative limits."
fi

run_stage() {
    local script="$1"
    log "================ stage $script ================"
    bash "$SCRIPT_DIR/$script"
}

run_stage 07_fetch_reference_proteomes.sh
run_stage 08_embed_cluster.sh
run_stage 09_fetch_alphafold.sh
run_stage 10_build_network.sh
run_stage 11_baseline_compare.sh

log "All cross-species stages finished."
log "  network tables/figures:  \$CROSS_SPECIES_ROOT/network/"
log "  benchmark tables/figures: \$CROSS_SPECIES_ROOT/benchmark/"
