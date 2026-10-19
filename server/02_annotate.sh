#!/usr/bin/env bash
#
# 02_annotate.sh — batch-annotate contig genomes with Bakta and/or Prokka.
#
# Each genome is a single (multi-contig) FASTA file in CONTIGS_DIR (set this
# variable in configs/server_config.yaml or via environment). Outputs are
# collected per strain and per annotator:
#   annotations/bakta/<sample>/<sample>.gff3 / .gbff / .faa
#   annotations/prokka/<sample>/<sample>.gff / .gbk / .faa / .ffn
#
# Existing bakta_env and checkm environments are never modified; both tools
# are invoked through `conda run -n ...`. Finished samples are skipped, and
# ANNOT_JOBS genomes are annotated in parallel.
#
#   MAX_SAMPLES=3 bash server/02_annotate.sh   # smoke: first 3 genomes
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="${PUFPROBE_SERVER_CONFIG:-$SCRIPT_DIR/configs/server_config.yaml}"
set -a
# shellcheck source=configs/server_config.yaml
source "$CONFIG_FILE"
set +a

mkdir -p "$LOG_DIR" "$ANNOT_ROOT/bakta" "$ANNOT_ROOT/prokka" "$ANNOT_ROOT/tmp"
LOG="$LOG_DIR/02_annotate_$(date +%Y%m%d_%H%M%S).log"

sample_id_of() {
    local name
    name="$(basename "$1")"
    # strip .gz, then the final extension
    name="${name%.gz}"
    name="${name%.*}"
    echo "$name"
}

# Decompress gz into a stable temp path when needed; otherwise pass through.
plain_input() {
    local file="$1" sample="$2"
    if [[ "$file" == *.gz ]]; then
        local out="$ANNOT_ROOT/tmp/$sample.fna"
        if [ ! -s "$out" ]; then gunzip -c "$file" > "$out"; fi
        echo "$out"
    else
        echo "$file"
    fi
}

run_bakta() {
    local file="$1" sample
    sample="$(sample_id_of "$file")"
    local out="$ANNOT_ROOT/bakta/$sample"
    if [ -f "$out/.bakta_done" ]; then
        echo "  [bakta] $sample already done; skip"; return 0
    fi
    mkdir -p "$out"
    local input; input="$(plain_input "$file" "$sample")"
    echo "  [bakta] $sample -> $out"
    conda run --no-capture-output -n "$BAKTA_ENV" \
        bakta --db "$BAKTA_DB" \
               --output "$out" --prefix "$sample" \
               --threads "$ANNOT_THREADS" --force \
               --keep-contig-headers \
               "$input"
    [ -s "$out/$sample.gff3" ] && [ -s "$out/$sample.faa" ] \
        || { echo "  [bakta] $sample FAILED: missing .gff3/.faa"; return 1; }
    touch "$out/.bakta_done"
}

run_prokka() {
    local file="$1" sample
    sample="$(sample_id_of "$file")"
    local out="$ANNOT_ROOT/prokka/$sample"
    if [ -f "$out/.prokka_done" ]; then
        echo "  [prokka] $sample already done; skip"; return 0
    fi
    mkdir -p "$out"
    local input; input="$(plain_input "$file" "$sample")"
    echo "  [prokka] $sample -> $out"
    conda run --no-capture-output -n "$PROKKA_ENV" \
        prokka --outdir "$out" --prefix "$sample" \
               --locustag "$sample" --kingdom "$PROKKA_KINGDOM" \
               --cpus "$ANNOT_THREADS" --force --addgenes \
               "$input"
    [ -s "$out/$sample.gff" ] && [ -s "$out/$sample.faa" ] \
        || { echo "  [prokka] $sample FAILED: missing .gff/.faa"; return 1; }
    touch "$out/.prokka_done"
}

# ---------------------------------------------------------------------------
# Self-dispatch mode for xargs parallelism: 02_annotate.sh --one <file>
# ---------------------------------------------------------------------------
if [ "${1:-}" = "--one" ]; then
    file="$2"
    exec > >(tee -a "$LOG") 2>&1
    [ "$RUN_BAKTA" = "1" ] && run_bakta "$file"
    [ "$RUN_PROKKA" = "1" ] && run_prokka "$file"
    exit 0
fi

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
exec > >(tee -a "$LOG") 2>&1
log() { echo "[02_annotate $(date +%H:%M:%S)] $*"; }
die() { echo "[02_annotate ERROR] $*" >&2; exit 1; }

[ -n "$CONTIGS_DIR" ] || die "CONTIGS_DIR is empty. Edit server/configs/server_config.yaml or run:
  CONTIGS_DIR=/path/to/your/contigs bash server/02_annotate.sh"
[ -d "$CONTIGS_DIR" ] || die "CONTIGS_DIR does not exist: $CONTIGS_DIR"

# Build the work list (one FASTA per genome).
LIST="$ANNOT_ROOT/tmp/genome_list.$(date +%s).txt"
find "$CONTIGS_DIR" -maxdepth 1 -type f \
    \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" -o -name "*.fas" \
       -o -name "*.fna.gz" -o -name "*.fa.gz" -o -name "*.fasta.gz" -o -name "*.fas.gz" \) \
    | sort > "$LIST"
[ -s "$LIST" ] || die "No FASTA files (*.fna/*.fa/*.fasta/*.fas) found in $CONTIGS_DIR"

if [ "$MAX_SAMPLES" -gt 0 ]; then
    head -n "$MAX_SAMPLES" "$LIST" > "$LIST.smoke"
    mv "$LIST.smoke" "$LIST"
fi

N=$(wc -l < "$LIST")
log "Found $N genome(s) to annotate (bakta=$RUN_BAKTA, prokka=$RUN_PROKKA, jobs=$ANNOT_JOBS, threads/job=$ANNOT_THREADS)"
log "Input list: $LIST"
cat "$LIST" | sed 's/^/  - /'

tr '\n' '\0' < "$LIST" | xargs -0 -P "$ANNOT_JOBS" -I{} \
    bash -c 'set -euo pipefail; "$0" --one "$1"' "$0" "{}"

log "Annotation complete. Outputs:"
echo "  bakta:  $ANNOT_ROOT/bakta"
echo "  prokka: $ANNOT_ROOT/prokka"
log "Next: FORCE_FALLBACK=1 bash server/03_run_pufprobe.sh   (smoke run)"
