#!/usr/bin/env bash
#
# 14_build_sequence_index.sh — export HPDB protein sequences and build an MMseqs2
# index so the database can be searched by an input sequence.
#
# The default mode exports the FASTA from the SQLite file, creates an MMseqs2
# sequence database and writes the prefetch index. Set QUERY_FASTA to run a
# sequence search against the finished index instead.
#
#   bash server/14_build_sequence_index.sh path/to/hpdb.sqlite
#   QUERY_FASTA=query.faa bash server/14_build_sequence_index.sh hpdb.sqlite
#
set -euo pipefail

DB_FILE="${1:?usage: $0 <hpdb.sqlite>}"
[ -f "$DB_FILE" ] || { echo "database not found: $DB_FILE" >&2; exit 1; }

OUT_DIR="${INDEX_DIR:-$(dirname "$DB_FILE")/sequence_index}"
THREADS="${MMSEQS_THREADS:-8}"
MMSEQS_BIN="${MMSEQS_BIN:-mmseqs}"

if ! command -v "$MMSEQS_BIN" >/dev/null 2>&1; then
    echo "mmseqs not found. Install MMseqs2 or set MMSEQS_BIN." >&2
    exit 1
fi

mkdir -p "$OUT_DIR"
FASTA_FILE="$OUT_DIR/hpdb_proteins.fasta"
SEQ_DB="$OUT_DIR/hpdb_proteins"
TMP_DIR="$OUT_DIR/tmp"

if [ -n "${QUERY_FASTA:-}" ]; then
    [ -f "$QUERY_FASTA" ] || { echo "query not found: $QUERY_FASTA" >&2; exit 1; }
    [ -f "$SEQ_DB" ] || { echo "index missing; run without QUERY_FASTA first." >&2; exit 1; }
    echo "Searching $QUERY_FASTA against the HPDB index."
    "$MMSEQS_BIN" easy-search "$QUERY_FASTA" "$SEQ_DB" "$OUT_DIR/sequence_hits.tsv" \
        "$TMP_DIR" --threads "$THREADS" -s 7.5 --format-mode 4
    echo "Hits written to $OUT_DIR/sequence_hits.tsv."
    exit 0
fi

echo "Exporting sequences to $FASTA_FILE."
pufprobe db export-fasta --db "$DB_FILE" --output "$FASTA_FILE" --force

echo "Creating MMseqs2 sequence database."
"$MMSEQS_BIN" createdb "$FASTA_FILE" "$SEQ_DB"

echo "Building prefetch index in $OUT_DIR."
"$MMSEQS_BIN" createindex "$SEQ_DB" "$TMP_DIR" --threads "$THREADS"

echo "Sequence index ready. Search with:"
echo "  QUERY_FASTA=query.faa bash server/14_build_sequence_index.sh $DB_FILE"
