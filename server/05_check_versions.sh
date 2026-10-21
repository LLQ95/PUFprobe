#!/usr/bin/env bash
#
# 05_check_versions.sh — read-only inventory of software and database versions.
#
# Reports the installed Python packages, external binaries and the build date
# or record count of every local database, so the environment can be compared
# against the latest upstream releases. This script changes nothing and is safe
# to run at any time.
#
# Usage:
#   bash server/05_check_versions.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="${PUFPROBE_SERVER_CONFIG:-$SCRIPT_DIR/configs/server_config.yaml}"
set -a
# shellcheck source=configs/server_config.yaml
source "$CONFIG_FILE"
set +a

PY="/data/liushiwei/anaconda3/envs/${PUFPROBE_ENV}/bin/python"
[ -x "$PY" ] || PY="$(command -v python)"

echo "############ ENVIRONMENT ############"
echo "date (UTC):      $(date -u '+%Y-%m-%d %H:%M:%S')"
echo "hostname:        $(hostname)"
echo "conda env:       ${PUFPROBE_ENV}"
"$PY" --version 2>&1 | sed 's/^/python:          /'
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=name,driver_version --format=csv,noheader \
    | sed 's/^/gpu:             /'
fi

echo
echo "############ PYTHON PACKAGES ############"
"$PY" - <<'PY'
import importlib.metadata as md

packages = [
    "torch", "fair-esm", "numpy", "scikit-learn", "click", "einops",
    "biopython", "igraph", "leidenalg", "pytest", "tqdm", "pandas",
]
for name in packages:
    try:
        print(f"  {name:14s} {md.version(name)}")
    except md.PackageNotFoundError:
        print(f"  {name:14s} (not installed)")
PY

echo
echo "############ EXTERNAL BINARIES ############"
report_bin() { # name command...
    local name="$1"; shift
    if command -v "$1" >/dev/null 2>&1; then
        ("$@" --version 2>&1 || true) | head -1 | sed "s/^/  $name: /"
    else
        echo "  $name: (not found)"
    fi
}

ENV_BIN="/data/liushiwei/anaconda3/envs/${PUFPROBE_ENV}/bin"
[ -x "$ENV_BIN/foldseek" ] && report_bin foldseek "$ENV_BIN/foldseek" || echo "  foldseek: (not found)"
[ -x "$ENV_BIN/mmseqs" ] && report_bin mmseqs "$ENV_BIN/mmseqs" || echo "  mmseqs: (not found)"
[ -x "$ENV_BIN/hhblits" ] && report_bin hhblits "$ENV_BIN/hhblits" || echo "  hhblits: (not found)"
[ -x "$ENV_BIN/pufprobe" ] && report_bin pufprobe "$ENV_BIN/pufprobe" || echo "  pufprobe: (not found)"
if [ -d "$DB_ROOT/p2rank" ]; then
  echo "  p2rank: $(find "$DB_ROOT/p2rank" -maxdepth 2 -name 'prank*' -type f 2>/dev/null | head -1 || echo 'present')"
fi

echo
echo "############ DATABASES (records / mtime) ############"
"$PY" - "$DBDIR" "$DB_ROOT" "$CARD_FA" "$VFDB_FA" <<'PY'
import gzip
import json
import os
import sys
from datetime import datetime, timezone

dbdir, db_root, card_fa, vfdb_fa = sys.argv[1:5]


def mtime(path):
    ts = os.path.getmtime(path)
    return datetime.fromtimestamp(ts, timezone.utc).strftime("%Y-%m-%d")


def count_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as handle:
        return sum(1 for line in handle if line.startswith(">"))


def report(label, path, detail=""):
    if os.path.exists(path):
        print(f"  {label:22s} {mtime(path):11s} {detail}")
    else:
        print(f"  {label:22s} (missing)")


pfam = os.path.join(dbdir, "pfam_names.json")
if os.path.exists(pfam):
    with open(pfam, encoding="utf-8") as handle:
        report("pfam_names", pfam, f"{len(json.load(handle))} families")
else:
    print("  pfam_names            (missing)")

xref = os.path.join(dbdir, "swissprot_xref.json")
if os.path.exists(xref):
    with open(xref, encoding="utf-8") as handle:
        report("swissprot_xref", xref, f"{len(json.load(handle))} accessions")
else:
    print("  swissprot_xref        (missing)")

sifts = os.path.join(dbdir, "pdb_chain_pfam.tsv.gz")
if os.path.exists(sifts):
    with gzip.open(sifts, "rt", errors="replace") as handle:
        rows = sum(1 for _ in handle) - 1
    report("sifts_pdb_chain_pfam", sifts, f"{rows} rows")
else:
    print("  sifts_pdb_chain_pfam  (missing)")

report("foldseek_db.index", os.path.join(dbdir, "foldseek_db.index"))
for name in ("uniref90", "uniref50", "swissprot"):
    report(f"mmseqs_{name}.index", os.path.join(dbdir, f"{name}.index"))
report("hhblits_hhm.ffdata", os.path.join(dbdir, "hhblits_db_hhm.ffdata"))

if os.path.exists(card_fa):
    report("CARD homologs", card_fa, f"{count_fasta(card_fa)} proteins")
else:
    print("  CARD homologs         (missing)")
if os.path.exists(vfdb_fa):
    report("VFDB setA", vfdb_fa, f"{count_fasta(vfdb_fa)} proteins")
else:
    print("  VFDB setA             (missing)")
PY

echo
echo "############ LATEST UPSTREAM (for reference) ############"
echo "  foldseek:  Release 10  (github.com/steineggerlab/foldseek/releases)"
echo "  mmseqs2:   Release 18  (github.com/soedinglab/MMseqs2/releases)"
echo "  p2rank:    2.5.1       (github.com/rdk/p2rank/releases)"
echo "  pfam:      current_release via InterPro (ftp.ebi.ac.uk Pfam)"
echo "  uniref30:  2023_02 (latest available on GWDG)"
echo "Done."
