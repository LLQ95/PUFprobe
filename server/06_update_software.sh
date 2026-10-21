#!/usr/bin/env bash
#
# 06_update_software.sh — update PUFprobe tools and databases to the latest
# releases, then rerun the affected stages.
#
#   1. Pin the package version string to the current release.
#   2. Upgrade Foldseek and MMseqs2 through conda (bioconda).
#   3. Refresh the small databases (SIFTS, CARD) and, when
#      UPDATE_FOLDSEEK_DB=1 (default), the Foldseek PDB structure database.
#   4. Rerun search, context, rank, report, annotate and the summary.
#   5. Verify versions and Pfam display names.
#
# The PyTorch build is left at the version already validated with ESMFold;
# fair-esm is no longer updated and a newer PyTorch can break the structure
# channel. Every step is idempotent and safe to re-run.
#
# Usage:
#   bash server/06_update_software.sh
#   UPDATE_FOLDSEEK_DB=0 bash server/06_update_software.sh   # skip 1.8 GB PDB
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="${PUFPROBE_SERVER_CONFIG:-$SCRIPT_DIR/configs/server_config.yaml}"
set -a
# shellcheck source=configs/server_config.yaml
source "$CONFIG_FILE"
set +a

ENV_PREFIX="/data/liushiwei/anaconda3/envs/${PUFPROBE_ENV}"
PY="$ENV_PREFIX/bin/python"
WD="$RUNS_ROOT/${BATCH_NAME}_${ANNOTATION_SOURCE}"
MARKERS="$DB_ROOT/.markers"
FOLDSEEK_TARGET="10.941cd33"
MMSEQS_TARGET="18.8cc5c"
PUF_VERSION="0.6.0"
UPDATE_FOLDSEEK_DB="${UPDATE_FOLDSEEK_DB:-1}"

log() { echo "[06_update $(date +%H:%M:%S)] $*"; }
die() { echo "[06_update ERROR] $*" >&2; exit 1; }

[ -d "$ENV_PREFIX" ] || die "conda env prefix $ENV_PREFIX not found."
[ -d "$WD" ] || die "working directory $WD not found; run 03_run_pufprobe.sh first."

# Pick conda or mamba for the install.
if command -v mamba >/dev/null 2>&1; then
    SOLVER=(mamba)
elif [ -x "$ENV_PREFIX/bin/mamba" ]; then
    SOLVER=("$ENV_PREFIX/bin/mamba")
else
    SOLVER=(conda)
fi
log "package solver: ${SOLVER[0]}"

echo "================ 1. pin package version ================"
"$PY" - "$REPO_ROOT/src/pufprobe/__init__.py" "$PUF_VERSION" <<'PY'
import pathlib
import re
import sys

path, version = pathlib.Path(sys.argv[1]), sys.argv[2]
text = path.read_text(encoding="utf-8")
updated = re.sub(r'__version__\s*=\s*"[^"]*"', f'__version__ = "{version}"', text)
if updated != text:
    path.write_text(updated, encoding="utf-8", newline="\n")
    print(f"set __version__ to {version}")
else:
    print(f"__version__ already {version}")
PY

echo "================ 2. install network analysis packages ================"
# igraph and leidenalg build the family-level network; pandas and tqdm support
# the data-processing and reporting scripts.
"$ENV_PREFIX/bin/pip" install --upgrade igraph leidenalg pandas tqdm
"$PY" - <<'PY'
import importlib.util
for name in ("igraph", "leidenalg", "pandas", "tqdm"):
    if importlib.util.find_spec(name) is None:
        raise SystemExit(f"missing required package: {name}")
print("network analysis packages available")
PY

echo "================ 3. upgrade Foldseek / MMseqs2 ================"
log "versions before:"
"$ENV_PREFIX/bin/foldseek" version 2>/dev/null | head -1 | sed 's/^/  foldseek: /' || true
"$ENV_PREFIX/bin/mmseqs" version 2>/dev/null | head -1 | sed 's/^/  mmseqs:  /' || true

"${SOLVER[@]}" install -y -n "$PUFPROBE_ENV" \
    -c bioconda -c conda-forge \
    "foldseek=${FOLDSEEK_TARGET}" "mmseqs2=${MMSEQS_TARGET}"

log "versions after:"
"$ENV_PREFIX/bin/foldseek" version 2>/dev/null | head -1 | sed 's/^/  foldseek: /' || true
"$ENV_PREFIX/bin/mmseqs" version 2>/dev/null | head -1 | sed 's/^/  mmseqs:  /' || true

echo "================ 4. refresh databases ================"
# SIFTS and CARD are small; remove them (and their markers) so the download
# script rebuilds them. The Foldseek PDB database is refreshed only when
# requested because the archive is about 1.8 GB.
log "removing SIFTS mapping for refresh"
rm -f "$DBDIR/pdb_chain_pfam.tsv.gz" "$MARKERS/sifts_pfam"
log "removing CARD bundle for refresh"
rm -f "$DB_ROOT/card/protein_fasta_protein_homolog_model.fasta" \
      "$DB_ROOT/card/card_latest.tar.bz2" "$MARKERS/card_vfdb"
if [ "$UPDATE_FOLDSEEK_DB" = "1" ]; then
    log "removing Foldseek PDB database for refresh"
    rm -f "$DBDIR/foldseek_db"* "$MARKERS/foldseek_smoke"
else
    log "keeping existing Foldseek PDB database (UPDATE_FOLDSEEK_DB=0)"
fi
log "running download script (rebuilds only the missing pieces)"
TIER=A SMOKE_FOLDSEEK="$UPDATE_FOLDSEEK_DB" \
    bash "$REPO_ROOT/server/01_download_dbs.sh"

echo "================ 5. rerun affected stages ================"
# New search databases invalidate the downstream results.
log "clearing search/context/rank/report/annotate outputs"
rm -rf "$WD/05_search" "$WD/06_context" "$WD/07_rank" \
       "$WD/08_report" "$WD/09_annotate"

log "rerunning pipeline stages"
TIER=A FORCE_FALLBACK=0 bash "$REPO_ROOT/server/03_run_pufprobe.sh"

log "rerunning summary"
bash "$REPO_ROOT/server/04_summary.sh"

echo "================ 6. verify ================"
"$PY" - "$DBDIR" "$REPO_ROOT" "$ENV_PREFIX" <<'PY'
import csv
import os
import re
import sys
from pathlib import Path

dbdir = Path(sys.argv[1])
repo = Path(sys.argv[2])
env = Path(sys.argv[3])

sys.path.insert(0, str(repo / "src"))
from pufprobe.ontology_names import load_pfam_names  # noqa: E402

result = repo / "results" / "salmonella_bakta" / "high_confidence_predictions.csv"
pdb_re = re.compile(r"^[0-9][A-Za-z0-9]{3}(?:[._-][A-Za-z0-9]+)*$")
names = load_pfam_names(dbdir / "pfam_names.json")

rows = pfam = 0
bad = []
with open(result, newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        rows += 1
        if row["term_type"] != "Pfam":
            continue
        pfam += 1
        label = row["term_name"].strip()
        if pdb_re.match(label):
            bad.append((row["term_id"], label, "bare PDB"))
        record = names.get(row["term_id"])
        if record is None:
            bad.append((row["term_id"], label, "not in index"))
        elif label != (record.get("name") or record.get("id")):
            bad.append((row["term_id"], label, "not official"))

print(f"high-confidence rows: {rows} pfam rows: {pfam}")
print(f"pfam name index families: {len(names)}")
if bad:
    for item in bad:
        print("PROBLEM:", item)
    sys.exit(1)
print("OK: every Pfam term shows its official Pfam name.")
PY

log "Update complete. Tool versions and databases are now current."
