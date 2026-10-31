#!/usr/bin/env bash
#
# verify_pocket_channel.sh - validate the pocket catalytic-template channel on
# the server end to end.
#
# The script syncs the latest code from GitHub main, checks the P2Rank and
# foldseek binaries, locates the M-CSA source data, builds the pocket template
# database when it is missing, and then exercises the online channel. A
# deterministic self-match runs a known template structure (1btl chain A)
# against its own template database and must recover its EC; the existing
# end-to-end proteins are searched as a supplementary batch.
#
# The build is idempotent and resumable; downloaded PDB files, extracted
# chains and the foldseek database are reused. Override the defaults with
# environment variables:
#   PUFPROBE_ROOT, E2E_WORKDIR, DBDIR, E2E_CONFIG, PUFPROBE_ENV,
#   MCSA_DIR, POCKET_WORK_DIR, FORCE_BUILD, THREADS
#
set -euo pipefail

REPO_ROOT="${PUFPROBE_ROOT:-/data/liushiwei/PUFprobe}"
WD="${E2E_WORKDIR:-$REPO_ROOT/test_e2e/run}"
DBDIR="${DBDIR:-$REPO_ROOT/dbs/dbdir}"
CONFIG="${E2E_CONFIG:-$REPO_ROOT/config_full.json}"
ENV_NAME="${PUFPROBE_ENV:-pufprobe}"
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 8)}"

cd "$REPO_ROOT"
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# 1. Sync the latest code from the GitHub main branch. --------------------
rm -rf /tmp/puf_extract && mkdir -p /tmp/puf_extract
curl -fL -o /tmp/puf.tar.gz \
  https://codeload.github.com/LLQ95/PUFprobe/tar.gz/refs/heads/main
tar -xzf /tmp/puf.tar.gz -C /tmp/puf_extract --strip-components=1
cp -a /tmp/puf_extract/src/pufprobe/. "$REPO_ROOT/src/pufprobe/"
cp -a /tmp/puf_extract/server/. "$REPO_ROOT/server/"
find "$REPO_ROOT/src" "$REPO_ROOT/server" -name __pycache__ -type d \
  -prune -exec rm -rf {} +
python -c "import pufprobe; print('synced pufprobe', pufprobe.__version__)"

# 2. Check the external binaries. -----------------------------------------
python - <<'PY'
import sys
from pufprobe import externals

checks = {"p2rank": externals.p2rank(), "foldseek": externals.foldseek()}
missing = []
for name, info in checks.items():
    print(f"{name}: available={info.available} path={info.path or '-'}")
    if not info.available:
        missing.append(name)
if missing:
    sys.exit("missing binary: " + ", ".join(missing) +
             "; see docs/server_setup.md")
PY

# 3. Locate the M-CSA source data. ----------------------------------------
MCSA_DIR="${MCSA_DIR:-}"
if [ -z "$MCSA_DIR" ]; then
    for candidate in \
        "${EVIDENCE_DB_ROOT:-/db/student/metagenome/evidence_dbs}/mcsa" \
        "$REPO_ROOT/dbs/evidence/mcsa"; do
        if [ -s "$candidate/mcsa_residues.json" ]; then
            MCSA_DIR="$candidate"
            break
        fi
    done
fi
if [ -z "$MCSA_DIR" ] || [ ! -s "$MCSA_DIR/mcsa_residues.json" ]; then
    echo "M-CSA JSON not found; run: bash server/12_download_evidence_dbs.sh mcsa"
    exit 1
fi
WORK_DIR="${POCKET_WORK_DIR:-$MCSA_DIR/build}"
echo "using M-CSA source: $MCSA_DIR"
echo "using build work dir: $WORK_DIR"

# 4. Build the template database when it is missing. ----------------------
if [ ! -s "$DBDIR/pocket_templates.json" ] || [ "${FORCE_BUILD:-0}" = "1" ]; then
    python server/build_pocket_templates.py \
        "$MCSA_DIR" "$WORK_DIR" "$DBDIR" --threads "$THREADS"
else
    echo "pocket templates present; skipping build (set FORCE_BUILD=1 to rebuild)"
fi

# 5. Verify the installed template artifacts. -----------------------------
python - <<PY
import json
from pathlib import Path

dbdir = Path("$DBDIR")
templates = json.load(open(dbdir / "pocket_templates.json"))
assert templates, "empty template index"
targets = {record["target"] for record in templates}
assert "1btl_A" in targets, "known template 1btl_A missing"
assert all(record["ec"] for record in templates), "a template lacks an EC"
assert Path(str(dbdir / "pocket_db") + ".index").exists(), "foldseek index missing"
print(f"OK: {len(templates)} templates, all with EC and a foldseek index")
PY

# 6. Exercise the online channel. -----------------------------------------
export POCKET_WD="$WD" POCKET_DBDIR="$DBDIR" POCKET_WORK="$WORK_DIR"
export POCKET_THREADS="$THREADS"
python - <<'PY'
import json
from pathlib import Path

from pufprobe.models import Family
from pufprobe.search import pocket

wd = Path(__import__("os").environ["POCKET_WD"])
dbdir = Path(__import__("os").environ["POCKET_DBDIR"])
work = Path(__import__("os").environ["POCKET_WORK"])
threads = int(__import__("os").environ["POCKET_THREADS"])

# 6a. Deterministic self-match: a template structure against its own db.
chain_pdb = work / "chain_pdb" / "1btl_A.pdb"
assert chain_pdb.exists(), (
    "1btl_A chain PDB not found under the build work dir; rerun with FORCE_BUILD=1"
)
self_base = dbdir.parent / "pocket_verify_self"
structures = self_base / "structures"
self_out = self_base / "search"
structures.mkdir(parents=True, exist_ok=True)
self_out.mkdir(parents=True, exist_ok=True)
(structures / "SELF1.pdb").write_text(chain_pdb.read_text(encoding="utf-8"),
                                      encoding="utf-8")
self_fam = Family(
    family_id="SELF1", representative_key="SELF1@g", members=["SELF1@g"],
    representative_seq="M", size=1, lineage_width=1,
    fold_status="folded", plddt=90.0,
)
self_ev, _ = pocket.run_pocket_search(
    [self_fam], structures, dbdir, self_out, min_plddt=70.0, threads=threads,
)
self_targets = [e.target_id for hits in self_ev.values() for e in hits]
print("self-match targets:", self_targets)
assert "EC:3.5.2.6" in self_targets, "self-match did not recover the 1btl EC"
print("OK: 1btl self-match recovered EC:3.5.2.6 through the pocket channel")

# 6b. Supplementary batch: the existing end-to-end folded proteins.
families_path = wd / "03_cluster" / "families.json"
if families_path.exists():
    families = [Family(**d) for d in json.load(open(families_path))]
    structures_dir = wd / "04_fold" / "structures"
    batch_out = dbdir.parent / "pocket_verify_e2e"
    batch_out.mkdir(parents=True, exist_ok=True)
    usable = pocket._select_query_structures(families, structures_dir, 70.0)
    print(f"end-to-end usable folded families: {len(usable)}")
    batch_ev, _ = pocket.run_pocket_search(
        families, structures_dir, dbdir, batch_out,
        min_plddt=70.0, threads=threads,
    )
    n_records = sum(len(v) for v in batch_ev.values())
    print(f"end-to-end representatives with pocket evidence: {len(batch_ev)}")
    print(f"end-to-end pocket evidence records: {n_records}")
    for rep, hits in batch_ev.items():
        for evidence in hits:
            assert evidence.channel.value == "pocket"
            assert evidence.target_id.startswith("EC:")
            assert evidence.details["db"] == "ec"
            assert evidence.details["catalytic_in_pocket"] >= 1
            assert evidence.details["pocket_consistency"] >= 0.5
            print(f"  {rep}: {evidence.target_id} score={evidence.score:.3f}")
    if n_records == 0:
        print("NOTE: no end-to-end protein matched a pocket template; the "
              "self-match above still confirms the channel works.")
else:
    print("end-to-end families not found; skipping the supplementary batch.")
PY

echo "=== pocket channel verification finished ==="
