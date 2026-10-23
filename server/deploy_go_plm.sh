#!/usr/bin/env bash
# Deploy and validate the GO true-path propagation and the pLM nearest-neighbour
# evidence channel on a server checkout.
#
# This script is stepwise and idempotent. Run each subcommand only when ready:
#
#   bash deploy_go_plm.sh patch <pufprobe_patch.tgz>   # apply code archive
#   bash deploy_go_plm.sh env                          # report model/data paths
#   bash deploy_go_plm.sh go                           # GO propagation rerun
#   bash deploy_go_plm.sh plm-smoke                    # tiny reference library
#   bash deploy_go_plm.sh plm-full                     # full library (GPU, hours)
#   bash deploy_go_plm.sh search-plm                   # search only plm_neighbor
#   bash deploy_go_plm.sh search-all                   # hhblits+eggnog+plm + rank
#
# The pLM library is built with the SAME checkpoint as the query embeddings;
# the model name is read from 02_embed/manifest.json, never assumed.
#
# Requires: conda env "pufprobe", an editable install of pufprobe, and (for the
# pLM stages) the Swiss-Prot flat file uniprot_sprot.dat[.gz].
set -euo pipefail

REPO="${REPO:-/data/liushiwei/PUFprobe}"
DBDIR="${DBDIR:-$REPO/dbs/dbdir}"
DBTMP="${DBTMP:-$REPO/dbs/tmp}"
WD="${WD:-$REPO/test_e2e/run}"
ENV_NAME="${ENV_NAME:-pufprobe}"
OBO_DEFAULT="/db/student/metagenome/eggnog_db_v7/go-basic.obo"
CRUN=(conda run --no-capture-output -n "$ENV_NAME")
SMOKE_LIMIT="${SMOKE_LIMIT:-1000}"
SHARD_SIZE="${SHARD_SIZE:-10000}"

log() { echo "[deploy $(date +%H:%M:%S)] $*"; }
die() { echo "[deploy ERROR] $*" >&2; exit 1; }

query_model() {
  "${CRUN[@]}" python -c \
    'import json,sys; print(json.load(open(sys.argv[1])).get("model",""))' \
    "$WD/02_embed/manifest.json"
}

# ---------------------------------------------------------------------------
# 1. Apply the patch archive over the checkout (one backup).
# ---------------------------------------------------------------------------
do_patch() {
  local tgz="${1:-}"
  [ -n "$tgz" ] || die "usage: deploy_go_plm.sh patch <pufprobe_patch.tgz>"
  [ -f "$tgz" ] || die "patch archive not found: $tgz"
  [ -d "$REPO/src/pufprobe" ] || die "package dir missing: $REPO/src/pufprobe"
  if [ ! -d "$REPO/src/pufprobe.pre_plm_bak" ]; then
    log "Backing up current package to src/pufprobe.pre_plm_bak"
    cp -r "$REPO/src/pufprobe" "$REPO/src/pufprobe.pre_plm_bak"
  fi
  log "Extracting $tgz over $REPO"
  tar xzf "$tgz" -C "$REPO"
  for f in pufprobe/ontology_graph.py pufprobe/search/plm_neighbor.py \
           pufprobe/embed/esm.py pufprobe/rank/pipeline.py; do
    [ -f "$REPO/src/$f" ] || die "expected file missing after patch: src/$f"
  done
  [ -f "$REPO/server/build_plm_reference.py" ] \
    || die "missing server/build_plm_reference.py after patch"
  log "Patch applied (editable install picks up the new code automatically)."
}

# ---------------------------------------------------------------------------
# 2. Report the query embedding model and locate the data files.
# ---------------------------------------------------------------------------
do_env() {
  log "Query workdir: $WD"
  if [ -f "$WD/02_embed/manifest.json" ]; then
    log "Query embedding model: $(query_model)"
  else
    log "WARNING: $WD/02_embed/manifest.json not found; run embed first."
  fi
  local obo=""
  for cand in "$DBDIR/go-basic.obo" "$OBO_DEFAULT"; do
    if [ -e "$cand" ]; then obo="$cand"; log "GO ontology found: $obo"; break; fi
  done
  [ -n "$obo" ] || log "WARNING: go-basic.obo not found."
  local dat=""
  for cand in "$DBTMP/uniprot_sprot.dat.gz" "$DBTMP/uniprot_sprot.dat" \
              "$DBDIR/uniprot_sprot.dat.gz" "$DBDIR/uniprot_sprot.dat"; do
    if [ -f "$cand" ]; then dat="$cand"; log "Swiss-Prot flat file found: $dat"; break; fi
  done
  [ -n "$dat" ] || log "WARNING: uniprot_sprot.dat[.gz] not found; re-run tier A download."
  if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv || true
  fi
}

# ---------------------------------------------------------------------------
# 3. GO true-path propagation: link ontology, rerun rank, verify.
# ---------------------------------------------------------------------------
do_go() {
  local obo=""
  for cand in "$DBDIR/go-basic.obo" "$OBO_DEFAULT"; do
    if [ -e "$cand" ]; then obo="$cand"; break; fi
  done
  [ -n "$obo" ] || die "go-basic.obo not found; supply it under $DBDIR or edit OBO_DEFAULT."
  if [ ! -e "$DBDIR/go-basic.obo" ]; then
    ln -s "$obo" "$DBDIR/go-basic.obo"
    log "Linked $obo -> $DBDIR/go-basic.obo"
  fi
  log "Direct GO evidence entries before propagation:"
  "${CRUN[@]}" python -c \
    'import json,sys; se=json.load(open(sys.argv[1])); print("direct GO evidence entries:", sum(1 for evs in se.values() for e in evs if e["target_id"].startswith("GO:")))' \
    "$WD/05_search/search_evidence.json"
  log "Rerunning rank with propagation"
  "${CRUN[@]}" pufprobe --workdir "$WD" --threads 96 rank --dbdir "$DBDIR" --force
  "${CRUN[@]}" python -c \
    'import json,sys; c=json.load(open(sys.argv[1])); print("go_true_path_propagation:", c.get("go_true_path_propagation")); print("n_go_propagated:", c.get("n_go_propagated"))' \
    "$WD/07_rank/calibration.json"
  log "GO propagation validation done."
}

# ---------------------------------------------------------------------------
# 4. Build the annotated pLM reference library.
# ---------------------------------------------------------------------------
find_dat() {
  local cand
  for cand in "$DBTMP/uniprot_sprot.dat.gz" "$DBTMP/uniprot_sprot.dat" \
              "$DBDIR/uniprot_sprot.dat.gz" "$DBDIR/uniprot_sprot.dat"; do
    if [ -f "$cand" ]; then echo "$cand"; return 0; fi
  done
  return 1
}

build_library() {
  local out="$1" limit="$2"
  local dat
  dat="$(find_dat)" || die "uniprot_sprot.dat[.gz] not found; rerun tier A download."
  local model
  model="$(query_model)"
  [ -n "$model" ] || die "query embedding model unknown; run embed first."
  log "Building reference library at $out with model $model (limit=${limit:-full})"
  local args=(--dat "$dat" --out "$out" --model "$model"
              --shard-size "$SHARD_SIZE" --batch-size 16)
  if [ "$limit" != "0" ]; then args+=(--limit "$limit"); fi
  if [ "${USE_CPU:-0}" = "1" ]; then args+=(--cpu); fi
  PYTHONUNBUFFERED=1 "${CRUN[@]}" python -u "$REPO/server/build_plm_reference.py" "${args[@]}"
  ls -la "$out"
}

do_plm_smoke() { build_library "$DBDIR/plm_ref_smoke" "$SMOKE_LIMIT"; }
do_plm_full()  { build_library "$DBDIR/plm_ref" 0; }

# ---------------------------------------------------------------------------
# 5. Search with plm_neighbor and rerank.
# ---------------------------------------------------------------------------
do_search() {
  local channels="$1"
  # Default to the full library; override with PLM_REF_DIR to smoke-test the
  # small library first (PLM_REF_DIR=$DBDIR/plm_ref_smoke).
  local refdir="${PLM_REF_DIR:-$DBDIR/plm_ref}"
  [ -f "$refdir/plm_ref_manifest.json" ] || \
    die "reference library missing at $refdir; run plm-smoke/plm-full first."
  log "Using reference library: $refdir"
  log "Search channels: $channels"
  local cfg="$WD/plm_${channels//[^a-z_]/_}_config.json"
  "${CRUN[@]}" python -c \
    'import json,os,sys; p,r=sys.argv[1],sys.argv[2]; cfg=json.load(open(p)) if os.path.exists(p) else {}; cfg.setdefault("search",{})["plm_ref_dir"]=r; json.dump(cfg,open(p,"w"),indent=2); print("wrote",p)' \
    "$cfg" "$refdir"
  "${CRUN[@]}" pufprobe --workdir "$WD" --threads 96 --config "$cfg" \
    search --channels "$channels" --dbdir "$DBDIR" --force
  "${CRUN[@]}" pufprobe --workdir "$WD" --threads 96 --config "$cfg" \
    rank --dbdir "$DBDIR" --force
  "${CRUN[@]}" python -c \
    'import json,sys; from collections import Counter; se=json.load(open(sys.argv[1])); c=Counter(e["channel"] for evs in se.values() for e in evs); print("channel counts:", dict(c))' \
    "$WD/05_search/search_evidence.json"
}

do_search_plm() { do_search "plm_neighbor"; }
do_search_all() { do_search "hhblits,eggnog,plm_neighbor"; }

case "${1:-}" in
  patch) shift; do_patch "$@" ;;
  env) do_env ;;
  go) do_go ;;
  plm-smoke) do_plm_smoke ;;
  plm-full) do_plm_full ;;
  search-plm) do_search_plm ;;
  search-all) do_search_all ;;
  *)
    echo "usage: $0 {patch <tgz>|env|go|plm-smoke|plm-full|search-plm|search-all}"
    exit 1
    ;;
esac
