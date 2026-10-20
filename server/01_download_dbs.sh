#!/usr/bin/env bash
#
# 01_download_dbs.sh — download and validate search databases in two tiers.
#
#   TIER=A (default, small): run this first to smoke-test the pipeline
#     - MMseqs2 Swiss-Prot sequence database   ($DBDIR/swissprot)
#     - Swiss-Prot cross-reference index       ($DBDIR/swissprot_xref.json)
#     - Foldseek PDB structure database        ($DBDIR/foldseek_db, optional)
#     - CARD protein homolog FASTA             (context --card)
#     - VFDB set A protein FASTA               (context --vfdb)
#     - ESM-2 35M weights (small, CPU friendly)
#
#   TIER=B (large, after the smoke run passes):
#     - MMseqs2 UniRef90 (UniRef50 optional)   ($DBDIR/uniref90)
#     - Foldseek AlphaFold/Proteome            ($DBDIR/foldseek_db)
#     - hhblits UniRef30 hhsuite database      ($DBDIR/hhblits_db_*)
#     - ESM-2 650M weights + ESMFold weights
#     - BFD off by default (DOWNLOAD_BFD=1, ~1.8 TB unpacked)
#
# Every step is idempotent (markers + file checks) and resumable.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="${PUFPROBE_SERVER_CONFIG:-$SCRIPT_DIR/configs/server_config.yaml}"
set -a
# shellcheck source=configs/server_config.yaml
source "$CONFIG_FILE"
set +a

mkdir -p "$LOG_DIR" "$DBDIR" "$MMSEQS_TMP" "$DB_ROOT/.markers" "$DB_ROOT/card" "$DB_ROOT/vfdb"
MARKERS="$DB_ROOT/.markers"
LOG="$LOG_DIR/01_download_dbs_${TIER}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log() { echo "[01_dbs:$TIER $(date +%H:%M:%S)] $*"; }
die() { echo "[01_dbs ERROR] $*" >&2; exit 1; }
done_step() { touch "$MARKERS/$1"; }
step_done() { [ -f "$MARKERS/$1" ]; }

export TORCH_HOME
export OMP_NUM_THREADS="$THREADS"
export MMSEQS_FORCE_CPU=1
# UniProt reference FASTA root (override to a mirror if the official host is slow)
UNIPROT_BASE="${UNIPROT_BASE:-https://ftp.uniprot.org/pub/databases/uniprot/current_release}"

# conda env tools (work without activating the shell)
CONDA_RUN=(conda run --no-capture-output -n "$PUFPROBE_ENV")
if ! conda env list | awk '{print $1}' | grep -qx "$PUFPROBE_ENV"; then
    die "conda environment '$PUFPROBE_ENV' not found; run 00_setup_env.sh first."
fi

download() { # url out
    local url="$1" out="$2"
    if [ -s "$out" ]; then log "  exists: $out"; return 0; fi
    log "  downloading $url"
    wget -c --tries "$WGET_TRIES" --timeout=120 --show-progress --progress=dot:giga -O "$out" "$url" \
        || die "download failed: $url
Re-run this script (resumes) or download on a networked machine and scp to $out."
}

# ---------------------------------------------------------------------------
# MMseqs2 databases
# ---------------------------------------------------------------------------
# Tries the built-in `mmseqs databases` first; if the bundled download URL is
# stale (the mmseqs.com tarballs were retired for Swiss-Prot), builds the
# database from the UniProt reference FASTA with createdb + createindex.
mmseqs_db() { # db_id slug prefix fasta_url
    local db_id="$1" slug="$2" prefix="$3" fasta_url="$4"
    local marker="mmseqs_${slug}"
    if step_done "$marker" && [ -f "${prefix}.index" ]; then
        log "MMseqs2 $db_id already present at $prefix"; return 0
    fi
    log "Building MMseqs2 database '$db_id' at $prefix"
    mkdir -p "$(dirname "$prefix")"
    rm -f "$MMSEQS_TMP/${slug}.tar.gz" 2>/dev/null || true
    if "${CONDA_RUN[@]}" mmseqs databases "$db_id" "$prefix" "$MMSEQS_TMP" \
            --threads "$THREADS" 2> "$MMSEQS_TMP/${slug}_db.err"; then
        log "  built-in downloader succeeded"
    else
        log "  built-in downloader unavailable; building from UniProt FASTA: $fasta_url"
        local fasta_gz="$MMSEQS_TMP/${slug}.fasta.gz"
        local fasta="$MMSEQS_TMP/${slug}.fasta"
        download "$fasta_url" "$fasta_gz"
        rm -f "$fasta"
        gunzip -f "$fasta_gz"
        [ -s "$fasta" ] || die "decompressed FASTA missing: $fasta"
        rm -f "${prefix}"* 2>/dev/null || true
        "${CONDA_RUN[@]}" mmseqs createdb "$fasta" "$prefix" --dbtype 1
        "${CONDA_RUN[@]}" mmseqs createindex "$prefix" "$MMSEQS_TMP" --threads "$THREADS"
    fi
    [ -f "${prefix}.index" ] || die "MMseqs2 $db_id missing ${prefix}.index after build"
    "${CONDA_RUN[@]}" mmseqs linfo "$prefix" | head -3
    done_step "$marker"
}

# Swiss-Prot cross-reference index for term normalization (Pfam/EC/GO).
sprot_xref() {
    local marker="sprot_xref"
    local xref="$DBDIR/swissprot_xref.json"
    if step_done "$marker" && [ -s "$xref" ]; then
        log "Swiss-Prot cross-reference index already present at $xref"; return 0
    fi
    local dat_gz="$MMSEQS_TMP/uniprot_sprot.dat.gz"
    download "${UNIPROT_BASE}/knowledgebase/complete/uniprot_sprot.dat.gz" "$dat_gz"
    log "Building Swiss-Prot cross-reference index (Pfam/EC/GO) ..."
    "${CONDA_RUN[@]}" python "$SCRIPT_DIR/build_xref.py" "$dat_gz" "$xref"
    [ -s "$xref" ] || die "Swiss-Prot cross-reference index missing: $xref"
    log "  index written: $xref"
    done_step "$marker"
}

# ---------------------------------------------------------------------------
# Foldseek databases
# ---------------------------------------------------------------------------
foldseek_db() { # db_id prefix
    local db_id="$1" prefix="$2" marker="$3"
    if step_done "$marker" && [ -f "${prefix}.index" ]; then
        log "Foldseek database already present at $prefix"; return 0
    fi
    log "Building Foldseek database '$db_id' at $prefix"
    # A smoke prefix at the same path must be removed before a full rebuild.
    rm -f "${prefix}"* 2>/dev/null || true
    "${CONDA_RUN[@]}" foldseek databases "$db_id" "$prefix" "$MMSEQS_TMP" \
        --threads "$THREADS" \
        || die "foldseek databases '$db_id' failed (network to foldseek.com is the common cause).
Download manually: foldseek databases '$db_id' <prefix> <tmp> on a networked machine, then copy."
    [ -f "${prefix}.index" ] || die "Foldseek $db_id missing ${prefix}.index after build"
    done_step "$marker"
}

# ---------------------------------------------------------------------------
# CARD / VFDB (context channel FASTA inputs)
# ---------------------------------------------------------------------------
card_vfdb() {
    if step_done card_vfdb; then log "CARD/VFDB already present"; return 0; fi
    log "Downloading CARD canonical bundle (latest/data) ..."
    local card_bz2="$DB_ROOT/card/card_latest.tar.bz2"
    download "https://card.mcmaster.ca/latest/data" "$card_bz2"
    # The legacy /download/6/...fasta link now returns an HTML page; the FASTA
    # ships inside the canonical tar.bz2 with a flat layout.
    if ! tar -tjf "$card_bz2" >/dev/null 2>&1; then
        log "  CARD bundle incomplete; re-downloading once"
        rm -f "$card_bz2"
        download "https://card.mcmaster.ca/latest/data" "$card_bz2"
    fi
    tar -xjf "$card_bz2" -C "$DB_ROOT/card" ./protein_fasta_protein_homolog_model.fasta
    [ -s "$CARD_FA" ] || die "CARD FASTA missing after extraction; check $card_bz2"
    log "Downloading VFDB set A protein sequences"
    download "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz" "$DB_ROOT/vfdb/VFDB_setA_pro.fas.gz"
    gunzip -f "$DB_ROOT/vfdb/VFDB_setA_pro.fas.gz"
    local n_card n_vf
    n_card=$(grep -c '^>' "$CARD_FA" || echo 0)
    n_vf=$(grep -c '^>' "$VFDB_FA" || echo 0)
    log "  CARD records: $n_card ; VFDB records: $n_vf"
    [ "$n_card" -gt 0 ] || die "CARD FASTA empty or blocked; get it from https://card.mcmaster.ca/download and place at $CARD_FA"
    [ "$n_vf" -gt 0 ] || die "VFDB FASTA empty or blocked; get it from http://www.mgc.ac.cn/VFs/download.htm and place at $VFDB_FA"
    done_step card_vfdb
}

# ---------------------------------------------------------------------------
# ESM / ESMFold weights (preload into $TORCH_HOME)
# ---------------------------------------------------------------------------
preload_esm() { # model_name marker
    local model="$1" marker="$2"
    if step_done "$marker"; then log "Weights for $model already cached"; return 0; fi
    log "Pre-loading $model weights into TORCH_HOME=$TORCH_HOME"
    if "${CONDA_RUN[@]}" python - "$model" <<'PY'
import sys
import esm
name = sys.argv[1]
esm.pretrained.load_model_and_alphabet(name)
print(f"preloaded {name}")
PY
    then
        done_step "$marker"
    else
        log "  WARN: could not preload $model (network to dl.fbaipublicfiles.com). The smoke run does not need it; rerun this script later or let embed fetch it on first use."
    fi
}

preload_esmfold() {
    if step_done esmfold; then log "ESMFold weights already cached"; return 0; fi
    log "Pre-loading ESMFold weights into TORCH_HOME=$TORCH_HOME"
    "${CONDA_RUN[@]}" python - <<'PY'
import esm
esm.pretrained.esmfold_v1()
print("preloaded esmfold_v1")
PY
    done_step esmfold
}

# ---------------------------------------------------------------------------
# hhblits UniRef30
# ---------------------------------------------------------------------------
uniref30() {
    local rel="$UNIREF30_RELEASE"
    local dir="$DB_ROOT/hhblits/UniRef30_$rel"
    local marker="hhblits_uniref30_$rel"
    if step_done "$marker" && [ -f "$DBDIR/hhblits_db_hhm.ffdata" ]; then
        log "hhblits UniRef30 already present"; return 0
    fi
    mkdir -p "$dir"
    local tar="$MMSEQS_TMP/UniRef30_${rel}_hhsuite.tar.gz"
    download "https://wwwuser.gwdg.de/~compbiol/uniclust/${rel}/UniRef30_${rel}_hhsuite.tar.gz" "$tar"
    log "extracting (this is tens of GB) ..."
    tar xzf "$tar" -C "$dir"
    local stem
    stem=$(find "$dir" -name "*_hhm.ffdata" -print -quit | sed 's/_hhm.ffdata$//')
    [ -n "$stem" ] || die "UniRef30 hhsuite extraction did not contain *_hhm.ffdata"
    for suffix in hhm.ffdata hhm.ffindex a3m.ffdata a3m.ffindex cs.ffdata cs.ffindex; do
        [ -e "${stem}_${suffix}" ] && ln -sf "${stem}_${suffix}" "$DBDIR/hhblits_db_${suffix}"
    done
    [ -f "$DBDIR/hhblits_db_hhm.ffdata" ] || die "hhblits_db_hhm.ffdata linking failed"
    done_step "$marker"
}

pdb70_smoke() {
    local marker="hhblits_pdb70"
    if step_done "$marker" && [ -f "$DBDIR/hhblits_db_hhm.ffdata" ]; then return 0; fi
    local dir="$DB_ROOT/hhblits/pdb70"
    mkdir -p "$dir"
    local tar="$MMSEQS_TMP/pdb70_from_mmcif_latest.tar.gz"
    download "https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_latest.tar.gz" "$tar"
    tar xzf "$tar" -C "$dir"
    local stem
    stem=$(find "$dir" -maxdepth 1 -name "pdb70*_hhm.ffdata" -print -quit | sed 's/_hhm.ffdata$//')
    [ -n "$stem" ] || die "PDB70 extraction did not contain *_hhm.ffdata"
    for suffix in hhm.ffdata hhm.ffindex a3m.ffdata a3m.ffindex cs.ffdata cs.ffindex; do
        [ -e "${stem}_${suffix}" ] && ln -sf "${stem}_${suffix}" "$DBDIR/hhblits_db_${suffix}"
    done
    done_step "$marker"
}

bfd() {
    if step_done bfd; then return 0; fi
    local dir="$DB_ROOT/bfd"
    mkdir -p "$dir"
    local tar="$MMSEQS_TMP/bfd.sorted_opt.tar.gz"
    download "https://bfd.mmseqs.com/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz" "$tar"
    log "extracting BFD (~1.8 TB unpacked) ..."
    tar xzf "$tar" -C "$dir"
    done_step bfd
}

# ---------------------------------------------------------------------------
# Tier A
# ---------------------------------------------------------------------------
if [ "$TIER" = "A" ]; then
    log "=== TIER A: smoke databases ==="
    log "Free disk at $DB_ROOT (must be >= ~40 GB for tier A):"
    df -h "$DB_ROOT" | tail -1
    mmseqs_db "UniProtKB/Swiss-Prot" swissprot "$DBDIR/swissprot" \
        "${UNIPROT_BASE}/knowledgebase/complete/uniprot_sprot.fasta.gz"
    sprot_xref
    if [ "$SMOKE_FOLDSEEK" = "1" ]; then
        foldseek_db "PDB" "$DBDIR/foldseek_db" foldseek_smoke
    fi
    card_vfdb
    preload_esm esm2_t12_35M_UR50D esm2_t12
    if [ "$SMOKE_HHBLITS" = "1" ]; then pdb70_smoke; fi

# ---------------------------------------------------------------------------
# Tier B
# ---------------------------------------------------------------------------
elif [ "$TIER" = "B" ]; then
    log "=== TIER B: full databases ==="
    log "Free disk at $DB_ROOT (need ~250 GB without Foldseek full, ~1.5 TB with it):"
    df -h "$DB_ROOT" | tail -1
    mmseqs_db "UniRef90" uniref90 "$DBDIR/uniref90" \
        "${UNIPROT_BASE}/uniref/uniref90/uniref90.fasta.gz"
    if [ "${DOWNLOAD_UNIREF50:-0}" = "1" ]; then
        mmseqs_db "UniRef50" uniref50 "$DBDIR/uniref50" \
            "${UNIPROT_BASE}/uniref/uniref50/uniref50.fasta.gz"
    fi
    foldseek_db "$FULL_FOLDSEEK_DB" "$DBDIR/foldseek_db" foldseek_full
    uniref30
    card_vfdb
    preload_esm esm2_t33_650M_UR50D esm2_t33
    if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
        preload_esmfold
    else
        log "No GPU detected; skipping ESMFold weights preload."
    fi
    if [ "$DOWNLOAD_BFD" = "1" ]; then bfd; fi
else
    die "TIER must be A or B (got '$TIER')"
fi

log "Database layout under $DB_ROOT:"
find "$DBDIR" -maxdepth 1 -type f -o -type l | sort | sed 's/^/  /'
log "Done. Next: bash server/02_annotate.sh   (after setting CONTIGS_DIR)"
