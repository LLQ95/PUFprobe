#!/usr/bin/env bash
#
# 12_download_evidence_dbs.sh — fetch auxiliary databases for the fusion
# evidence channels that are not covered by 01_download_dbs.sh.
#
#   string   STRING v12.5 DATABASE-corroboration channel
#            - COG-level association links (tiny; pair with eggNOG COG IDs)
#            - species list and species tree
#            - per-species protein networks for the target taxid panel
#            - optional global protein.aliases table (STRING_WITH_ALIASES=1)
#   go       GO ontology (go-basic.obo) for true-path propagation; installs
#            goatools. Set WITH_GOA=1 to also fetch the GOA reference-proteome
#            GAF (~23 GB gz) used for BPO model training.
#   mcsa     M-CSA catalytic residues / entries (JSON API) for the pocket
#            catalytic-template channel
#   hhsuite  UniRef30 HHblits database (large; stored on the data volume and
#            linked into $DBDIR as hhblits_db_*, matching the 01 script)
#   small    string small files + go + mcsa (no UniRef30, no aliases)
#   all      every module
#
# Everything is idempotent and resumable: existing files are skipped, curl
# resumes partial transfers, and completed modules are marked under
# $EVIDENCE_DB_ROOT/.markers.
#
# Usage:
#   bash server/12_download_evidence_dbs.sh small      # quick, ~300 MB
#   bash server/12_download_evidence_dbs.sh hhsuite    # big, in background
#   bash server/12_download_evidence_dbs.sh all
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="${PUFPROBE_SERVER_CONFIG:-$SCRIPT_DIR/configs/server_config.yaml}"
set -a
# shellcheck source=configs/server_config.yaml
source "$CONFIG_FILE"
set +a

EVID="${EVIDENCE_DB_ROOT:-/db/student/metagenome/evidence_dbs}"
MARKERS="$EVID/.markers"
mkdir -p "$MARKERS" "$LOG_DIR"
LOG="$LOG_DIR/12_evidence_dbs_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log() { echo "[12_evid $(date +%H:%M:%S)] $*"; }
die() { echo "[12_evid ERROR] $*" >&2; exit 1; }
done_step() { touch "$MARKERS/$1"; }
step_done() { [ -f "$MARKERS/$1" ]; }

STRING_VER="${STRING_VER:-12.5}"
STRING_BASE="https://stringdb-downloads.org/download"
GO_OBO_URL="http://purl.obolibrary.org/obo/go/go-basic.obo"
GOA_BASE="https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT"
MCSA_BASE="https://www.ebi.ac.uk/thornton-srv/m-csa"
# UniRef30 hhsuite build. The Steinegger lab S3 bucket (served through
# opendata.mmseqs.org) is fast and multi-connection friendly; GWDG throttles
# per-IP and repeatedly drops connections, so it is the last-resort mirror.
HH_S3_BASE="https://opendata.mmseqs.org/uniclust"
HH_GWDG_BASE="https://wwwuser.gwdg.de/~compbiol/uniclust"

# Required artifact: abort on failure after removing the partial file.
download() { # url out
    local url="$1" out="$2"
    if [ -s "$out" ]; then log "  exists: $out"; return 0; fi
    log "  downloading $url"
    mkdir -p "$(dirname "$out")"
    if command -v curl >/dev/null 2>&1; then
        curl -fL -C - --retry 5 --retry-delay 5 --connect-timeout 60 -o "$out" "$url" \
            || { rm -f "$out"; die "download failed: $url (re-run to resume)"; }
    else
        wget -c --tries 5 --timeout=120 --progress=dot:giga -O "$out" "$url" \
            || { rm -f "$out"; die "download failed: $url (re-run to resume)"; }
    fi
}

# Optional artifact: warn and continue when the remote is missing.
try_download() { # url out
    local url="$1" out="$2"
    if [ -s "$out" ]; then log "  exists: $out"; return 0; fi
    log "  downloading (optional) $url"
    mkdir -p "$(dirname "$out")"
    if command -v curl >/dev/null 2>&1; then
        curl -fL -C - --retry 3 --connect-timeout 60 -o "$out" "$url" \
            || { rm -f "$out"; log "  WARN: could not fetch $url; continuing"; return 1; }
    else
        wget -c --tries 3 --timeout=120 --progress=dot:giga -O "$out" "$url" \
            || { rm -f "$out"; log "  WARN: could not fetch $url; continuing"; return 1; }
    fi
}

# Large artifact: try each mirror in order. Prefer aria2c multi-connection
# segmented download; fall back to resumable single-connection curl/wget.
# Usage: fetch_parallel <out> <url> [mirror_url ...]
fetch_parallel() {
    local out="$1"; shift
    local urls=("$@")
    [ "${#urls[@]}" -ge 1 ] || die "fetch_parallel: no URL given for $out"
    if [ -s "$out" ]; then log "  exists: $out"; return 0; fi
    mkdir -p "$(dirname "$out")"
    local url
    for url in "${urls[@]}"; do
        log "  downloading $url"
        if command -v aria2c >/dev/null 2>&1; then
            # Do NOT set --lowest-speed-limit: when a host throttles every
            # connection below the threshold, aria2 aborts them all and fails
            # the whole job. --timeout only drops stalled (no-data) connections;
            # --max-tries=0 keeps reconnecting and -c resumes from the .aria2
            # control file when the job is restarted.
            if (cd "$(dirname "$out")" && aria2c -x 16 -s 16 -k 10M -c \
                --max-tries=0 --retry-wait=5 --connect-timeout=60 \
                --timeout=60 --file-allocation=none \
                -o "$(basename "$out")" "$url") && [ -s "$out" ]; then
                return 0
            fi
            log "  mirror failed: $url (trying next mirror)"
        else
            log "  aria2c not found; using single-connection curl (install aria2 for speed)"
            if command -v curl >/dev/null 2>&1; then
                curl -fL -C - --retry 5 --retry-delay 5 --connect-timeout 60 -o "$out" "$url" \
                    && [ -s "$out" ] && return 0
            elif command -v wget >/dev/null 2>&1; then
                wget -c --tries 5 --timeout=120 --progress=dot:giga -O "$out" "$url" \
                    && [ -s "$out" ] && return 0
            fi
            log "  mirror failed: $url (trying next mirror)"
        fi
    done
    rm -f "$out" "$out.aria2"
    die "download failed from every mirror: $out"
}

# ---------------------------------------------------------------------------
# STRING
# ---------------------------------------------------------------------------
string_db() {
    local v="$STRING_VER"
    local d="$EVID/string/v$v"
    mkdir -p "$d/per_taxa"
    log "=== STRING v$v: COG-level association links ==="
    download "$STRING_BASE/COG.links.v$v.txt.gz"          "$d/COG.links.v$v.txt.gz"
    download "$STRING_BASE/COG.links.detailed.v$v.txt.gz" "$d/COG.links.detailed.v$v.txt.gz"
    log "=== STRING v$v: species list and tree ==="
    download "$STRING_BASE/species.v$v.txt"      "$d/species.v$v.txt"
    download "$STRING_BASE/species.tree.v$v.txt" "$d/species.tree.v$v.txt"

    log "=== STRING v$v: per-species networks ==="
    local taxids=""
    if [ -n "${STRING_TAXIDS_FILE:-}" ] && [ -f "${STRING_TAXIDS_FILE}" ]; then
        taxids="$(grep -vE '^\s*($|#)' "${STRING_TAXIDS_FILE}" || true)"
    fi
    taxids="${taxids:-${STRING_TAXIDS:-99287 511145 386585 272620 584 35703 300267 169963 214092 243277 192222 226185 93061 208964 400668}}"
    local t
    for t in $taxids; do
        try_download \
            "$STRING_BASE/protein.links.v$v/${t}.protein.links.v$v.txt.gz" \
            "$d/per_taxa/${t}.protein.links.v$v.txt.gz" || true
    done

    if [ "${STRING_WITH_ALIASES:-0}" = "1" ]; then
        log "=== STRING v$v: global protein.aliases (~3.2 GB) ==="
        download "$STRING_BASE/protein.aliases.v$v.txt.gz" "$d/protein.aliases.v$v.txt.gz"
    else
        log "skipping global protein.aliases (set STRING_WITH_ALIASES=1 to fetch ~3.2 GB)"
    fi
    done_step "string_$v"
}

# ---------------------------------------------------------------------------
# GO ontology (+ optional GOA reference-proteome GAF)
# ---------------------------------------------------------------------------
go_db() {
    local d="$EVID/go"
    mkdir -p "$d"
    log "=== GO: go-basic.obo (true-path rule) ==="
    download "$GO_OBO_URL" "$d/go-basic.obo"
    # GOA no longer ships a prokaryotic-only aggregate (goa_pdr was retired).
    # The ontology is all that is needed for true-path propagation; the GO
    # labels themselves come from eggNOG 7 and the Swiss-Prot cross-references.
    # For BPO model training, set WITH_GOA=1 to fetch the GO reference-proteome
    # GAF (~23 GB gz, bacterial plus eukaryotic reference proteomes).
    if [ "${WITH_GOA:-0}" = "1" ]; then
        log "=== GOA: reference-proteome GAF (gcrp, ~23 GB gz) ==="
        download "$GOA_BASE/goa_uniprot_gcrp.gaf.gz" "$d/goa_uniprot_gcrp.gaf.gz"
    else
        log "skipping GOA GAF (set WITH_GOA=1 to fetch the ~23 GB reference-proteome GAF)"
    fi
    local env="${PUFPROBE_ENV:-pufprobe}"
    if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$env"; then
        log "ensuring goatools is available in conda env '$env'"
        conda run --no-capture-output -n "$env" \
            python -m pip install --quiet --upgrade goatools \
            || log "  WARN: goatools install failed; install it manually"
    else
        log "conda env '$env' not found; skip goatools install (pip install goatools later)"
    fi
    done_step "go"
}

# ---------------------------------------------------------------------------
# M-CSA catalytic residues
# ---------------------------------------------------------------------------
mcsa_db() {
    local d="$EVID/mcsa"
    mkdir -p "$d"
    log "=== M-CSA: entries and catalytic residues (JSON API) ==="
    download "$MCSA_BASE/api/entries/?format=json"  "$d/mcsa_entries.json"
    download "$MCSA_BASE/api/residues/?format=json" "$d/mcsa_residues.json"
    done_step "mcsa"
}

# ---------------------------------------------------------------------------
# UniRef30 hhsuite (large)
# ---------------------------------------------------------------------------
hhsuite_db() {
    local rel="${EVIDENCE_UNIREF30_RELEASE:-2023_02}"
    local dir="$EVID/hhsuite/UniRef30_$rel"
    local marker="hhblits_uniref30_$rel"
    mkdir -p "$dir"
    if step_done "$marker" && [ -f "$DBDIR/hhblits_db_hhm.ffdata" ] && [ -f "$DBDIR/hhblits_db_cs219.ffdata" ]; then
        log "UniRef30 already linked into $DBDIR"; return 0
    fi
    log "free disk at $EVID (need ~250 GB for the tar plus the extracted build):"
    df -h "$EVID" | tail -1
    local tar="$dir/UniRef30_${rel}_hhsuite.tar.gz"
    # S3 places the file directly under /uniclust/; GWDG nests it under /<rel>/.
    fetch_parallel "$tar" \
        "$HH_S3_BASE/UniRef30_${rel}_hhsuite.tar.gz" \
        "$HH_GWDG_BASE/${rel}/UniRef30_${rel}_hhsuite.tar.gz"
    log "extracting UniRef30 (tens of GB) ..."
    tar xzf "$tar" -C "$dir"
    local stem
    stem="$(find "$dir" -name "*_hhm.ffdata" -print -quit | sed 's/_hhm.ffdata$//')"
    [ -n "$stem" ] || die "UniRef30 extraction did not contain *_hhm.ffdata"
    mkdir -p "$DBDIR"
    # HH-suite3 ships column-state files as cs219 (legacy HH-suite2 used "cs").
    local suffix
    for suffix in hhm.ffdata hhm.ffindex a3m.ffdata a3m.ffindex cs219.ffdata cs219.ffindex; do
        [ -e "${stem}_${suffix}" ] || die "UniRef30 extraction missing ${stem}_${suffix}"
        ln -sfn "${stem}_${suffix}" "$DBDIR/hhblits_db_${suffix}"
        [ -e "$DBDIR/hhblits_db_${suffix}" ] || die "hhblits_db_${suffix} link is dangling"
    done
    rm -f "$tar"
    log "UniRef30 linked: $DBDIR/hhblits_db_* -> $(basename "$stem")_*"
    done_step "$marker"
}

case "${1:-small}" in
    small)   string_db; go_db; mcsa_db ;;
    string)  string_db ;;
    go)      go_db ;;
    mcsa)    mcsa_db ;;
    hhsuite) hhsuite_db ;;
    all)     string_db; go_db; mcsa_db; hhsuite_db ;;
    *) die "usage: $0 {small|string|go|mcsa|hhsuite|all}" ;;
esac

log "Evidence database layout under $EVID:"
find "$EVID" -maxdepth 3 \( -type f -o -type l \) -not -path "*/.markers/*" 2>/dev/null | sort | sed 's/^/  /'
log "Done. HHblits prefix is linked into $DBDIR/hhblits_db_*."
