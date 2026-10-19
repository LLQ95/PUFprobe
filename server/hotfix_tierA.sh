#!/usr/bin/env bash
# =============================================================================
# hotfix_tierA.sh - 修复 v0.2.1 及更早版本 01_download_dbs.sh 在 Tier A 的下载
#
# 问题（已在 v0.2.2 修复）：
#   1) MMseqs2 v18 中 Swiss-Prot 的库名是 "UniProtKB/Swiss-Prot"，旧脚本的
#      "Swiss-Prot" 报 "Selected database ... not found"；
#   2) 旧回退 mmseqs.com/.../swissprot.tar.gz 已 302 到网页，下载到 HTML；
#   3) CARD 旧直链 /download/6/*.fasta 现在返回 HTML。
#
# 本脚本直接用已验证可用的官方源把 Tier A 冒烟数据库建好，并写入与
# 01_download_dbs.sh 一致的完成标记，之后可直接进入第 3 步（02_annotate.sh）。
# 幂等：已存在且完整的库会跳过，可中断重跑。
#
# 用法：
#   bash server/hotfix_tierA.sh                 # 建库（含 Foldseek PDB）
#   SKIP_FOLDSEEK=1 bash server/hotfix_tierA.sh # 不下载 Foldseek PDB（约省 10-15GB）
# =============================================================================
set -euo pipefail

# ---- 路径（按实际服务器环境；与 server/configs/server_config.yaml 一致） ----
CONDA_BASE="${CONDA_BASE:-/data/liushiwei/anaconda3}"
ENV_NAME="pufprobe"
DB_ROOT="${DB_ROOT:-/data/liushiwei/PUFprobe/dbs}"
DBDIR="$DB_ROOT/dbdir"
TMP="$DB_ROOT/tmp"
DONE="$DB_ROOT/done"
THREADS="${THREADS:-112}"
UNIPROT_BASE="${UNIPROT_BASE:-https://ftp.uniprot.org/pub/databases/uniprot/current_release}"

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

mkdir -p "$DBDIR" "$TMP" "$DONE" "$DB_ROOT/card" "$DB_ROOT/vfdb"
log() { echo "[hotfix_tierA $(date +%H:%M:%S)] $*"; }

# ---- 0. 清理旧版失败留下的 HTML 假文件 ---------------------------------------
rm -f "$TMP/swissprot.tar.gz" "$TMP/Swiss-Prot.tar.gz" 2>/dev/null || true
log "cleaned stale partial files under $TMP"

# ---- 1. Swiss-Prot：UniProt FASTA 自建 MMseqs2 库 ----------------------------
if [ -f "$DBDIR/swissprot.index" ] && [ -f "$DONE/mmseqs_swissprot" ]; then
    log "Swiss-Prot MMseqs2 DB already present, skip"
else
    log "Downloading Swiss-Prot FASTA from UniProt (~90 MB gz) ..."
    cd "$TMP"
    wget -c -q --show-progress "${UNIPROT_BASE}/knowledgebase/complete/uniprot_sprot.fasta.gz" \
        -O uniprot_sprot.fasta.gz || wget -c "${UNIPROT_BASE}/knowledgebase/complete/uniprot_sprot.fasta.gz"
    gunzip -f uniprot_sprot.fasta.gz
    head -1 uniprot_sprot.fasta | grep -q '^>' || { echo "ERROR: Swiss-Prot FASTA invalid"; exit 1; }
    log "Creating MMseqs2 database (createdb + createindex) ..."
    rm -f "$DBDIR/swissprot"* 2>/dev/null || true
    mmseqs createdb uniprot_sprot.fasta "$DBDIR/swissprot" --dbtype 1
    mmseqs createindex "$DBDIR/swissprot" "$TMP" --threads "$THREADS"
    mmseqs linfo "$DBDIR/swissprot" | head -3
    touch "$DONE/mmseqs_swissprot"
    log "Swiss-Prot done"
fi

# ---- 2. CARD：官方 latest/data 数据包 ----------------------------------------
CARD_FA="$DB_ROOT/card/protein_fasta_protein_homolog_model.fasta"
if [ -s "$CARD_FA" ] && [ -f "$DONE/card_vfdb" ]; then
    log "CARD/VFDB already present, skip"
else
    log "Downloading CARD canonical bundle (latest/data, ~4 MB) ..."
    CARD_BZ2="$DB_ROOT/card/card_latest.tar.bz2"
    wget -c -q --show-progress "https://card.mcmaster.ca/latest/data" -O "$CARD_BZ2" \
        || wget -c "https://card.mcmaster.ca/latest/data" -O "$CARD_BZ2"
    if ! tar -tjf "$CARD_BZ2" >/dev/null 2>&1; then
        log "bundle incomplete, retrying once"
        rm -f "$CARD_BZ2"; wget "https://card.mcmaster.ca/latest/data" -O "$CARD_BZ2"
    fi
    tar -xjf "$CARD_BZ2" -C "$DB_ROOT/card" ./protein_fasta_protein_homolog_model.fasta
    n_card=$(grep -c '^>' "$CARD_FA" || echo 0)
    log "CARD records: $n_card"
    [ "$n_card" -gt 0 ] || { echo "ERROR: CARD FASTA empty"; exit 1; }

    # ---- 3. VFDB：set A 蛋白集 ----
    log "Downloading VFDB set A protein (~1.3 MB) ..."
    wget -c -q --show-progress "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz" \
        -O "$DB_ROOT/vfdb/VFDB_setA_pro.fas.gz" || \
        wget -c "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz" -O "$DB_ROOT/vfdb/VFDB_setA_pro.fas.gz"
    gunzip -f "$DB_ROOT/vfdb/VFDB_setA_pro.fas.gz"
    n_vf=$(grep -c '^>' "$DB_ROOT/vfdb/VFDB_setA_pro.fas" || echo 0)
    log "VFDB records: $n_vf"
    [ "$n_vf" -gt 0 ] || { echo "ERROR: VFDB FASTA empty"; exit 1; }
    touch "$DONE/card_vfdb"
fi

# ---- 4. Foldseek PDB（可选，失败不阻塞） -------------------------------------
if [ "${SKIP_FOLDSEEK:-0}" = "1" ]; then
    log "SKIP_FOLDSEEK=1, skipping structure database (structure channel will be off)"
elif [ -f "$DBDIR/foldseek_db.index" ] && [ -f "$DONE/foldseek_smoke" ]; then
    log "Foldseek PDB already present, skip"
else
    log "Building Foldseek PDB via official command (reads search.foldseek.com, ~10-15 GB) ..."
    if foldseek databases PDB "$DBDIR/foldseek_db" "$TMP" --threads "$THREADS"; then
        touch "$DONE/foldseek_smoke"
        log "Foldseek PDB done"
    else
        log "WARN: Foldseek PDB failed (network). This only disables the structure channel."
        log "      Re-run later with: foldseek databases PDB $DBDIR/foldseek_db $TMP"
        log "      or set SKIP_FOLDSEEK=1. The remaining pipeline is unaffected."
    fi
fi

# ---- 5. ESM-2 35M 权重（可选，冒烟回退特征不需要，失败不阻塞） ---------------
if [ "${SKIP_ESM:-0}" = "1" ]; then
    log "SKIP_ESM=1, skipping ESM-2 weight preload"
elif [ -f "$DONE/esm2_t12" ]; then
    log "ESM-2 35M weights already cached, skip"
else
    log "Pre-loading ESM-2 35M weights (optional) ..."
    if python - <<'PY'
import esm
esm.pretrained.load_model_and_alphabet("esm2_t12_35M_UR50D")
print("preloaded esm2_t12_35M")
PY
    then touch "$DONE/esm2_t12"; else
        log "WARN: ESM preload failed (network to dl.fbaipublicfiles.com); harmless for smoke."
    fi
fi

echo
log "==== Tier A databases ready under $DBDIR ===="
ls -lh "$DBDIR"/swissprot.index 2>/dev/null || true
ls -lh "$DBDIR"/foldseek_db.index 2>/dev/null || echo "  (foldseek PDB not built)"
echo "  CARD records: $(grep -c '^>' "$DB_ROOT/card/protein_fasta_protein_homolog_model.fasta" || echo 0)"
echo "  VFDB records: $(grep -c '^>' "$DB_ROOT/vfdb/VFDB_setA_pro.fas" || echo 0)"
log "Next: edit CONTIGS_DIR in server/configs/server_config.yaml, then"
log "      MAX_SAMPLES=3 bash server/02_annotate.sh"
