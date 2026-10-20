#!/usr/bin/env bash
#
# 00_setup_env.sh — create the standalone pufprobe conda environment.
#
# This script never touches the existing bakta_env / checkm environments. It
# creates a new environment (default name: pufprobe) with MMseqs2, Foldseek,
# HH-suite, scikit-learn, torch (GPU or CPU auto-detected) and the PUFprobe
# package itself installed editable. P2Rank is not on bioconda; it is fetched
# as a portable GitHub release into thirdparty/ (optional, non-blocking).
#
#
# Re-runnable: an existing environment at the correct location is reused and
# only missing packages are installed.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="${PUFPROBE_SERVER_CONFIG:-$SCRIPT_DIR/configs/server_config.yaml}"
set -a
# shellcheck source=configs/server_config.yaml
source "$CONFIG_FILE"
set +a

mkdir -p "$LOG_DIR" "$MMSEQS_TMP"
LOG="$LOG_DIR/00_setup_env_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log() { echo "[00_setup_env $(date +%H:%M:%S)] $*"; }
warn() { echo "[00_setup_env WARN] $*"; }
die() { echo "[00_setup_env ERROR] $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 1. Locate conda (or mamba), and PIN the root prefix
# ---------------------------------------------------------------------------
CONDA_BIN=""
for candidate in "$(command -v conda 2>/dev/null || true)" \
                 "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
                 "/opt/conda/bin/conda" "/usr/local/miniconda3/bin/conda"; do
    if [ -n "$candidate" ] && [ -x "$candidate" ]; then
        CONDA_BIN="$candidate"
        break
    fi
done
[ -n "$CONDA_BIN" ] || die "conda not found. Install Miniconda first:
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh -b
then re-run this script."

# shellcheck disable=SC1091
source "$("$CONDA_BIN" info --base)/etc/profile.d/conda.sh"

CONDA_BASE="$("$CONDA_BIN" info --base)"
# Pin the root prefix so libmamba never mistakes the currently activated
# environment (e.g. checkm) for the root and nests envs under it.
export MAMBA_ROOT_PREFIX="$CONDA_BASE"
log "Conda base: $CONDA_BASE (MAMBA_ROOT_PREFIX pinned)"

# Prefer an absolute mamba path when present; otherwise use conda directly.
# Never rely on a bare "mamba" name resolved from inside another env.
RUNNER="$CONDA_BIN"
if [ -x "$CONDA_BASE/bin/mamba" ]; then
    RUNNER="$CONDA_BASE/bin/mamba"
elif command -v mamba >/dev/null 2>&1; then
    RUNNER="$(command -v mamba)"
fi
log "Using $RUNNER"

# ---------------------------------------------------------------------------
# 2. Detect GPU and select a CUDA tag
# ---------------------------------------------------------------------------
GPU=0
CUDA_TAG=""
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    GPU=1
    CUDA_VERSION="$(nvidia-smi | awk '/CUDA Version/ {for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+\.[0-9]+/) {print $i; exit}}')"
    CUDA_MAJOR="$(echo "${CUDA_VERSION:-0.0}" | cut -d. -f1)"
    CUDA_MINOR="$(echo "${CUDA_VERSION:-0.0}" | cut -d. -f2)"
    log "GPU detected (driver CUDA runtime ${CUDA_VERSION:-unknown})"
    if [ "$CUDA_MAJOR" -ge 12 ] || { [ "$CUDA_MAJOR" -eq 11 ] && [ "$CUDA_MINOR" -ge 8 ]; }; then
        if [ "$CUDA_MAJOR" -ge 12 ]; then CUDA_TAG="cu121"; else CUDA_TAG="cu118"; fi
    else
        log "Driver CUDA version older than 11.8; falling back to CPU torch. Update the NVIDIA driver for GPU support."
    fi
else
    log "No nvidia-smi / GPU available; installing CPU-only torch."
fi

# ---------------------------------------------------------------------------
# 3. Remove any incorrectly nested environment left by older script versions
# ---------------------------------------------------------------------------
PUF_ENV_PATH="$CONDA_BASE/envs/$PUFPROBE_ENV"
nested_envs="$("$CONDA_BIN" env list | awk '{print $NF}' | grep -E "/envs/[^/]+/envs/$PUFPROBE_ENV\$" || true)"
if [ -n "$nested_envs" ]; then
    while read -r bad_env; do
        if [ "$bad_env" != "$PUF_ENV_PATH" ]; then
            log "Removing incorrectly nested environment created by an older script: $bad_env"
            rm -rf "$bad_env"
        fi
    done <<< "$nested_envs"
fi

# ---------------------------------------------------------------------------
# 4. Create the environment with China-accessible mirrors
# ---------------------------------------------------------------------------
CF="$CONDA_MIRROR_HOST/cloud/conda-forge"
BIO="$CONDA_MIRROR_HOST/cloud/bioconda"
PKG="$CONDA_MIRROR_HOST/pkgs/main"

if [ ! -d "$PUF_ENV_PATH" ]; then
    log "Creating conda environment '$PUFPROBE_ENV' at $PUF_ENV_PATH (python $PYTHON_VERSION) ..."
    "$RUNNER" create -y -p "$PUF_ENV_PATH" \
        --override-channels \
        -c "$CF" -c "$BIO" -c "$PKG" \
        python="$PYTHON_VERSION"
else
    log "Environment '$PUFPROBE_ENV' already exists at $PUF_ENV_PATH; reusing it."
fi

# conda-forge first, then bioconda. P2Rank is NOT on bioconda and is installed
# separately in step 5b; do not add it here or the solver aborts.
log "Installing bioinformatics tools (mmseqs2, foldseek, hhsuite, scikit-learn) ..."
"$RUNNER" install -y -p "$PUF_ENV_PATH" \
    --override-channels -c "$CF" -c "$BIO" -c "$PKG" \
    mmseqs2 foldseek hhsuite scikit-learn numpy click pip openjdk=17

conda activate "$PUFPROBE_ENV"
export JAVA_HOME="$PUF_ENV_PATH"

# ---------------------------------------------------------------------------
# 5a. torch (GPU or CPU) with China mirrors first, official index last
# ---------------------------------------------------------------------------
log "Installing torch ..."
if ! python -c "import torch" >/dev/null 2>&1; then
    if [ "$GPU" -eq 1 ] && [ -n "$CUDA_TAG" ]; then
        TORCH_INDEX_URLS=(
            "https://mirrors.aliyun.com/pytorch-wheels/$CUDA_TAG/"
            "https://mirror.sjtu.edu.cn/pytorch-wheels/$CUDA_TAG/"
            "https://download.pytorch.org/whl/$CUDA_TAG"
        )
    else
        TORCH_INDEX_URLS=(
            "https://mirrors.aliyun.com/pytorch-wheels/cpu/"
            "https://mirror.sjtu.edu.cn/pytorch-wheels/cpu/"
            "https://download.pytorch.org/whl/cpu"
        )
    fi
    torch_ok=0
    for idx_url in "${TORCH_INDEX_URLS[@]}"; do
        log "Trying torch index: $idx_url"
        if pip install --timeout 120 torch --index-url "$idx_url"; then
            torch_ok=1
            break
        fi
        warn "torch install from $idx_url failed; trying next index."
    done
    [ "$torch_ok" -eq 1 ] || die "torch installation failed from all indexes. Retry later or install manually, then re-run this script."
else
    log "torch already importable ($(python -c 'import torch; print(torch.__version__)' 2>/dev/null || echo unknown)); skipped."
fi

log "Installing fair-esm and ESMFold runtime dependencies ..."
# einops is required by ESMFold but is not pulled in by fair-esm.
pip install -i "$PIP_INDEX_URL" fair-esm einops || pip install fair-esm einops

# ---------------------------------------------------------------------------
# 5b. P2Rank portable release (optional; pocket channel is non-essential)
# ---------------------------------------------------------------------------
P2RANK_DIR="$REPO_ROOT/thirdparty/p2rank"
P2RANK_VERSION="${P2RANK_VERSION:-2.5.1}"
if [ -x "$P2RANK_DIR/distro/prank" ]; then
    log "P2Rank already present at $P2RANK_DIR; skipped."
else
    log "Installing P2Rank $P2RANK_VERSION portable release (optional) ..."
    p2rank_ok=0
    P2RANK_TARBALL="p2rank_${P2RANK_VERSION}.tar.gz"
    P2RANK_URLS=(
        "https://github.com/rdk/p2rank/releases/download/${P2RANK_VERSION}/${P2RANK_TARBALL}"
        "https://ghfast.top/https://github.com/rdk/p2rank/releases/download/${P2RANK_VERSION}/${P2RANK_TARBALL}"
        "https://gh-proxy.com/https://github.com/rdk/p2rank/releases/download/${P2RANK_VERSION}/${P2RANK_TARBALL}"
    )
    mkdir -p "$REPO_ROOT/thirdparty"
    for u in "${P2RANK_URLS[@]}"; do
        log "Trying P2Rank download: $u"
        if wget -q --timeout=120 --tries=2 -O "$REPO_ROOT/thirdparty/$P2RANK_TARBALL" "$u"; then
            p2rank_ok=1
            break
        fi
        warn "P2Rank download failed from $u; trying next mirror."
    done
    if [ "$p2rank_ok" -eq 1 ]; then
        tar -xzf "$REPO_ROOT/thirdparty/$P2RANK_TARBALL" -C "$REPO_ROOT/thirdparty"
        rm -rf "$P2RANK_DIR"
        mv "$REPO_ROOT/thirdparty/p2rank_${P2RANK_VERSION}" "$P2RANK_DIR"
        rm -f "$REPO_ROOT/thirdparty/$P2RANK_TARBALL"
        chmod +x "$P2RANK_DIR/distro/prank" 2>/dev/null || true
        log "P2Rank installed at $P2RANK_DIR (launcher: distro/prank)."
    else
        warn "P2Rank download failed from all mirrors. The pocket evidence channel will be skipped; other channels are unaffected."
        warn "You can install it later by downloading $P2RANK_TARBALL and extracting it to $P2RANK_DIR."
    fi
fi
export P2RANK_HOME="$P2RANK_DIR"

# ---------------------------------------------------------------------------
# 6. Install PUFprobe editable
# ---------------------------------------------------------------------------
log "Installing PUFprobe from $REPO_ROOT ..."
pip install -i "$PIP_INDEX_URL" -e "$REPO_ROOT" || pip install -e "$REPO_ROOT"

# ---------------------------------------------------------------------------
# 7. Verification
# ---------------------------------------------------------------------------
log "Tool versions:"
mmseqs version 2>/dev/null | head -1 | sed 's/^/  mmseqs: /'  || true
foldseek version 2>/dev/null | head -1 | sed 's/^/  foldseek: /' || true
hhblits 2>&1 | head -2 | grep -i version | sed 's/^/  hhblits: /' || true
java -version 2>&1 | head -1 | sed 's/^/  java: /' || true
[ -x "$P2RANK_DIR/distro/prank" ] && echo "  p2rank: $P2RANK_DIR/distro/prank" || echo "  p2rank: not installed (optional)"
python - <<'PY'
try:
    import torch
    print(f"  torch: {torch.__version__} cuda_available={torch.cuda.is_available()}")
except Exception as exc:
    print(f"  torch: NOT IMPORTABLE ({exc})")
try:
    import esm
    print("  fair-esm: importable")
except Exception as exc:
    print(f"  fair-esm: NOT IMPORTABLE ({exc})")
PY

DOCTOR_LOG="$LOG_DIR/doctor_$(date +%Y%m%d_%H%M%S).json"
pufprobe doctor | tee "$DOCTOR_LOG"

log "Environment ready. Doctor report saved to $DOCTOR_LOG"
log "Next: bash server/01_download_dbs.sh   (TIER=A smoke databases first)"
