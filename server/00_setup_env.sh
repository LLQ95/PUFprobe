#!/usr/bin/env bash
#
# 00_setup_env.sh — create the standalone pufprobe conda environment.
#
# This script never touches the existing bakta_env / checkm environments. It
# creates a new environment (default name: pufprobe) with MMseqs2, Foldseek,
# HH-suite, P2Rank, scikit-learn, torch (GPU or CPU auto-detected) and the
# PUFprobe package itself installed editable.
#
# Re-runnable: an existing environment is reused and only missing packages are
# installed.
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
die() { echo "[00_setup_env ERROR] $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 1. Locate conda (or mamba)
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

RUNNER="conda"
if command -v mamba >/dev/null 2>&1; then
    RUNNER="mamba"
fi
log "Using $RUNNER ($CONDA_BIN)"

# ---------------------------------------------------------------------------
# 2. Detect GPU and select a CUDA tag
# ---------------------------------------------------------------------------
GPU=0
TORCH_INDEX_URL="$TORCH_INDEX_URL_CPU"
CUDA_TAG=""
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    GPU=1
    CUDA_VERSION="$(nvidia-smi | awk '/CUDA Version/ {for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+\.[0-9]+/) {print $i; exit}}')"
    CUDA_MAJOR="$(echo "${CUDA_VERSION:-0.0}" | cut -d. -f1)"
    CUDA_MINOR="$(echo "${CUDA_VERSION:-0.0}" | cut -d. -f2)"
    log "GPU detected (driver CUDA runtime ${CUDA_VERSION:-unknown})"
    if [ "$CUDA_MAJOR" -ge 12 ] || { [ "$CUDA_MAJOR" -eq 11 ] && [ "$CUDA_MINOR" -ge 8 ]; }; then
        if [ "$CUDA_MAJOR" -ge 12 ]; then CUDA_TAG="cu121"; else CUDA_TAG="cu118"; fi
        TORCH_INDEX_URL="$(echo "$TORCH_INDEX_URL_GPU" | sed "s/cu121/$CUDA_TAG/")"
    else
        log "Driver CUDA version older than 11.8; falling back to CPU torch. Update the NVIDIA driver for GPU support."
    fi
else
    log "No nvidia-smi / GPU available; installing CPU-only torch."
fi

# ---------------------------------------------------------------------------
# 3. Create the environment with China-accessible mirrors
# ---------------------------------------------------------------------------
CF="$CONDA_MIRROR_HOST/cloud/conda-forge"
BIO="$CONDA_MIRROR_HOST/cloud/bioconda"
PKG="$CONDA_MIRROR_HOST/pkgs/main"

if ! conda env list | awk '{print $1}' | grep -qx "$PUFPROBE_ENV"; then
    log "Creating conda environment '$PUFPROBE_ENV' (python $PYTHON_VERSION) ..."
    "$RUNNER" create -y -n "$PUFPROBE_ENV" \
        --override-channels \
        -c "$CF" -c "$BIO" -c "$PKG" \
        python="$PYTHON_VERSION"
else
    log "Environment '$PUFPROBE_ENV' already exists; reusing it."
fi

# conda-forge first, then bioconda; bioconda tools must never pull a different
# python from a lower-priority channel.
log "Installing bioinformatics tools (mmseqs2, foldseek, hhsuite, p2rank, scikit-learn) ..."
"$RUNNER" install -y -n "$PUFPROBE_ENV" \
    --override-channels -c "$CF" -c "$BIO" -c "$PKG" \
    mmseqs2 foldseek hhsuite p2rank scikit-learn numpy click pip

conda activate "$PUFPROBE_ENV"

# ---------------------------------------------------------------------------
# 4. torch (GPU or CPU) + fair-esm
# ---------------------------------------------------------------------------
log "Installing torch from $TORCH_INDEX_URL"
if ! python -c "import torch" >/dev/null 2>&1; then
    pip install -i "$PIP_INDEX_URL" torch --index-url "$TORCH_INDEX_URL" \
        || {
            echo "Mirror torch install failed; retrying with official PyTorch index."
            if [ "$GPU" -eq 1 ]; then
                pip install torch --index-url "https://download.pytorch.org/whl/$CUDA_TAG"
            else
                pip install torch --index-url "https://download.pytorch.org/whl/cpu"
            fi
        }
else
    log "torch already importable ($(python -c 'import torch; print(torch.__version__)' 2>/dev/null || echo unknown)); skipped."
fi

log "Installing fair-esm ..."
pip install -i "$PIP_INDEX_URL" fair-esm || pip install fair-esm

# ---------------------------------------------------------------------------
# 5. Install PUFprobe editable
# ---------------------------------------------------------------------------
log "Installing PUFprobe from $REPO_ROOT ..."
pip install -i "$PIP_INDEX_URL" -e "$REPO_ROOT" || pip install -e "$REPO_ROOT"

# ---------------------------------------------------------------------------
# 6. Verification
# ---------------------------------------------------------------------------
log "Tool versions:"
mmseqs version 2>/dev/null | head -1 | sed 's/^/  mmseqs: /'  || true
foldseek version 2>/dev/null | head -1 | sed 's/^/  foldseek: /' || true
hhblits 2>&1 | head -2 | grep -i version | sed 's/^/  hhblits: /' || true
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
