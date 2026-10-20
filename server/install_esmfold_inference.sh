#!/usr/bin/env bash
#
# install_esmfold_inference.sh — prepare fair-esm ESMFold for inference.
#
# ESMFold needs OpenFold, whose official installer requires nvcc and pulls in
# heavy training-only packages (DeepSpeed, pytorch-lightning). Those are not
# needed for inference: this script installs only the pure-Python
# dependencies, places a vendored, compiler-free OpenFold into the active
# environment, and applies the fair-esm Python 3.11+ dataclass compatibility
# patch. It is re-runnable and installs nothing for a component that already
# imports.
#
# Run with the target environment active, e.g.
#     conda run -n pufprobe bash server/install_esmfold_inference.sh
# or
#     conda activate pufprobe && bash server/install_esmfold_inference.sh
#
# ESMFold weights (~10 GB) are downloaded separately on first use (or run the
# pre-download command printed at the end).
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VENDOR_DIR="$SCRIPT_DIR/vendor/esmfold_inference"

PYTHON_BIN="${PYTHON_BIN:-python}"

log()  { echo "[install_esmfold $(date +%H:%M:%S)] $*"; }
warn() { echo "[install_esmfold WARN] $*"; }
die()  { echo "[install_esmfold ERROR] $*" >&2; exit 1; }

PIP_INDEX_URL="${PIP_INDEX_URL:-https://pypi.tuna.tsinghua.edu.cn/simple}"

# ---------------------------------------------------------------------------
# 0. Pre-flight checks
# ---------------------------------------------------------------------------
"$PYTHON_BIN" -c 'import torch' 2>/dev/null \
  || die "torch is not installed in this environment. Run server/00_setup_env.sh first."
"$PYTHON_BIN" -c 'import esm' 2>/dev/null \
  || die "fair-esm is not installed. Run server/00_setup_env.sh (or: pip install fair-esm)."

[ -f "$VENDOR_DIR/install_vendor.py" ] \
  || die "Vendor directory missing: $VENDOR_DIR"

# ---------------------------------------------------------------------------
# 1. Pure-Python dependencies (no DeepSpeed / pytorch-lightning; the vendor
#    supplies inference-only stubs for those).
# ---------------------------------------------------------------------------
PURE_PKGS=(omegaconf biopython ml_collections scipy dm-tree einops)

need=()
for pkg in "${PURE_PKGS[@]}"; do
  mod="$pkg"
  case "$pkg" in
    ml_collections) mod="ml_collections" ;;
    dm-tree) mod="tree" ;;
    biopython) mod="Bio" ;;
  esac
  if ! "$PYTHON_BIN" -c "import ${mod}" 2>/dev/null; then
    need+=("$pkg")
  fi
done

if [ "${#need[@]}" -gt 0 ]; then
  log "Installing pure-Python ESMFold dependencies: ${need[*]}"
  if ! "$PYTHON_BIN" -m pip install -i "$PIP_INDEX_URL" "${need[@]}"; then
    warn "Install from configured index failed; falling back to PyPI."
    "$PYTHON_BIN" -m pip install "${need[@]}"
  fi
else
  log "All pure-Python ESMFold dependencies already present."
fi

# ---------------------------------------------------------------------------
# 2. Install the compiler-free OpenFold vendor and training stubs.
# ---------------------------------------------------------------------------
log "Installing inference-only OpenFold vendor ..."
"$PYTHON_BIN" "$VENDOR_DIR/install_vendor.py"

# ---------------------------------------------------------------------------
# 3. Apply the fair-esm Python 3.11+ nested-dataclass compatibility patch.
# ---------------------------------------------------------------------------
log "Applying fair-esm Python 3.11+ compatibility patch (idempotent) ..."
"$PYTHON_BIN" "$SCRIPT_DIR/patch_esmfold_py311.py"

# ---------------------------------------------------------------------------
# 4. Verify the ESMFold class imports without downloading weights.
# ---------------------------------------------------------------------------
log "Verifying ESMFold import ..."
"$PYTHON_BIN" - <<'PYEOF'
import warnings
warnings.filterwarnings("ignore")
from esm.esmfold.v1.esmfold import ESMFold  # noqa: F401
print("ESMFold import OK (openfold + dependencies resolved).")
PYEOF

log "Setup complete. ESMFold weights (~10 GB) download on first use."
cat <<'EOF'
Pre-download weights into the PUFprobe cache (optional but recommended):
    TORCH_HOME=/data/liushiwei/PUFprobe/dbs/torch_cache \
      conda run -n pufprobe python -c "import esm; esm.pretrained.esmfold_v1(); print('weights ready')"
EOF
