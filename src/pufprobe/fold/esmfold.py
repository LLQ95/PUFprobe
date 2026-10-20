"""ESMFold structure prediction with fully lazy torch imports.

Torch and the ESMFold weights are loaded once through :func:`load_esmfold`; the
fold pipeline loads the model a single time before iterating over families and
calls :func:`fold_with_model` for each representative. :func:`fold_sequence` is
a convenience that combines both steps, used by the standalone call and the
tests. Without torch or weights, :class:`ModelUnavailable` is raised with
install guidance, and the fold pipeline marks the families skipped rather than
crashing.

ESMFold convention: the B-factor column of the emitted PDB holds the
per-residue pLDDT on a 0-100 scale.
"""

from __future__ import annotations

from typing import Any, Optional, Tuple

import numpy as np

from pufprobe.logging_utils import get_logger

logger = get_logger(__name__)

INSTALL_HINT = (
    "Install PyTorch and the fair-esm ESMFold weights to enable end-to-end "
    "folding: `pip install fair-esm torch` (GPU recommended; the ~3B-parameter "
    "ESMFold model downloads weights on first use)."
)


class ModelUnavailable(RuntimeError):
    """Raised when torch / ESMFold weights cannot be loaded."""


class LoadedESMFold:
    """A loaded ESMFold model bound to its inference device."""

    def __init__(self, model: Any, device: str, torch_module: Any) -> None:
        self.model = model
        self.device = device
        self.torch = torch_module


def load_esmfold() -> LoadedESMFold:
    """Load torch and the ESMFold model once, bound to the available device."""
    try:
        import torch  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise ModelUnavailable(
            f"PyTorch is not available ({exc}). {INSTALL_HINT}"
        ) from exc

    try:
        model = torch.hub.load("facebookresearch/esm:main", "esmfold_v1")
    except Exception as exc:  # pragma: no cover
        raise ModelUnavailable(
            f"Could not load ESMFold weights ({exc}). {INSTALL_HINT}"
        ) from exc

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device).eval()
    logger.info("Loaded ESMFold on device '%s'", device)
    return LoadedESMFold(model, device, torch)


def fold_with_model(
    loaded: LoadedESMFold, sequence: str
) -> Tuple[str, float, Optional[float]]:
    """Run a single sequence through a loaded model.

    Returns ``(pdb_str, median_plddt, mean_pae)``. ``median_plddt`` is on the
    0-100 scale; ``mean_pae`` may be ``None`` when the backend does not expose it.
    """
    try:
        with loaded.torch.no_grad():
            output = loaded.model.infer(sequence)
        pdb_str = loaded.model.output_to_pdb(output)[0]
        plddt = float(np.median(output["plddt"].cpu().numpy()))
        mean_pae = None
        if "predicted_aligned_error" in output:
            pae = output["predicted_aligned_error"].cpu().numpy()
            mean_pae = float(np.mean(pae))
    except ModelUnavailable:
        raise
    except Exception as exc:  # pragma: no cover
        raise ModelUnavailable(
            f"ESMFold forward pass failed ({exc}). {INSTALL_HINT}"
        ) from exc
    return pdb_str, plddt, mean_pae


def fold_sequence(sequence: str) -> Tuple[str, float, Optional[float]]:
    """Fold one sequence, loading the model on demand and releasing it after."""
    loaded = load_esmfold()
    try:
        return fold_with_model(loaded, sequence)
    finally:
        if loaded.device == "cuda":
            loaded.torch.cuda.empty_cache()


# ---------------------------------------------------------------------------
# Shared PDB B-factor parsing (used by precomputed and esmfold paths)
# ---------------------------------------------------------------------------

def parse_plddt_from_pdb(pdb_text: str) -> Optional[float]:
    """Return median per-residue pLDDT from a PDB's B-factor column.

    Auto-detects whether B-factors are on a 0-1 or 0-100 scale and normalizes
    the result to 0-100. Prefers CA atoms; falls back to all ATOM records.
    """
    ca_values: list = []
    all_values: list = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            bfactor = float(line[60:66])
            atom = line[12:16].strip()
        except (ValueError, IndexError):
            continue
        all_values.append(bfactor)
        if atom == "CA":
            ca_values.append(bfactor)

    values = ca_values or all_values
    if not values:
        return None
    arr = np.asarray(values, dtype=float)
    if arr.max() <= 1.0 + 1e-6:
        arr = arr * 100.0
    return float(np.median(arr))
