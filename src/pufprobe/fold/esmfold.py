"""ESMFold structure prediction with fully lazy torch imports.

``torch`` and the ESMFold weights are loaded only inside :func:`fold_sequence`.
Without torch / weights a module-local :class:`ModelUnavailable` is raised with
install guidance; the fold pipeline catches it and marks families failed/low
confidence rather than crashing.

ESMFold convention: the B-factor column of the emitted PDB holds the
per-residue pLDDT on a 0-100 scale.
"""

from __future__ import annotations

from typing import Optional, Tuple

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


def fold_sequence(sequence: str) -> Tuple[str, float, Optional[float]]:
    """Fold a single protein sequence.

    Returns ``(pdb_str, median_plddt, mean_pae)``. ``median_plddt`` is on the
    0-100 scale; ``mean_pae`` may be ``None`` if the backend does not expose it.

    Raises
    ------
    ModelUnavailable
        If torch, the ESMFold weights or the forward pass are unavailable.
    """

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

    device = "cpu"
    if torch.cuda.is_available():
        device = "cuda"
    logger.info("Loaded ESMFold on device '%s' (seq len=%d)", device, len(sequence))

    try:
        model = model.to(device).eval()
        with torch.no_grad():
            output = model.infer(sequence)
        pdb_str = model.output_to_pdb(output)[0]
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
    finally:
        del model
        if device == "cuda":
            torch.cuda.empty_cache()

    return pdb_str, plddt, mean_pae


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
