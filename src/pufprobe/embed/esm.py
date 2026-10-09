"""ESM-2 sequence embeddings with fully lazy heavy imports.

torch and fair-esm are imported only inside :func:`compute_embeddings`, so the
rest of PUFprobe (including the pure-numpy fallback path) stays importable on a
machine without them. When torch, the ``esm`` package or model weights are not
available, a module-local :class:`ModelUnavailable` exception is raised carrying
a one-line install instruction; callers (the embed pipeline) catch it and fall
back to hand-crafted sequence features.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

from pufprobe.logging_utils import get_logger

logger = get_logger(__name__)

INSTALL_HINT = (
    "Install PyTorch and fair-esm to enable ESM-2 embeddings: "
    "`pip install fair-esm torch` (for GPU acceleration install a CUDA build of "
    "torch, see https://pytorch.org/get-started/)."
)


class ModelUnavailable(RuntimeError):
    """Raised when torch / fair-esm / model weights cannot be loaded.

    The message always ends with the install guidance so an operator reading
    the log knows exactly what to do.
    """


def compute_embeddings(
    records: Sequence,
    model_name: str,
    cpu: bool = False,
    batch_size: int = 8,
) -> Dict[str, np.ndarray]:
    """Return ``{global_key: L2-ready residue-mean vector}`` for every record.

    Parameters
    ----------
    records:
        Iterable of objects exposing ``.key`` (``sample_id@locus_tag``) and
        ``.sequence_aa``.
    model_name:
        A fair-esm checkpoint name such as ``esm2_t33_650M_UR50D``.
    cpu:
        Force CPU inference even when CUDA is available.
    batch_size:
        Number of sequences per forward pass.

    Raises
    ------
    ModelUnavailable
        If torch, fair-esm or the checkpoint weights cannot be loaded (no
        install, no GPU/CUDA runtime, or no network to download weights).
    """

    # --- lazy heavy imports -------------------------------------------------
    try:
        import torch  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on host install
        raise ModelUnavailable(
            f"PyTorch is not available ({exc}). {INSTALL_HINT}"
        ) from exc

    try:
        import esm  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on host install
        raise ModelUnavailable(
            f"The fair-esm package is not importable ({exc}). {INSTALL_HINT}"
        ) from exc

    # --- load model + alphabet (downloads weights on first use) --------------
    try:
        model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    except Exception as exc:  # pragma: no cover - network / weights
        raise ModelUnavailable(
            f"Could not load ESM checkpoint '{model_name}' ({exc}). {INSTALL_HINT}"
        ) from exc

    device = "cpu"
    if not cpu and torch.cuda.is_available():
        device = "cuda"
    logger.info("Loaded ESM model '%s' on device '%s' (%d records queued)",
                model_name, device, len(records))

    try:
        model = model.to(device)
        model.eval()
        batch_converter = alphabet.get_batch_converter()
        num_layers = model.num_layers

        result: Dict[str, np.ndarray] = {}
        items = [(rec.key, rec.sequence_aa) for rec in records]

        for start in range(0, len(items), max(1, batch_size)):
            batch = items[start : start + max(1, batch_size)]
            labels = [label for label, _ in batch]
            _, _, batch_tokens = batch_converter(batch)
            batch_tokens = batch_tokens.to(device)

            with torch.no_grad():
                out = model(batch_tokens, repr_layers=[num_layers],
                            return_contacts=False)
            reps = out["representations"][num_layers]  # (B, L, D)

            for i, (label, seq) in enumerate(batch):
                # Mean over real residue tokens only: drop the prepended <cls>
                # (index 0), trailing <eos> and any padding. Lengths from the
                # alphabet include the two special tokens.
                seq_len = len(seq)
                token_len = min(seq_len, reps.shape[1] - 2)
                residue_tokens = reps[i, 1 : 1 + token_len]
                vector = residue_tokens.mean(dim=0)
                result[label] = vector.detach().cpu().numpy().astype(np.float32)
    except ModelUnavailable:
        raise
    except Exception as exc:  # pragma: no cover - forward-pass failure
        raise ModelUnavailable(
            f"ESM forward pass failed for model '{model_name}' ({exc}). "
            f"{INSTALL_HINT}"
        ) from exc
    finally:
        del model
        if device == "cuda":
            torch.cuda.empty_cache()

    return result
