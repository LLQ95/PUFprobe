"""ESM-2 sequence embeddings with fully lazy heavy imports.

torch and fair-esm are imported only inside :func:`load_esm_model`, so the rest
of PUFprobe (including the pure-numpy fallback path) stays importable on a
machine without them. When torch, the ``esm`` package or model weights are not
available, a module-local :class:`ModelUnavailable` exception is raised carrying
a one-line install instruction; callers (the embed pipeline) catch it and fall
back to hand-crafted sequence features.

The model loading and the residue-mean pooling are split into
:func:`load_esm_model` and :func:`embed_items` so that bulk jobs that embed
many independent batches (for example the annotated reference library built by
``server/build_plm_reference.py``) load the checkpoint once and reuse the exact
same pooling as the query stage.
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence, Tuple

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


def load_esm_model(model_name: str, cpu: bool = False):
    """Load an ESM checkpoint and return the model plus batching handles.

    Returns ``(model, alphabet, batch_converter, device, num_layers)``.
    Raises :class:`ModelUnavailable` when torch, fair-esm or the checkpoint
    weights cannot be loaded.
    """
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

    try:
        model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    except Exception as exc:  # pragma: no cover - network / weights
        raise ModelUnavailable(
            f"Could not load ESM checkpoint '{model_name}' ({exc}). {INSTALL_HINT}"
        ) from exc

    device = "cpu"
    if not cpu and torch.cuda.is_available():
        device = "cuda"
    model = model.to(device)
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    num_layers = model.num_layers
    logger.info("Loaded ESM model '%s' on device '%s'.", model_name, device)
    return model, alphabet, batch_converter, device, num_layers


def embed_items(
    model: Any,
    alphabet: Any,
    batch_converter: Any,
    num_layers: int,
    device: str,
    items: Sequence[Tuple[str, str]],
    batch_size: int = 8,
) -> Dict[str, np.ndarray]:
    """Embed ``(key, sequence)`` items, returning ``{key: mean-pooled vector}``.

    Vectors are the mean of the final-layer residue representations over the
    real residue tokens only (the prepended ``<cls>``, trailing ``<eos>`` and
    any padding are excluded), in float32.
    """
    import torch  # type: ignore

    result: Dict[str, np.ndarray] = {}
    batch_size = max(1, int(batch_size))
    for start in range(0, len(items), batch_size):
        batch = list(items[start : start + batch_size])
        labels = [label for label, _ in batch]
        _, _, batch_tokens = batch_converter(batch)
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            out = model(batch_tokens, repr_layers=[num_layers], return_contacts=False)
        reps = out["representations"][num_layers]  # (B, L, D)

        for i, (label, seq) in enumerate(batch):
            token_len = min(len(seq), reps.shape[1] - 2)
            residue_tokens = reps[i, 1 : 1 + token_len]
            vector = residue_tokens.mean(dim=0)
            result[label] = vector.detach().cpu().numpy().astype(np.float32)
    return result


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
    handles = load_esm_model(model_name, cpu=cpu)
    model, alphabet, batch_converter, device, num_layers = handles
    logger.info("Embedding %d records (batch_size=%d).", len(records), batch_size)
    try:
        items = [(rec.key, rec.sequence_aa) for rec in records]
        return embed_items(
            model, alphabet, batch_converter, num_layers, device,
            items, batch_size=batch_size,
        )
    except Exception as exc:  # pragma: no cover - forward-pass failure
        if isinstance(exc, ModelUnavailable):
            raise
        raise ModelUnavailable(
            f"ESM forward pass failed for model '{model_name}' ({exc}). "
            f"{INSTALL_HINT}"
        ) from exc
    finally:
        del model
        if device == "cuda":
            import torch  # type: ignore
            torch.cuda.empty_cache()
