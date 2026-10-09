"""Dependency-free fallback embeddings built from :mod:`pufprobe.features`.

These vectors keep the pipeline runnable without torch / ESM weights. The raw
descriptor is 30 floats (20 amino-acid frequencies + 10 physicochemical). It is
padded with zeros up to ``config["embed"]["fallback_feature_dim"]`` so that the
embedding dimension stays consistent regardless of which source produced the
vectors; the manifest records how many trailing dimensions were zero-padded.
"""

from __future__ import annotations

from typing import Dict, Sequence, Tuple

import numpy as np

from pufprobe.features import fallback_features
from pufprobe.logging_utils import get_logger

logger = get_logger(__name__)


def fallback_vector(sequence: str, dim: int) -> Tuple[np.ndarray, int]:
    """Return a length-``dim`` float32 vector for ``sequence``.

    Returns the vector and the number of components actually filled by the
    raw descriptor (the rest are zeros).
    """
    raw = fallback_features(sequence)
    filled = min(len(raw), dim)
    vector = np.zeros(dim, dtype=np.float32)
    vector[:filled] = np.asarray(raw[:filled], dtype=np.float32)
    if len(raw) < dim:
        logger.debug(
            "Fallback descriptor (%d-d) shorter than configured dim (%d); "
            "zero-padding %d trailing dimensions.",
            len(raw), dim, dim - len(raw),
        )
    elif len(raw) > dim:
        logger.warning(
            "Fallback descriptor (%d-d) longer than configured dim (%d); "
            "truncating the trailing %d dimensions.",
            len(raw), dim, len(raw) - dim,
        )
    return vector, filled


def compute_fallback_embeddings(
    records: Sequence, dim: int
) -> Tuple[Dict[str, np.ndarray], int]:
    """Compute fallback vectors for every record.

    Returns ``({key: vector}, filled_dimension)``. ``filled_dimension`` is the
    number of leading components backed by real descriptors; the manifest
    records ``dim`` (padded size) and ``padded`` (whether trailing zeros were
    added).
    """
    vectors: Dict[str, np.ndarray] = {}
    filled = 0
    for rec in records:
        vector, filled = fallback_vector(rec.sequence_aa, dim)
        vectors[rec.key] = vector
    return vectors, filled
