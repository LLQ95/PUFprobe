"""Whitening transform for protein-language-model embeddings.

Residue-mean protein embeddings are anisotropic: unrelated globular proteins
end up with a high cosine similarity, which compresses the geometry and makes
an absolute cosine threshold uninformative. BERT-whitening removes the mean
and decorrelates the embedding dimensions so the resulting vectors have zero
mean and identity covariance, spreading the cosine distribution around zero.

The transform is estimated once on the reference library and stored with it.
The same mean and matrix are applied to query embeddings at search time,
keeping the query and reference spaces aligned.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

# Eigenvalues below this floor are treated as zero directions, preventing
# numerical blow-up when inverting the covariance.
_EIGEN_FLOOR = 1e-8


def whitening_matrix(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return ``(mean, W)`` where W whitens the centered embeddings.

    Parameters
    ----------
    matrix:
        Reference embeddings, shape (n, d). The same input used to estimate
        the transform should be L2-normalized consistently for queries and
        references.
    """
    ref = np.asarray(matrix, dtype=np.float64)
    mean = ref.mean(axis=0)
    centered = ref - mean
    covariance = (centered.T @ centered) / max(1, ref.shape[0] - 1)
    eigvals, eigvecs = np.linalg.eigh(covariance)
    eigvals = np.maximum(eigvals, _EIGEN_FLOOR)
    inv_sqrt = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T
    return mean.astype(np.float32), inv_sqrt.astype(np.float32)


def l2_normalize(matrix: np.ndarray) -> np.ndarray:
    """L2-normalize rows of a matrix."""
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    return np.asarray(matrix, dtype=np.float32) / np.maximum(norms, 1e-8)


def apply_whitening(
    matrix: np.ndarray, mean: np.ndarray, whiten: np.ndarray,
) -> np.ndarray:
    """Whiten one or more vectors and return L2-normalized results."""
    arr = np.asarray(matrix, dtype=np.float32)
    single = arr.ndim == 1
    if single:
        arr = arr.reshape(1, -1)
    transformed = (arr - mean) @ whiten
    normalized = l2_normalize(transformed)
    return normalized[0] if single else normalized
