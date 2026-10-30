#!/usr/bin/env python
"""Verify whether whitening the pLM embeddings restores family separation.

The residue-mean ESM protein embeddings are anisotropic: the cosine between
unrelated globular proteins is compressed near 0.93, so an absolute cosine
threshold and a top-k neighbour list lose functional specificity. This script
estimates the embedding mean and covariance on the reference library, applies
BERT-whitening (z = (x - mu) @ W with W = U S^-1/2 U^T), and reports, for a
few query proteins, the cosine distribution and the separation between
same-Pfam and other proteins before and after the transform.

The transform is diagnostic only; no files are overwritten. If whitening
spreads the cosine distribution and lifts same-Pfam members clearly above
other proteins, the same transform can be applied when building the reference
library and the query embeddings.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pufprobe.logging_utils import setup_logging, get_logger  # noqa: E402
from pufprobe.search.plm_neighbor import (  # noqa: E402
    load_reference,
    load_query_vectors,
)

logger = get_logger("verify_whitening")


def whitening_matrix(matrix: np.ndarray):
    """Return (mean, W) where W whitens the centered embeddings."""
    mean = matrix.mean(axis=0)
    centered = matrix - mean
    # The covariance is a 1280x1280 matrix regardless of the library size.
    cov = (centered.T @ centered) / max(1, matrix.shape[0] - 1)
    eigvals, eigvecs = np.linalg.eigh(cov)
    # Clamp tiny eigenvalues for numerical stability, then build S^-1/2.
    eigvals = np.maximum(eigvals, 1e-8)
    inv_sqrt = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T
    return mean, inv_sqrt.astype(np.float32)


def l2_normalize(matrix: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    return matrix / np.maximum(norms, 1e-8)


def report(lib_vectors: np.ndarray, query: np.ndarray, lib, label: str) -> None:
    sims = np.asarray(lib_vectors @ query, dtype=np.float64)
    print("[%s] cosine p50=%.4f p90=%.4f p99=%.4f" % (
        label, np.percentile(sims, 50), np.percentile(sims, 90),
        np.percentile(sims, 99)))
    same_idx = [i for i, acc in enumerate(lib.ids)
                if any(p in (lib.terms.get(acc) or {}).get("Pfam", [])
                       for p in _query_pfams)]
    mask = np.ones(len(lib.ids), bool)
    mask[same_idx] = False
    if same_idx:
        print("    same-Pfam median=%.4f | other median=%.4f other p99=%.4f" % (
            np.median(sims[same_idx]), np.median(sims[mask]),
            np.percentile(sims[mask], 99)))


def main() -> int:
    setup_logging()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-dir", type=Path, required=True)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--queries", type=str, nargs="+", required=True)
    args = parser.parse_args()

    lib = load_reference(args.ref_dir)
    raw = np.asarray(lib.vectors, dtype=np.float32)
    mean, whiten = whitening_matrix(raw)
    white = l2_normalize((raw - mean) @ whiten)

    global _query_pfams
    for key in args.queries:
        accession = key.split("@")[-1]
        _query_pfams = (lib.terms.get(accession) or {}).get("Pfam", [])
        query = load_query_vectors(args.embeddings, [key])[key]
        print("=== %s (Pfam %s) ===" % (key, ",".join(_query_pfams) or "-"))
        report(raw, query, lib, "raw")
        query_white = l2_normalize(((query - mean) @ whiten).reshape(1, -1))[0]
        report(white, query_white, lib, "whitened")
        print()
    return 0


_query_pfams: list = []


if __name__ == "__main__":
    sys.exit(main())
