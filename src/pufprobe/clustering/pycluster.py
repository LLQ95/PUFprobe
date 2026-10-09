"""Pure-numpy vector clustering (no mmseqs binary required).

Steps: L2-normalize every embedding, connect pairs whose cosine similarity is
at or above ``config.cluster.cosine_threshold``, then read out connected
components with union-find. Pairwise similarities are computed in row batches
so ten-thousand-sequence inputs never materialize a dense N*N matrix.

If ``igraph`` and ``leidenalg`` are both importable the thresholded graph is
additionally refined with the Leiden algorithm; this is purely optional and
never a hard dependency. The deterministic union-find partition is always
available and is what the tests exercise.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

from pufprobe.externals import has_module
from pufprobe.logging_utils import get_logger

logger = get_logger(__name__)


def cluster_vectors(
    vectors: Dict[str, np.ndarray],
    cosine_threshold: float,
    batch_size: int = 256,
    use_leiden: bool = True,
    sequences: Optional[Dict[str, str]] = None,
    kmer_gate: float = 0.0,
    kmer_size: int = 3,
) -> Dict[str, int]:
    """Return ``{global_key: cluster_index}`` from cosine-connected components.

    When ``sequences`` and a positive ``kmer_gate`` are supplied, an edge also
    requires Jaccard overlap of the proteins' k-mer sets at or above the gate.
    This stops coarse composition descriptors from merging proteins that have
    similar amino-acid counts but unrelated sequence order.
    """
    keys: List[str] = list(vectors.keys())
    n = len(keys)
    if n == 0:
        return {}
    matrix = np.vstack([np.asarray(vectors[k], dtype=np.float32) for k in keys])

    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    normalized = matrix / norms

    kmer_sets: Optional[List[set]] = None
    if sequences is not None and kmer_gate > 0:
        kmer_sets = [_kmer_set(sequences.get(k, ""), kmer_size) for k in keys]
    edges = _cosine_edges(normalized, cosine_threshold, batch_size, kmer_sets, kmer_gate)

    if use_leiden and has_module("igraph") and has_module("leidenalg"):
        partition = _leiden_partition(keys, normalized, edges, cosine_threshold)
        if partition is not None:
            return partition

    return _union_find(keys, edges)


def _kmer_set(sequence: str, k: int) -> set:
    seq = sequence.upper().replace("*", "")
    if len(seq) < k:
        return {seq} if seq else set()
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def _cosine_edges(
    normalized: np.ndarray,
    threshold: float,
    batch_size: int,
    kmer_sets: Optional[List[set]] = None,
    kmer_gate: float = 0.0,
) -> List[tuple]:
    """Collect (i, j) pairs with cosine >= threshold, i < j, in batches."""
    n = normalized.shape[0]
    edges: List[tuple] = []
    for start in range(0, n, batch_size):
        stop = min(start + batch_size, n)
        block = normalized[start:stop]  # (b, d)
        sim = block @ normalized.T  # (b, n)
        rows, cols = np.where(sim >= threshold)
        for r, c in zip(rows.tolist(), cols.tolist()):
            j_global = start + r
            if j_global < c:
                if kmer_sets is not None:
                    left, right = kmer_sets[j_global], kmer_sets[c]
                    union = left | right
                    overlap = len(left & right) / len(union) if union else 0.0
                    if overlap < kmer_gate:
                        continue
                edges.append((j_global, int(c)))
    return edges


def _union_find(keys: List[str], edges: List[tuple]) -> Dict[str, int]:
    parent = list(range(len(keys)))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i, j in edges:
        union(i, j)

    root_to_index: Dict[int, int] = {}
    result: Dict[str, int] = {}
    next_index = 0
    for idx, key in enumerate(keys):
        root = find(idx)
        if root not in root_to_index:
            root_to_index[root] = next_index
            next_index += 1
        result[key] = root_to_index[root]
    return result


def _leiden_partition(
    keys: List[str],
    normalized: np.ndarray,
    edges: List[tuple],
    resolution: float,
) -> Optional[Dict[str, int]]:
    """Optional Leiden refinement; returns None on any failure."""
    try:  # pragma: no cover - optional path, depends on host modules
        import igraph  # type: ignore
        import leidenalg  # type: ignore

        graph = igraph.Graph(n=len(keys), edges=edges, directed=False)
        graph.simplify()
        membership = leidenalg.find_partition(
            graph, leidenalg.RBConfigurationVertexPartition,
            resolution_parameter=resolution,
        ).membership
        return {key: int(membership[i]) for i, key in enumerate(keys)}
    except Exception as exc:  # pragma: no cover
        logger.warning("Leiden refinement failed (%s); using union-find clusters.", exc)
        return None
