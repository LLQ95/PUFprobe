"""Protein-language-model nearest-neighbour evidence channel.

This channel transfers functional annotation from an *annotated* reference
library to each representative protein by nearest-neighbour search in protein
language model (pLM) embedding space. The query embeddings produced in the
embed stage (``02_embed/embeddings.npz``) are compared against a library of
reference proteins that carry GO, EC and Pfam annotations; a reference library
is built offline by ``server/build_plm_reference.py`` from Swiss-Prot with the
same ESM-2 checkpoint.

The query and reference embeddings MUST come from the same checkpoint and the
same residue-mean pooling, otherwise distances are not comparable. The channel
verifies this from the two manifest files and disables itself on a mismatch.

For each representative the cosine similarity to every reference vector is
computed, the top-k neighbours are kept, and neighbour weights are derived as
a linear ramp of cosine above ``cos_floor``. Each candidate term is scored by a
noisy-OR over the weights of neighbours that carry that term, so agreement
across several neighbours accumulates while a single very close annotated
neighbour can still be decisive. The resulting evidence is recorded on
``EvidenceChannel.PLM_NEIGHBOR`` (the independent "sequence" evidence group).

Reference library layout (``<ref_dir>``, default ``<dbdir>/plm_ref``):

    plm_ref_vectors.npy    (N, D) float16, L2-normalized reference vectors
    plm_ref_ids.json       [accession, ...] in matrix row order
    plm_ref_terms.json     {accession: {"GO": [...], "EC": [...], "Pfam": [...]}}
    plm_ref_manifest.json  {"model", "dim", "n", "normalized", ...}
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from pufprobe.io_utils import read_fasta
from pufprobe.logging_utils import get_logger
from pufprobe.models import Evidence, EvidenceChannel

logger = get_logger(__name__)

REF_VECTORS = "plm_ref_vectors.npy"
REF_IDS = "plm_ref_ids.json"
REF_TERMS = "plm_ref_terms.json"
REF_MANIFEST = "plm_ref_manifest.json"

DEFAULT_TOP_K = 20
DEFAULT_COS_FLOOR = 0.50
DEFAULT_COS_MIN = 0.50
DEFAULT_MIN_TERM_SCORE = 0.50
DEFAULT_MIN_VOTERS = 1
_QUERY_BATCH = 64

# Term lists are emitted in this order.
_TERM_KINDS = ("GO", "EC", "Pfam")


@dataclass
class ReferenceLibrary:
    """In-memory reference library used for the nearest-neighbour search."""

    vectors: np.ndarray  # (N, D) float32, L2-normalized
    ids: List[str]
    terms: Dict[str, Dict[str, List[str]]]
    manifest: Dict[str, Any]

    @property
    def model(self) -> str:
        return str(self.manifest.get("model", ""))


def _l2_normalize(vectors: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    if vectors.ndim == 1:
        norm = float(np.linalg.norm(vectors))
        return vectors / (norm if norm > eps else 1.0)
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    norms = np.maximum(norms, eps)
    return vectors / norms


def load_reference(ref_dir: Optional[Path]) -> Optional[ReferenceLibrary]:
    """Load a reference library; return None when the directory is incomplete."""
    if ref_dir is None:
        return None
    vec_path = ref_dir / REF_VECTORS
    ids_path = ref_dir / REF_IDS
    terms_path = ref_dir / REF_TERMS
    manifest_path = ref_dir / REF_MANIFEST
    missing = [p.name for p in (vec_path, ids_path, terms_path, manifest_path) if not p.exists()]
    if missing:
        logger.info(
            "pLM reference library incomplete in %s (missing %s); skipping plm_neighbor channel.",
            ref_dir, ", ".join(missing),
        )
        return None

    vectors = np.load(vec_path, mmap_mode="r").astype(np.float32, copy=False)
    with ids_path.open("r", encoding="utf-8") as handle:
        ids = json.load(handle)
    with terms_path.open("r", encoding="utf-8") as handle:
        terms = json.load(handle)
    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    if len(ids) != vectors.shape[0]:
        raise ValueError(
            f"pLM reference library mismatch: {len(ids)} ids but "
            f"{vectors.shape[0]} vectors in {ref_dir}"
        )
    vectors = _l2_normalize(np.asarray(vectors, dtype=np.float32))
    logger.info(
        "pLM reference library loaded: %d proteins x %d dims (model=%s).",
        vectors.shape[0], vectors.shape[1], manifest.get("model", "?"),
    )
    return ReferenceLibrary(vectors=vectors, ids=ids, terms=terms, manifest=manifest)


def load_query_vectors(
    embeddings_path: Path, keys: List[str]
) -> Dict[str, np.ndarray]:
    """Read representative embeddings from the embed stage and normalise them."""
    if not embeddings_path.exists():
        raise FileNotFoundError(
            f"Missing {embeddings_path}; run `pufprobe embed` before plm_neighbor."
        )
    out: Dict[str, np.ndarray] = {}
    with np.load(embeddings_path, allow_pickle=True) as data:
        for key in keys:
            if key not in data.files:
                continue
            vector = np.asarray(data[key], dtype=np.float32)
            norm = float(np.linalg.norm(vector))
            if norm > 1e-8:
                out[key] = vector / norm
    return out


def _batch_topk(
    ref: np.ndarray,
    queries: np.ndarray,
    query_keys: List[str],
    k: int,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Return ``{query_key: (neighbour_indices, cosine)}`` for the top-k."""
    n_ref = ref.shape[0]
    k = min(k, n_ref)
    results: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for start in range(0, len(query_keys), _QUERY_BATCH):
        batch = queries[start : start + _QUERY_BATCH]
        sims = ref @ batch.T  # (N, batch)
        for j in range(batch.shape[0]):
            col = sims[:, j]
            if n_ref <= k:
                order = np.argsort(-col)
            else:
                pooled = np.argpartition(col, n_ref - k)[-k:]
                order = pooled[np.argsort(-col[pooled])]
            results[query_keys[start + j]] = (order, col[order])
        del sims
    return results


def _cosine_weight(cosine: float, cos_floor: float) -> float:
    if cosine <= cos_floor:
        return 0.0
    return float(min(1.0, max(0.0, (cosine - cos_floor) / (1.0 - cos_floor))))


def _normalise_term(kind: str, raw: str) -> str:
    if kind == "EC" and not raw.startswith("EC:"):
        return f"EC:{raw}"
    return raw


def vote_terms(
    neighbour_idx: np.ndarray,
    neighbour_cos: np.ndarray,
    ref_ids: List[str],
    ref_terms: Dict[str, Dict[str, List[str]]],
    cos_floor: float = DEFAULT_COS_FLOOR,
    cos_min: float = DEFAULT_COS_MIN,
    min_term_score: float = DEFAULT_MIN_TERM_SCORE,
    min_voters: int = DEFAULT_MIN_VOTERS,
    self_accession: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Aggregate neighbour annotations into scored term candidates.

    Terms are scored by noisy-OR over the ramped cosine weights of neighbours
    that carry them. Only neighbours at or above ``cos_min`` count, and a term
    needs at least ``min_voters`` supporting neighbours and a score at or above
    ``min_term_score``. When ``self_accession`` is given, the reference entry
    with that accession is skipped so a query that is itself in the reference
    library cannot transfer its own annotation.
    """
    # Per-term accumulation of (1 - weight) products plus voter records.
    failure: Dict[Tuple[str, str], float] = {}
    voters: Dict[Tuple[str, str], List[Tuple[str, float]]] = {}

    for idx, cos in zip(neighbour_idx, neighbour_cos):
        cos_f = float(cos)
        if cos_f < cos_min:
            continue
        weight = _cosine_weight(cos_f, cos_floor)
        if weight <= 0.0:
            continue
        accession = ref_ids[int(idx)]
        if self_accession is not None and accession == self_accession:
            continue
        entry = ref_terms.get(accession)
        if not entry:
            continue
        for kind in _TERM_KINDS:
            for raw_term in entry.get(kind, []):
                term = _normalise_term(kind, raw_term)
                slot = (kind, term)
                failure[slot] = failure.get(slot, 1.0) * (1.0 - weight)
                voters.setdefault(slot, []).append((accession, cos_f))

    candidates: List[Dict[str, Any]] = []
    for (kind, term), prod in failure.items():
        term_voters = voters[(kind, term)]
        n_voters = len(term_voters)
        if n_voters < min_voters:
            continue
        score = 1.0 - prod
        if score < min_term_score:
            continue
        term_voters.sort(key=lambda item: -item[1])
        candidates.append({
            "kind": kind.lower(),
            "term_id": term,
            "score": float(min(1.0, score)),
            "n_voters": n_voters,
            "top_cos": term_voters[0][1],
            "neighbours": term_voters[:5],
        })
    candidates.sort(key=lambda c: (-c["score"], c["term_id"]))
    return candidates


def run_plm_neighbor_search(
    reps_faa: Path,
    embeddings_path: Path,
    dbdir: Optional[Path],
    out_dir: Path,
    ref_dir: Optional[Path] = None,
    query_manifest_path: Optional[Path] = None,
    top_k: int = DEFAULT_TOP_K,
    cos_floor: float = DEFAULT_COS_FLOOR,
    cos_min: float = DEFAULT_COS_MIN,
    min_term_score: float = DEFAULT_MIN_TERM_SCORE,
    min_voters: int = DEFAULT_MIN_VOTERS,
) -> Tuple[Dict[str, List[Evidence]], Optional[Path]]:
    """Find annotated pLM neighbours for each representative.

    Returns a map ``global key -> [Evidence]`` and the raw TSV path, following
    the same contract as the other search channels. A missing library or a
    model mismatch disables the channel and returns an empty map with None.
    """
    ref_dir = ref_dir or (dbdir / "plm_ref" if dbdir is not None else None)
    library = load_reference(ref_dir)
    if library is None:
        return {}, None

    # The query and reference embeddings must share a checkpoint.
    if query_manifest_path is not None and query_manifest_path.exists():
        with query_manifest_path.open("r", encoding="utf-8") as handle:
            query_manifest = json.load(handle)
        query_model = str(query_manifest.get("model", ""))
        if query_model and library.model and query_model != library.model:
            logger.warning(
                "pLM model mismatch: query embeddings use '%s' but the reference "
                "library uses '%s'; rebuild the library with the same checkpoint "
                "or re-embed. Skipping plm_neighbor channel.",
                query_model, library.model,
            )
            return {}, None

    rep_keys = [header.split()[0] for header, _ in read_fasta(reps_faa)]
    query_vectors = load_query_vectors(embeddings_path, rep_keys)
    missing = [key for key in rep_keys if key not in query_vectors]
    if missing:
        logger.info("pLM neighbours: %d representatives lack embeddings.", len(missing))
    keys = [key for key in rep_keys if key in query_vectors]
    if not keys:
        logger.info("pLM neighbours: no representative embeddings available; nothing to search.")
        return {}, None

    query_matrix = np.stack([query_vectors[key] for key in keys])
    top_hits = _batch_topk(library.vectors, query_matrix, keys, top_k)

    evidence: Dict[str, List[Evidence]] = {}
    rows = ["query\tkind\ttarget\tscore\tn_voters\ttop_cos\tneighbours"]
    n_terms = 0
    for key in keys:
        idx, cos = top_hits[key]
        # A global key is "<sample>@<accession>". When the accession is itself
        # a reference entry (curated training proteins), exclude it so the
        # query cannot label itself; ordinary genome loci do not match.
        self_accession = key.split("@", 1)[1] if "@" in key else None
        candidates = vote_terms(
            idx, cos, library.ids, library.terms,
            cos_floor=cos_floor, cos_min=cos_min,
            min_term_score=min_term_score, min_voters=min_voters,
            self_accession=self_accession,
        )
        for cand in candidates:
            neighbour_text = ",".join(f"{acc}:{c:.3f}" for acc, c in cand["neighbours"])
            rows.append(
                f"{key}\t{cand['kind']}\t{cand['term_id']}\t{cand['score']:.3f}\t"
                f"{cand['n_voters']}\t{cand['top_cos']:.3f}\t{neighbour_text}"
            )
            evidence.setdefault(key, []).append(Evidence(
                channel=EvidenceChannel.PLM_NEIGHBOR,
                target_id=cand["term_id"],
                target_label=cand["term_id"],
                score=cand["score"],
                details={
                    "db": cand["kind"],
                    "n_voters": cand["n_voters"],
                    "n_neighbors": int(top_k),
                    "top_cos": round(cand["top_cos"], 4),
                    "neighbours": [
                        {"accession": acc, "cos": round(c, 4)}
                        for acc, c in cand["neighbours"]
                    ],
                },
            ))
            n_terms += 1

    raw_tsv = out_dir / "plm_neighbor_raw.tsv"
    raw_tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")
    logger.info(
        "plm_neighbor transferred %d term annotations to %d / %d representatives "
        "(top_k=%d, cos_floor=%.2f).",
        n_terms, len(evidence), len(rep_keys), top_k, cos_floor,
    )
    return evidence, raw_tsv
