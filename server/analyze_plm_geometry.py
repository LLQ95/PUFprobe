#!/usr/bin/env python
"""Quantify how functional the raw pLM nearest-neighbour geometry is.

The plm_neighbor channel transfers annotations from the closest proteins in
ESM-2 residue-mean embedding space. Residue-mean protein vectors are known to
be anisotropic: most globular proteins sit at high cosine regardless of
function, so an absolute cosine threshold may not separate homologues from
unrelated proteins. This script measures the symptom across every hold-out
representative without rerunning the embeddings.

For each query it takes the top-k neighbours, excludes the query's own
reference entry, and reports how often the neighbours share a Pfam family with
the query, both for the single closest neighbour and for any of the top-k. It
also summarises the cosine distribution so the overlap between same-family
and different-family proteins can be judged directly.

Example
-------
::

    python -u analyze_plm_geometry.py \
        --reps ltr_work/03_cluster/representatives.faa \
        --embeddings ltr_work/02_embed/embeddings.npz \
        --ref-dir dbs/dbdir/plm_ref \
        --top-k 20
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pufprobe.io_utils import read_fasta  # noqa: E402
from pufprobe.logging_utils import get_logger, setup_logging  # noqa: E402
from pufprobe.search.plm_neighbor import (  # noqa: E402
    _batch_topk,
    load_query_vectors,
    load_reference,
)

logger = get_logger("analyze_plm_geometry")


def _pfam(lib: Any, index: int) -> set:
    return set((lib.terms.get(lib.ids[int(index)]) or {}).get("Pfam", []))


def main() -> int:
    setup_logging()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=Path, required=True)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--ref-dir", type=Path, required=True)
    parser.add_argument("--top-k", type=int, default=20)
    args = parser.parse_args()

    lib = load_reference(args.ref_dir)
    rep_keys = [header.split()[0] for header, _ in read_fasta(args.reps)]
    query_vectors = load_query_vectors(args.embeddings, rep_keys)
    keys = [key for key in rep_keys if key in query_vectors]
    query_matrix = np.stack([query_vectors[key] for key in keys])
    top_hits = _batch_topk(lib.vectors, query_matrix, keys, args.top_k)

    n = 0
    top1_match = 0
    any_match = 0
    self_in_topk = 0
    same_cos: List[float] = []
    diff_cos: List[float] = []

    for key in keys:
        self_accession = key.split("@", 1)[1] if "@" in key else None
        idx, cos = top_hits[key]
        # Locate the query's own Pfam set from its reference entry, then drop
        # the self neighbour so it cannot inflate the agreement.
        query_pfam: set = set()
        kept_idx: List[int] = []
        kept_cos: List[float] = []
        for j, i in enumerate(idx):
            acc = lib.ids[int(i)]
            if acc == self_accession:
                self_in_topk += 1
                query_pfam = _pfam(lib, int(i))
                continue
            kept_idx.append(int(i))
            kept_cos.append(float(cos[j]))
        if not kept_idx or not query_pfam:
            continue
        n += 1
        neighbour_pfam = _pfam(lib, kept_idx[0])
        if neighbour_pfam & query_pfam:
            top1_match += 1
        if any(_pfam(lib, i) & query_pfam for i in kept_idx):
            any_match += 1
        for i, c in zip(kept_idx, kept_cos):
            if _pfam(lib, i) & query_pfam:
                same_cos.append(c)
            else:
                diff_cos.append(c)

    summary: Dict[str, Any] = {
        "queries": n,
        "top1_pfam_match": top1_match,
        "top1_pfam_match_rate": round(top1_match / n, 4) if n else 0.0,
        "any_topk_pfam_match": any_match,
        "any_topk_pfam_match_rate": round(any_match / n, 4) if n else 0.0,
        "self_present_in_topk": self_in_topk,
        "same_family_neighbour_cos_median": round(float(np.median(same_cos)), 4) if same_cos else None,
        "different_family_neighbour_cos_median": round(float(np.median(diff_cos)), 4) if diff_cos else None,
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
