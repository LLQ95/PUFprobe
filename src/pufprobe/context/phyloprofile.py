"""Phylogenetic presence/absence profiles and their correlations.

Each family gets a binary vector across samples.  PUF families are correlated
against the presence/absence profiles of annotated functional genes, with
emphasis on genes carrying AMR / mobilisation / virulence keywords.  Similarity
is reported as a Jaccard index (the phi coefficient is undefined on constant
two-sample vectors, so Jaccard is used as the robust default).  A family must
occur in at least ``min_genomes_for_profile`` samples before it is scored.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Set, Tuple

from pufprobe.models import CDSRecord, GeneStatus
from pufprobe.context.operons import candidate_pufs

__all__ = ["build_profiles", "AMR_KEYWORDS", "jaccard"]

AMR_KEYWORDS = (
    "resistance",
    "efflux",
    "transposase",
    "integrase",
    "toxin",
    "virulence",
)


def jaccard(a: Set[str], b: Set[str]) -> float:
    union = len(a | b)
    if union == 0:
        return 0.0
    return len(a & b) / union


def build_profiles(
    records: List[CDSRecord], config: Dict[str, Any]
) -> Tuple[List[Dict[str, Any]], Dict[str, Set[str]], List[str]]:
    """Return (correlation rows, family->samples, ordered sample list)."""
    min_genomes = int(config["context"]["min_genomes_for_profile"])
    samples = sorted({r.sample_id for r in records})

    family_samples: Dict[str, Set[str]] = defaultdict(set)
    annotated_samples: Dict[str, Set[str]] = defaultdict(set)
    for record in records:
        if record.family_id:
            family_samples[record.family_id].add(record.sample_id)
        if record.status == GeneStatus.ANNOTATED:
            key = record.product.lower().strip()
            annotated_samples[key].add(record.sample_id)

    puf_families = {r.family_id for r in candidate_pufs(records) if r.family_id}

    rows: List[Dict[str, Any]] = []
    for family_id in sorted(puf_families):
        fset = family_samples[family_id]
        if len(fset) < min_genomes:
            continue
        vector = "".join("1" if s in fset else "0" for s in samples)
        for product, aset in sorted(annotated_samples.items()):
            score = jaccard(fset, aset)
            emphasis = any(keyword in product for keyword in AMR_KEYWORDS)
            rows.append(
                {
                    "family_id": family_id,
                    "family_profile": vector,
                    "annotated_gene": product,
                    "annotated_profile": "".join("1" if s in aset else "0" for s in samples),
                    "jaccard": round(score, 4),
                    "emphasis": emphasis,
                }
            )
    return rows, dict(family_samples), samples
