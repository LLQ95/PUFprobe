"""Cross-strain micro-synteny and core/accessory classification.

For every PUF family we look at the neighbour annotations of its members across
strains and measure how often a given neighbour is conserved.  A family present
in every sample is ``core``; one restricted to a single sample is
``clonal_specific``; everything else is ``accessory``.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any, Dict, List, Set, Tuple

from pufprobe.context.operons import candidate_pufs
from pufprobe.models import CDSRecord

__all__ = ["build_microsynteny"]


def build_microsynteny(
    records: List[CDSRecord],
    neighbour_rows: List[Dict[str, Any]],
    family_samples: Dict[str, Set[str]],
    samples: List[str],
) -> List[Dict[str, Any]]:
    n_samples = len(samples)

    # Collect the set of neighbour products seen next to each family member.
    family_neighbours: Dict[str, Counter] = defaultdict(Counter)
    for row in neighbour_rows:
        family_id = row["family_id"]
        if family_id:
            family_neighbours[family_id][row["neighbour_product"].lower().strip()] += 1

    puf_families = {r.family_id for r in candidate_pufs(records) if r.family_id}
    rows: List[Dict[str, Any]] = []
    for family_id in sorted(puf_families):
        present = family_samples.get(family_id, set())
        width = len(present)
        if width == n_samples and n_samples > 0:
            classification = "core"
        elif width <= 1:
            classification = "clonal_specific"
        else:
            classification = "accessory"

        counter = family_neighbours.get(family_id, Counter())
        if counter:
            conserved_neighbour, total = counter.most_common(1)[0]
            conservation = round(total / max(1, sum(counter.values())), 4)
        else:
            conserved_neighbour, conservation = "", 0.0

        rows.append(
            {
                "family_id": family_id,
                "lineage_width": width,
                "n_samples": n_samples,
                "classification": classification,
                "conserved_neighbour": conserved_neighbour,
                "conservation_fraction": conservation,
                "samples": sorted(present),
            }
        )
    return rows
