"""Wet-lab priority scoring for families that received no anchored term.

The priority is a transparent weighted sum of five interpretable components
(config ``rank.priority_weights``):

* structural novelty  — ``novel_score`` when the representative was folded,
  neutral 0.5 when no usable structure exists;
* pocket evidence     — best score from POCKET-channel evidence;
* clonal specificity  — rare families (few samples) score high;
* AMR/virulence neighbourhood — context evidence hitting AMR/virulence keywords;
* small protein       — representatives shorter than 50 aa.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

from pufprobe.logging_utils import get_logger
from pufprobe.models import Evidence, EvidenceChannel, Family
from pufprobe.rank.fuse import AMR_VIRULENCE_KEYWORDS

log = get_logger("rank.priority")

_SMALL_PROTEIN_AA = 50


def _structural_novelty(family: Family) -> Tuple[float, str]:
    if family.fold_status == "folded":
        if family.novel_score is not None:
            return float(family.novel_score), f"structural novelty {family.novel_score:.2f}"
        return 0.5, "folded but no novelty score; neutral 0.5"
    return 0.5, f"no usable structure (fold_status={family.fold_status}); neutral 0.5"


def _pocket_score(rep_evidence: Sequence[Evidence]) -> Tuple[float, Optional[str]]:
    pockets = [ev for ev in rep_evidence if ev.channel == EvidenceChannel.POCKET]
    if not pockets:
        return 0.0, None
    best = max(pockets, key=lambda e: e.score)
    return float(best.score), f"predicted pocket (score {best.score:.2f}, {best.target_label})"


def _clonal_specificity(family: Family, total_samples: int) -> Tuple[float, str]:
    metrics = family.metrics or {}
    n_samples = metrics.get("sample_count", metrics.get("n_samples"))
    if n_samples is None:
        n_samples = metrics.get("lineage_width", family.lineage_width)
    if total_samples <= 0 or n_samples is None:
        return 0.5, "sample distribution unknown; neutral 0.5"
    prevalence = min(max(n_samples / total_samples, 0.0), 1.0)
    score = 1.0 - prevalence
    return score, (
        f"present in {n_samples}/{total_samples} samples "
        f"({'clonal/accessory' if prevalence <= 0.25 else 'broadly distributed'})"
    )


def _amr_virulence(rep_evidence: Sequence[Evidence]) -> Tuple[float, Optional[str]]:
    hits: List[Evidence] = []
    for ev in rep_evidence:
        haystack = f"{ev.target_label} {ev.target_id} {' '.join(str(v) for v in (ev.details or {}).values())}".lower()
        if any(kw.lower() in haystack for kw in AMR_VIRULENCE_KEYWORDS):
            hits.append(ev)
    if not hits:
        return 0.0, None
    best = max(hits, key=lambda e: e.score)
    return float(best.score), f"AMR/virulence neighbour hit: {best.target_label} ({best.score:.2f})"


def _small_protein(family: Family) -> Tuple[float, Optional[str]]:
    length = len(family.representative_seq or "")
    if length and length < _SMALL_PROTEIN_AA:
        return 1.0, f"small protein ({length} aa < {_SMALL_PROTEIN_AA})"
    return 0.0, None


def compute_family_priority(
    family: Family,
    rep_evidence: Sequence[Evidence],
    total_samples: int,
    config: Dict[str, Any],
) -> Tuple[float, List[str]]:
    """Return ``(priority_score, reasons)`` for one unanchored family."""
    weights: Dict[str, float] = config["rank"]["priority_weights"]
    components = [
        ("structural_novelty", *_structural_novelty(family)),
        ("pocket", *_pocket_score(rep_evidence)),
        ("clonal_specificity", *_clonal_specificity(family, total_samples)),
        ("amr_virulence_neighbourhood", *_amr_virulence(rep_evidence)),
        ("small_protein", *_small_protein(family)),
    ]
    score = 0.0
    reasons: List[str] = []
    for name, value, reason in components:
        score += weights.get(name, 0.0) * value
        if reason:
            reasons.append(f"{name}: {reason}")
    return min(max(score, 0.0), 1.0), reasons
