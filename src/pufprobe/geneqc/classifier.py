"""True-gene quality control: score each candidate CDS as a likely false ORF.

The classifier is deliberately transparent: every verdict carries a list of
human-readable reasons in the form ``"check: value (threshold)"`` so a reviewer
can see exactly which signals drove the decision.  All thresholds are read from
the ``qc`` block of the configuration; nothing here hard-codes a cutoff.

A record receives a false-ORF score in ``[0, 1]`` that accumulates over the
independent signals below:

* coding length (with the permissive smORF path when ``include_smorf`` is set),
* low-complexity fraction of the amino-acid sequence,
* longest single-residue homopolymer run,
* coordinate overlap with other CDS on the same contig,
* and, when a nucleotide sequence is available, a robust GC-content z-score,
  GC skew and internal stop codons.

The ribosome-binding-site check requires an upstream intergenic sequence that isnormally absent from GFF3/GBK-only inputs; rather than faking it, the reason
list always records ``"RBS: not evaluated"``.
"""

from __future__ import annotations

import statistics
from typing import Any, Dict, List, Optional, Sequence

from pufprobe.features import (
    gc_content,
    gc_skew,
    has_internal_stop,
    longest_homopolymer,
    low_complexity_fraction,
)
from pufprobe.models import CDSRecord, GeneVerdict

__all__ = ["classify_all", "classify_record"]


def _overlap_fraction(record: CDSRecord, siblings: Sequence[CDSRecord]) -> float:
    """Largest fraction of ``record`` overlapped by another CDS on the same contig."""
    own_len = max(record.end - record.start + 1, 1)
    best = 0.0
    for other in siblings:
        if other is record or other.contig != record.contig:
            continue
        overlap = min(record.end, other.end) - max(record.start, other.start) + 1
        if overlap > 0:
            best = max(best, overlap / own_len)
    return best


def _gc_stats(records: Sequence[CDSRecord]) -> Optional[tuple[float, float]]:
    values = [gc_content(r.sequence_nt) for r in records if r.sequence_nt]
    if len(values) < 2:
        return None
    median = statistics.median(values)
    mad = statistics.median([abs(v - median) for v in values])
    return median, mad


def classify_record(
    record: CDSRecord,
    qc: Dict[str, Any],
    include_smorf: bool,
    gc_stats: Optional[tuple[float, float]],
    overlap_fraction: float,
) -> float:
    """Populate ``verdict``/``qc_score``/``qc_reasons`` for one record; return score."""
    seq = record.sequence_aa.upper().replace("*", "")
    length = len(seq)
    score = 0.0
    reasons: List[str] = []

    # -- length -----------------------------------------------------------
    smorf_length = int(qc["smorf_length_aa"])
    if include_smorf and length < smorf_length:
        min_len = int(qc["min_length_smorf_aa"])
        if length < min_len:
            score += 0.20
            reasons.append(f"length: {length} aa below smORF minimum (min {min_len})")
    else:
        min_len = int(qc["min_length_aa"])
        if length < min_len:
            score += 0.30
            reasons.append(f"length: {length} aa (min {min_len})")

    # -- low complexity ---------------------------------------------------
    lc = low_complexity_fraction(seq)
    max_lc = float(qc["max_low_complexity_fraction"])
    if lc > max_lc:
        score += 0.25
        reasons.append(f"low_complexity: {lc:.2f} (max {max_lc:.2f})")

    # -- homopolymer run --------------------------------------------------
    homo = longest_homopolymer(seq)
    max_homo = int(qc["max_homopolymer_run"])
    if homo > max_homo:
        score += 0.15
        reasons.append(f"homopolymer: {homo} aa (max {max_homo})")

    # -- coordinate overlap ------------------------------------------------
    max_overlap = float(qc["max_overlap_fraction"])
    if overlap_fraction > max_overlap:
        score += 0.20
        reasons.append(
            f"overlap: {overlap_fraction:.2f} (max {max_overlap:.2f})"
        )

    # -- nucleotide-level signals (only when a CDS translation exists) -----
    nt = record.sequence_nt
    if nt:
        if has_internal_stop(nt):
            score += 0.30
            reasons.append("internal_stop: present")
        if gc_stats is not None:
            median, mad = gc_stats
            gc = gc_content(nt)
            z = 0.0 if mad == 0 else 0.6745 * (gc - median) / mad
            if abs(z) > float(qc["gc_zscore_cutoff"]):
                score += 0.15
                reasons.append(f"gc_zscore: {z:+.2f} (cutoff {qc['gc_zscore_cutoff']})")
        skew = gc_skew(nt)
        if abs(skew) > float(qc["max_gc_skew_abs"]):
            score += 0.10
            reasons.append(f"gc_skew: {skew:+.2f} (max {qc['max_gc_skew_abs']:.2f})")
    else:
        reasons.append("gc: not evaluated (no nucleotide sequence)")

    # -- RBS requires an upstream intergenic window we usually do not have -
    reasons.append("RBS: not evaluated")

    score = round(min(1.0, score), 4)
    if score >= float(qc["false_orf_score"]):
        verdict = GeneVerdict.LIKELY_FALSE_ORF
    elif score >= float(qc["borderline_score"]):
        verdict = GeneVerdict.BORDERLINE
    else:
        verdict = GeneVerdict.TRUE_PROTEIN_UNKNOWN

    record.verdict = verdict
    record.qc_score = score
    record.qc_reasons = reasons
    return score


def classify_all(
    records: List[CDSRecord],
    config: Dict[str, Any],
    include_smorf: bool = False,
) -> Dict[str, Any]:
    """Classify every record in place and return a small summary dictionary."""
    qc = config["qc"]
    gc_stats = _gc_stats(records)
    overlap = {r.key: _overlap_fraction(r, records) for r in records}

    counts = {v.value: 0 for v in GeneVerdict}
    for record in records:
        classify_record(record, qc, include_smorf, gc_stats, overlap[record.key])
        counts[record.verdict.value] += 1  # type: ignore[union-attr]

    return {"total": len(records), "by_verdict": counts}
