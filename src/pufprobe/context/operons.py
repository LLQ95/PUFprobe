"""Genomic neighbourhood / operon-coupling signal.

Records are grouped by sample and contig and ordered by start coordinate.  For
each candidate PUF the ``neighbour_window`` genes on either side are collected;
same-strand genes separated by fewer than ``operon_intergenic_bp`` base pairs
are treated as operon-coupled neighbours.
"""

from __future__ import annotations

from typing import Any, Dict, List

from pufprobe.models import CDSRecord, GeneStatus, GeneVerdict

__all__ = ["candidate_pufs", "build_neighbours"]


def candidate_pufs(records: List[CDSRecord]) -> List[CDSRecord]:
    """Hypothetical / conserved-hypothetical records that passed QC."""
    return [
        r
        for r in records
        if r.status in (GeneStatus.HYPOTHETICAL, GeneStatus.CONSERVED_HYPOTHETICAL)
        and r.verdict != GeneVerdict.LIKELY_FALSE_ORF
    ]


def _intergenic(a: CDSRecord, b: CDSRecord) -> int:
    gap = max(a.start, b.start) - min(a.end, b.end) - 1
    return max(0, gap)


def build_neighbours(records: List[CDSRecord], config: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Return one row per (candidate, neighbour) pair."""
    ctx = config["context"]
    window = int(ctx["neighbour_window"])
    intergenic_bp = int(ctx["operon_intergenic_bp"])
    same_strand_only = bool(ctx.get("same_strand_operon", True))

    ordered: Dict[tuple, List[CDSRecord]] = {}
    for record in records:
        ordered.setdefault((record.sample_id, record.contig), []).append(record)
    for rows in ordered.values():
        rows.sort(key=lambda r: r.start)

    rows: List[Dict[str, Any]] = []
    for candidate in candidate_pufs(records):
        contig_rows = ordered[(candidate.sample_id, candidate.contig)]
        try:
            idx = contig_rows.index(candidate)
        except ValueError:
            continue
        for j in range(max(0, idx - window), min(len(contig_rows), idx + window + 1)):
            if j == idx:
                continue
            nb = contig_rows[j]
            gap = _intergenic(candidate, nb)
            strand_match = candidate.strand == nb.strand
            coupled = (strand_match or not same_strand_only) and gap < intergenic_bp
            rows.append(
                {
                    "sample_id": candidate.sample_id,
                    "locus_tag": candidate.locus_tag,
                    "family_id": candidate.family_id or "",
                    "neighbour_locus": nb.locus_tag,
                    "neighbour_product": nb.product,
                    "neighbour_status": nb.status.value,
                    "position": "upstream" if j < idx else "downstream",
                    "intergenic_bp": gap,
                    "strand_match": strand_match,
                    "coupled": coupled,
                }
            )
    return rows
