"""context stage pipeline.

Pure-Python aggregation of operon, phylogenetic-profile and micro-synteny
signals into unified Evidence records and tabular outputs.  See the frozen
contract docstring for the input/output specification.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from pufprobe.io_utils import read_fasta, read_json, write_json
from pufprobe.logging_utils import get_logger
from pufprobe.models import (
    CDSRecord,
    Evidence,
    EvidenceChannel,
    GeneStatus,
    GeneVerdict,
    to_dict,
)
from pufprobe.ontology import parse_dbxrefs
from pufprobe.workdir import WorkDir
from pufprobe.context import microsynteny, operons, phyloprofile

log = get_logger("context")


def _records_from_dicts(rows: List[dict]) -> List[CDSRecord]:
    out: List[CDSRecord] = []
    for row in rows:
        data = dict(row)
        data["status"] = GeneStatus(data["status"])
        if data.get("verdict") is not None:
            data["verdict"] = GeneVerdict(data["verdict"])
        out.append(CDSRecord(**data))
    return out


def _coerce_enums(records: List[CDSRecord]) -> None:
    for record in records:
        if isinstance(record.status, str):
            record.status = GeneStatus(record.status)
        if record.verdict is not None and isinstance(record.verdict, str):
            record.verdict = GeneVerdict(record.verdict)


def _assign_families(records: List[CDSRecord]) -> None:
    """Greedy sequence-identity family assignment when clustering has not run."""
    if records and all(r.family_id for r in records):
        return
    reps: Dict[str, str] = {}
    counter = 0
    for record in records:
        seq = record.sequence_aa.upper().replace("*", "")
        assigned = None
        for family_id, rep in reps.items():
            if seq == rep:
                assigned = family_id
                break
        if assigned is None:
            counter += 1
            assigned = f"FAM_{counter:04d}"
            reps[assigned] = seq
        record.family_id = assigned


def _write_tsv(path: Path, header: List[str], rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(col, "")) for col in header) + "\n")


def _reference_terms(path: Path) -> Set[str]:
    terms: Set[str] = set()
    for header, _ in read_fasta(path):
        for token in re.split(r"[^A-Za-z0-9]+", header.lower()):
            if len(token) >= 4:
                terms.add(token)
    return terms


def _context_evidence(
    product: str,
    ontologies: Dict[str, List[str]],
    score: float,
    details: Dict[str, Any],
    fallback_id: str,
) -> List[Evidence]:
    """Build ontology Evidence for one context signal.

    When the associated product has Pfam/EC/GO cross-references the evidence
    is emitted on those ontology terms so it can converge with direct
    evidence. Without an ontology id a free-text label is kept (it cannot
    converge across channels and stays in the indirect context group).
    """
    out: List[Evidence] = []
    for term_type in ("Pfam", "EC", "GO"):
        for term_id in ontologies.get(term_type, []):
            d = dict(details)
            d["db"] = term_type.lower()
            out.append(
                Evidence(
                    channel=EvidenceChannel.CONTEXT,
                    target_id=term_id,
                    target_label=product,
                    score=score,
                    details=d,
                )
            )
    if not out:
        out.append(
            Evidence(
                channel=EvidenceChannel.CONTEXT,
                target_id=fallback_id,
                target_label=product,
                score=score,
                details=dict(details),
            )
        )
    return out


def run(
    workdir: WorkDir,
    config: Dict[str, Any],
    card: Optional[Path] = None,
    vfdb: Optional[Path] = None,
    reference_fna: Optional[Path] = None,
) -> None:
    records_path = workdir.path("prepare", "cds_records.json")
    loaded = read_json(records_path)  # type: ignore[assignment]
    records = _records_from_dicts(loaded)  # type: ignore[arg-type]
    _coerce_enums(records)
    _assign_families(records)
    log.info("Loaded %d CDS records for context analysis", len(records))

    card_terms = _reference_terms(Path(card)) if card else None
    vfdb_terms = _reference_terms(Path(vfdb)) if vfdb else None
    if card is None:
        log.info("CARD reference not supplied; resistance-neighbourhood check skipped.")
    if vfdb is None:
        log.info("VFDB reference not supplied; virulence-neighbourhood check skipped.")

    # ---- signals ---------------------------------------------------------
    neighbour_rows = operons.build_neighbours(records, config)
    profile_rows, family_samples, samples = phyloprofile.build_profiles(records, config)
    micro_rows = microsynteny.build_microsynteny(records, neighbour_rows, family_samples, samples)

    # Map every observed product to its ontology cross-references.
    product_ont: Dict[str, Dict[str, List[str]]] = {}
    for r in records:
        key = r.product.lower().strip()
        if key and key not in product_ont:
            product_ont[key] = parse_dbxrefs(r.dbxrefs)

    # ---- evidence --------------------------------------------------------
    evidence: Dict[str, List[Evidence]] = {}
    by_key = {r.key: r for r in records}

    for row in neighbour_rows:
        if not row["coupled"]:
            continue
        product = row["neighbour_product"]
        nb_key = f"{row['sample_id']}@{row['neighbour_locus']}"
        nb_rec = by_key.get(nb_key)
        if nb_rec is not None:
            ont = parse_dbxrefs(nb_rec.dbxrefs)
        else:
            ont = product_ont.get(product.lower().strip(), {})
        target_label = product
        pname = product.lower()
        if card_terms and any(term in pname for term in card_terms):
            target_label = f"[CARD] {product}"
        if vfdb_terms and any(term in pname for term in vfdb_terms):
            target_label = f"[VFDB] {product}"
        key = f"{row['sample_id']}@{row['locus_tag']}"
        details = {"kind": "operon_coupling", "intergenic_bp": row["intergenic_bp"]}
        evidence.setdefault(key, []).extend(
            _context_evidence(
                target_label, ont, 1.0, details, row["neighbour_locus"]
            )
        )

    # Profile correlation -> one best evidence entry per family.
    best: Dict[str, Dict[str, Any]] = {}
    for prow in profile_rows:
        if not prow["emphasis"]:
            continue
        cur = best.get(prow["family_id"])
        if cur is None or prow["jaccard"] > cur["jaccard"]:
            best[prow["family_id"]] = prow
    for family_id, prow in best.items():
        product = prow["annotated_gene"]
        ont = product_ont.get(product, {})
        for record in records:
            if record.family_id != family_id:
                continue
            evidence.setdefault(record.key, []).extend(
                _context_evidence(
                    product,
                    ont,
                    prow["jaccard"],
                    {"kind": "phylogenetic_profile"},
                    product,
                )
            )

    # Micro-synteny conservation -> evidence per family member.
    for mrow in micro_rows:
        if mrow["conservation_fraction"] <= 0:
            continue
        product = mrow["conserved_neighbour"]
        ont = product_ont.get(product, {})
        for record in records:
            if record.family_id != mrow["family_id"]:
                continue
            details = {
                "kind": "microsynteny",
                "classification": mrow["classification"],
            }
            evidence.setdefault(record.key, []).extend(
                _context_evidence(
                    product,
                    ont,
                    mrow["conservation_fraction"],
                    details,
                    product,
                )
            )

    # ---- outputs ----------------------------------------------------------
    _write_tsv(
        workdir.path("context", "neighbours.tsv"),
        ["sample_id", "locus_tag", "family_id", "neighbour_locus", "neighbour_product",
         "neighbour_status", "position", "intergenic_bp", "strand_match", "coupled"],
        neighbour_rows,
    )
    _write_tsv(
        workdir.path("context", "phyloprofiles.tsv"),
        ["family_id", "family_profile", "annotated_gene", "annotated_profile", "jaccard", "emphasis"],
        profile_rows,
    )
    _write_tsv(
        workdir.path("context", "microsynteny.tsv"),
        ["family_id", "lineage_width", "n_samples", "classification",
         "conserved_neighbour", "conservation_fraction", "samples"],
        micro_rows,
    )
    write_json(
        workdir.path("context", "context_evidence.json"),
        {key: [to_dict(e) for e in items] for key, items in evidence.items()},
    )

    workdir.mark_complete(
        "context",
        {"neighbours": len(neighbour_rows), "profiles": len(profile_rows),
         "families": len(micro_rows), "evidence_keys": len(evidence)},
    )
    log.info("context stage done: %d evidence-bearing records", len(evidence))
