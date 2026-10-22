"""Term normalization through UniProt cross-references.

Sequence search returns hits in the target database namespace: an MMseqs2 hit
against Swiss-Prot carries an accession such as ``sp|P0ACF0|DBHA_ECOLI``, a
Foldseek hit carries a PDB identifier, and context carries a neighbour gene
name. Independent evidence groups can only reinforce the same functional term
when those names are mapped onto a shared ontology (Pfam, EC, GO).

This module parses the Swiss-Prot flat text file (``uniprot_sprot.dat.gz``),
which holds a curated cross-reference block for every entry, into a compact
``accession -> {Pfam, EC, GO}`` index. Homology hits are then supplemented with
ontology evidence. The evidence channel is unchanged (the hit still originates
from sequence homology), so the same cross-reference from several Swiss-Prot
accessions remains one evidence group; a structurally derived hit mapped to the
same Pfam later supplies a genuinely independent group.

Gene Ontology terms are emitted as identifiers without true-path propagation,
because no ontology graph is loaded in this release.
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path
from typing import Dict, IO, List, Optional, Tuple

from pufprobe.logging_utils import get_logger
from pufprobe.models import Evidence, EvidenceChannel

logger = get_logger(__name__)

# Pfam short identifier for an uncharacterised family, e.g. "DUF874". Profile
# hits frequently carry such a token in the description ("DUF874 family").
_DUF_NAME_RE = re.compile(r"\bDUF\d+\b")

# DE line carrying an EC number, e.g. "DE            EC=3.4.24.-;".
_EC_DE_RE = re.compile(r"EC=([0-9]+\.[0-9]+\.[0-9]+\.[0-9A-Za-z\-]+)")
_DR_PREFIX = "DR   "
_DE_PREFIX = "DE"


def _open_text(path: Path) -> IO[str]:
    """Open a flat text file, transparently decoding gzip."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_sprot_dat(handle: IO[str]) -> Dict[str, Dict[str, List[str]]]:
    """Parse a Swiss-Prot flat text stream into an accession cross-reference index.

    Every primary and secondary accession of an entry is mapped to the same
    record. Entries with no Pfam, EC, or GO cross-reference are omitted.
    """
    xref: Dict[str, Dict[str, List[str]]] = {}

    accessions: List[str] = []
    pfam: List[str] = []
    ec: List[str] = []
    go: List[str] = []

    def flush() -> None:
        if not accessions:
            return
        entry: Dict[str, List[str]] = {}
        if pfam:
            entry["Pfam"] = sorted(set(pfam))
        if ec:
            entry["EC"] = sorted(set(ec))
        if go:
            entry["GO"] = sorted(set(go))
        if entry:
            for accession in accessions:
                xref[accession] = entry

    for raw_line in handle:
        line = raw_line.rstrip("\n")

        if line.startswith("AC"):
            # Accession lines may list several accessions, ;-separated.
            payload = line[5:].lstrip()
            for token in payload.split(";"):
                token = token.strip()
                if token:
                    accessions.append(token)
        elif line.startswith(_DR_PREFIX + "Pfam;"):
            accession = line[len(_DR_PREFIX + "Pfam;"):].split(";", 1)[0].strip()
            if accession:
                pfam.append(accession)
        elif line.startswith(_DR_PREFIX + "GO;"):
            accession = line[len(_DR_PREFIX + "GO;"):].split(";", 1)[0].strip()
            if accession:
                go.append(accession)
        elif line.startswith(_DE_PREFIX):
            for match in _EC_DE_RE.finditer(line):
                ec.append(match.group(1))
        elif line.strip() == "//":
            flush()
            accessions = []
            pfam = []
            ec = []
            go = []

    flush()
    return xref


def build_xref(dat_path: Path, out_path: Path) -> Dict[str, Dict[str, List[str]]]:
    """Parse ``uniprot_sprot.dat[.gz]`` and write a compact JSON index."""
    with _open_text(dat_path) as handle:
        xref = parse_sprot_dat(handle)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(xref, handle, separators=(",", ":"))
    logger.info(
        "Swiss-Prot cross-reference index: %d accessions with ontology terms -> %s",
        len(xref), out_path,
    )
    return xref


def load_xref(path: Optional[Path]) -> Dict[str, Dict[str, List[str]]]:
    """Load a cross-reference index, returning an empty dict when unavailable."""
    if path is None or not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        return {}
    return data


def accessions_from_target(target_id: str) -> List[str]:
    """Extract UniProt accession(s) from an MMseqs2 target identifier.

    Swiss-Prot identifiers use the form ``sp|ACCESSION|ENTRY_NAME``; targets
    built from other sources are returned as a single bare identifier.
    """
    tid = (target_id or "").strip()
    if not tid:
        return []
    parts = tid.split("|")
    if len(parts) >= 2 and parts[1].strip():
        return [parts[1].strip()]
    return [tid]


def _ontology_id(term_type: str, raw_id: str) -> str:
    if term_type == "EC" and not raw_id.startswith("EC:"):
        return f"EC:{raw_id}"
    return raw_id


def expand_evidence(
    evidence: Dict[str, List[Evidence]],
    xref: Dict[str, Dict[str, List[str]]],
) -> Tuple[Dict[str, List[Evidence]], int]:
    """Supplement homology hits with Pfam/EC/GO evidence from the index.

    The original hit is preserved; each cross-reference becomes an additional
    Evidence object in the same channel so that ``rank.classify_term`` assigns
    it to the ontology term. Returns the expanded evidence and the number of
    added objects.
    """
    if not xref:
        return evidence, 0

    expanded: Dict[str, List[Evidence]] = {}
    added = 0

    for key, hits in evidence.items():
        out = list(hits)
        for hit in hits:
            if hit.channel != EvidenceChannel.REMOTE_HOMOLOGY:
                continue
            for accession in accessions_from_target(hit.target_id):
                entry = xref.get(accession)
                if not entry:
                    continue
                for term_type in ("Pfam", "EC", "GO"):
                    for raw_id in entry.get(term_type, []):
                        details = dict(hit.details or {})
                        details["db"] = term_type.lower()
                        details["source_accession"] = accession
                        details["xref"] = True
                        out.append(
                            Evidence(
                                channel=hit.channel,
                                target_id=_ontology_id(term_type, raw_id),
                                target_label=hit.target_label,
                                score=hit.score,
                                evalue=hit.evalue,
                                coverage=hit.coverage,
                                details=details,
                            )
                        )
                        added += 1
        expanded[key] = out
    return expanded, added


def expand_homology_pfam(
    evidence: Dict[str, List[Evidence]],
    pfam_id_index: Dict[str, str],
) -> Tuple[Dict[str, List[Evidence]], int]:
    """Supplement REMOTE_HOMOLOGY hits with Pfam accessions.

    Profile hits carry a description that often contains a Pfam short
    identifier for an uncharacterised family (e.g. "DUF874 family"). Mapping
    that token to the Pfam accession lets the sequence-profile group converge
    on the same term as a database channel that reports the same Pfam. The
    channel stays REMOTE_HOMOLOGY (the hit still originates from the profile
    search), so the added evidence remains one group; a database channel on
    the same accession is the genuinely independent second group.
    """
    if not pfam_id_index:
        return evidence, 0

    expanded: Dict[str, List[Evidence]] = {}
    seen: set = set()
    added = 0

    for key, hits in evidence.items():
        out = list(hits)
        for hit in hits:
            if hit.channel != EvidenceChannel.REMOTE_HOMOLOGY:
                continue
            for short in _DUF_NAME_RE.findall(hit.target_label or ""):
                accession = pfam_id_index.get(short)
                if not accession or (key, accession) in seen:
                    continue
                seen.add((key, accession))
                details = dict(hit.details or {})
                details["db"] = "pfam"
                details["pfam_short"] = short
                details["xref"] = True
                out.append(
                    Evidence(
                        channel=hit.channel,
                        target_id=accession,
                        target_label=accession,
                        score=hit.score,
                        evalue=hit.evalue,
                        coverage=hit.coverage,
                        details=details,
                    )
                )
                added += 1
        expanded[key] = out
    return expanded, added
