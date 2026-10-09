"""A small, robust GFF3 parser.

The parser is deliberately dependency-free and tolerant of real-world
annotations from Bakta, Prokka and PGAP.  It understands the structural comment
lines (``##gff-version``, ``##sequence-region``, ``##FASTA``) and turns CDS
rows into feature dictionaries.  ``gene`` rows are only used to anchor the
contig/coordinate neighbourhood; the identity of a coding sequence is taken
from its own CDS row attributes (``ID``, ``Name``, ``locus_tag``, ``product``,
``Dbxref``).

Coordinates stay 1-based as written in the file and the strand character is
preserved verbatim.  Multi-value attributes are split on ``;`` into ``key=value``
pairs and percent-encoded values are URL-decoded.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import unquote

__all__ = ["parse_gff3", "detect_source", "GFFFeature"]


class GFFFeature:
    """A single parsed feature row (CDS or gene)."""

    __slots__ = (
        "contig",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    )

    def __init__(
        self,
        contig: str,
        feature_source: str,
        feature_type: str,
        start: int,
        end: int,
        score: str,
        strand: str,
        phase: str,
        attributes: Dict[str, str],
    ) -> None:
        self.contig = contig
        self.source = feature_source
        self.type = feature_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

    @property
    def locus_tag(self) -> Optional[str]:
        return self.attributes.get("locus_tag")

    @property
    def product(self) -> str:
        return self.attributes.get("product", "")

    @property
    def protein_id(self) -> Optional[str]:
        return (
            self.attributes.get("protein_id")
            or self.attributes.get("ID")
            or self.attributes.get("Name")
        )

    def dbxrefs(self) -> List[str]:
        raw = self.attributes.get("Dbxref") or self.attributes.get("dbxref") or ""
        return [item for item in (part.strip() for part in raw.split(",")) if item]

    def to_cds_dict(self) -> Dict[str, Any]:
        return {
            "contig": self.contig,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "locus_tag": self.locus_tag,
            "protein_id": self.protein_id,
            "product": self.product,
            "dbxrefs": self.dbxrefs(),
            "attributes": dict(self.attributes),
        }


_ATTR_SPLIT = re.compile(r";\s*")


def _parse_attributes(column: str) -> Dict[str, str]:
    """Parse the ninth GFF3 column into a URL-decoded key/value dict."""
    attributes: Dict[str, str] = {}
    if not column or column == ".":
        return attributes
    for pair in _ATTR_SPLIT.split(column.strip()):
        if not pair:
            continue
        # GFF3 allows key "=" value; the value itself may contain "=" (e.g. Dbxref URLs).
        key, sep, value = pair.partition("=")
        key = key.strip()
        if not key:
            continue
        value = unquote(value.strip()) if sep else ""
        # A repeated key keeps the first occurrence (Prokka repeats locus_tag).
        attributes.setdefault(key, value)
    return attributes


def parse_gff3(path: Path) -> List[GFFFeature]:
    """Parse all CDS/gene feature rows from a GFF3 file.

    Comment lines (starting with ``#``) are skipped structurally; the parser
    stops consuming feature rows when it hits the ``##FASTA`` delimiter.
    """
    path = Path(path)
    features: List[GFFFeature] = []
    in_fasta = False
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if not line.strip():
                continue
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if in_fasta:
                continue
            if line.startswith("#"):
                # Structural metadata: ##gff-version, ##sequence-region, etc.
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                # Tolerate whitespace-aligned (non-tab) rows minimally.
                fields = line.split()
                if len(fields) < 9:
                    continue
            contig, feature_source, ftype, start, end, score, strand, phase, attr_col = fields[:9]
            ftype = ftype.strip()
            if ftype not in ("CDS", "gene"):
                continue
            try:
                start_i = int(float(start))
                end_i = int(float(end))
            except ValueError:
                continue
            features.append(
                GFFFeature(
                    contig=contig.strip(),
                    feature_source=feature_source.strip(),
                    feature_type=ftype,
                    start=start_i,
                    end=end_i,
                    score=score.strip(),
                    strand=strand.strip(),
                    phase=phase.strip(),
                    attributes=_parse_attributes(attr_col),
                )
            )
    return features


def parse_cds(path: Path) -> List[Dict[str, Any]]:
    """Convenience wrapper returning only CDS rows as plain dictionaries."""
    return [f.to_cds_dict() for f in parse_gff3(path) if f.type == "CDS"]


def detect_source(path: Path) -> str:
    """Heuristically identify the annotation engine from file contents.

    Returns ``"bakta"``, ``"prokka"`` or ``"generic"``.  Bakta GFF files carry
    explanatory comment lines (e.g. ``##Annotation ...`` / "bakta" in the
    header), whereas Prokka features use ``PROKKA_``-prefixed identifiers.
    """
    path = Path(path)
    header_comments: List[str] = []
    prokka_id = False
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n\r")
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#"):
                    header_comments.append(line.lower())
                    continue
                fields = line.split("\t")
                if len(fields) >= 9:
                    attr = fields[8]
                    if "PROKKA_" in attr:
                        prokka_id = True
    except OSError:
        return "generic"

    blob = "\n".join(header_comments)
    if "bakta" in blob:
        return "bakta"
    if prokka_id or "prokka" in blob:
        return "prokka"
    return "generic"
