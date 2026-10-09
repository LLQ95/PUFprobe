"""A minimal GenBank flat-file parser.

This is not a full GenBank parser: it only extracts what PUFprobe needs, namely
the ``CDS`` features with their ``/locus_tag``, ``/product``, ``/protein_id``,
``/db_xref``, ``/transl_table`` and optional ``/translation`` qualifiers, plus
the nucleotide sequence from ``ORIGIN`` when present.  Feature coordinates
understand ``complement(...)`` and ``join(...)`` location expressions.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional

__all__ = ["parse_gbk", "parse_gbk_cds"]

_INT = re.compile(r"(\d+)")
_QUAL = re.compile(r'^\s*/([^=]+)(?:=(.*))?$')


def _parse_location(location: str) -> Dict[str, Any]:
    """Return start/end/strand from a GenBank location expression."""
    location = location.strip()
    strand = "-" if "complement" in location else "+"
    numbers = [int(n) for n in _INT.findall(location)]
    if not numbers:
        return {"start": 0, "end": 0, "strand": strand, "raw": location}
    return {
        "start": min(numbers),
        "end": max(numbers),
        "strand": strand,
        "raw": location,
    }


def _clean_qualifier_value(value: str) -> str:
    value = value.strip()
    # Quoted strings may span lines; the closing quote may be on the last line.
    if value.startswith('"'):
        value = value[1:]
    if value.endswith('"'):
        value = value[:-1]
    return value.strip()


def parse_gbk(path: Path) -> Dict[str, Any]:
    """Parse a GenBank file into ``{locus, accession, cds: [...], sequence_nt}``."""
    path = Path(path)
    locus = ""
    accession = ""
    features: List[Dict[str, Any]] = []

    section: Optional[str] = None
    current: Optional[Dict[str, Any]] = None
    current_qual: Optional[str] = None
    origin_lines: List[str] = []
    in_origin = False

    def _flush() -> None:
        nonlocal current
        if current is not None:
            features.append(current)
        current = None

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")

            if line.startswith("//"):
                _flush()
                break
            if line.startswith("ORIGIN"):
                _flush()
                in_origin = True
                continue
            if in_origin:
                # Strip the leading offset and whitespace; keep nucleotide letters.
                origin_lines.append(re.sub(r"[^a-zA-Z]", "", line))
                continue

            key = line[:12].strip()
            if key == "LOCUS":
                locus = line[12:].strip().split()[0]
                continue
            if key == "ACCESSION":
                accession = line[12:].strip().split()[0]
                continue
            if key == "FEATURES":
                section = "FEATURES"
                continue

            if section == "FEATURES":
                # Feature name starts in column ~6; qualifiers start with '/'.
                field_start = line[:21]
                body = line[21:] if len(line) > 21 else ""

                qual_match = _QUAL.match(line)
                if qual_match and line.lstrip().startswith("/"):
                    qname = qual_match.group(1).strip()
                    qvalue = qual_match.group(2) or ""
                    current_qual = qname
                    if current is not None:
                        current["qualifiers"][qname] = _clean_qualifier_value(qvalue)
                    continue

                # A new feature begins when columns 6-20 hold a non-empty token.
                feature_token = field_start[5:21].strip()
                if feature_token and not body.startswith(" "):
                    _flush()
                    if feature_token == "CDS":
                        loc = _parse_location(body)
                        current = {
                            "type": "CDS",
                            "contig": locus or accession or "chr",
                            "start": loc["start"],
                            "end": loc["end"],
                            "strand": loc["strand"],
                            "location": loc["raw"],
                            "qualifiers": {},
                        }
                        current_qual = None
                    else:
                        current = None
                        current_qual = None
                    continue

                # Continuation of the current feature (location or qualifier value).
                if current is not None:
                    if current_qual is not None:
                        prev = current["qualifiers"].get(current_qual, "")
                        current["qualifiers"][current_qual] = _clean_qualifier_value(
                            (prev + " " + body.strip()).strip()
                        )
                    else:
                        current["location"] += " " + body.strip()
                        loc = _parse_location(current["location"])
                        current["start"] = loc["start"]
                        current["end"] = loc["end"]
                        current["strand"] = loc["strand"]
                continue

    cds_records: List[Dict[str, Any]] = []
    for feat in features:
        if feat["type"] != "CDS":
            continue
        q = feat["qualifiers"]
        dbxrefs_raw = q.get("db_xref", "")
        dbxrefs = [d.strip() for d in dbxrefs_raw.split(",") if d.strip()]
        cds_records.append(
            {
                "contig": feat["contig"],
                "start": feat["start"],
                "end": feat["end"],
                "strand": feat["strand"],
                "locus_tag": q.get("locus_tag"),
                "protein_id": q.get("protein_id"),
                "product": q.get("product", ""),
                "dbxrefs": dbxrefs,
                "transl_table": int(q["transl_table"]) if q.get("transl_table", "").isdigit() else 11,
                "translation": q.get("translation", "").replace(" ", ""),
                "attributes": dict(q),
            }
        )

    return {
        "locus": locus,
        "accession": accession,
        "cds": cds_records,
        "sequence_nt": "".join(origin_lines).upper(),
    }


def parse_gbk_cds(path: Path) -> List[Dict[str, Any]]:
    """Return only the CDS dictionaries from a GenBank file."""
    return parse_gbk(path)["cds"]
