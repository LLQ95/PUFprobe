"""Ontology identifier parsing shared by search and context.

Functional annotations from GFF Dbxref attributes and structure-database
cross-references are written as "database:identifier" tokens (for example
GO:0003735, Pfam:PF00042, EC:3.4.24.77). These helpers convert such tokens
into normalized identifiers so that evidence from different channels can
converge on the same functional term.
"""

from __future__ import annotations

import re
from typing import Dict, List

_GO_BARE_RE = re.compile(r"^\d{7}$")
_GO_RE = re.compile(r"^GO:\d{7}$")
_EC_RE = re.compile(r"^\d+\.\d+\.\d+\.[\dR\-]+(?:\.[A-Z])?$")
_PFAM_RE = re.compile(r"^PF\d{5}")


def normalize_go(token: str) -> str:
    token = token.strip()
    return token if token.startswith("GO:") else f"GO:{token}"


def normalize_ec(token: str) -> str:
    token = token.strip()
    return token if token.startswith("EC:") else f"EC:{token}"


def parse_dbxrefs(dbxrefs: List[str]) -> Dict[str, List[str]]:
    """Split GFF Dbxref tokens into normalized Pfam, EC and GO identifiers.

    Tokens that use another database prefix (KEGG, RefSeq, COG) or that do not
    match an ontology pattern are ignored.
    """
    pfam: List[str] = []
    ec: List[str] = []
    go: List[str] = []

    for raw in dbxrefs:
        token = raw.strip()
        if not token:
            continue
        head, sep, value = token.partition(":")
        key = head.strip().lower()
        val = value.strip()

        if key == "pfam" and _PFAM_RE.match(val):
            pfam.append(val)
        elif key == "go" and _GO_BARE_RE.match(val):
            go.append(normalize_go(val))
        elif key == "ec" and _EC_RE.match(val):
            ec.append(normalize_ec(val))
        elif _PFAM_RE.match(token):
            pfam.append(token)
        elif _GO_RE.match(token):
            go.append(token)
        elif _EC_RE.match(token):
            ec.append(normalize_ec(token))

    return {
        "Pfam": sorted(set(pfam)),
        "EC": sorted(set(ec)),
        "GO": sorted(set(go)),
    }
