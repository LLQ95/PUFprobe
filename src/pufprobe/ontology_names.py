"""Official family names and descriptions for ontology terms.

The sequence and structure channels converge on Pfam family accessions. A
Foldseek hit carries the PDB template name, and the best sequence hit may
carry a distant eukaryotic homolog; neither is the family's own name. Pfam
ships a small STOCKHOLM metadata file, ``Pfam-A.hmm.dat``, in which every
family record carries its accession (``#=GF AC``), short identifier
(``#=GF ID``), official description (``#=GF DE``) and entry type
(``#=GF TP``). These helpers parse that file into a compact
``accession -> name`` index and resolve the display name of a Pfam term.
Accessions are stored without the trailing version number, so PF00005.34 and
PF00005 converge on the same key.
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path
from typing import Dict, IO, Optional

from pufprobe.logging_utils import get_logger

logger = get_logger(__name__)

_GF_PREFIX = "#=GF "
_VERSION_RE = re.compile(r"^(PF\d{5})\.\d+$")


def _open_text(path: Path) -> IO[str]:
    """Open a flat text file, transparently decoding gzip."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _strip_version(accession: str) -> str:
    """Drop the trailing version number from a Pfam accession (PF00005.34)."""
    match = _VERSION_RE.match(accession.strip())
    if match:
        return match.group(1)
    return accession.strip()


def parse_pfam_dat(handle: IO[str]) -> Dict[str, Dict[str, str]]:
    """Parse a Pfam-A.hmm.dat stream into an accession-to-name index.

    Family records are separated by ``//`` lines. Accessions are stored
    without their version number. Records that carry neither an identifier
    nor a description are omitted.
    """
    names: Dict[str, Dict[str, str]] = {}

    acc = ""
    fam_id = ""
    desc = ""
    entry_type = ""

    def flush() -> None:
        if not acc:
            return
        record: Dict[str, str] = {}
        if fam_id:
            record["id"] = fam_id
        if desc:
            record["name"] = desc
        if entry_type:
            record["type"] = entry_type
        if record:
            names[acc] = record

    for raw_line in handle:
        line = raw_line.rstrip("\n")
        if line.startswith(_GF_PREFIX):
            key, _, value = line[len(_GF_PREFIX):].partition(" ")
            value = value.strip()
            if key == "AC":
                acc = _strip_version(value)
            elif key == "ID":
                fam_id = value
            elif key == "DE":
                desc = value
            elif key == "TP":
                entry_type = value
        elif line.strip() == "//":
            flush()
            acc = ""
            fam_id = ""
            desc = ""
            entry_type = ""

    flush()
    return names


def build_pfam_names(
    dat_path: Path, out_path: Path
) -> Dict[str, Dict[str, str]]:
    """Parse ``Pfam-A.hmm.dat[.gz]`` and write a compact JSON index."""
    with _open_text(dat_path) as handle:
        names = parse_pfam_dat(handle)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(names, handle, separators=(",", ":"))
    logger.info("Pfam family name index: %d families -> %s", len(names), out_path)
    return names


def load_pfam_names(path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    """Load a Pfam name index, returning an empty dict when unavailable."""
    if path is None or not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        return {}
    return data


def pfam_term_name(
    names: Dict[str, Dict[str, str]], accession: str
) -> Optional[str]:
    """Resolve the official family description for a Pfam accession.

    The official description is preferred, then the short family identifier.
    ``None`` is returned when the accession is absent so the caller can fall
    back to its own label logic.
    """
    if not names:
        return None
    record = names.get(_strip_version(accession))
    if not record:
        return None
    return record.get("name") or record.get("id") or None
