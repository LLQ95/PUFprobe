"""SIFTS structure-to-Pfam mapping.

The SIFTS project (PDBe) provides a weekly mapping of PDB chains to UniProt
entries and Pfam families. PUFprobe uses the small chain-level summary
``pdb_chain_pfam.tsv.gz`` to turn Foldseek PDB-chain hits into Pfam family
identifiers so that structure evidence converges on the same functional terms
as sequence evidence.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, List

_PDB_PREFIX = "pdb"


def _open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def pdb_chain_key(target: str) -> str:
    """Normalize a Foldseek PDB target to the SIFTS key ``pdb_chain``.

    Foldseek PDB targets are written as ``1abc_A`` for a chain and may carry
    an additional domain suffix; only the PDB id and chain are used.
    """
    parts = target.split("_")
    if len(parts) >= 2:
        return f"{parts[0].lower()}_{parts[1]}"
    return target.lower()


def load_pdb_chain_pfam(path: Path) -> Dict[str, List[str]]:
    """Load the SIFTS summary into a ``pdb_chain -> [Pfam id]`` mapping."""
    mapping: Dict[str, List[str]] = {}
    with _open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            pdb, chain, _sp, pfam = cols[:4]
            if pdb == "PDB" or not pfam:
                continue
            key = f"{pdb.lower()}_{chain}"
            bucket = mapping.setdefault(key, [])
            if pfam not in bucket:
                bucket.append(pfam)
    return mapping
