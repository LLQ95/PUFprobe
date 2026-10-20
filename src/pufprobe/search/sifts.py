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
    an additional domain suffix; only the PDB id and chain are used. A target
    without a chain suffix (single-chain-type entries such as ``101m``) is
    returned as the bare, lowercased PDB id.
    """
    parts = target.split("_")
    if len(parts) >= 2:
        return f"{parts[0].lower()}_{parts[1]}"
    return target.lower()


def load_pdb_chain_pfam(path: Path) -> Dict[str, List[str]]:
    """Load the SIFTS summary into a PDB-to-Pfam mapping.

    The returned dictionary carries both chain-level keys (``101m_A``) and
    bare PDB-id keys (``101m``). Foldseek writes a bare PDB id for structures
    that contain a single protein chain type (homomeric assemblies), so the
    bare-id value is the union of the Pfam families across every chain of
    that PDB. Heteromeric structures, which Foldseek indexes as explicit
    ``pdb_chain`` targets, are only queried through their chain-level keys.
    """
    chain_map: Dict[str, List[str]] = {}
    pdb_map: Dict[str, List[str]] = {}
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
            pdb = pdb.lower()
            chain_key = f"{pdb}_{chain}"
            bucket = chain_map.setdefault(chain_key, [])
            if pfam not in bucket:
                bucket.append(pfam)
            pdb_bucket = pdb_map.setdefault(pdb, [])
            if pfam not in pdb_bucket:
                pdb_bucket.append(pfam)
    mapping = dict(chain_map)
    mapping.update(pdb_map)
    return mapping
