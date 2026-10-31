"""Pure data logic for building the M-CSA pocket template index.

The offline build (driven by ``server/build_pocket_templates.py``) downloads
reference structures and runs foldseek, but the transformations of M-CSA JSON
into template records live here so they can be unit-tested without external
programs or network access. A template is one reference ``(pdb, chain)`` with
its catalytic residues and the M-CSA entry metadata (EC, enzyme name,
reference UniProt accession).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple


def collect_reference_sites(
    residues: List[dict],
) -> Dict[Tuple[str, str], dict]:
    """Group reference catalytic residues by ``(pdb, chain)``.

    Reads the flat M-CSA residues list and keeps only reference chain rows
    (``is_reference`` true). Each site stores its M-CSA id and a map of
    catalytic residue sequence positions (the M-CSA ``resid``) to three-letter
    residue codes. Residue positions are deduplicated within a site.
    """
    sites: Dict[Tuple[str, str], dict] = {}
    for residue in residues:
        mcsa_id = residue.get("mcsa_id")
        for chain_row in residue.get("residue_chains", []) or []:
            if not chain_row.get("is_reference"):
                continue
            pdb = str(chain_row.get("pdb_id", "")).strip().lower()
            chain = str(chain_row.get("chain_name", "")).strip()
            resid = chain_row.get("resid")
            code = str(chain_row.get("code", "")).strip()
            if not pdb or not chain or resid is None:
                continue
            key = (pdb, chain)
            site = sites.setdefault(
                key, {"mcsa_id": mcsa_id, "residues": {}}
            )
            site["residues"][int(resid)] = code
    return sites


def entries_by_id(entries: object) -> Dict[int, dict]:
    """Normalize an M-CSA entries payload into a map keyed by M-CSA id.

    Accepts either a full list of entries or a paginated dict with a
    ``results`` list (the first API page).
    """
    if isinstance(entries, dict):
        rows = entries.get("results", [])
    elif isinstance(entries, list):
        rows = entries
    else:
        rows = []
    out: Dict[int, dict] = {}
    for entry in rows:
        if isinstance(entry, dict) and entry.get("mcsa_id") is not None:
            out[int(entry["mcsa_id"])] = entry
    return out


def make_template_record(
    target: str,
    site: dict,
    entry: Optional[dict],
    catalytic_residues: List[dict],
) -> dict:
    """Assemble one template record for a reference ``(pdb, chain)``.

    ``catalytic_residues`` contains only the residues verified against the
    downloaded structure (``[{"resid": int, "code": str}, ...]``). The EC is
    taken from the entry's single ``all_ecs`` value; entries without an EC
    yield an empty string rather than a guessed value.
    """
    entry = entry or {}
    mcsa_id = int(site["mcsa_id"])
    pdb, chain = _split_target(target)
    ecs = entry.get("all_ecs") or []
    ec = str(ecs[0]) if ecs else ""
    return {
        "target": target,
        "mcsa_id": mcsa_id,
        "pdb": pdb,
        "chain": chain,
        "ec": ec,
        "enzyme_name": str(entry.get("enzyme_name", "")),
        "reference_uniprot": str(entry.get("reference_uniprot_id", "")),
        "catalytic_residues": [
            {"resid": int(r["resid"]), "code": str(r["code"])}
            for r in catalytic_residues
        ],
    }


def _split_target(target: str) -> Tuple[str, str]:
    parts = target.rsplit("_", 1)
    if len(parts) != 2:
        return target, ""
    return parts[0], parts[1]
