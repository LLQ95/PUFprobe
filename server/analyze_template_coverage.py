#!/usr/bin/env python
"""Diagnose why some M-CSA reference sites do not become pocket templates.

The build keeps a site only when at least one catalytic residue is found in the
downloaded reference chain with the same amino-acid code recorded by M-CSA.
This script re-runs that verification without foldseek, classifies every
catalytic residue (verified, amino-acid mismatch, or position out of range),
and groups the sites so the coverage loss can be attributed to a concrete
cause rather than a single count.

Usage:
    python analyze_template_coverage.py <mcsa_dir> <work_dir>
        [--report PATH]

Inputs and caches follow build_pocket_templates.py: the flat M-CSA files in
<mcsa_dir>, downloaded structures under <work_dir>/raw_pdb, and the merged
entry list under <work_dir>/entries_all.json.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

from pufprobe.search.pocket import (
    extract_chain_pdb,
    parse_pdb_residues,
)
from pufprobe.search.pocket_templates import collect_reference_sites

# Reuse the build orchestration helpers so download and parsing stay identical.
import build_pocket_templates as builder

_MCSA_RESIDUES_NAME = "mcsa_residues.json"
_REPORT_NAME = "coverage_report.json"


def classify_site(pdb_text: str, chain: str, expected: Dict[int, str]) -> dict:
    """Classify one reference chain's catalytic residues against its PDB."""
    chain_text = extract_chain_pdb(pdb_text, chain)
    residues = parse_pdb_residues(chain_text)
    by_index = {r.seq_index: r for r in residues}
    chain_len = len(residues)

    residue_results: List[dict] = []
    for resid, expected_code in sorted(expected.items()):
        residue = by_index.get(resid)
        if residue is None:
            status, actual = "out_of_range", None
        elif residue.resname.upper() != expected_code.upper():
            status, actual = "mismatch", residue.resname
        else:
            status, actual = "verified", residue.resname
        residue_results.append(
            {
                "resid": resid,
                "expected": expected_code,
                "actual": actual,
                "status": status,
            }
        )

    verified = [r for r in residue_results if r["status"] == "verified"]
    mismatches = [r for r in residue_results if r["status"] == "mismatch"]
    out_of_range = [r for r in residue_results if r["status"] == "out_of_range"]
    if len(verified) == len(residue_results):
        site_class = "full"
    elif verified:
        site_class = "partial"
    else:
        site_class = "none"
    return {
        "chain_length": chain_len,
        "site_class": site_class,
        "residues": residue_results,
        "counts": {
            "verified": len(verified),
            "mismatch": len(mismatches),
            "out_of_range": len(out_of_range),
        },
    }


def offset_guess(pdb_text: str, chain: str, expected: Dict[int, str]) -> Optional[int]:
    """Find a constant sequence-position offset that makes every residue match.

    A single offset that repairs an otherwise failing site points to a numbering
    convention mismatch rather than unrelated structural differences.
    """
    chain_text = extract_chain_pdb(pdb_text, chain)
    by_index = {r.seq_index: r for r in parse_pdb_residues(chain_text)}
    offsets = {
        candidate
        for residue in by_index.values()
        for resid, code in expected.items()
        if residue.resname.upper() == code.upper()
        for candidate in (residue.seq_index - resid,)
    }
    for offset in offsets:
        if all(
            (r := by_index.get(resid + offset)) is not None
            and r.resname.upper() == code.upper()
            for resid, code in expected.items()
        ):
            return offset
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mcsa_dir")
    parser.add_argument("work_dir")
    parser.add_argument("--report", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mcsa_dir = Path(args.mcsa_dir)
    work_dir = Path(args.work_dir)
    raw_dir = work_dir / "raw_pdb"
    report_path = Path(args.report) if args.report else work_dir / _REPORT_NAME

    residues_path = mcsa_dir / _MCSA_RESIDUES_NAME
    with residues_path.open("r", encoding="utf-8") as handle:
        residues = json.load(handle)
    sites = collect_reference_sites(residues)

    detailed: Dict[str, dict] = {}
    for index, ((pdb_id, chain), site) in enumerate(sorted(sites.items()), 1):
        target = f"{pdb_id}_{chain}"
        try:
            pdb_path = builder.download_pdb(pdb_id, raw_dir)
            pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
        except RuntimeError as exc:
            detailed[target] = {
                "mcsa_id": site["mcsa_id"],
                "site_class": "download_failed",
                "error": str(exc),
                "chain_length": 0,
                "residues": [],
                "counts": {
                    "verified": 0,
                    "mismatch": 0,
                    "out_of_range": 0,
                },
            }
            print(f"  [warn] {target}: {exc}")
            continue
        result = classify_site(pdb_text, chain, site["residues"])
        result["mcsa_id"] = site["mcsa_id"]
        if result["site_class"] == "none":
            result["offset_guess"] = offset_guess(
                pdb_text, chain, site["residues"]
            )
        detailed[target] = result
        if index % 100 == 0:
            print(f"processed {index}/{len(sites)} sites")

    site_classes = {
        label: sum(1 for r in detailed.values() if r["site_class"] == label)
        for label in ("full", "partial", "none", "download_failed")
    }
    residue_status = {label: 0 for label in ("verified", "mismatch", "out_of_range")}
    for result in detailed.values():
        for label, count in result["counts"].items():
            residue_status[label] += count
    offset_fixed = sorted(
        (target, result["offset_guess"])
        for target, result in detailed.items()
        if result["site_class"] == "none" and result["offset_guess"] is not None
    )

    summary = {
        "total_sites": len(detailed),
        "site_classes": site_classes,
        "residue_status": residue_status,
        "none_sites_with_constant_offset": len(offset_fixed),
    }
    payload = {"summary": summary, "sites": detailed}
    report_path.write_text(
        json.dumps(payload, indent=2), encoding="utf-8", newline="\n"
    )

    print("=== coverage summary ===")
    print(json.dumps(summary, indent=2))
    if offset_fixed:
        print("constant-offset repairs (target, offset):", offset_fixed[:20])
    print("detailed report written to", report_path)


if __name__ == "__main__":
    main()
