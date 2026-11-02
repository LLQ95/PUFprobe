#!/usr/bin/env python
"""Diagnose why some M-CSA reference sites do not become pocket templates.

The build keeps a site only when at least one catalytic residue is found in the
downloaded reference chain with the same amino-acid code recorded by M-CSA.
The current build treats the M-CSA ``resid`` as a foldseek 1-based in-chain
sequence index. M-CSA also provides ``auth_resid``, the PDB author numbering,
which is the numbering carried by the parsed PDB residues. This script runs
both numbering schemes without foldseek and reports how many sites each one
recovers, so the coverage loss can be attributed to a concrete cause.

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

# Reuse the build orchestration helpers so download and parsing stay identical.
import build_pocket_templates as builder

_MCSA_RESIDUES_NAME = "mcsa_residues.json"
_REPORT_NAME = "coverage_report.json"


def collect_sites(residues: List[dict]) -> Dict[tuple, dict]:
    """Group reference catalytic residues, keeping both numbering schemes."""
    sites: Dict[tuple, dict] = {}
    for residue in residues:
        mcsa_id = residue.get("mcsa_id")
        for chain_row in residue.get("residue_chains", []) or []:
            if not chain_row.get("is_reference"):
                continue
            pdb = str(chain_row.get("pdb_id", "")).strip().lower()
            chain = str(chain_row.get("chain_name", "")).strip()
            resid = chain_row.get("resid")
            if not pdb or not chain or resid is None:
                continue
            key = (pdb, chain)
            site = sites.setdefault(
                key, {"mcsa_id": mcsa_id, "residues": {}}
            )
            site["residues"][int(resid)] = {
                "resid": int(resid),
                "auth_resid": chain_row.get("auth_resid"),
                "code": str(chain_row.get("code", "")).strip(),
            }
    return sites


def _classify(residues: list, expected: Dict[int, dict], mode: str) -> dict:
    """Classify one chain under a selected numbering scheme."""
    by_index = {r.seq_index: r for r in residues}
    by_auth = {r.auth_resid: r for r in residues}
    auth_values = [r.auth_resid for r in residues]
    auth_min = min(auth_values) if auth_values else 0
    auth_max = max(auth_values) if auth_values else 0

    residue_results: List[dict] = []
    for resid, item in sorted(expected.items()):
        use_auth = mode == "auth" and item["auth_resid"] is not None
        if use_auth:
            residue = by_auth.get(item["auth_resid"])
            anchor = item["auth_resid"]
        else:
            residue = by_index.get(resid)
            anchor = resid
        expected_code = item["code"]
        missing_reason = None
        if residue is None:
            status, actual = "missing", None
            if use_auth:
                within = auth_min <= anchor <= auth_max
            else:
                within = 1 <= anchor <= len(residues)
            missing_reason = "internal_gap" if within else "out_of_range"
        elif residue.resname.upper() != expected_code.upper():
            status, actual = "mismatch", residue.resname
        else:
            status, actual = "verified", residue.resname
        residue_results.append(
            {
                "anchor": anchor,
                "resid": resid,
                "auth_resid": item["auth_resid"],
                "expected": expected_code,
                "actual": actual,
                "status": status,
                "missing_reason": missing_reason,
            }
        )

    verified = [r for r in residue_results if r["status"] == "verified"]
    mismatch = [r for r in residue_results if r["status"] == "mismatch"]
    missing = [r for r in residue_results if r["status"] == "missing"]
    if len(verified) == len(residue_results):
        site_class = "full"
    elif verified:
        site_class = "partial"
    else:
        site_class = "none"
    return {
        "chain_length": len(residues),
        "site_class": site_class,
        "residues": residue_results,
        "counts": {
            "verified": len(verified),
            "mismatch": len(mismatch),
            "missing": len(missing),
        },
    }


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
    sites = collect_sites(residues)

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
            }
            print(f"  [warn] {target}: {exc}")
            continue
        chain_residues = parse_pdb_residues(extract_chain_pdb(pdb_text, chain))
        detailed[target] = {
            "mcsa_id": site["mcsa_id"],
            "seq_based": _classify(chain_residues, site["residues"], "seq"),
            "auth_based": _classify(chain_residues, site["residues"], "auth"),
        }
        if index % 100 == 0:
            print(f"processed {index}/{len(sites)} sites")

    def class_counts(mode: str) -> Dict[str, int]:
        key = "seq_based" if mode == "seq" else "auth_based"
        return {
            label: sum(
                1 for r in detailed.values()
                if r.get("site_class") == "download_failed"
                and label == "download_failed"
                or r.get(key, {}).get("site_class") == label
            )
            for label in ("full", "partial", "none", "download_failed")
        }

    def residue_counts(mode: str) -> Dict[str, int]:
        key = "seq_based" if mode == "seq" else "auth_based"
        totals = {"verified": 0, "mismatch": 0, "missing": 0}
        for r in detailed.values():
            if key in r:
                for label, count in r[key]["counts"].items():
                    totals[label] += count
        return totals

    summary = {
        "total_sites": len(detailed),
        "seq_based_sites": class_counts("seq"),
        "auth_based_sites": class_counts("auth"),
        "seq_based_residues": residue_counts("seq"),
        "auth_based_residues": residue_counts("auth"),
    }
    report_path.write_text(
        json.dumps({"summary": summary, "sites": detailed}, indent=2),
        encoding="utf-8",
        newline="\n",
    )

    print("=== coverage summary ===")
    print(json.dumps(summary, indent=2))
    print("detailed report written to", report_path)


if __name__ == "__main__":
    main()
