"""Static analysis of the M-CSA pocket template library.

This script works from the M-CSA source JSON files and does not download any
structure or run foldseek. It reports the functional coverage of the template
library (EC classes, distinct EC numbers, templates per EC), the composition of
the catalytic residues, the numbering conventions that distinguish the M-CSA
``resid`` from the PDB author number, and the small set of sites that cannot be
retained. When a coverage report from ``analyze_template_coverage.py`` is given,
the retained-site statistics use the verified set instead of the full set.

Usage:
    python server/analyze_template_library.py <mcsa_dir> [--report PATH]
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional

from pufprobe.search.pocket_templates import (
    collect_reference_sites,
    entries_by_id,
)

_RESIDUES_NAME = "mcsa_residues.json"
_ENTRIES_NAME = "entries_all.json"
# The four sites the verified build does not retain. Two have no catalytic
# residue that maps to a structure, and two reference PDB entries that have
# been withdrawn.
_DEFAULT_EXCLUDED = {"1g64_B", "1g9z_B", "4v40_A", "4v4e_A"}


def _load_json(path: Path) -> object:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _ec_class(ec: str) -> str:
    return ec.split(".", 1)[0] if ec else "unknown"


def _percentiles(values: List[int]) -> Dict[str, float]:
    ordered = sorted(values)
    if not ordered:
        return {}

    def pick(fraction: float) -> int:
        index = round(fraction * (len(ordered) - 1))
        return ordered[index]

    return {
        "min": ordered[0],
        "p25": pick(0.25),
        "median": pick(0.5),
        "p75": pick(0.75),
        "max": ordered[-1],
    }


def _retained_targets(sites: dict, coverage: Optional[dict]) -> set:
    all_targets = {f"{pdb}_{chain}" for (pdb, chain) in sites}
    if not coverage:
        return all_targets - set(_DEFAULT_EXCLUDED)
    retained = set()
    for target, record in coverage.get("sites", {}).items():
        auth = record.get("auth_based", {})
        if auth.get("site_class") in {"full", "partial"}:
            retained.add(target)
    return retained


def analyze(
    residues: List[dict],
    entries: List[dict],
    coverage: Optional[dict],
) -> dict:
    sites = collect_reference_sites(residues)
    entry_map = entries_by_id(entries)
    retained = _retained_targets(sites, coverage)

    distinct_pdb = {pdb for (pdb, _chain) in sites}
    distinct_mcsa = {site["mcsa_id"] for site in sites.values()}

    # Functional coverage. Each entry has exactly one EC.
    ec_by_target = {}
    for (pdb, chain), site in sites.items():
        target = f"{pdb}_{chain}"
        entry = entry_map.get(int(site["mcsa_id"]), {})
        ecs = entry.get("all_ecs") or []
        ec_by_target[target] = str(ecs[0]) if ecs else ""

    retained_ec = [ec_by_target[t] for t in sorted(retained) if ec_by_target[t]]
    ec_class_counts = Counter(_ec_class(ec) for ec in retained_ec)
    templates_per_ec = Counter(retained_ec)
    singleton_ec = sorted(ec for ec, count in templates_per_ec.items() if count == 1)

    # Catalytic residue composition over the reference chain rows.
    residue_code_counts = Counter()
    residues_per_site = {}
    single_residue_sites = []
    for (pdb, chain), site in sorted(sites.items()):
        target = f"{pdb}_{chain}"
        codes = [info["code"] for info in site["residues"].values() if info["code"]]
        residues_per_site[target] = len(codes)
        residue_code_counts.update(codes)
        if len(codes) == 1:
            single_residue_sites.append(target)

    # Numbering conventions, examined on the raw chain rows.
    row_equal = row_different = row_auth_none = row_non_int_auth = 0
    per_site_numbering = Counter()
    non_int_auth_examples = []
    for residue in residues:
        for chain_row in residue.get("residue_chains", []) or []:
            if not chain_row.get("is_reference"):
                continue
            resid = chain_row.get("resid")
            auth = chain_row.get("auth_resid")
            if resid is None:
                continue
            if auth is None:
                row_auth_none += 1
            elif not isinstance(auth, int):
                row_non_int_auth += 1
                if len(non_int_auth_examples) < 10:
                    non_int_auth_examples.append(
                        {
                            "pdb": str(chain_row.get("pdb_id", "")).lower(),
                            "chain": chain_row.get("chain_name"),
                            "resid": resid,
                            "auth_resid": auth,
                        }
                    )
            elif int(resid) == int(auth):
                row_equal += 1
            else:
                row_different += 1

    for (pdb, chain), site in sites.items():
        has_missing = has_non_int = has_shifted = False
        for resid_key, info in site["residues"].items():
            auth = info["auth_resid"]
            if auth is None:
                has_missing = True
            elif not isinstance(auth, int):
                has_non_int = True
            elif int(auth) != int(resid_key):
                has_shifted = True
        if has_missing:
            per_site_numbering["contains_missing_auth"] += 1
        if has_non_int:
            per_site_numbering["contains_non_integer_auth"] += 1
        if has_shifted:
            per_site_numbering["contains_shifted_auth"] += 1
        if not has_missing and not has_non_int and not has_shifted:
            per_site_numbering["all_resid_equals_auth"] += 1

    excluded = sorted(
        f"{pdb}_{chain}"
        for (pdb, chain) in sites
        if f"{pdb}_{chain}" not in retained
    )

    return {
        "scale": {
            "reference_sites": len(sites),
            "retained_templates": len(retained),
            "distinct_reference_pdb": len(distinct_pdb),
            "distinct_mcsa_entries": len(distinct_mcsa),
        },
        "ec_coverage": {
            "distinct_ec_numbers": len(templates_per_ec),
            "ec_class_counts": dict(sorted(ec_class_counts.items())),
            "templates_per_ec_distribution": _percentiles(
                list(templates_per_ec.values())
            ),
            "single_template_ec_count": len(singleton_ec),
            "single_template_ec_examples": singleton_ec[:15],
        },
        "residue_composition": {
            "residue_code_counts": dict(
                residue_code_counts.most_common()
            ),
            "residues_per_site_distribution": _percentiles(
                list(residues_per_site.values())
            ),
            "single_residue_site_count": len(single_residue_sites),
            "single_residue_site_examples": single_residue_sites[:15],
        },
        "numbering_conventions": {
            "reference_chain_rows": {
                "resid_equals_auth": row_equal,
                "resid_differs_from_auth": row_different,
                "missing_auth": row_auth_none,
                "non_integer_auth": row_non_int_auth,
            },
            "per_site": dict(per_site_numbering),
            "non_integer_auth_examples": non_int_auth_examples,
        },
        "excluded_sites": excluded,
    }


def _print_summary(report: dict) -> None:
    scale = report["scale"]
    print("=== template library scale ===")
    for key, value in scale.items():
        print(f"  {key}: {value}")

    ec = report["ec_coverage"]
    print("=== EC coverage ===")
    print(f"  distinct EC numbers: {ec['distinct_ec_numbers']}")
    print(f"  EC class counts: {ec['ec_class_counts']}")
    print(
        "  templates per EC: "
        f"{ec['templates_per_ec_distribution']}"
    )
    print(
        f"  single-template EC numbers: {ec['single_template_ec_count']}"
    )

    composition = report["residue_composition"]
    print("=== catalytic residue composition ===")
    print(f"  residue codes: {composition['residue_code_counts']}")
    print(
        "  residues per site: "
        f"{composition['residues_per_site_distribution']}"
    )
    print(
        f"  single-residue sites: {composition['single_residue_site_count']}"
    )

    numbering = report["numbering_conventions"]
    print("=== numbering conventions ===")
    print(f"  chain rows: {numbering['reference_chain_rows']}")
    print(f"  per site: {numbering['per_site']}")
    if numbering["non_integer_auth_examples"]:
        print(
            "  non-integer auth examples: "
            f"{numbering['non_integer_auth_examples']}"
        )

    print("=== excluded sites ===")
    print(f"  {report['excluded_sites']}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mcsa_dir", help="directory with the M-CSA JSON files")
    parser.add_argument(
        "--report",
        help="optional path for the JSON report (default: stdout summary only)",
    )
    parser.add_argument(
        "--coverage",
        help="optional coverage_report.json from analyze_template_coverage.py",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mcsa_dir = Path(args.mcsa_dir)
    residues_path = mcsa_dir / _RESIDUES_NAME
    if not residues_path.exists():
        # The local probe directory keeps the un-prefixed name.
        residues_path = mcsa_dir / "residues.json"
    entries_path = mcsa_dir / _ENTRIES_NAME
    if not residues_path.exists() or not entries_path.exists():
        raise SystemExit(
            f"missing M-CSA source files in {mcsa_dir}; run the mcsa download"
        )

    residues = _load_json(residues_path)
    entries = _load_json(entries_path)
    coverage = _load_json(Path(args.coverage)) if args.coverage else None
    report = analyze(residues, entries, coverage)
    _print_summary(report)

    if args.report:
        report_path = Path(args.report)
        with report_path.open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(report, handle, ensure_ascii=False, indent=2)
        print(f"JSON report written to {report_path}")


if __name__ == "__main__":
    main()
