#!/usr/bin/env python3
"""Summarize downloaded pathogen genomes.

For every sample directory under the genomes root this reports the number of
sequence contigs, protein records, GFF3 CDS features, and proteins whose RefSeq
product description marks them as hypothetical. The totals quantify the scale
of the hypothetical-protein set that feeds the HPDB build.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

HYPOTHETICAL_MARKERS = (
    "hypothetical protein",
    "uncharacterized protein",
    "putative protein",
)


def count_fasta_records(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def count_cds(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 3 and fields[2] == "CDS":
                count += 1
    return count


def count_hypothetical(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line.lstrip(">").lower()
            if any(marker in header for marker in HYPOTHETICAL_MARKERS):
                count += 1
    return count


def inspect_sample(sample_dir: Path) -> Dict[str, object]:
    return {
        "sample": sample_dir.name,
        "contigs": count_fasta_records(sample_dir / "genomic.fna"),
        "proteins": count_fasta_records(sample_dir / "protein.faa"),
        "cds_features": count_cds(sample_dir / "genomic.gff"),
        "hypothetical": count_hypothetical(sample_dir / "protein.faa"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--genomes",
        required=True,
        help="Path to the genomes directory holding sample folders.",
    )
    parser.add_argument(
        "--json",
        dest="json_path",
        default="",
        help="Optional path to write the per-sample summary as JSON.",
    )
    args = parser.parse_args()

    root = Path(args.genomes)
    sample_dirs = sorted(path for path in root.iterdir() if path.is_dir())
    rows: List[Dict[str, object]] = [
        inspect_sample(sample_dir) for sample_dir in sample_dirs
    ]

    header = (
        f"{'sample':<48}{'contigs':>8}{'proteins':>10}"
        f"{'cds':>8}{'hypo':>8}{'hypo%':>8}"
    )
    print(header)
    for row in rows:
        proteins = int(row["proteins"])
        fraction = (int(row["hypothetical"]) / proteins * 100.0) if proteins else 0.0
        print(
            f"{str(row['sample']):<48}{int(row['contigs']):>8}"
            f"{proteins:>10}{int(row['cds_features']):>8}"
            f"{int(row['hypothetical']):>8}{fraction:>7.1f}%"
        )

    totals = {
        "samples": len(rows),
        "contigs": sum(int(row["contigs"]) for row in rows),
        "proteins": sum(int(row["proteins"]) for row in rows),
        "cds_features": sum(int(row["cds_features"]) for row in rows),
        "hypothetical": sum(int(row["hypothetical"]) for row in rows),
    }
    print("-" * len(header))
    print(
        f"{'TOTAL':<48}{totals['contigs']:>8}{totals['proteins']:>10}"
        f"{totals['cds_features']:>8}{totals['hypothetical']:>8}"
    )

    if args.json_path:
        payload = {"per_sample": rows, "totals": totals}
        out = Path(args.json_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
