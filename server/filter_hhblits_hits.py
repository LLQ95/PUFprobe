#!/usr/bin/env python
"""Filter hhblits evidence by posterior probability without rerunning hhblits.

Running hhblits with the ``-p`` gate drops marginal hits at the source. When a
search_evidence.json was produced before that gate existed, this script applies
the same filter to the stored evidence: every hhblits entry (the raw profile
hit and the ontology terms normalised from it) whose posterior probability is
below the threshold is removed. Other channels, including mmseqs entries that
share the REMOTE_HOMOLOGY channel, are left untouched.

The hhblits entries are identified by ``details.database == "hhblits"``; the
normalised copies inherit that marker and the original ``details.prob``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple


def filter_evidence(
    data: Dict[str, List[Dict[str, Any]]],
    min_prob: float,
) -> Tuple[Dict[str, List[Dict[str, Any]]], int]:
    filtered: Dict[str, List[Dict[str, Any]]] = {}
    removed = 0
    for key, rows in data.items():
        kept: List[Dict[str, Any]] = []
        for entry in rows:
            details = entry.get("details") or {}
            if details.get("database") == "hhblits":
                prob = details.get("prob")
                if prob is not None and float(prob) < min_prob:
                    removed += 1
                    continue
            kept.append(entry)
        filtered[key] = kept
    return filtered, removed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True,
                        help="Input search_evidence.json.")
    parser.add_argument("--out", type=Path, required=True,
                        help="Filtered output path.")
    parser.add_argument("--min-prob", type=float, default=90.0,
                        help="Minimum hhblits posterior probability (percent).")
    args = parser.parse_args()

    with args.evidence.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    filtered, removed = filter_evidence(data, args.min_prob)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(filtered, handle)

    before = sum(len(r) for r in data.values())
    after = sum(len(r) for r in filtered.values())
    print(json.dumps({
        "input": str(args.evidence),
        "output": str(args.out),
        "min_prob": args.min_prob,
        "rows_before": before,
        "rows_removed": removed,
        "rows_after": after,
        "proteins": len(filtered),
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
