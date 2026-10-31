#!/usr/bin/env python
"""Evaluate the term-level precision of the pLM-neighbour channel.

The neighbour search can produce many term candidates, but only the terms
that match the hidden annotations are useful. This script compares every
pLM-neighbour term against the labels (with GO true-path closure) and
reports the overall precision, the precision of pLM-only terms (terms not
supported by any other channel), and a breakdown by the number of voting
neighbours.

The label processing mirrors train_ltr.py so the comparison uses the same
aligned keys and closed GO labels. Nothing is written.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from train_ltr import (  # noqa: E402
    align_to_keys,
    normalize_labels,
    propagate_go_labels,
)
from pufprobe.logging_utils import setup_logging, get_logger  # noqa: E402

logger = get_logger("eval_plm_precision")


def kind_of(target_id: str) -> Optional[str]:
    if target_id.startswith("GO:"):
        return "GO"
    if target_id.startswith("EC:"):
        return "EC"
    if target_id.startswith("PF"):
        return "Pfam"
    return None


def voter_bucket(n_voters: int) -> str:
    if n_voters <= 1:
        return "1"
    if n_voters <= 3:
        return "2-3"
    return "4+"


def main() -> int:
    setup_logging()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--go-obo", type=Path, default=None)
    args = parser.parse_args()

    with args.evidence.open("r", encoding="utf-8") as handle:
        evidence = json.load(handle)
    with args.labels.open("r", encoding="utf-8") as handle:
        raw_labels = json.load(handle)

    labels = normalize_labels(raw_labels)
    labels = align_to_keys(labels, sorted(evidence))
    labels = {k: v for k, v in labels.items() if v}
    propagate_go_labels(labels, args.go_obo)

    # Counters: (subset, voter bucket) -> [correct, total]
    counters: Dict[Tuple[str, str], List[int]] = {}

    def record(subset: str, n_voters: int, correct: bool) -> None:
        bucket = voter_bucket(n_voters)
        slot = counters.setdefault((subset, bucket), [0, 0])
        slot[1] += 1
        if correct:
            slot[0] += 1

    for key, rows in evidence.items():
        truth = labels.get(key)
        if not truth:
            continue
        true_terms: Dict[str, Set[str]] = {
            kind: set(term_list) for kind, term_list in truth.items()}
        other_targets = {
            row.get("target_id") for row in rows
            if row.get("channel") != "plm_neighbor"}

        for row in rows:
            if row.get("channel") != "plm_neighbor":
                continue
            target = str(row.get("target_id", ""))
            kind = kind_of(target)
            if kind is None:
                continue
            n_voters = int((row.get("details") or {}).get("n_voters", 1))
            correct = target in true_terms.get(kind, set())
            record("all_plm", n_voters, correct)
            if target not in other_targets:
                record("plm_only", n_voters, correct)

    for subset in ("all_plm", "plm_only"):
        print("=== %s ===" % subset)
        subset_total = [0, 0]
        for bucket in ("1", "2-3", "4+"):
            correct, total = counters.get((subset, bucket), [0, 0])
            precision = correct / total if total else 0.0
            print("voters %-3s: precision %.4f (%d/%d)" % (
                bucket, precision, correct, total))
            subset_total[0] += correct
            subset_total[1] += total
        overall = subset_total[0] / subset_total[1] if subset_total[1] else 0.0
        print("overall: precision %.4f (%d/%d)" % (
            overall, subset_total[0], subset_total[1]))
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
