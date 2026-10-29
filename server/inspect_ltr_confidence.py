#!/usr/bin/env python
"""Inspect the evidence strength behind LTR confidence tiers.

The LTR ranker replaces the rule-based single-group cap with a continuous
score. To confirm that a high LTR confidence is backed by genuinely strong
evidence rather than a model artifact, this script joins the LTR predictions
with the raw search evidence and reports the e-value, coverage and channel
score distribution for each confidence tier.

GO ancestor terms have no direct evidence (they are produced by true-path
propagation), so they are reported separately rather than counted as weak.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple


def classify_term(term_id: str) -> str:
    if term_id.startswith("GO:"):
        return "GO"
    if term_id.startswith("PF"):
        return "Pfam"
    if term_id.startswith("COG") or term_id.startswith("KOG"):
        return "COG"
    parts = term_id.split(".")
    if len(parts) == 2 and all(p.isdigit() for p in parts):
        return "EC"
    return "COG"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True,
                        help="Raw search_evidence.json before GO propagation.")
    parser.add_argument("--threshold", type=float, default=0.80)
    args = parser.parse_args()

    with args.predictions.open("r", encoding="utf-8") as handle:
        predictions = json.load(handle)
    with args.evidence.open("r", encoding="utf-8") as handle:
        evidence: Dict[str, List[dict]] = json.load(handle)

    # Tier -> list of (protein, term, confidence)
    tiers: Dict[str, List[Tuple[str, str, str, float]]] = {
        "high": [], "mid": [], "low": []
    }
    for prot in predictions:
        pkey = f"{prot.get('sample_id', '')}@{prot.get('locus_tag', '')}"
        for term in prot.get("terms", []):
            conf = float(term["confidence"])
            tid = term["term_id"]
            ttype = classify_term(tid)
            if conf >= args.threshold:
                tiers["high"].append((pkey, ttype, tid, conf))
            elif conf >= 0.60:
                tiers["mid"].append((pkey, ttype, tid, conf))
            else:
                tiers["low"].append((pkey, ttype, tid, conf))

    def summarize(pkey: str, ttype: str, tid: str) -> dict:
        """Return the strongest direct evidence for one term."""
        rows = evidence.get(pkey, [])
        direct = [e for e in rows if e.get("target_id") == tid]
        if not direct:
            return {"kind": "indirect"}
        best = max(direct, key=lambda e: e.get("score", 0.0))
        return {
            "kind": "direct",
            "evalue": best.get("evalue"),
            "coverage": best.get("coverage"),
            "score": best.get("score", 0.0),
        }

    for tier_name in ("high", "mid", "low"):
        items = tiers[tier_name]
        direct = [summarize(*item[:3]) for item in items]
        direct = [d for d in direct if d["kind"] == "direct"]
        indirect = len(items) - len(direct)
        print(f"--- {tier_name} tier: {len(items)} terms "
              f"({indirect} GO-ancestor/indirect) ---")
        if not direct:
            continue
        scores = [d["score"] for d in direct]
        coverages = [d["coverage"] for d in direct if d["coverage"] is not None]
        evalues = [d["evalue"] for d in direct if d["evalue"] is not None]
        print(f"  channel score: median {statistics.median(scores):.2f}, "
              f"min {min(scores):.2f}")
        if coverages:
            print(f"  coverage: median {statistics.median(coverages):.2f}, "
                  f"min {min(coverages):.2f}")
        if evalues:
            logev = [math.log10(max(e, 1e-300)) for e in evalues]
            print(f"  log10 evalue: median {statistics.median(logev):.1f}, "
                  f"weakest {max(logev):.1f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
