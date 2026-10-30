#!/usr/bin/env python
"""Evaluate candidate cosine floors on the whitened top-k neighbours.

The cosine distributions alone do not show how a floor behaves inside the
actual top-k neighbour list. For each sampled query this script takes the
whitened top-k neighbours and, for every candidate floor, reports the mean
number of retained neighbours, the fraction of retained neighbours that
share a Pfam family with the query (neighbour precision), and the fraction
of queries that keep at least one same-Pfam neighbour (family coverage).

A usable floor keeps neighbour precision high without emptying the result
for small or weakly separated families. Nothing is written.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pufprobe.embed.whitening import apply_whitening  # noqa: E402
from pufprobe.io_utils import read_fasta  # noqa: E402
from pufprobe.logging_utils import setup_logging, get_logger  # noqa: E402
from pufprobe.search.plm_neighbor import (  # noqa: E402
    load_reference,
    load_query_vectors,
)

logger = get_logger("eval_whitening_floor")


def main() -> int:
    setup_logging()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-dir", type=Path, required=True)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--reps", type=Path, required=True)
    parser.add_argument("--sample", type=int, default=300)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--floors", type=float, nargs="+",
                        default=[0.05, 0.08, 0.10, 0.12])
    args = parser.parse_args()

    lib = load_reference(args.ref_dir)
    if not lib.is_whitened:
        raise SystemExit("The reference library is not whitened.")

    pfam_to_idx: dict = {}
    for i, ref_acc in enumerate(lib.ids):
        for pf in (lib.terms.get(ref_acc) or {}).get("Pfam", []):
            pfam_to_idx.setdefault(pf, []).append(i)

    rep_keys = [header.split()[0] for header, _ in read_fasta(args.reps)]
    raw_queries = load_query_vectors(args.embeddings, rep_keys)
    eligible = []
    for key in rep_keys:
        if key not in raw_queries:
            continue
        acc = key.split("@")[-1]
        if (lib.terms.get(acc) or {}).get("Pfam", []):
            eligible.append(key)

    rng = np.random.default_rng(42)
    rng.shuffle(eligible)
    chosen = eligible[: args.sample]
    logger.info("Sampling %d queries.", len(chosen))

    stats = {floor: {"n": [], "correct": 0, "total": 0, "covered": 0}
             for floor in args.floors}
    for start in range(0, len(chosen), 50):
        keys = chosen[start:start + 50]
        query_matrix = apply_whitening(
            np.stack([raw_queries[key] for key in keys]),
            lib.whiten_mean, lib.whiten_w)
        sims = query_matrix @ lib.vectors.T
        top_idx = np.argpartition(-sims, args.top_k, axis=1)[:, : args.top_k]
        for key, row_idx, row_sim in zip(keys, top_idx, sims):
            acc = key.split("@")[-1]
            pfams = set((lib.terms.get(acc) or {}).get("Pfam", []))
            same_set = {i for pf in pfams for i in pfam_to_idx.get(pf, [])}
            for floor in args.floors:
                kept = [int(i) for i in row_idx if row_sim[int(i)] >= floor]
                bucket = stats[floor]
                bucket["n"].append(len(kept))
                correct = sum(1 for i in kept if i in same_set)
                bucket["correct"] += correct
                bucket["total"] += len(kept)
                if any(i in same_set for i in kept):
                    bucket["covered"] += 1

    print("floor\tmean_n\tneigh_precision\tfamily_coverage")
    for floor in args.floors:
        bucket = stats[floor]
        precision = bucket["correct"] / bucket["total"] if bucket["total"] else 0.0
        coverage = bucket["covered"] / len(chosen)
        print("%.2f\t%.2f\t%.4f\t\t%.4f" % (
            floor, np.mean(bucket["n"]), precision, coverage))
    return 0


if __name__ == "__main__":
    sys.exit(main())
