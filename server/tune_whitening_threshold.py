#!/usr/bin/env python
"""Choose a cosine floor for the whitened pLM neighbour search.

After whitening, unrelated proteins have a cosine near zero while same-Pfam
members have a positive cosine. This script samples query proteins, computes
their whitened cosine against the whole reference library, and reports the
cosine distribution separately for same-Pfam and other proteins. A good
floor sits just above the extreme tail of the unrelated distribution but
below the low tail of the same-Pfam distribution, so unrelated neighbours are
removed without discarding genuine family members.

The script is diagnostic and writes nothing.
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

logger = get_logger("tune_whitening_threshold")


def main() -> int:
    setup_logging()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-dir", type=Path, required=True)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--reps", type=Path, required=True)
    parser.add_argument("--sample", type=int, default=200,
                        help="Number of query proteins to sample.")
    parser.add_argument("--batch", type=int, default=50)
    args = parser.parse_args()

    lib = load_reference(args.ref_dir)
    if not lib.is_whitened:
        raise SystemExit("The reference library is not whitened.")

    rep_keys = [header.split()[0] for header, _ in read_fasta(args.reps)]
    raw_queries = load_query_vectors(args.embeddings, rep_keys)
    rep_keys = [key for key in rep_keys if key in raw_queries]

    # Keep only queries whose Pfam is represented in the library, otherwise
    # there is no positive same-Pfam set to compare.
    eligible = []
    for key in rep_keys:
        acc = key.split("@")[-1]
        pfams = (lib.terms.get(acc) or {}).get("Pfam", [])
        if pfams:
            eligible.append(key)
    rng = np.random.default_rng(42)
    rng.shuffle(eligible)
    chosen = eligible[: args.sample]
    logger.info("Sampling %d of %d eligible queries.", len(chosen), len(eligible))

    same_values = []
    other_values = []
    for start in range(0, len(chosen), args.batch):
        keys = chosen[start:start + args.batch]
        query_matrix = np.stack([raw_queries[key] for key in keys])
        query_matrix = apply_whitening(
            query_matrix, lib.whiten_mean, lib.whiten_w)
        sims = query_matrix @ lib.vectors.T
        for row, key in zip(sims, keys):
            acc = key.split("@")[-1]
            pfams = set((lib.terms.get(acc) or {}).get("Pfam", []))
            same_idx = [i for i, ref_acc in enumerate(lib.ids)
                        if pfams & set((lib.terms.get(ref_acc) or {}).get("Pfam", []))]
            mask = np.ones(len(lib.ids), bool)
            mask[same_idx] = False
            same_values.append(row[same_idx])
            other_values.append(row[mask])

    same = np.concatenate(same_values)
    other = np.concatenate(other_values)
    print("pairs: same-Pfam %d, other %d" % (same.size, other.size))
    for p in (50, 90, 99, 99.9):
        print("other cosine p%-5s: %.4f" % (p, np.percentile(other, p)))
    for p in (1, 5, 10, 50):
        print("same-Pfam cosine p%-5s: %.4f" % (p, np.percentile(same, p)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
