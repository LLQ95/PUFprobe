#!/usr/bin/env python
"""Add the pLM-neighbour channel to an existing LTR dataset.

The LTR training set is built by ``build_ltr_trainset.py`` with the sequence
and database channels. When the pLM reference library was not available at
build time, the trained model never sees the ``plm_neighbor`` feature and
cannot use it at inference. Rerunning the whole build would also rerun the
expensive hhblits searches. This script runs only the pLM-neighbour channel
against the already-clustered representatives and merges the new evidence
into the existing ``search_evidence.json``, leaving the other channels
untouched.

The hold-out proteins are themselves part of the reviewed Swiss-Prot library,
so the full reference library contains each query's own entry. The exact
self-hit is removed after the search (a query cannot label itself); other
close neighbours are kept, matching the leave-one-out logic of the training
set builder.

Example
-------
::

    python -u add_plm_channel.py \
        --reps       ltr_work/03_cluster/representatives.faa \
        --embeddings ltr_work/02_embed/embeddings.npz \
        --dbdir      dbs/dbdir \
        --ref-dir    dbs/dbdir/plm_ref \
        --evidence   ltr_dataset/search_evidence.json \
        --out        ltr_dataset/search_evidence_plm.json
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pufprobe.logging_utils import setup_logging, get_logger  # noqa: E402
from pufprobe.models import to_dict  # noqa: E402
from pufprobe.search.idmapping import accessions_from_target  # noqa: E402
from pufprobe.search import plm_neighbor  # noqa: E402

logger = get_logger("add_plm_channel")


def remove_self_hits(
    evidence: Dict[str, List[Dict[str, Any]]],
) -> int:
    """Drop entries whose target resolves to the query's own accession."""
    removed = 0
    for key, rows in evidence.items():
        query_acc = key.split("@")[-1]
        kept: List[Dict[str, Any]] = []
        for row in rows:
            targets = accessions_from_target(str(row.get("target_id", "")))
            if query_acc in targets:
                removed += 1
                continue
            kept.append(row)
        evidence[key] = kept
    return removed


def main() -> int:
    setup_logging()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=Path, required=True,
                        help="Clustered representatives fasta.")
    parser.add_argument("--embeddings", type=Path, required=True,
                        help="Query embeddings from the embed stage.")
    parser.add_argument("--dbdir", type=Path, required=True)
    parser.add_argument("--ref-dir", type=Path, required=True,
                        help="Annotated pLM reference library directory.")
    parser.add_argument("--evidence", type=Path, required=True,
                        help="Existing search_evidence.json (other channels).")
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--cos-floor", type=float, default=0.50)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    if not args.evidence.exists():
        raise SystemExit(f"Evidence file not found: {args.evidence}")
    if not args.ref_dir.exists():
        raise SystemExit(f"Reference library not found: {args.ref_dir}")

    tmp_out = Path(tempfile.mkdtemp(
        prefix="add_plm_", dir=str(args.evidence.parent)))

    logger.info("Running pLM-neighbour search on %s ...", args.reps)
    plm_ev, _ = plm_neighbor.run_plm_neighbor_search(
        args.reps,
        args.embeddings,
        args.dbdir,
        tmp_out,
        ref_dir=args.ref_dir,
        top_k=args.top_k,
        cos_floor=args.cos_floor,
    )

    plm_dict = {key: [to_dict(ev) for ev in evs] for key, evs in plm_ev.items()}
    n_self = remove_self_hits(plm_dict)
    plm_dict = {key: rows for key, rows in plm_dict.items() if rows}

    with args.evidence.open("r", encoding="utf-8") as handle:
        merged: Dict[str, Any] = json.load(handle)

    n_added = 0
    for key, rows in plm_dict.items():
        merged.setdefault(key, []).extend(rows)
        n_added += len(rows)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(merged, handle)

    summary = {
        "out": str(args.out),
        "plm_proteins": len(plm_dict),
        "self_hits_removed": n_self,
        "plm_evidence_added": n_added,
        "total_proteins": len(merged),
    }
    print(json.dumps(summary, indent=2))
    logger.info("plm merged: %d evidence entries across %d proteins (%d self removed).",
                n_added, len(plm_dict), n_self)
    return 0


if __name__ == "__main__":
    sys.exit(main())
