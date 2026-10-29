#!/usr/bin/env python
"""Add the eggNOG channel to an existing LTR dataset.

The LTR training set is built by ``build_ltr_trainset.py``, which runs every
search channel in one pass. When the eggNOG channel was skipped because
``emapper`` lives in a different conda environment, rerunning the whole build
would also rerun the expensive hhblits searches. This script runs only the
eggnog channel against the already-clustered representatives and merges the
new evidence into the existing ``search_evidence.json``, leaving the other
channels untouched.

The eggnog evidence is already expressed as ontology terms (GO/EC/Pfam/COG), so
it needs no further accession-to-term normalization; the per-term feature
builder groups it with the other channels at training time.

Example
-------
::

    python -u add_eggnog_channel.py \
        --reps     ltr_work/03_cluster/representatives.faa \
        --dbdir    dbs/dbdir \
        --evidence ltr_dataset/search_evidence.json \
        --emapper  /data/liushiwei/anaconda3/envs/eggnog/bin/emapper.py \
        --data-dir /db/student/metagenome/eggnog_db_v7 \
        --threads  96 \
        --out      ltr_dataset/search_evidence_eggnog.json
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pufprobe.logging_utils import setup_logging, get_logger  # noqa: E402
from pufprobe.models import to_dict  # noqa: E402
from pufprobe.ontology_names import load_pfam_names  # noqa: E402
from pufprobe.search import eggnog as eggnog_channel  # noqa: E402

logger = get_logger("add_eggnog_channel")


def main() -> int:
    setup_logging()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=Path, required=True,
                        help="Clustered representatives fasta.")
    parser.add_argument("--dbdir", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True,
                        help="Existing search_evidence.json (other channels).")
    parser.add_argument("--emapper", type=Path, required=True,
                        help="Absolute path to emapper.py in its conda env.")
    parser.add_argument("--data-dir", type=str, default="",
                        help="eggNOG database directory (emapper --data_dir).")
    parser.add_argument("--search-mode", type=str, default="diamond",
                        choices=["diamond", "mmseqs", "hmmer"])
    parser.add_argument("--threads", type=int, default=96)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    if not args.evidence.exists():
        raise SystemExit(f"Evidence file not found: {args.evidence}")
    if not args.emapper.exists():
        raise SystemExit(f"emapper not found: {args.emapper}")

    pfam_names = load_pfam_names(args.dbdir / "pfam_names.json")

    # Keep the eggnog raw files in a throwaway directory next to the dataset so
    # the working directory is not polluted.
    tmp_out = Path(tempfile.mkdtemp(
        prefix="add_eggnog_", dir=str(args.evidence.parent)))

    logger.info("Running eggnog on %s ...", args.reps)
    egg_ev, _ = eggnog_channel.run_eggnog_search(
        args.reps,
        args.dbdir,
        tmp_out,
        threads=args.threads,
        data_dir=args.data_dir or None,
        search_mode=args.search_mode,
        pfam_names=pfam_names,
        emapper_path=str(args.emapper),
    )

    with args.evidence.open("r", encoding="utf-8") as handle:
        merged: Dict[str, Any] = json.load(handle)

    n_added = 0
    for key, evs in egg_ev.items():
        merged.setdefault(key, []).extend(to_dict(ev) for ev in evs)
        n_added += len(evs)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(merged, handle)

    summary = {
        "out": str(args.out),
        "eggnog_proteins": len(egg_ev),
        "eggnog_evidence_added": n_added,
        "total_proteins": len(merged),
    }
    print(json.dumps(summary, indent=2))
    logger.info("eggnog merged: %d evidence entries across %d proteins.",
                n_added, len(egg_ev))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
