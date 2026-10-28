#!/usr/bin/env python
"""Build an annotated training set for the learning-to-rank fusion model.

A collection of proteins whose annotations are already known (the reviewed
Swiss-Prot hold-out produced by ``baseline_holdout.py`` or any FASTA of
curated proteins) is treated as if its function were unknown. The proteins are
embedded and clustered with the normal PUFprobe stages, the selected evidence
channels are run against the reference databases, and the self-hit of every
query is removed so a protein cannot label itself. The resulting evidence is
then paired with the curated truth, giving the exact feature distribution the
ranker sees at inference time.

The truth is read from a Swiss-Prot cross-reference index
(``swissprot_xref.json``, built from the same reviewed release) and matched
by accession. A temporal split is downstream: pass the dates file (accession
-> ISO date) through to ``train_ltr.py``.

Example
-------
Build a training set from a 5,000-protein reviewed hold-out using the
sequence-profile and database channels, then train LightGBM::

    python -u build_ltr_trainset.py \
        --fasta  holdout/holdout.faa \
        --xref   dbdir/swissprot_xref.json \
        --dbdir  dbdir --channels hhblits,eggnog,mmseqs \
        --model  esm2_t33_650M_UR50D --threads 96 \
        --workdir ltr_work --out ltr_dataset

    python -u train_ltr.py \
        --evidence ltr_dataset/search_evidence.json \
        --labels   ltr_dataset/labels.json \
        --dates    ltr_dataset/dates.json \
        --go-obo   dbdir/go-basic.obo \
        --model-out dbdir/ltr_fusion.pkl --backend auto
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pufprobe.config import load_config  # noqa: E402
from pufprobe.io_utils import read_fasta, write_json  # noqa: E402
from pufprobe.logging_utils import setup_logging, get_logger  # noqa: E402
from pufprobe.models import (  # noqa: E402
    CDSRecord,
    GeneStatus,
    to_dict,
)
from pufprobe.search.idmapping import (  # noqa: E402
    accessions_from_target,
    load_xref,
)
from pufprobe.workdir import WorkDir  # noqa: E402

logger = get_logger("build_ltr_trainset")


def parse_accession(header: str) -> str:
    """Extract the primary accession from a UniProt or plain FASTA header."""
    header = header.strip()
    if header.startswith("sp|") or header.startswith("tr|"):
        parts = header.split("|")
        if len(parts) >= 2 and parts[1]:
            return parts[1].strip()
    return header.split()[0].strip()


def read_holdout(fasta: Path) -> List[Tuple[str, str]]:
    proteins: List[Tuple[str, str]] = []
    seen = set()
    for header, sequence in read_fasta(fasta):
        accession = parse_accession(header)
        if not accession or accession in seen or not sequence:
            continue
        seen.add(accession)
        proteins.append((accession, sequence))
    return proteins


def write_prepare(workdir: WorkDir, proteins: List[Tuple[str, str]]) -> None:
    """Create one conserved-hypothetical CDS record per hold-out protein."""
    records: List[Dict[str, Any]] = []
    for accession, sequence in proteins:
        record = CDSRecord(
            sample_id="holdout",
            locus_tag=accession,
            protein_id=accession,
            sequence_aa=sequence,
            start=1,
            end=len(sequence) * 3,
            strand="+",
            product="hypothetical protein",
            status=GeneStatus.CONSERVED_HYPOTHETICAL,
        )
        records.append(to_dict(record))
    write_json(workdir.path("prepare", "cds_records.json"), records)
    workdir.mark_complete("prepare", payload={"n_records": len(records)})


def remove_self_hits(
    evidence: Dict[str, List[Dict[str, Any]]],
) -> Tuple[Dict[str, List[Dict[str, Any]]], int]:
    """Drop any evidence whose target resolves to the query's own accession."""
    cleaned: Dict[str, List[Dict[str, Any]]] = {}
    removed = 0
    for key, rows in evidence.items():
        query_acc = key.split("@")[-1]
        kept = []
        for row in rows:
            target_accs = accessions_from_target(str(row.get("target_id", "")))
            if query_acc in target_accs:
                removed += 1
                continue
            kept.append(row)
        cleaned[key] = kept
    return cleaned, removed


def load_raw_evidence(path: Path) -> Dict[str, List[Dict[str, Any]]]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"expected a keyed evidence file: {path}")
    return data


def build_labels(
    proteins: List[Tuple[str, str]],
    xref: Dict[str, Dict[str, List[str]]],
) -> Dict[str, Dict[str, List[str]]]:
    """Map every protein key to its curated term lists."""
    labels: Dict[str, Dict[str, List[str]]] = {}
    for accession, _ in proteins:
        terms = xref.get(accession)
        if not terms:
            continue
        labels[f"holdout@{accession}"] = {
            "GO": sorted(terms.get("GO", [])),
            "EC": sorted(terms.get("EC", [])),
            "Pfam": sorted(terms.get("Pfam", [])),
            "COG": sorted(terms.get("COG", [])),
        }
    return labels


def align_dates(
    proteins: List[Tuple[str, str]], dates_in: Dict[str, str]
) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for accession, _ in proteins:
        date = dates_in.get(accession)
        if date:
            out[f"holdout@{accession}"] = date
    return out


def main() -> int:
    setup_logging()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True,
                        help="Hold-out proteins (UniProt or plain FASTA).")
    parser.add_argument("--xref", type=Path, required=True,
                        help="swissprot_xref.json with curated GO/EC/Pfam.")
    parser.add_argument("--dbdir", type=Path, required=True,
                        help="Database directory for the search channels.")
    parser.add_argument("--workdir", type=Path, required=True,
                        help="Scratch workdir for embed/cluster/search.")
    parser.add_argument("--out", type=Path, required=True,
                        help="Output dataset directory.")
    parser.add_argument("--channels", default="hhblits,eggnog,mmseqs",
                        help="Comma-separated search channels to run.")
    parser.add_argument("--model", default="esm2_t33_650M_UR50D",
                        help="Embedding model; should match the deployed model.")
    parser.add_argument("--threads", type=int, default=96)
    parser.add_argument("--hhblits-jobs", type=int, default=0,
                        help="Concurrent independent hhblits searches. 0 runs "
                             "one job per thread (recommended on a many-core host).")
    parser.add_argument("--hhblits-timeout", type=int, default=300,
                        help="Per-query hhblits wall-clock limit in seconds.")
    parser.add_argument("--dates", type=Path, default=None,
                        help="accession -> ISO date, copied for temporal split.")
    parser.add_argument("--cpu", action="store_true",
                        help="Force CPU embedding (slow).")
    args = parser.parse_args()

    proteins = read_holdout(args.fasta)
    if not proteins:
        raise SystemExit(f"No proteins parsed from {args.fasta}")
    logger.info("Hold-out proteins: %d.", len(proteins))

    xref = load_xref(args.xref)
    if not xref:
        raise SystemExit(f"Could not load cross-reference index from {args.xref}")
    labels = build_labels(proteins, xref)
    if not labels:
        raise SystemExit("No hold-out protein has curated terms in the xref.")
    logger.info("Proteins with curated labels: %d.", len(labels))

    workdir = WorkDir(args.workdir, cpu=args.cpu, threads=args.threads)
    write_prepare(workdir, proteins)

    config = load_config(None)
    # Keep families close to singletons so each curated protein contributes a
    # ranking group; 0.95 identity only merges near-identical entries.
    config["cluster"]["identity"] = 0.95

    # Parallelize the independent hhblits queries. Default to one single-thread
    # search per thread, which uses a many-core host far better than one
    # -cpu=$threads invocation (a single query scales sublinearly).
    hh_jobs = args.hhblits_jobs if args.hhblits_jobs > 0 else max(1, args.threads)
    config["search"]["hhblits_jobs"] = hh_jobs
    config["search"]["hhblits_query_timeout"] = args.hhblits_timeout

    from pufprobe.embed.pipeline import run as embed_run
    if workdir.is_complete("embed") and workdir.path("cluster", "families.json").exists():
        logger.info("embed/cluster already complete in %s; reusing it.", workdir.root)
    else:
        embed_run(workdir, config, model=args.model)

    channels = [c.strip().lower() for c in args.channels.split(",") if c.strip()]
    from pufprobe.search.pipeline import run as search_run
    if workdir.is_complete("search") and workdir.path("search", "search_evidence.json").exists():
        logger.info("search already complete; reusing its evidence.")
    else:
        search_run(workdir, config, channels=channels, dbdir=args.dbdir)

    raw_path = workdir.path("search", "search_evidence.json")
    evidence = load_raw_evidence(raw_path)
    evidence, n_removed = remove_self_hits(evidence)
    logger.info("Self-hit removal dropped %d evidence entries.", n_removed)

    # Restrict the dataset to proteins that have both labels and (any)
    # evidence after the self-hit filter.
    keys_with_evidence = {k for k, rows in evidence.items() if rows}
    labels = {k: v for k, v in labels.items() if k in keys_with_evidence}
    evidence = {k: evidence[k] for k in labels}
    if not labels:
        raise SystemExit("No usable training rows after self-hit removal.")

    args.out.mkdir(parents=True, exist_ok=True)
    write_json(args.out / "search_evidence.json", evidence)
    write_json(args.out / "labels.json", labels)
    if args.dates:
        with args.dates.open("r", encoding="utf-8") as handle:
            dates_in = json.load(handle)
        dates = align_dates(proteins, dates_in)
        write_json(args.out / "dates.json", dates)
        logger.info("Dates exported for %d proteins.", len(dates))

    n_terms = sum(len(r) for r in evidence.values())
    logger.info("Training set written to %s: %d proteins, %d evidence rows.",
                args.out, len(labels), n_terms)
    print(json.dumps({
        "dataset": str(args.out),
        "n_proteins": len(labels),
        "n_evidence_rows": n_terms,
        "self_hits_removed": n_removed,
        "channels": channels,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
