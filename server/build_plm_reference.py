#!/usr/bin/env python
"""Build the annotated pLM reference library for the plm_neighbor channel.

The reference library lets the search stage transfer GO/EC/Pfam annotations to
representative proteins by nearest-neighbour search in ESM-2 embedding space.
Reference proteins are taken from the reviewed Swiss-Prot database and filtered
to entries that carry at least one ontology term. The embeddings are produced
with the same ESM-2 checkpoint and residue-mean pooling as the query embed
stage, so the model name MUST match ``embed.default_model`` (or the model the
query workdir actually used).

Preferred input is the Swiss-Prot flat file, which carries the sequence and
its Pfam/EC/GO cross-references in one place:

    python build_plm_reference.py \
        --dat  /db/.../uniprot_sprot.dat.gz \
        --out  /data/.../dbdir/plm_ref \
        --model esm2_t33_650M_UR50D --batch-size 16

A FASTA plus a prebuilt cross-reference index is also accepted:

    python build_plm_reference.py \
        --fasta /db/.../uniprot_sprot.fasta.gz \
        --xref  /data/.../dbdir/swissprot_xref.json \
        --out   /data/.../dbdir/plm_ref

The run is sharded and resumable: each shard is written under
``<out>/_shards/`` and re-running skips completed shards. After all shards
finish the final library files are written and the shards are removed unless
``--keep-shards`` is given.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np

from pufprobe.embed.esm import ModelUnavailable, embed_items, load_esm_model
from pufprobe.io_utils import read_fasta
from pufprobe.logging_utils import get_logger
from pufprobe.search.idmapping import (
    _open_text,
    iter_sprot_records,
    load_xref,
)

logger = get_logger("build_plm_reference")

VECTORS = "plm_ref_vectors.npy"
IDS = "plm_ref_ids.json"
TERMS = "plm_ref_terms.json"
MANIFEST = "plm_ref_manifest.json"


def _iter_fasta_accessions(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (primary accession, sequence) pairs from a UniProt FASTA."""
    for header, sequence in read_fasta(path):
        if header.startswith("sp|") or header.startswith("tr|"):
            fields = header.split("|")
            if len(fields) >= 3:
                yield fields[1].strip(), sequence
                continue
        yield header.split()[0], sequence


def _filter(
    records: Iterator[Tuple[str, str, Dict[str, List[str]]]],
    min_len: int,
    max_len: int,
    limit: int,
) -> List[Tuple[str, str, Dict[str, List[str]]]]:
    """Sort and de-duplicate annotated records within the length window."""
    selected: List[Tuple[str, str, Dict[str, List[str]]]] = []
    seen = set()
    for accession, sequence, terms in records:
        if accession in seen or not terms:
            continue
        if not (min_len <= len(sequence) <= max_len):
            continue
        seen.add(accession)
        selected.append((accession, sequence, terms))
    selected.sort(key=lambda item: item[0])
    if limit and limit > 0:
        selected = selected[:limit]
    return selected


def records_from_dat(dat_path: Path) -> Iterator[Tuple[str, str, Dict[str, List[str]]]]:
    with _open_text(dat_path) as handle:
        for record in iter_sprot_records(handle):
            yield record["accession"], record["sequence"], record["terms"]


def records_from_fasta(
    fasta: Path, xref: Dict[str, Dict[str, List[str]]]
) -> Iterator[Tuple[str, str, Dict[str, List[str]]]]:
    for accession, sequence in _iter_fasta_accessions(fasta):
        terms = xref.get(accession)
        if terms:
            yield accession, sequence, terms


def _normalise(matrix: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-8)
    return matrix / norms


def build(
    out_dir: Path,
    selected: List[Tuple[str, str, Dict[str, List[str]]]],
    source: str,
    model_name: str,
    batch_size: int,
    shard_size: int,
    cpu: bool,
    keep_shards: bool,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    shard_dir = out_dir / "_shards"
    shard_dir.mkdir(parents=True, exist_ok=True)

    n_total = len(selected)
    n_shards = math.ceil(n_total / shard_size)
    logger.info(
        "Selected %d annotated reference proteins -> %d shard(s) of up to %d (model=%s).",
        n_total, n_shards, shard_size, model_name,
    )

    # Load the checkpoint once and reuse it across shards.
    model, alphabet, batch_converter, device, num_layers = load_esm_model(
        model_name, cpu=cpu
    )

    for shard_index in range(n_shards):
        chunk = selected[shard_index * shard_size : (shard_index + 1) * shard_size]
        shard_path = shard_dir / f"shard_{shard_index:05d}.npz"
        if shard_path.exists():
            logger.info("Shard %d/%d already present; skipping.", shard_index + 1, n_shards)
            continue
        items = [(acc, seq) for acc, seq, _ in chunk]
        try:
            vectors = embed_items(
                model, alphabet, batch_converter, num_layers, device,
                items, batch_size=batch_size,
            )
        except ModelUnavailable:
            raise
        ids = [acc for acc, _, _ in chunk if acc in vectors]
        if not ids:
            raise SystemExit(f"Shard {shard_index} produced no embeddings; aborting.")
        matrix = np.stack([vectors[acc] for acc in ids]).astype(np.float32)
        matrix = _normalise(matrix)
        np.savez(
            shard_path,
            vectors=matrix.astype(np.float16),
            ids=np.asarray(ids, dtype=object),
        )
        logger.info("Shard %d/%d written (%d proteins).", shard_index + 1, n_shards, len(ids))

    # Merge shards into the final library.
    all_ids: List[str] = []
    all_vectors: List[np.ndarray] = []
    for shard_index in range(n_shards):
        shard_path = shard_dir / f"shard_{shard_index:05d}.npz"
        with np.load(shard_path, allow_pickle=True) as data:
            shard_ids = [str(x) for x in data["ids"].tolist()]
            shard_vectors = data["vectors"].astype(np.float32)
        all_ids.extend(shard_ids)
        all_vectors.append(shard_vectors)

    terms_by_acc = {acc: terms for acc, _, terms in selected}
    vectors = np.concatenate(all_vectors, axis=0).astype(np.float16)
    terms = {acc: terms_by_acc[acc] for acc in all_ids if acc in terms_by_acc}

    np.save(out_dir / VECTORS, vectors)
    with (out_dir / IDS).open("w", encoding="utf-8") as handle:
        json.dump(all_ids, handle, separators=(",", ":"))
    with (out_dir / TERMS).open("w", encoding="utf-8") as handle:
        json.dump(terms, handle, separators=(",", ":"))
    manifest = {
        "model": model_name,
        "pooling": "residue_mean",
        "normalized": True,
        "dim": int(vectors.shape[1]),
        "n": int(vectors.shape[0]),
        "n_with_terms": int(len(terms)),
        "source": source,
    }
    with (out_dir / MANIFEST).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    if not keep_shards:
        for shard_index in range(n_shards):
            (shard_dir / f"shard_{shard_index:05d}.npz").unlink(missing_ok=True)
        shard_dir.rmdir()

    logger.info(
        "Reference library written to %s: %d proteins x %d dims (%.1f GB float16).",
        out_dir, vectors.shape[0], vectors.shape[1],
        vectors.nbytes / 1e9,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dat", type=Path, default=None,
                        help="UniProt Swiss-Prot flat file (.dat or .dat.gz). Preferred.")
    parser.add_argument("--fasta", type=Path, default=None,
                        help="UniProt Swiss-Prot FASTA (used with --xref).")
    parser.add_argument("--xref", type=Path, default=None,
                        help="swissprot_xref.json (used with --fasta).")
    parser.add_argument("--out", type=Path, required=True,
                        help="Output directory, e.g. <dbdir>/plm_ref.")
    parser.add_argument("--model", default="esm2_t33_650M_UR50D",
                        help="ESM checkpoint; must match the query embed model.")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--shard-size", type=int, default=20000,
                        help="Proteins per resumable shard.")
    parser.add_argument("--cpu", action="store_true",
                        help="Force CPU inference (slow).")
    parser.add_argument("--min-len", type=int, default=30)
    parser.add_argument("--max-len", type=int, default=1000,
                        help="Skip proteins longer than the ESM context window.")
    parser.add_argument("--limit", type=int, default=0,
                        help="Use only the first N selected proteins (for testing).")
    parser.add_argument("--keep-shards", action="store_true",
                        help="Keep intermediate shard files after merging.")
    args = parser.parse_args()

    if args.dat is not None:
        if not args.dat.exists():
            raise SystemExit(f"flat file not found: {args.dat}")
        records_iter = records_from_dat(args.dat)
        source = str(args.dat)
    elif args.fasta is not None and args.xref is not None:
        if not args.fasta.exists():
            raise SystemExit(f"FASTA not found: {args.fasta}")
        if not args.xref.exists():
            raise SystemExit(f"xref index not found: {args.xref}")
        xref = load_xref(args.xref)
        if not xref:
            raise SystemExit(f"No cross-reference index loaded from {args.xref}")
        records_iter = records_from_fasta(args.fasta, xref)
        source = f"{args.fasta} + {args.xref}"
    else:
        raise SystemExit("Provide either --dat, or both --fasta and --xref.")

    selected = _filter(records_iter, args.min_len, args.max_len, args.limit)
    if not selected:
        raise SystemExit("No annotated reference proteins selected; check inputs.")
    build(
        out_dir=args.out,
        selected=selected,
        source=source,
        model_name=args.model,
        batch_size=args.batch_size,
        shard_size=args.shard_size,
        cpu=args.cpu,
        keep_shards=args.keep_shards,
    )


if __name__ == "__main__":
    main()
