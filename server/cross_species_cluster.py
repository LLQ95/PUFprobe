#!/usr/bin/env python3
"""Cluster cross-species hypothetical proteins into families and embed representatives.

The cross-species dataset produced by fetch_reference_proteomes.py is a single
TSV of uncharacterized proteins. This script turns it into a family-level view:

1. Write every protein to a FASTA file keyed by UniProt accession.
2. Cluster the sequences with MMseqs2 ``linclust`` (linear time, suitable for
   hundreds of thousands of proteins), using identity and coverage cutoffs.
3. Collapse each cluster into a family and choose the longest member as the
   representative. Family spread across phyla and organisms is recorded.
4. Optionally compute ESM-2 mean-pooled embeddings for the representatives in
   shards, so a long GPU run can be interrupted and resumed without redoing the
   finished shards.

The MMseqs2 step and every embedding shard are idempotent; existing outputs are
reused on re-run.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

TSV_FIELDS = ["Entry", "Protein names", "Sequence", "Length", "Annotation",
              "Protein existence", "Organism (ID)", "Organism", "Proteome",
              "Phylum"]

EMBED_SHARD = 500


@dataclass
class Protein:
    """Minimal record exposing the attributes ESM embedding expects."""

    key: str
    sequence_aa: str
    phylum: str
    organism: str
    proteome: str
    name: str
    length: int


def write_fasta(path: Path, items: List[tuple]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for key, sequence in items:
            handle.write(f">{key}\n{sequence}\n")


def load_proteins(input_tsv: Path) -> Dict[str, Protein]:
    proteins: Dict[str, Protein] = {}
    with input_tsv.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row.get("Entry", "").strip()
            sequence = row.get("Sequence", "").strip()
            if not accession or not sequence:
                continue
            proteins[accession] = Protein(
                key=accession,
                sequence_aa=sequence,
                phylum=row.get("Phylum", "").strip(),
                organism=row.get("Organism", "").strip(),
                proteome=row.get("Proteome", "").strip(),
                name=row.get("Protein names", "").strip(),
                length=int(row.get("Length") or len(sequence)),
            )
    return proteins


def run_mmseqs(cmd: List[str]) -> None:
    completed = subprocess.run(cmd, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(cmd)}\n"
            f"{completed.stdout}\n{completed.stderr}"
        )


def linclust(
    fasta: Path,
    work_dir: Path,
    min_seq_id: float,
    coverage: float,
    cov_mode: int,
    cluster_mode: int,
    threads: int,
) -> Path:
    """Run MMseqs2 linclust and return the representative/member TSV."""
    mmseqs = shutil.which("mmseqs")
    if mmseqs is None:
        raise RuntimeError(
            "mmseqs binary not found on PATH; activate the pufprobe conda "
            "environment or install MMseqs2 (conda install -c conda-forge mmseqs2)."
        )
    work_dir.mkdir(parents=True, exist_ok=True)
    tmp = work_dir / "tmp"
    seq_db = work_dir / "seqdb"
    clu_db = work_dir / "clu"
    tsv = work_dir / "clusters.tsv"

    if tsv.exists() and tsv.stat().st_size > 0:
        return tsv

    run_mmseqs([mmseqs, "createdb", str(fasta), str(seq_db)])
    run_mmseqs(
        [
            mmseqs, "linclust", str(seq_db), str(clu_db), str(tmp),
            "--min-seq-id", str(min_seq_id),
            "-c", str(coverage),
            "--cov-mode", str(cov_mode),
            "--cluster-mode", str(cluster_mode),
            "--threads", str(threads),
        ]
    )
    run_mmseqs([mmseqs, "createtsv", str(seq_db), str(seq_db), str(clu_db), str(tsv)])
    return tsv


def parse_clusters(tsv: Path) -> Dict[str, List[str]]:
    """Return {representative: [members]} from createtsv output."""
    groups: Dict[str, List[str]] = {}
    with tsv.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            groups.setdefault(rep, []).append(member)
    return groups


def build_families(
    groups: Dict[str, List[str]], proteins: Dict[str, Protein]
) -> List[dict]:
    families: List[dict] = []
    ordered = sorted(
        groups.values(),
        key=lambda members: (-len(members), min(members)),
    )
    for order, members in enumerate(ordered, start=1):
        valid = [m for m in members if m in proteins]
        if not valid:
            continue
        rep_key = max(valid, key=lambda m: (proteins[m].length, m))
        phyla = sorted({proteins[m].phylum for m in valid if proteins[m].phylum})
        organisms = sorted({proteins[m].organism for m in valid if proteins[m].organism})
        proteomes = sorted({proteins[m].proteome for m in valid if proteins[m].proteome})
        lengths = [proteins[m].length for m in valid]
        families.append({
            "family_id": f"X-FAM{order:06d}",
            "representative": rep_key,
            "size": len(valid),
            "n_phyla": len(phyla),
            "n_organisms": len(organisms),
            "n_proteomes": len(proteomes),
            "mean_length": round(float(np.mean(lengths)), 1),
            "phyla": ";".join(phyla),
            "members": sorted(valid),
        })
    return families


def embed_representatives(
    reps: List[Protein],
    output_dir: Path,
    model_name: str,
    batch_size: int,
    threads_gpu: int = 1,
) -> Path:
    """Embed representatives in resumable shards, then merge into one npz."""
    shard_dir = output_dir / "embedding_shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    merged = output_dir / "representative_embeddings.npz"

    shard_paths = sorted(shard_dir.glob("shard_*.npz"))
    done = 0
    for shard_path in shard_paths:
        with np.load(shard_path) as data:
            done += len(data.files)
    if done < len(reps) or not merged.exists():
        from pufprobe.embed.esm import ModelUnavailable, compute_embeddings

        for start in range(0, len(reps), EMBED_SHARD):
            shard_index = start // EMBED_SHARD
            shard_path = shard_dir / f"shard_{shard_index:05d}.npz"
            if shard_path.exists():
                continue
            batch = reps[start:start + EMBED_SHARD]
            vectors = compute_embeddings(
                batch, model_name=model_name, cpu=False, batch_size=batch_size
            )
            if not vectors:
                raise ModelUnavailable("embedding shard returned no vectors")
            np.savez(shard_path, **{k: np.asarray(v, dtype=np.float32)
                                    for k, v in vectors.items()})
            print(f"    embedded {min(start + EMBED_SHARD, len(reps))}/{len(reps)} "
                  f"representatives")

    merged_vectors: Dict[str, np.ndarray] = {}
    for shard_path in sorted(shard_dir.glob("shard_*.npz")):
        with np.load(shard_path) as data:
            for key in data.files:
                merged_vectors[key] = np.asarray(data[key], dtype=np.float32)
    if merged_vectors:
        np.savez(merged, **merged_vectors)
    return merged


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="all_hypothetical_proteins.tsv from stage 07")
    parser.add_argument("--output", type=Path, required=True,
                        help="directory for families and representative embeddings")
    parser.add_argument("--min-seq-id", type=float, default=0.30)
    parser.add_argument("--coverage", type=float, default=0.80)
    parser.add_argument("--cov-mode", type=int, default=0)
    parser.add_argument("--cluster-mode", type=int, default=1)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--embed-model", default="esm2_t33_650M_UR50D")
    parser.add_argument("--embed-batch", type=int, default=16)
    parser.add_argument("--skip-embedding", action="store_true",
                        help="cluster only; do not compute ESM-2 embeddings")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    proteins = load_proteins(args.input)
    if not proteins:
        print("No proteins loaded from input TSV.", file=sys.stderr)
        return 1
    print(f"Loaded {len(proteins)} hypothetical proteins.")

    fasta = args.output / "all_hypothetical.faa"
    if not fasta.exists():
        write_fasta(fasta, [(p.key, p.sequence_aa) for p in
                            sorted(proteins.values(), key=lambda x: x.key)])

    print("Running MMseqs2 linclust ...")
    cluster_tsv = linclust(
        fasta, args.output / "linclust",
        args.min_seq_id, args.coverage, args.cov_mode, args.cluster_mode,
        args.threads,
    )
    groups = parse_clusters(cluster_tsv)
    families = build_families(groups, proteins)

    # Family table
    family_path = args.output / "families.tsv"
    with family_path.open("w", encoding="utf-8", newline="\n") as handle:
        header = ["family_id", "representative", "size", "n_phyla", "n_organisms",
                  "n_proteomes", "mean_length", "phyla"]
        handle.write("\t".join(header) + "\n")
        for fam in families:
            handle.write("\t".join(str(fam[k]) for k in header) + "\n")

    # Member-level mapping (accession -> family) with taxonomy
    member_path = args.output / "family_members.tsv"
    with member_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("family_id\taccession\tphylum\torganism\tproteome\tlength\n")
        for fam in families:
            for member in fam["members"]:
                p = proteins[member]
                handle.write(f"{fam['family_id']}\t{member}\t{p.phylum}\t"
                             f"{p.organism}\t{p.proteome}\t{p.length}\n")

    # Representative FASTA
    rep_fasta = args.output / "representatives.faa"
    reps = [proteins[fam["representative"]] for fam in families]
    write_fasta(rep_fasta, [(p.key, p.sequence_aa) for p in reps])

    summary = {
        "n_proteins": len(proteins),
        "n_families": len(families),
        "singletons": sum(1 for fam in families if fam["size"] == 1),
        "multi_phylum_families": sum(1 for fam in families if fam["n_phyla"] > 1),
        "n_representatives": len(reps),
        "linclust": {"min_seq_id": args.min_seq_id, "coverage": args.coverage,
                     "cov_mode": args.cov_mode, "cluster_mode": args.cluster_mode},
    }
    (args.output / "cluster_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print(f"Clustering done: {len(proteins)} proteins -> {len(families)} families "
          f"({summary['singletons']} singletons, "
          f"{summary['multi_phylum_families']} cross-phylum families).")

    if not args.skip_embedding:
        print(f"Computing ESM-2 embeddings for {len(reps)} representatives ...")
        try:
            embed_path = embed_representatives(
                reps, args.output, args.embed_model, args.embed_batch
            )
            manifest = {
                "model": args.embed_model,
                "n_vectors": len(reps),
                "path": str(embed_path),
                "pooling": "residue-mean",
            }
            (args.output / "embedding_manifest.json").write_text(
                json.dumps(manifest, indent=2), encoding="utf-8")
            print(f"Embeddings written to {embed_path}.")
        except Exception as exc:  # noqa: BLE001 - clustering is already complete
            print(f"Embedding step skipped: {exc}", file=sys.stderr)
            print("Clustering outputs are complete; re-run after fixing the GPU "
                  "environment to generate embeddings.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
