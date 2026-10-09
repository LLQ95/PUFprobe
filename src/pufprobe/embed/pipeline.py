"""embed stage pipeline.

Contract implemented here:

1. Load CDSRecord list from 00_prepare/cds_records.json and keep records that
   are candidates for inference (verdict != likely_false_orf; status
   hypothetical or conserved_hypothetical).
2. Produce one vector per candidate, in this order of preference:
   - --precomputed: a .npz mapping global key to vector, or a directory of
     <key>.npy files;
   - ESM-2 through pufprobe.embed.esm when torch and weights are available
     (use the config cpu_model when workdir.cpu);
   - fallback features from pufprobe.features when torch/weights are missing,
     with an explicit log message.
3. Cluster candidates into families:
   - MMseqs2 through pufprobe.clustering.mmseqs when the binary is present;
   - pure-Python graph clustering through pufprobe.clustering.pycluster using
     cosine similarity of the vectors.
4. Select one representative per family (longest sequence, ties by key).

Outputs:
- 02_embed/embeddings.npz    key -> vector
- 02_embed/manifest.json     vector dimension and embedding source
- 03_cluster/families.json   list of Family dicts
- 03_cluster/clusters.tsv    family_id, member global keys
- 03_cluster/representatives.faa
The family_id and is_representative fields are added back to the records and
rewritten to 00_prepare/cds_records.json (additive update only).
Completion markers are written for embed and cluster.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from pufprobe.clustering import mmseqs as mmseqs_clustering
from pufprobe.clustering import pycluster
from pufprobe.embed.esm import ModelUnavailable, compute_embeddings
from pufprobe.embed.fallback import compute_fallback_embeddings
from pufprobe.io_utils import read_json, write_fasta, write_json
from pufprobe.logging_utils import get_logger
from pufprobe.models import (
    CDSRecord,
    Family,
    GeneStatus,
    GeneVerdict,
    to_dict,
)
from pufprobe.workdir import WorkDir

logger = get_logger(__name__)

_CANDIDATE_STATUS = {GeneStatus.HYPOTHETICAL.value, GeneStatus.CONSERVED_HYPOTHETICAL.value}


def run(
    workdir: WorkDir,
    config: Dict[str, Any],
    model: Optional[str] = None,
    precomputed: Optional[Path] = None,
    fallback_only: bool = False,
) -> None:
    records_path = workdir.path("prepare", "cds_records.json")
    if not records_path.exists():
        raise FileNotFoundError(f"Missing {records_path}; run `pufprobe prepare` first.")

    raw_rows: List[dict] = read_json(records_path)  # type: ignore[assignment]
    all_records = [_row_to_record(row) for row in raw_rows]
    candidates = [rec for rec in all_records if _is_candidate(rec)]
    logger.info(
        "Loaded %d CDS records; %d candidates for embedding/clustering.",
        len(all_records), len(candidates),
    )
    if not candidates:
        logger.warning("No candidate records; writing empty embedding/cluster outputs.")

    # --- 1. compute embeddings ---------------------------------------------
    embed_cfg = config["embed"]
    fb_dim = int(embed_cfg["fallback_feature_dim"])

    vectors: Dict[str, np.ndarray] = {}
    source: str
    manifest: Dict[str, Any] = {"n_records": len(candidates)}

    if precomputed is not None:
        vectors, source = _load_precomputed(precomputed, candidates)
        manifest["precomputed_path"] = str(precomputed)
    elif fallback_only or workdir.cpu:
        vectors, filled = compute_fallback_embeddings(candidates, fb_dim)
        source = "fallback"
        manifest["filled_dim"] = filled
        logger.warning(
            "Using sequence-feature fallback embeddings (fallback_only=%s, cpu=%s).",
            fallback_only, workdir.cpu,
        )
    else:
        model_name = model or (
            embed_cfg["cpu_model"] if workdir.cpu else embed_cfg["default_model"]
        )
        try:
            vectors = compute_embeddings(
                candidates,
                model_name=model_name,
                cpu=workdir.cpu,
                batch_size=int(embed_cfg["batch_size"]),
            )
            source = "esm"
            manifest["model"] = model_name
        except ModelUnavailable as exc:
            logger.warning("ESM unavailable (%s); falling back to sequence features.", exc)
            vectors, filled = compute_fallback_embeddings(candidates, fb_dim)
            source = "fallback"
            manifest["filled_dim"] = filled
            manifest["esm_error"] = str(exc)

    # keep only candidates that actually produced a vector
    candidates = [rec for rec in candidates if rec.key in vectors]
    dims = {int(v.shape[0]) for v in vectors.values()}
    if not dims:
        dim = fb_dim
    else:
        dim = max(dims)
        if len(dims) > 1:
            logger.warning("Inconsistent embedding dimensions %s; padding all to %d.", dims, dim)
            padded: Dict[str, np.ndarray] = {}
            for key, vec in vectors.items():
                if vec.shape[0] < dim:
                    fixed = np.zeros(dim, dtype=np.float32)
                    fixed[: vec.shape[0]] = vec
                    padded[key] = fixed
                else:
                    padded[key] = vec
            vectors = padded

    # --- 2. write embeddings.npz + manifest --------------------------------
    embed_dir = workdir.stage_dir("embed")
    np.savez(embed_dir / "embeddings.npz", **vectors)
    manifest.update({"source": source, "dim": dim, "n_records": len(candidates)})
    if source == "fallback":
        manifest["padded"] = dim > manifest.get("filled_dim", dim)
    write_json(embed_dir / "manifest.json", manifest)

    # --- 3. cluster --------------------------------------------------------
    families = _cluster(workdir, config, candidates, vectors, source)

    # --- 4. write cluster outputs -----------------------------------------
    cluster_dir = workdir.stage_dir("cluster")
    write_json(cluster_dir / "families.json", [to_dict(f) for f in families])

    fam_by_id = {f.family_id: f for f in families}
    member_to_family: Dict[str, str] = {}
    with (cluster_dir / "clusters.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("family_id\tmember_key\n")
        for fam in families:
            for member in fam.members:
                member_to_family[member] = fam.family_id
                handle.write(f"{fam.family_id}\t{member}\n")

    write_fasta(
        cluster_dir / "representatives.faa",
        [(fam.representative_key, fam.representative_seq) for fam in families],
    )

    # --- 5. back-write family_id / is_representative into cds_records.json --
    _backwrite_records(records_path, raw_rows, member_to_family, fam_by_id)

    workdir.mark_complete("embed", payload={"source": source, "dim": dim,
                                            "n_records": len(candidates)})
    workdir.mark_complete("cluster", payload={"n_families": len(families)})
    logger.info("embed+cluster done: %d candidates -> %d families (source=%s, dim=%d).",
                len(candidates), len(families), source, dim)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _row_to_record(row: dict) -> CDSRecord:
    """Rebuild a CDSRecord from its to_dict row, resolving enum fields.

    Local copy because the shared ``models.from_dict`` relies on
    ``annotation.startswith`` which fails on PEP 585 / runtime typing generics.
    """
    row = dict(row)
    row["status"] = GeneStatus(row["status"])
    verdict = row.get("verdict")
    row["verdict"] = GeneVerdict(verdict) if verdict else None
    return CDSRecord(**row)


def _is_candidate(rec: CDSRecord) -> bool:
    if rec.status.value not in _CANDIDATE_STATUS:
        return False
    if rec.verdict is not None and rec.verdict.value == GeneVerdict.LIKELY_FALSE_ORF.value:
        return False
    return True


def _load_precomputed(
    precomputed: Path, candidates: List[CDSRecord]
) -> tuple:
    vectors: Dict[str, np.ndarray] = {}
    cand_keys = {rec.key for rec in candidates}
    if precomputed.is_file() and precomputed.suffix.lower() == ".npz":
        with np.load(precomputed, allow_pickle=True) as data:
            for key in data.files:
                if key in cand_keys:
                    vectors[key] = np.asarray(data[key], dtype=np.float32)
    elif precomputed.is_dir():
        for rec in candidates:
            npy = precomputed / f"{rec.key}.npy"
            if npy.exists():
                vectors[rec.key] = np.load(npy).astype(np.float32)
    else:
        raise FileNotFoundError(
            f"--precomputed must be a .npz file or a directory of <key>.npy: {precomputed}"
        )
    missing = cand_keys - set(vectors)
    if missing:
        logger.warning("%d candidate(s) lacked precomputed vectors; excluded: %s",
                       len(missing), sorted(missing)[:5])
    return vectors, "precomputed"


def _cluster(
    workdir: WorkDir,
    config: Dict[str, Any],
    candidates: List[CDSRecord],
    vectors: Dict[str, np.ndarray],
    source: str,
) -> List[Family]:
    seq_by_key = {rec.key: rec.sequence_aa for rec in candidates}
    sample_by_key = {rec.key: rec.sample_id for rec in candidates}

    cluster_dir = workdir.stage_dir("cluster")
    cand_fasta = cluster_dir / "_candidates.faa"
    write_fasta(cand_fasta, [(k, seq_by_key[k]) for k in sorted(seq_by_key)])

    cl_cfg = config["cluster"]
    mapping = mmseqs_clustering.run_mmseqs_cluster(
        cand_fasta,
        cluster_dir / "_mmseqs",
        identity=float(cl_cfg["identity"]),
        coverage=float(cl_cfg["coverage"]),
        threads=int(workdir.threads),
    )
    if mapping is None:
        logger.info("Clustering candidates with pure-numpy cosine/pycluster.")
        threshold = float(cl_cfg["cosine_threshold"])
        kmer_gate = 0.0
        if source != "esm":
            # Coarse composition descriptors need a stricter merge threshold
            # than protein language model embeddings, plus a k-mer order gate.
            threshold = max(threshold, 0.95)
            kmer_gate = 0.5
        mapping = pycluster.cluster_vectors(
            vectors,
            cosine_threshold=threshold,
            sequences=seq_by_key if kmer_gate else None,
            kmer_gate=kmer_gate,
        )

    # group members by cluster index
    groups: Dict[int, List[str]] = {}
    for key, cidx in mapping.items():
        groups.setdefault(cidx, []).append(key)

    families: List[Family] = []
    ordered = sorted(
        groups.values(),
        key=lambda members: (-len(members), min(members)),
    )
    for order, members in enumerate(ordered, start=1):
        rep_key = _pick_representative(members, seq_by_key)
        samples = [sample_by_key[m] for m in members]
        sample_counts: Dict[str, int] = {}
        for s in samples:
            sample_counts[s] = sample_counts.get(s, 0) + 1
        families.append(
            Family(
                family_id=f"FAM{order:05d}",
                representative_key=rep_key,
                members=sorted(members),
                representative_seq=seq_by_key[rep_key],
                size=len(members),
                lineage_width=len(set(samples)),
                metrics={"sample_distribution": sample_counts},
            )
        )
    return families


def _pick_representative(members: List[str], seq_by_key: Dict[str, str]) -> str:
    # longest sequence; ties broken by lexicographically smallest key
    return sorted(members, key=lambda k: (-len(seq_by_key[k]), k))[0]


def _backwrite_records(
    records_path: Path,
    raw_rows: List[dict],
    member_to_family: Dict[str, str],
    fam_by_id: Dict[str, Family],
) -> None:
    for row in raw_rows:
        key = f"{row.get('sample_id')}@{row.get('locus_tag')}"
        fam_id = member_to_family.get(key)
        row["family_id"] = fam_id
        row["is_representative"] = (
            fam_id is not None and key == fam_by_id[fam_id].representative_key
        )
    write_json(records_path, raw_rows)
