#!/usr/bin/env python3
"""ProstT5 structure-aware embeddings for family representatives (optional).

ProstT5 is a sequence-to-structure language model that translates amino-acid
sequences into the Foldseek 3Di structural alphabet. Its encoder provides a
fold-aware embedding that captures relationships invisible to a sequence-only
model. The embeddings are produced in resumable shards and written in the same
``.npz`` layout as the ESM-2 embeddings, so the network stage can combine them.

This channel is optional and is skipped cleanly when transformers or the
ProstT5 checkpoint cannot be loaded (for example without GPU or offline). It
does not replace the AlphaFold/foldseek step; it adds a second, database-free
structural representation that also covers representatives with no foldseek
hit.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np

EMBED_SHARD = 500
MAX_LEN = 1000


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    key = None
    chunks: List[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if key is not None:
                records[key] = "".join(chunks)
            key = line[1:].strip().split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if key is not None:
        records[key] = "".join(chunks)
    return records


def load_model(model_name: str):
    try:
        import torch  # noqa: F401
        from transformers import AutoModelForSeq2SeqLM, T5Tokenizer
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            "ProstT5 needs transformers and torch; install with "
            "`pip install transformers sentencepiece torch`."
        ) from exc
    try:
        tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
        model = AutoModelForSeq2SeqLM.from_pretrained(model_name)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"Could not load ProstT5 checkpoint {model_name}.") from exc

    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda" and model_name.endswith("fp16"):
        model = model.half()
    model = model.to(device).eval()
    return model, tokenizer, device, torch


def encode_batch(model, tokenizer, torch_module, sequences: List[str], device: str
                 ) -> np.ndarray:
    prepared = ["<AA2fold> " + " ".join(seq[:MAX_LEN]) for seq in sequences]
    tokens = tokenizer.batch_encode_plus(
        prepared, add_special_tokens=True, padding="longest",
        truncation=True, max_length=MAX_LEN + 2, return_tensors="pt",
    ).to(device)
    with torch_module.no_grad():
        output = model.encoder(
            input_ids=tokens.input_ids,
            attention_mask=tokens.attention_mask,
            return_dict=True,
        )
    hidden = output.last_hidden_state  # (B, L, D)
    mask = tokens.attention_mask.unsqueeze(-1).float()
    mean = (hidden * mask).sum(1) / mask.sum(1).clamp(min=1)
    return mean.detach().cpu().float().numpy().astype(np.float32)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model", default="Rostlab/ProstT5_fp16")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--max-reps", type=int, default=0)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    shard_dir = args.output / "prostt5_shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    merged = args.output / "prostt5_embeddings.npz"

    records = read_fasta(args.reps)
    keys = sorted(records)
    if args.max_reps > 0:
        keys = keys[:args.max_reps]
    print(f"ProstT5 representatives: {len(keys)}")

    try:
        model, tokenizer, device, torch_module = load_model(args.model)
    except RuntimeError as exc:
        print(f"ProstT5 channel skipped: {exc}", file=sys.stderr)
        return 0

    for start in range(0, len(keys), EMBED_SHARD):
        shard_index = start // EMBED_SHARD
        shard_path = shard_dir / f"shard_{shard_index:05d}.npz"
        if shard_path.exists():
            continue
        shard_keys = keys[start:start + EMBED_SHARD]
        vectors: Dict[str, np.ndarray] = {}
        for inner in range(0, len(shard_keys), args.batch_size):
            batch_keys = shard_keys[inner:inner + args.batch_size]
            batch_seqs = [records[k] for k in batch_keys]
            emb = encode_batch(model, tokenizer, torch_module, batch_seqs, device)
            for key, vector in zip(batch_keys, emb):
                vectors[key] = vector
        np.savez(shard_path, **{k: np.asarray(v, dtype=np.float32)
                                for k, v in vectors.items()})
        print(f"    ProstT5: {min(start + EMBED_SHARD, len(keys))}/{len(keys)}")

    merged_vectors: Dict[str, np.ndarray] = {}
    for shard_path in sorted(shard_dir.glob("shard_*.npz")):
        with np.load(shard_path) as data:
            for key in data.files:
                merged_vectors[key] = np.asarray(data[key], dtype=np.float32)
    if merged_vectors:
        np.savez(merged, **merged_vectors)
        dim = next(iter(merged_vectors.values())).shape[0]
        (args.output / "prostt5_manifest.json").write_text(json.dumps({
            "model": args.model, "n_vectors": len(merged_vectors),
            "dim": dim, "pooling": "encoder-mean",
            "path": str(merged),
        }, indent=2), encoding="utf-8")
        print(f"ProstT5 embeddings written to {merged} ({len(merged_vectors)} x {dim}).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
