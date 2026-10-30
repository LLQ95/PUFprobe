#!/usr/bin/env python
"""Whiten an already-built pLM reference library without re-embedding.

Building the full reference library takes many hours of GPU inference. This
script reuses the stored residue-mean embeddings, estimates the whitening
transform on them, and writes a new library whose vectors are whitened and
L2-normalized together with the transform parameters. The identifier and
term files are copied unchanged, so the expensive embedding step is not
repeated.

The output directory is created or overwritten independently; the source
library is left untouched. Re-running the script is safe and simply
regenerates the whitened output.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pufprobe.embed.whitening import (  # noqa: E402
    apply_whitening,
    whitening_matrix,
)
from pufprobe.logging_utils import setup_logging, get_logger  # noqa: E402
from pufprobe.search.plm_neighbor import (  # noqa: E402
    REF_IDS,
    REF_MANIFEST,
    REF_TERMS,
    REF_VECTORS,
    REF_WHITEN_MEAN,
    REF_WHITEN_W,
)

logger = get_logger("whiten_existing_library")


def main() -> int:
    setup_logging()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in", dest="src", type=Path, required=True,
                        help="Existing (unwhitened) reference library directory.")
    parser.add_argument("--out", dest="out", type=Path, required=True,
                        help="Output directory for the whitened library.")
    args = parser.parse_args()

    src = args.src
    out = args.out
    for name in (REF_VECTORS, REF_IDS, REF_TERMS, REF_MANIFEST):
        if not (src / name).exists():
            raise SystemExit(f"source library missing {name} in {src}")

    with (src / REF_MANIFEST).open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if manifest.get("whitened", False):
        logger.warning("Source library is already whitened; recomputing anyway.")

    out.mkdir(parents=True, exist_ok=True)
    raw = np.load(src / REF_VECTORS, mmap_mode="r").astype(np.float32)
    logger.info("Estimating whitening transform on %d x %d embeddings ...", *raw.shape)
    mean, whiten = whitening_matrix(raw)
    whitened = apply_whitening(raw, mean, whiten).astype(np.float16)

    np.save(out / REF_VECTORS, whitened)
    np.save(out / REF_WHITEN_MEAN, mean)
    np.save(out / REF_WHITEN_W, whiten)
    shutil.copyfile(src / REF_IDS, out / REF_IDS)
    shutil.copyfile(src / REF_TERMS, out / REF_TERMS)

    manifest["whitened"] = True
    manifest["normalized"] = True
    with (out / REF_MANIFEST).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    logger.info("Whitened library written to %s: %d proteins.", out, whitened.shape[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
