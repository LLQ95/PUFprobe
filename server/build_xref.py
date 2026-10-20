#!/usr/bin/env python
"""Build the Swiss-Prot cross-reference index used for term normalization.

Usage:
    python build_xref.py <uniprot_sprot.dat.gz> <swissprot_xref.json>

The flat text file is parsed into a compact accession -> {Pfam, EC, GO} JSON
index consumed by the search stage.
"""

from __future__ import annotations

import sys
from pathlib import Path

from pufprobe.search.idmapping import build_xref


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: build_xref.py <uniprot_sprot.dat.gz> <swissprot_xref.json>"
        )
    dat_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
    if not dat_path.exists():
        raise SystemExit(f"input not found: {dat_path}")
    build_xref(dat_path, out_path)


if __name__ == "__main__":
    main()
