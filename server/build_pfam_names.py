#!/usr/bin/env python
"""Build the Pfam family name index used for term display names.

Usage:
    python build_pfam_names.py <Pfam-A.hmm.dat.gz> <pfam_names.json>

The Pfam metadata file is parsed into a compact accession -> official
description JSON index consumed by the rank stage so that Pfam terms show
the family's own name rather than a PDB id or a distant homolog's name.
"""

from __future__ import annotations

import sys
from pathlib import Path

from pufprobe.ontology_names import build_pfam_names


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: build_pfam_names.py <Pfam-A.hmm.dat.gz> <pfam_names.json>"
        )
    dat_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
    if not dat_path.exists():
        raise SystemExit(f"input not found: {dat_path}")
    build_pfam_names(dat_path, out_path)


if __name__ == "__main__":
    main()
