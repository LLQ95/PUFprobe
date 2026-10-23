#!/usr/bin/env python
"""Extract a per-accession annotation date from a Swiss-Prot flat file.

The CAFA-style temporal split needs a time anchor for every curated protein.
For each entry the "sequence version" date from the DT lines is used (the date
at which the sequence was finalized), falling back to the date the entry was
integrated into UniProtKB when no sequence-version line is present. The output
is a JSON object ``{accession: "YYYY-MM-DD"}`` that can be passed to
``build_ltr_trainset.py --dates`` and then to ``train_ltr.py``.

The reader is streaming and accepts a plain or gzip-compressed flat file.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path
from typing import Dict, Iterator, Optional

MONTHS = {
    "JAN": "01", "FEB": "02", "MAR": "03", "APR": "04",
    "MAY": "05", "JUN": "06", "JUL": "07", "AUG": "08",
    "SEP": "09", "OCT": "10", "NOV": "11", "DEC": "12",
}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_dt_date(token: str) -> Optional[str]:
    """Convert a UniProt date token '01-JUN-1999' to '1999-06-01'."""
    parts = token.strip().split("-")
    if len(parts) != 3:
        return None
    day, month, year = parts
    if month not in MONTHS:
        return None
    try:
        return f"{int(year):04d}-{MONTHS[month]}-{int(day):02d}"
    except ValueError:
        return None


def iter_dates(path: Path) -> Iterator[tuple]:
    """Yield (primary_accession, ISO date) for every flat-file entry."""
    primary: Optional[str] = None
    integrated: Optional[str] = None
    sequence_date: Optional[str] = None

    with open_text(path) as handle:
        for line in handle:
            code = line[:2]
            if code == "DT":
                body = line[5:].strip()
                date_token = body.split(",", 1)[0]
                iso = parse_dt_date(date_token)
                if iso is None:
                    continue
                if "integrated into UniProtKB" in body and integrated is None:
                    integrated = iso
                elif "sequence version" in body:
                    sequence_date = iso
            elif code == "AC":
                if primary is None:
                    first = line[5:].split(";", 1)[0].strip()
                    if first:
                        primary = first
            elif line.startswith("//"):
                if primary is not None:
                    date = sequence_date or integrated
                    if date:
                        yield primary, date
                primary = None
                integrated = None
                sequence_date = None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sprot", type=Path, required=True,
                        help="uniprot_sprot.dat or .dat.gz")
    parser.add_argument("--out", type=Path, required=True,
                        help="Output JSON {accession: YYYY-MM-DD}")
    args = parser.parse_args()

    dates: Dict[str, str] = {}
    for accession, iso in iter_dates(args.sprot):
        dates[accession] = iso

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(dates, handle)
    print(json.dumps({"out": str(args.out), "n_dates": len(dates)}, indent=2))
    return 0 if dates else 1


if __name__ == "__main__":
    raise SystemExit(main())
