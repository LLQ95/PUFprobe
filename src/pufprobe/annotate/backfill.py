"""Backfill high-confidence predictions into new annotation files.

Original GFF3 / GenBank / FASTA files are only ever opened for reading.  For
every sample the stage writes parallel files under ``09_annotate``:

* ``<sample>.pufprobe.gff3`` — every line copied verbatim except matched CDS
  lines, whose ``product`` attribute is replaced, ``Dbxref`` entries appended
  and a ``PUFprobe=...`` evidence note added;
* ``<sample>.pufprobe.gbk``  — minimal line-based GenBank rewrite when a
  GenBank source exists (``/product`` replaced, ``/db_xref`` and ``/note``
  appended);
* ``<sample>.pufprobe.faa``  — protein records renamed via
  :func:`pufprobe.io_utils.write_fasta`.

A ``changes.tsv`` ledger records every edit.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from pufprobe.io_utils import read_fasta, write_fasta
from pufprobe.logging_utils import get_logger
from pufprobe.models import ProteinPrediction, Sample, TermPrediction

log = get_logger("annotate.backfill")


def select_terms(
    predictions: List[ProteinPrediction], min_confidence: float
) -> Dict[str, Dict[str, Any]]:
    """Map global key -> accepted edit info for terms >= min_confidence."""
    accepted: Dict[str, Dict[str, Any]] = {}
    for pred in predictions:
        terms = sorted(
            (t for t in pred.terms if t.confidence >= min_confidence),
            key=lambda t: (-t.confidence, t.term_type, t.term_id),
        )
        if not terms:
            continue
        groups = sorted({g for t in terms for g in t.evidence_groups})
        accepted[pred.key] = {
            "new_product": " / ".join(t.term_name for t in terms),
            "confidence": max(t.confidence for t in terms),
            "evidence_groups": groups,
            "term_ids": [f"{t.term_type}:{t.term_id}" for t in terms],
            "terms": terms,
        }
    return accepted


# ---------------------------------------------------------------------------
# GFF3
# ---------------------------------------------------------------------------
def parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in attr_string.split(";"):
        if not part:
            continue
        key, _, value = part.partition("=")
        out[key] = value
    return out


def rewrite_gff3(
    src: Path,
    out: Path,
    sample_id: str,
    accepted: Dict[str, Dict[str, Any]],
) -> List[Dict[str, str]]:
    """Copy ``src`` to ``out``, editing matched CDS lines. Returns change rows."""
    changes: List[Dict[str, str]] = []
    out.parent.mkdir(parents=True, exist_ok=True)
    with src.open("r", encoding="utf-8", errors="replace") as handle:
        lines = handle.read().splitlines()

    rewritten: List[str] = []
    for line in lines:
        if not line or line.startswith("#"):
            rewritten.append(line)
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            rewritten.append(line)
            continue
        attrs = parse_gff_attributes(fields[8])
        locus = attrs.get("locus_tag") or attrs.get("ID", "")
        key = f"{sample_id}@{locus}"
        if key not in accepted:
            rewritten.append(line)
            continue

        edit = accepted[key]
        old_product = attrs.get("product", "")
        attrs["product"] = edit["new_product"]

        dbxref = [x for x in attrs.get("Dbxref", "").split(",") if x]
        for term_id in edit["term_ids"]:
            if term_id not in dbxref:
                dbxref.append(term_id)
        if dbxref:
            attrs["Dbxref"] = ",".join(dbxref)

        note = f"PUFprobe:confidence={edit['confidence']:.3f};evidence_groups={'+'.join(edit['evidence_groups'])}"
        if attrs.get("note"):
            attrs["note"] = f"{attrs['note']},{note}"
        else:
            attrs["note"] = note

        fields[8] = ";".join(f"{k}={v}" for k, v in attrs.items())
        rewritten.append("\t".join(fields))
        changes.append(
            {
                "sample": sample_id,
                "locus": locus,
                "old_product": old_product,
                "new_product": edit["new_product"],
                "confidence": f"{edit['confidence']:.3f}",
                "supporting_terms": ";".join(edit["term_ids"]),
            }
        )

    with out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(rewritten) + "\n")
    return changes


# ---------------------------------------------------------------------------
# GenBank (minimal line-based rewrite)
# ---------------------------------------------------------------------------
_PRODUCT_RE = re.compile(r"^(\s+)/product=")
_LOCUS_RE = re.compile(r'^\s+/locus_tag="([^"]+)"')


def rewrite_gbk(
    src: Path,
    out: Path,
    sample_id: str,
    accepted: Dict[str, Dict[str, Any]],
) -> int:
    """Minimal GenBank rewrite; returns the number of edited feature blocks."""
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = src.read_text(encoding="utf-8", errors="replace").splitlines()
    rewritten: List[str] = []
    i = 0
    edited = 0
    while i < len(lines):
        line = lines[i]
        if not line.startswith(" ") or line.startswith("FEATURES"):
            rewritten.append(line)
            i += 1
            continue
        # Collect one feature block: a non-qualifier line starts a new feature.
        block = [line]
        j = i + 1
        while j < len(lines):
            nxt = lines[j]
            if nxt.startswith("ORIGIN") or nxt.startswith("//"):
                break
            if not nxt.startswith(" "):
                break
            stripped = nxt.lstrip()
            if not stripped.startswith("/"):
                break
            block.append(nxt)
            j += 1

        locus = None
        for bl in block:
            m = _LOCUS_RE.match(bl)
            if m:
                locus = m.group(1)
                break
        key = f"{sample_id}@{locus}" if locus else None
        if key is None or key not in accepted:
            rewritten.extend(block)
        else:
            edit = accepted[key]
            for bl in block:
                rewritten.append(bl)
                if _PRODUCT_RE.match(bl):
                    indent = _PRODUCT_RE.match(bl).group(1)
                    rewritten.append(f'{indent}/db_xref="{";".join(edit["term_ids"])}"')
                    rewritten.append(
                        f'{indent}/note="PUFprobe:confidence={edit["confidence"]:.3f};'
                        f'evidence_groups={"+".join(edit["evidence_groups"])}"'
                    )
            edited += 1
        i = j
    out.write_text("\n".join(rewritten) + "\n", encoding="utf-8", newline="\n")
    return edited


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------
def rewrite_faa(
    src: Path,
    out: Path,
    sample_id: str,
    accepted: Dict[str, Dict[str, Any]],
) -> int:
    """Rename matched protein records; returns the number of renamed records."""
    renamed = 0
    records = []
    for header, seq in read_fasta(src):
        first = header.split()[0] if header else ""
        new_header = header
        for key, edit in accepted.items():
            _, _, locus = key.partition("@")
            if locus and locus in header:
                new_header = f"{first} {edit['new_product']}"
                renamed += 1
                break
        records.append((new_header, seq))
    write_fasta(out, records)
    return renamed


def run_backfill(
    out_dir: Path,
    samples: List[Sample],
    accepted: Dict[str, Dict[str, Any]],
) -> List[Dict[str, str]]:
    """Backfill every sample; returns the combined change ledger."""
    all_changes: List[Dict[str, str]] = []
    for sample in samples:
        sid = sample.sample_id
        if sample.gff_path:
            changes = rewrite_gff3(
                Path(sample.gff_path), out_dir / f"{sid}.pufprobe.gff3", sid, accepted
            )
            all_changes.extend(changes)
            log.info("%s: %d CDS rows backfilled in GFF3.", sid, len(changes))
        if sample.gbk_path:
            n = rewrite_gbk(
                Path(sample.gbk_path), out_dir / f"{sid}.pufprobe.gbk", sid, accepted
            )
            log.info("%s: %d feature blocks edited in GenBank.", sid, n)
        if sample.faa_path:
            n = rewrite_faa(
                Path(sample.faa_path), out_dir / f"{sid}.pufprobe.faa", sid, accepted
            )
            log.info("%s: %d protein records renamed in FAA.", sid, n)
    return all_changes
