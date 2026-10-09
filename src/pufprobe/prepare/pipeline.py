"""prepare stage pipeline.

Grouping, source detection, CDS parsing, FASTA alignment, GeneStatus
assignment and true-gene quality control.  See the module docstring of the
frozen contract for the exact input/output contract.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from pufprobe.geneqc.classifier import classify_all
from pufprobe.io_utils import read_fasta, write_json
from pufprobe.logging_utils import get_logger
from pufprobe.models import CDSRecord, GeneStatus, Sample, to_dict
from pufprobe.prepare import gbk, gff
from pufprobe.workdir import WorkDir

log = get_logger("prepare")


def _sample_id_for(path: Path) -> str:
    """Default sample id: the containing directory name, else the file stem."""
    parent = path.parent.name
    return parent if parent not in ("", ".") else path.stem


def _kind_of(path: Path) -> Optional[str]:
    suffix = path.suffix.lower()
    if suffix in (".gff", ".gff3"):
        return "gff"
    if suffix in (".gbk", ".gbff", ".genbank"):
        return "gbk"
    if suffix in (".faa", ".fasta", ".fa"):
        return "faa"
    return None


def _assign_status(product: str, patterns: List[str]) -> GeneStatus:
    lowered = product.lower()
    if "duf" in lowered or "domain of unknown function" in lowered:
        return GeneStatus.CONSERVED_HYPOTHETICAL
    if any(pattern.lower() in lowered for pattern in patterns):
        return GeneStatus.HYPOTHETICAL
    return GeneStatus.ANNOTATED


def _group_samples(paths: List[Path]) -> Dict[str, Sample]:
    """Group annotation/protein files into samples by containing directory."""
    groups: Dict[str, Sample] = {}
    for path in paths:
        kind = _kind_of(path)
        if kind is None:
            log.warning("Skipping unrecognised input file: %s", path)
            continue
        sid = _sample_id_for(path)
        key = str(path.parent.resolve())
        sample = groups.get(key)
        if sample is None:
            sample = Sample(sample_id=sid)
            groups[key] = sample
        if kind == "gff":
            sample.gff_path = str(path)
        elif kind == "gbk":
            sample.gbk_path = str(path)
        elif kind == "faa":
            sample.faa_path = str(path)
    return groups


def _load_faa_index(faa_path: Path) -> Dict[str, str]:
    """Map the first whitespace-delimited FASTA header token to the sequence."""
    index: Dict[str, str] = {}
    for header, sequence in read_fasta(faa_path):
        token = header.split()[0] if header else ""
        if token:
            index[token] = sequence.upper().replace("*", "")
    return index


def _parse_sample_cds(sample: Sample) -> List[Dict[str, Any]]:
    if sample.gff_path is not None:
        return gff.parse_cds(Path(sample.gff_path))
    if sample.gbk_path is not None:
        return gbk.parse_gbk_cds(Path(sample.gbk_path))
    return []


def run(
    workdir: WorkDir,
    config: Dict[str, Any],
    paths: List[Path],
    include_smorf: bool = False,
) -> None:
    patterns = config["hypothetical_patterns"]
    groups = _group_samples([Path(p) for p in paths])
    if not groups:
        raise RuntimeError("No usable GFF3/GBK/FAA inputs were found.")

    records: List[CDSRecord] = []
    unaligned: List[Dict[str, str]] = []

    for key, sample in sorted(groups.items(), key=lambda kv: kv[1].sample_id):
        annotation_path = sample.gff_path or sample.gbk_path
        if annotation_path is None:
            log.warning("Sample %s has no annotation file; skipped.", sample.sample_id)
            continue

        if sample.gff_path:
            sample.source = gff.detect_source(Path(sample.gff_path))
        log.info(
            "Parsing sample '%s' (%s) from %s",
            sample.sample_id,
            sample.source,
            annotation_path,
        )

        cds_rows = _parse_sample_cds(sample)
        faa_index: Dict[str, str] = {}
        if sample.faa_path:
            faa_index = _load_faa_index(Path(sample.faa_path))

        for row in cds_rows:
            locus = row.get("locus_tag") or row.get("protein_id") or ""
            product = row.get("product", "")
            sequence_aa = ""
            # Align by locus tag first, then by protein id / FASTA header token.
            if locus in faa_index:
                sequence_aa = faa_index[locus]
            elif row.get("protein_id") and row["protein_id"] in faa_index:
                sequence_aa = faa_index[row["protein_id"]]
            elif faa_index:
                # No direct match: record it explicitly rather than dropping.
                unaligned.append({"sample": sample.sample_id, "locus_tag": locus, "product": product})
                log.warning(
                    "CDS %s in %s could not be aligned to any .faa sequence",
                    locus,
                    sample.sample_id,
                )
            # GenBank /translation is a fallback protein sequence with no .faa.
            if not sequence_aa and row.get("translation"):
                sequence_aa = row["translation"]

            record = CDSRecord(
                sample_id=sample.sample_id,
                locus_tag=locus,
                protein_id=row.get("protein_id"),
                sequence_aa=sequence_aa,
                start=int(row.get("start", 0)),
                end=int(row.get("end", 0)),
                strand=row.get("strand", ""),
                product=product,
                status=_assign_status(product, patterns),
                contig=row.get("contig", ""),
                dbxrefs=list(row.get("dbxrefs", [])),
                attributes=dict(row.get("attributes", {})),
                sequence_nt=None,
            )
            records.append(record)

    # True-gene quality control across the whole batch.
    summary = classify_all(records, config, include_smorf=include_smorf)

    # ---- write prepare outputs ------------------------------------------
    write_json(workdir.path("prepare", "samples.json"), [to_dict(s) for s in groups.values()])
    write_json(workdir.path("prepare", "cds_records.json"), [to_dict(r) for r in records])

    by_status: Dict[str, int] = {}
    for record in records:
        by_status[record.status.value] = by_status.get(record.status.value, 0) + 1

    qc_summary = {
        "total_cds": len(records),
        "include_smorf": include_smorf,
        "by_status": by_status,
        "by_verdict": summary["by_verdict"],
        "unaligned": unaligned,
    }
    write_json(workdir.path("qc", "qc_summary.json"), qc_summary)

    workdir.mark_complete("prepare", {"records": len(records), "samples": len(groups)})
    workdir.mark_complete("qc", qc_summary)
    log.info(
        "prepare/qc done: %d CDS across %d samples, %d unaligned",
        len(records),
        len(groups),
        len(unaligned),
    )
