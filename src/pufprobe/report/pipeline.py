"""report stage pipeline.

Contract implemented here:

1. Load predictions, CDS records and families.
2. Always write machine-readable 08_report/evidence.json containing every
   prediction with its full evidence list, raw and calibrated scores.
3. Write report.md and, when requested, a self-contained report.html with
   inline CSS. Each protein section lists the proposed terms and an auditable
   per-channel evidence table (target id, label, score, e-value, coverage,
   group). Unanchored families and their priority scores are listed separately.

Completion marker written for report.
"""

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple

from pufprobe.io_utils import read_json, write_json
from pufprobe.logging_utils import get_logger
from pufprobe.models import (
    CDSRecord,
    families_from_dicts,
    to_dict,
)
from pufprobe.rank.fuse import cds_record_from_row, load_predictions
from pufprobe.report.html_report import compute_overview, render_html, render_markdown
from pufprobe.workdir import WorkDir

log = get_logger("report.pipeline")


def _load_priorities(path: Path) -> List[Tuple[str, float, List[str]]]:
    rows: List[Tuple[str, float, List[str]]] = []
    if not path.exists():
        return rows
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("family_id"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        fid, score, reasons = parts[0], float(parts[1]), parts[2]
        rows.append((fid, score, reasons.split("; ") if reasons else []))
    return rows


def run(workdir: WorkDir, config: Dict[str, Any], formats: List[str]) -> None:
    predictions = load_predictions(workdir.path("rank", "predictions.json"))
    records: List[CDSRecord] = [
        cds_record_from_row(row) for row in read_json(workdir.path("prepare", "cds_records.json"))  # type: ignore[arg-type]
    ]
    families_path = workdir.path("cluster", "families.json")
    families = (
        families_from_dicts(read_json(families_path))  # type: ignore[arg-type]
        if families_path.exists()
        else []
    )
    priority_rows = _load_priorities(workdir.path("rank", "priorities.tsv"))
    high_conf = float(config["rank"]["high_confidence"])
    generated_at = datetime.now().isoformat(timespec="seconds")

    overview = compute_overview(records, predictions, high_conf)

    # Record index: one row per CDS, used by downstream consumers.
    family_of = {}
    for fam in families:
        for member in fam.members:
            family_of[member] = fam.family_id
    index = []
    for rec in records:
        pred = next((p for p in predictions if p.key == rec.key), None)
        accepted = [t for t in (pred.terms if pred else [])]
        index.append(
            {
                "key": rec.key,
                "sample_id": rec.sample_id,
                "locus_tag": rec.locus_tag,
                "product": rec.product,
                "status": rec.status.value,
                "verdict": rec.verdict.value if rec.verdict else None,
                "family_id": family_of.get(rec.key),
                "n_terms": len(accepted),
                "max_confidence": max((t.confidence for t in accepted), default=None),
            }
        )

    evidence_payload = {
        "generated": generated_at,
        "overview": overview,
        "predictions": [to_dict(p) for p in predictions],
        "index": index,
    }
    write_json(workdir.path("report", "evidence.json"), evidence_payload)

    if "md" in formats:
        md = render_markdown(overview, predictions, records, priority_rows, high_conf, generated_at)
        workdir.path("report", "report.md").write_text(md, encoding="utf-8", newline="\n")
    if "html" in formats:
        html_doc = render_html(overview, predictions, records, priority_rows, high_conf, generated_at)
        workdir.path("report", "report.html").write_text(html_doc, encoding="utf-8", newline="\n")

    workdir.mark_complete("report")
    log.info(
        "report done: %d predictions, %d CDS indexed; formats=%s.",
        len(predictions),
        len(index),
        ",".join(formats),
    )
