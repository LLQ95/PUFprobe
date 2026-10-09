"""annotate stage pipeline.

Contract implemented here:

1. Load predictions and select terms whose confidence is at least
   --min-confidence (default config rank.min_confidence_annotate).
2. For each sample, produce new annotation files under 09_annotate/:
   - <sample>.pufprobe.gff3: original GFF rows with the product and Dbxref
     fields updated for matched CDS rows, plus a PUFprobe evidence note;
   - <sample>.pufprobe.gbk when a GBK source exists;
   - <sample>.pufprobe.faa: protein records renamed with the new product.
3. Original files are never modified or overwritten; a mapping table
   09_annotate/changes.tsv records every edit with the supporting term and
   confidence.

Completion marker written for annotate.
"""

from pathlib import Path
from typing import Any, Dict, Optional

from pufprobe.io_utils import read_json
from pufprobe.logging_utils import get_logger
from pufprobe.models import Sample
from pufprobe.annotate.backfill import run_backfill, select_terms
from pufprobe.rank.fuse import load_predictions
from pufprobe.workdir import WorkDir

log = get_logger("annotate.pipeline")


def run(
    workdir: WorkDir,
    config: Dict[str, Any],
    min_confidence: Optional[float] = None,
) -> None:
    threshold = (
        float(min_confidence)
        if min_confidence is not None
        else float(config["rank"]["min_confidence_annotate"])
    )

    predictions = load_predictions(workdir.path("rank", "predictions.json"))
    samples_path = workdir.path("prepare", "samples.json")
    if not samples_path.exists():
        raise FileNotFoundError(f"missing {samples_path}; run prepare first")
    samples = [Sample(**row) for row in read_json(samples_path)]  # type: ignore[arg-type]

    accepted = select_terms(predictions, threshold)
    log.info(
        "annotate: %d terms pass confidence threshold %.2f across %d proteins.",
        sum(len(edit["terms"]) for edit in accepted.values()),
        threshold,
        len(accepted),
    )

    out_dir = workdir.stage_dir("annotate")
    changes = run_backfill(out_dir, samples, accepted)

    changes_path = out_dir / "changes.tsv"
    with changes_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "sample\tlocus\told_product\tnew_product\tconfidence\tsupporting terms\n"
        )
        for row in changes:
            handle.write(
                f"{row['sample']}\t{row['locus']}\t{row['old_product']}\t"
                f"{row['new_product']}\t{row['confidence']}\t{row['supporting_terms']}\n"
            )

    workdir.mark_complete("annotate")
    log.info("annotate done: %d edits recorded in changes.tsv.", len(changes))
