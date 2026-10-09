"""prepare stage: parse the demo genome and check record/status/QC bookkeeping."""

from pathlib import Path

from pufprobe.config import load_config
from pufprobe.io_utils import read_json
from pufprobe.prepare.pipeline import run as prepare_run
from pufprobe.workdir import WorkDir

EXAMPLES = Path(__file__).resolve().parent.parent / "examples" / "demo_genome"


def _inputs() -> list[Path]:
    return sorted(list(EXAMPLES.rglob("*.gff3")) + list(EXAMPLES.rglob("*.faa")))


def test_prepare_parses_both_samples(tmp_path: Path) -> None:
    wd = WorkDir(tmp_path / "run")
    config = load_config(None)
    prepare_run(wd, config, _inputs())

    records = read_json(wd.path("prepare", "cds_records.json"))  # type: ignore[arg-type]
    assert len(records) == 25  # sampleA 13 + sampleB 12

    # Every CDS must align to a .faa sequence; none silently dropped.
    assert all(r["sequence_aa"] for r in records)

    samples = read_json(wd.path("prepare", "samples.json"))  # type: ignore[arg-type]
    sample_ids = {s["sample_id"] for s in samples}
    assert sample_ids == {"sampleA", "sampleB"}

    by_status = {}
    for r in records:
        by_status[r["status"]] = by_status.get(r["status"], 0) + 1
    assert by_status["annotated"] == 18
    assert by_status["hypothetical"] == 5
    assert by_status["conserved_hypothetical"] == 2

    qc = read_json(wd.path("qc", "qc_summary.json"))  # type: ignore[arg-type]
    assert qc["total_cds"] == 25
    assert qc["by_verdict"]["likely_false_orf"] == 1
    assert qc["unaligned"] == []

    # Completion markers written for both stages.
    assert wd.marker("prepare").is_file()
    assert wd.marker("qc").is_file()


def test_false_orf_locus_is_flagged(tmp_path: Path) -> None:
    wd = WorkDir(tmp_path / "run")
    prepare_run(wd, load_config(None), _inputs())
    records = read_json(wd.path("prepare", "cds_records.json"))  # type: ignore[arg-type]
    flagged = [r for r in records if r["verdict"] == "likely_false_orf"]
    assert len(flagged) == 1
    assert flagged[0]["locus_tag"] == "SA_0012"
