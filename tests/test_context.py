"""context stage: operon neighbours, phylogenetic profiles and core/accessory."""

from pathlib import Path

from pufprobe.config import load_config
from pufprobe.context.pipeline import run as context_run
from pufprobe.io_utils import read_json
from pufprobe.prepare.pipeline import run as prepare_run
from pufprobe.workdir import WorkDir

EXAMPLES = Path(__file__).resolve().parent.parent / "examples" / "demo_genome"


def _inputs() -> list[Path]:
    return sorted(list(EXAMPLES.rglob("*.gff3")) + list(EXAMPLES.rglob("*.faa")))


def _run(tmp_path: Path) -> WorkDir:
    wd = WorkDir(tmp_path / "run")
    config = load_config(None)
    config["context"]["min_genomes_for_profile"] = 2
    prepare_run(wd, config, _inputs())
    context_run(wd, config)
    return wd


def test_context_outputs_exist(tmp_path: Path) -> None:
    wd = _run(tmp_path)
    for name in ("neighbours.tsv", "phyloprofiles.tsv", "microsynteny.tsv"):
        assert wd.path("context", name).is_file(), name
    assert wd.path("context", "context_evidence.json").is_file()
    assert wd.marker("context").is_file()


def test_operon_neighbours_extracted(tmp_path: Path) -> None:
    wd = _run(tmp_path)
    rows = wd.path("context", "neighbours.tsv").read_text(encoding="utf-8").splitlines()
    # header + at least one neighbour row
    assert len(rows) >= 2
    coupled = [line for line in rows[1:] if line.endswith("True")]
    assert coupled, "expected at least one operon-coupled neighbour"


def test_shared_family_is_core_orphan_is_clonal(tmp_path: Path) -> None:
    wd = _run(tmp_path)
    rows = wd.path("context", "microsynteny.tsv").read_text(encoding="utf-8").splitlines()
    header = rows[0].split("\t")
    fid_i = header.index("family_id")
    cls_i = header.index("classification")

    classifications = {}
    for line in rows[1:]:
        cols = line.split("\t")
        classifications[cols[fid_i]] = cols[cls_i]

    # The shared conserved-hypothetical family appears in both samples -> core.
    core = [c for c in classifications.values() if c == "core"]
    assert core, "expected at least one core family"
    # Sample-restricted families are flagged clonal_specific.
    assert "clonal_specific" in classifications.values()

    evidence = read_json(wd.path("context", "context_evidence.json"))  # type: ignore[arg-type]
    assert evidence, "expected context Evidence entries"
