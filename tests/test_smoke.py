"""End-to-end smoke test over the frozen examples/demo_genome layout.

Runs the full stage chain with CliRunner:

    prepare -> embed (fallback-only) -> fold -> search -> context -> rank
    -> report -> annotate

Each stage must exit with code 0 and leave its marker + key artefacts behind.
If the example genome has not been generated yet, the whole test is skipped;
the integration stage runs it once prepare/embed/fold/search/context are
available.
"""

from pathlib import Path

import pytest
from click.testing import CliRunner

from pufprobe.cli import cli

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = REPO_ROOT / "examples" / "demo_genome"

pytestmark = pytest.mark.skipif(
    not (EXAMPLE / "sampleA" / "sampleA.gff3").exists(),
    reason="examples/demo_genome not generated yet; run during integration",
)


def _run(runner: CliRunner, workdir: Path, *args: str):
    result = runner.invoke(cli, ["--workdir", str(workdir), "--threads", "2", *args])
    assert result.exit_code == 0, f"{args[0]} failed: {result.output}\n{result.exception}"
    return result


def test_full_chain_smoke(tmp_path: Path):
    workdir = tmp_path / "smoke_wd"
    runner = CliRunner()

    _run(runner, workdir, "prepare", "--input", str(EXAMPLE))
    _run(runner, workdir, "embed", "--fallback-only")
    _run(runner, workdir, "fold")
    _run(runner, workdir, "search")
    _run(runner, workdir, "context")
    _run(runner, workdir, "rank")
    _run(runner, workdir, "report", "--formats", "json,html,md")
    _run(runner, workdir, "annotate")

    for stage in ("prepare", "embed", "fold", "search", "context", "rank", "report", "annotate"):
        assert (workdir / stage_dir(stage) / "_DONE").exists(), f"{stage} marker missing"

    assert (workdir / "00_prepare" / "cds_records.json").exists()
    assert (workdir / "03_cluster" / "families.json").exists()
    rank_dir = workdir / "07_rank"
    for name in ("predictions.json", "calibration.json", "priorities.tsv"):
        assert (rank_dir / name).exists(), name
    report_dir = workdir / "08_report"
    for name in ("evidence.json", "report.md", "report.html"):
        assert (report_dir / name).exists(), name
    assert (workdir / "09_annotate" / "changes.tsv").exists()


def stage_dir(stage: str) -> str:
    mapping = {
        "prepare": "00_prepare",
        "embed": "02_embed",
        "fold": "04_fold",
        "search": "05_search",
        "context": "06_context",
        "rank": "07_rank",
        "report": "08_report",
        "annotate": "09_annotate",
    }
    return mapping[stage]
