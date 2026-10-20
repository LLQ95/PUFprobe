"""Regression: Bakta's extra <prefix>.hypothetical.faa subset must not be indexed.

Bakta writes both a complete ``<prefix>.faa`` and a ``<prefix>.hypothetical.faa``
holding only the hypothetical proteins. The prepare stage indexes the complete
set; if the subset overrides it, every annotated protein is reported as
unaligned.
"""

from pathlib import Path

from pufprobe.cli import _collect_inputs
from pufprobe.config import load_config
from pufprobe.io_utils import read_json
from pufprobe.prepare.pipeline import run as prepare_run
from pufprobe.workdir import WorkDir


def _make_sample(root: Path) -> Path:
    sample = root / "GEN"
    sample.mkdir(parents=True)
    (sample / "GEN.gff3").write_text(
        "##gff-version 3\n"
        "contig1\tbakta\tCDS\t1\t90\t.\t+\t0\t"
        "ID=GEN_00001;locus_tag=GEN_00001;product=DNA polymerase\n"
        "contig1\tbakta\tCDS\t100\t160\t.\t+\t0\t"
        "ID=GEN_00002;locus_tag=GEN_00002;product=hypothetical protein\n"
        "contig1\tbakta\tCDS\t200\t290\t.\t+\t0\t"
        "ID=GEN_00003;locus_tag=GEN_00003;product=ribosomal protein\n",
        encoding="utf-8",
    )
    (sample / "GEN.faa").write_text(
        ">GEN_00001\n" + "MA" + "AAA" * 12 + "\n"
        ">GEN_00002\n" + "ML" + "LLL" * 12 + "\n"
        ">GEN_00003\n" + "MR" + "RRR" * 12 + "\n",
        encoding="utf-8",
    )
    # The subset contains only the hypothetical protein.
    (sample / "GEN.hypothetical.faa").write_text(
        ">GEN_00002\n" + "ML" + "LLL" * 12 + "\n",
        encoding="utf-8",
    )
    return sample


def test_collect_inputs_excludes_hypothetical_subset(tmp_path: Path) -> None:
    sample = _make_sample(tmp_path)
    paths = _collect_inputs([sample], None)
    names = [p.name for p in paths]
    assert "GEN.faa" in names
    assert "GEN.hypothetical.faa" not in names
    assert "GEN.gff3" in names


def test_prepare_aligns_all_proteins_with_subset_present(tmp_path: Path) -> None:
    sample = _make_sample(tmp_path)
    paths = _collect_inputs([sample], None)

    wd = WorkDir(tmp_path / "run")
    prepare_run(wd, load_config(None), paths)

    records = read_json(wd.path("prepare", "cds_records.json"))  # type: ignore[arg-type]
    assert len(records) == 3
    # All CDS, including the two annotated proteins, must carry a sequence.
    assert all(r["sequence_aa"] for r in records)
    qc = read_json(wd.path("qc", "qc_summary.json"))  # type: ignore[arg-type]
    assert qc["unaligned"] == []
