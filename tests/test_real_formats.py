"""Real-world format tests for Bakta .gff3 and Prokka .gff outputs.

These fixtures imitate the structural quirks of the two annotators:
Prokka uses a ``.gff`` extension, ``PROKKA_``-prefixed ids, an embedded
``##FASTA`` section and percent-encoded spaces in product strings. Bakta uses
``.gff3`` with ``locus_tag`` ids. The prepare stage must parse both, align
every CDS to its .faa sequence and classify the product strings.
"""

from pathlib import Path

from pufprobe.cli import _collect_inputs
from pufprobe.config import load_config
from pufprobe.io_utils import read_json
from pufprobe.models import GeneStatus, cds_from_dict
from pufprobe.prepare import gff
from pufprobe.workdir import WorkDir
from pufprobe.prepare.pipeline import run as prepare_run

_PROKKA_GFF = """##gff-version 3
##sequence-region contig001 1 5000
contig001\tProkka\tgene\t1\t300\t.\t+\t.\tID=PROKKA_00001;locus_tag=PROKKA_00001
contig001\tProkka\tCDS\t1\t300\t.\t+\t0\tID=PROKKA_00001;Parent=PROKKA_00001;Name=dnaN;gene=dnaN;locus_tag=PROKKA_00001;product=DNA%20polymerase%20III%20subunit%20beta
contig001\tProkka\tgene\t350\t500\t.\t+\t.\tID=PROKKA_00002;locus_tag=PROKKA_00002
contig001\tProkka\tCDS\t350\t500\t.\t+\t0\tID=PROKKA_00002;Parent=PROKKA_00002;locus_tag=PROKKA_00002;product=hypothetical%20protein
contig001\tProkka\tgene\t550\t750\t.\t-\t.\tID=PROKKA_00003;locus_tag=PROKKA_00003
contig001\tProkka\tCDS\t550\t750\t.\t-\t0\tID=PROKKA_00003;Parent=PROKKA_00003;locus_tag=PROKKA_00003;product=DUF1234%20domain-containing%20protein
##FASTA
>contig001
ATGCCCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""

_PROKKA_FAA = """\
>PROKKA_00001 DNA polymerase III subunit beta
MPDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR
>PROKKA_00002 hypothetical protein
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKISVGDEVEVKPEM
>PROKKA_00003 DUF1234 domain-containing protein
MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDAALNLA
"""

_BAKTA_GFF3 = """##gff-version 3
##sequence-region contig1 1 4000
contig1\tbakta\tCDS\t1\t400\t.\t+\t0\tID=SAL001_00001;locus_tag=SAL001_00001;Name=gyrA;product=DNA%20gyrase%20subunit%20A;Dbxref=COG:COG0188
contig1\tbakta\tCDS\t500\t700\t.\t-\t0\tID=SAL001_00002;locus_tag=SAL001_00002;Name=hypothetical protein;product=hypothetical protein
contig1\tbakta\tCDS\t800\t990\t.\t+\t0\tID=SAL001_00003;locus_tag=SAL001_00003;product=uncharacterized%20protein
"""

_BAKTA_FAA = """\
>SAL001_00001 DNA gyrase subunit A
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKISVGDEVEVKPEM
>SAL001_00002 hypothetical protein
MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDAALNLA
>SAL001_00003 uncharacterized protein
MARGKKIGYSGPRKGTGRKGTAAKGTQVEVKPEMGQDAALNLAKQRQISFVKSHFSRQLEERLGLI
"""


def _write_fixture(root: Path, sample: str, gff_name: str, gff: str, faa: str) -> Path:
    directory = root / sample
    directory.mkdir(parents=True)
    stem = gff_name.split(".")[0]
    (directory / gff_name).write_text(gff, encoding="utf-8", newline="\n")
    (directory / f"{stem}.faa").write_text(faa, encoding="utf-8", newline="\n")
    return directory


def test_prokka_gff_parses_and_aligns(tmp_path: Path) -> None:
    directory = _write_fixture(tmp_path, "prokka_sample", "prokka_sample.gff",
                               _PROKKA_GFF, _PROKKA_FAA)
    gff_path = directory / "prokka_sample.gff"
    assert gff.detect_source(gff_path) == "prokka"

    wd = WorkDir(tmp_path / "run")
    config = load_config(None)
    prepare_run(wd, config, [gff_path, directory / "prokka_sample.faa"])

    rows = read_json(wd.path("prepare", "cds_records.json"))
    records = [cds_from_dict(row) for row in rows]
    assert len(records) == 3
    assert all(record.sequence_aa for record in records)
    assert records[0].status == GeneStatus.ANNOTATED
    assert records[1].status == GeneStatus.HYPOTHETICAL
    assert records[2].status == GeneStatus.CONSERVED_HYPOTHETICAL

    qc = read_json(wd.path("qc", "qc_summary.json"))
    assert qc["unaligned"] == []


def test_bakta_gff3_parses_and_aligns(tmp_path: Path) -> None:
    directory = _write_fixture(tmp_path, "bakta_sample", "bakta_sample.gff3",
                               _BAKTA_GFF3, _BAKTA_FAA)
    gff_path = directory / "bakta_sample.gff3"
    assert gff.detect_source(gff_path) == "bakta"

    wd = WorkDir(tmp_path / "run")
    config = load_config(None)
    prepare_run(wd, config, [gff_path, directory / "bakta_sample.faa"])

    rows = read_json(wd.path("prepare", "cds_records.json"))
    records = [cds_from_dict(row) for row in rows]
    assert len(records) == 3
    assert all(record.sequence_aa for record in records)
    assert records[0].status == GeneStatus.ANNOTATED
    assert records[1].status == GeneStatus.HYPOTHETICAL
    assert records[2].status == GeneStatus.HYPOTHETICAL

    qc = read_json(wd.path("qc", "qc_summary.json"))
    assert qc["unaligned"] == []


def test_gbff_files_are_discovered(tmp_path: Path) -> None:
    root = tmp_path / "tree"
    sample = root / "genome1"
    sample.mkdir(parents=True)
    (sample / "genome1.gff3").write_text("##gff-version 3\n", encoding="utf-8")
    (sample / "genome1.gbff").write_text("LOCUS genome1\n", encoding="utf-8")
    (sample / "genome1.faa").write_text(">a\nMKK\n", encoding="utf-8")

    paths = {p.name for p in _collect_inputs([root], None)}
    assert "genome1.gff3" in paths
    assert "genome1.gbff" in paths
    assert "genome1.faa" in paths
