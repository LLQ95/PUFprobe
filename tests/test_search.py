"""search: evidence parsing from raw TSV/HHR, missing-tool skip policy."""

import json
from pathlib import Path

import pytest

from pufprobe.config import load_config
from pufprobe.externals import ToolInfo
from pufprobe.io_utils import read_json, write_json
from pufprobe.models import Family
from pufprobe.search import foldseek as foldseek_mod
from pufprobe.search import hhblits as hhblits_mod
from pufprobe.search import mmseqs_search as mmseqs_mod
from pufprobe.search.pipeline import run as search_run
from pufprobe.workdir import WorkDir

COLS = "query,target,pident,evalue,qstart,qend,tstart,tend,qcov,tcov,desc"


def test_parse_convertalis_evidence() -> None:
    rows = [
        "sA@g1\tsp|X1|\t45.2\t1e-20\t1\t50\t1\t50\t0.9\t0.85\tsome protein desc",
        "sA@g1\tsp|X2|\t20.0\t1.0\t1\t10\t1\t10\t0.2\t0.2\tbad evalue",
    ]
    tsv = Path("dummy.tsv")
    tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")
    try:
        ev = mmseqs_mod.parse_convertalis(tsv, report_evalue=1e-3)
    finally:
        tsv.unlink()
    assert set(ev.keys()) == {"sA@g1"}
    hits = ev["sA@g1"]
    assert len(hits) == 1  # second row filtered by e-value
    h = hits[0]
    assert h.channel.value == "remote_homology"
    assert h.target_id == "sp|X1|"
    assert h.target_label == "some protein desc"
    assert h.evalue == pytest.approx(1e-20)
    assert h.coverage == pytest.approx(0.9)
    assert 0.0 <= h.score <= 1.0
    assert h.score == pytest.approx(1.0)  # -log10(1e-20)=20 -> capped at 1


def test_parse_foldseek_tsv() -> None:
    fam_to_rep = {"FAM00001": "sA@g1"}
    rows = ["FAM00001\tAF-A0A-1\t0.85\t98.0\t1e-15\t0.9\t0.88"]
    tsv = Path("fs.tsv")
    tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")
    try:
        ev = foldseek_mod.parse_foldseek_tsv(tsv, fam_to_rep)
    finally:
        tsv.unlink()
    assert "sA@g1" in ev
    h = ev["sA@g1"][0]
    assert h.channel.value == "foldseek"
    assert h.score == pytest.approx(0.98)  # prob 98 -> /100


def test_parse_hhr_table() -> None:
    hhr = [
        "Query:      q",
        "Match_columns 50 No_seqs_used 1",
        "No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM",
        "  1 sp|P1|name_one                95.0 1.2E-45  1.2E-45  450.0  0.00   45   1-45     1-45",
        "  2 sp|P2|name_two                 12.0 5.0E+01  5.0E+01   50.0  0.00   10   1-10     1-10",
        ">sp|P3|not_a_table_row",
    ]
    p = Path("x.hhr")
    p.write_text("\n".join(hhr) + "\n", encoding="utf-8")
    try:
        hits = hhblits_mod.parse_hhr(p)
    finally:
        p.unlink()
    assert [h[0] for h in hits] == ["sp|P1|name_one", "sp|P2|name_two"]
    assert hits[0][1] == pytest.approx(95.0)
    assert hits[0][2] == pytest.approx(1.2e-45)


def _minimal_workdir(tmp_path: Path) -> WorkDir:
    wd = WorkDir(tmp_path / "run")
    families = [
        Family(family_id="FAM00001", representative_key="sA@g1",
               members=["sA@g1"], representative_seq="M" * 50,
               size=1, lineage_width=1, fold_status="folded", plddt=90.0),
    ]
    write_json(wd.path("cluster", "families.json"),
               [json.loads(json.dumps(f.__dict__)) for f in families])
    wd.path("cluster", "representatives.faa").write_text(
        ">sA@g1\n" + "M" * 50 + "\n", encoding="utf-8"
    )
    return wd


def test_search_pipeline_skips_all_missing_tools(tmp_path: Path, monkeypatch) -> None:
    from pufprobe import externals
    monkeypatch.setattr(externals, "mmseqs",
                        lambda: ToolInfo(name="mmseqs", available=False))
    monkeypatch.setattr(externals, "foldseek",
                        lambda: ToolInfo(name="foldseek", available=False))
    monkeypatch.setattr(externals, "hhblits",
                        lambda: ToolInfo(name="hhblits", available=False))
    monkeypatch.setattr(externals, "tmvec",
                        lambda: ToolInfo(name="tmvec", available=False))

    wd = _minimal_workdir(tmp_path)
    search_run(wd, load_config(None), dbdir=None)

    out = read_json(wd.path("search", "search_evidence.json"))  # type: ignore[arg-type]
    assert out == {}
    assert wd.marker("search").is_file()
