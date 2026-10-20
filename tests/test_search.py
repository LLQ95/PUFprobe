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

COLS = "query,target,pident,evalue,qstart,qend,tstart,tend,qcov,tcov,theader"


def test_parse_convertalis_evidence() -> None:
    rows = [
        "sA@g1\tX1\t45.2\t1e-20\t1\t50\t1\t50\t0.9\t0.85\tsp|X1|NAME_ONE some protein desc",
        "sA@g1\tX2\t20.0\t1.0\t1\t10\t1\t10\t0.2\t0.2\tsp|X2|NAME_TWO bad evalue",
        "sA@g1\tX3\t55.0\t1e-25\t1\t50\t1\t50\t0.9\t0.9\tsp|X3|NAME_THR hypothetical protein",
    ]
    tsv = Path("dummy.tsv")
    tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")
    try:
        ev = mmseqs_mod.parse_convertalis(tsv, report_evalue=1e-3)
    finally:
        tsv.unlink()
    assert set(ev.keys()) == {"sA@g1"}
    hits = ev["sA@g1"]
    # X2 filtered by e-value; X3 is a strong hit but the target itself has no
    # function (hypothetical protein), so it cannot anchor a function term.
    assert len(hits) == 1
    h = hits[0]
    assert h.channel.value == "remote_homology"
    assert h.target_id == "X1"
    assert h.target_label == "some protein desc"
    assert h.evalue == pytest.approx(1e-20)
    assert h.coverage == pytest.approx(0.9)
    assert 0.0 <= h.score <= 1.0
    # identity strength 0.84 (pident 45.2) x coverage 0.9 x significance 1.0
    assert h.score == pytest.approx(0.756, abs=1e-3)


def test_uninformative_labels() -> None:
    assert mmseqs_mod.is_uninformative("hypothetical protein") is True
    assert mmseqs_mod.is_uninformative("Uncharacterized protein OS=Escherichia coli") is True
    assert mmseqs_mod.is_uninformative("protein of unknown function") is True
    # A trailing ORF/locus name does not add a function.
    assert mmseqs_mod.is_uninformative(
        "Hypothetical protein At1g04090 OS=Arabidopsis thaliana") is True
    assert mmseqs_mod.is_uninformative(
        "Uncharacterized protein 002L OS=Frog virus 3") is True
    assert mmseqs_mod.is_uninformative("Putative protein") is True
    # A specific name (even after a qualifier) still carries function info.
    assert mmseqs_mod.is_uninformative("Putative protein IntB") is False
    assert mmseqs_mod.is_uninformative("Tyrosine recombinase XerC") is False
    assert mmseqs_mod.is_uninformative("Probable protease SohB") is False


def test_description_extraction_from_header() -> None:
    # Full UniProt header as emitted by the theader column.
    assert mmseqs_mod._description(
        "P0AG14", "sp|P0AG14|SOHB_ECOLI Probable protease SohB OS=Escherichia coli"
    ) == "Probable protease SohB OS=Escherichia coli"
    # Header with no description: the whole header is returned.
    assert mmseqs_mod._description("X1", "sp|X1|NAME") == "sp|X1|NAME"
    # Empty header: fall back to the target identifier.
    assert mmseqs_mod._description("X1", "") == "X1"


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
    # Score is TM-score x query coverage x significance (prob is retained only
    # in details): 0.85 x 0.90 x 1.0 (E=1e-15) = 0.765.
    assert h.score == pytest.approx(0.765, abs=1e-3)
    assert h.details.get("prob") == pytest.approx(98.0)


def test_parse_foldseek_tsv_with_sifts() -> None:
    fam_to_rep = {"FAM00001": "sA@g1"}
    rows = ["FAM00001\t101m_A\t0.90\t95.0\t1e-10\t0.9\t0.9"]
    tsv = Path("fs2.tsv")
    tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")
    sifts = {"101m_A": ["PF00042"]}
    try:
        ev = foldseek_mod.parse_foldseek_tsv(tsv, fam_to_rep, sifts)
    finally:
        tsv.unlink()
    h = ev["sA@g1"][0]
    assert h.target_id == "PF00042"
    assert h.details.get("db") == "pfam"
    assert h.details.get("pdb_target") == "101m_A"


def test_parse_foldseek_tsv_bare_pdb_id_with_sifts() -> None:
    # Foldseek writes a bare PDB id for single-chain-type structures.
    fam_to_rep = {"FAM00001": "sA@g1"}
    rows = ["FAM00001\t101m\t0.90\t95.0\t1e-10\t0.9\t0.9"]
    tsv = Path("fs3.tsv")
    tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")
    sifts = {"101m_A": ["PF00042"], "101m": ["PF00042"]}
    try:
        ev = foldseek_mod.parse_foldseek_tsv(tsv, fam_to_rep, sifts)
    finally:
        tsv.unlink()
    h = ev["sA@g1"][0]
    assert h.target_id == "PF00042"
    assert h.details.get("pdb_target") == "101m"


def test_foldseek_structure_quality_gates() -> None:
    fam_to_rep = {"FAM00001": "sA@g1"}
    # Columns: query, target, alntmscore, prob, evalue, qcov, tcov.
    rows = [
        # Below the same-fold TM-score boundary: dropped even with prob=1.
        "FAM00001\t1aaa_A\t0.40\t1.0\t1e-10\t0.9\t0.9",
        # Non-significant E-value: dropped despite a high TM-score.
        "FAM00001\t1bbb_A\t0.90\t1.0\t5.0\t0.9\t0.9",
        # Local fragment: high TM-score over only ~35% of the query.
        "FAM00001\t7ax1_A\t0.956\t1.0\t8.32e-05\t0.352\t0.197",
        # Strong domain-level hit: high TM-score, high query coverage, significant.
        "FAM00001\t9kb6\t0.814\t1.0\t1.99e-09\t0.962\t0.416",
    ]
    tsv = Path("fsg.tsv")
    tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")
    try:
        ev = foldseek_mod.parse_foldseek_tsv(tsv, fam_to_rep)
    finally:
        tsv.unlink()
    by_target = {e.details.get("pdb_target"): e for e in ev["sA@g1"]}
    assert set(by_target) == {"7ax1_A", "9kb6"}
    # Local fragment kept but down-weighted: 0.956 x 0.352 x 0.408 = 0.138.
    assert by_target["7ax1_A"].score == pytest.approx(0.138, abs=2e-3)
    # Strong domain-level hit: 0.814 x 0.962 x 0.870 = 0.682.
    assert by_target["9kb6"].score == pytest.approx(0.682, abs=2e-3)
    # prob is no longer the score but stays available for auditing.
    assert by_target["9kb6"].details.get("prob") == pytest.approx(1.0)


def test_structure_score_formula() -> None:
    # fold x coverage x significance, clamped to [0, 1].
    assert foldseek_mod._structure_score(0.9, 0.9, 1e-15) == pytest.approx(0.81)
    assert foldseek_mod._structure_score(0.0, 0.9, 1e-15) == 0.0
    assert foldseek_mod._structure_score(1.5, 1.0, 0.0) == 1.0


def test_load_sifts_pdb_chain_pfam(tmp_path: Path) -> None:
    from pufprobe.search import sifts as sifts_mod

    raw = (
        "# 2026/10/11 - 19:53\n"
        "PDB\tCHAIN\tSP_PRIMARY\tPFAM_ID\tCOVERAGE\n"
        "101m\tA\tP02185\tPF00042\t1\n"
        "101m\tB\tP02185\tPF00042\t1\n"
        "102m\tA\tP02185\tPF00042\t1\n"
        "1abc\tA\tP12345\tPF00001\t1\n"
        "201l\tA\tP11111\tPF00001\t1\n"
        "201l\tB\tP22222\tPF00002\t1\n"
    )
    gz = tmp_path / "pdb_chain_pfam.tsv.gz"
    import gzip as _gzip
    with _gzip.open(gz, "wt", encoding="utf-8") as handle:
        handle.write(raw)
    mapping = sifts_mod.load_pdb_chain_pfam(gz)
    assert mapping["101m_A"] == ["PF00042"]
    assert mapping["1abc_A"] == ["PF00001"]
    # Bare PDB id exposes the union of Pfam families across chains.
    assert mapping["101m"] == ["PF00042"]
    assert mapping["201l"] == ["PF00001", "PF00002"]
    assert sifts_mod.pdb_chain_key("101m") == "101m"
    assert sifts_mod.pdb_chain_key("1abc_A_2") == "1abc_A"


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


def test_foldseek_finds_models_from_absolute_structures_dir(tmp_path: Path, monkeypatch) -> None:
    # The fold stage stores Family.structure_path relative to the run workdir,
    # while the search process runs from the repository root. Structures must
    # still be resolved against the absolute structures directory.
    from pufprobe import externals

    wd = WorkDir(tmp_path / "run")
    structures_dir = wd.stage_dir("fold") / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)
    (structures_dir / "FAM00001.pdb").write_text("ATOM\n", encoding="utf-8")

    dbdir = tmp_path / "dbdir"
    dbdir.mkdir(parents=True, exist_ok=True)
    (dbdir / "foldseek_db.index").write_text("idx\n", encoding="utf-8")

    out_dir = wd.stage_dir("search")
    fam = Family(
        family_id="FAM00001", representative_key="sA@g1",
        members=["sA@g1"], representative_seq="M" * 50,
        size=1, lineage_width=1,
        fold_status="folded", plddt=90.0,
        structure_path="04_fold/structures/FAM00001.pdb",
    )

    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    monkeypatch.setattr(
        externals, "foldseek",
        lambda: ToolInfo(name="foldseek", available=True, path="foldseek"),
    )
    calls: list = []
    monkeypatch.setattr(foldseek_mod, "_run", lambda cmd: calls.append(cmd))

    foldseek_mod.run_foldseek_search(
        [fam], structures_dir, dbdir, out_dir, min_plddt=70.0,
    )

    assert calls, "foldseek must be invoked, not skipped on a relative-path check"
    assert calls[0][1] == "createdb"
    assert (out_dir / "_foldseek_query" / "FAM00001.pdb").exists()
    search_cmd = next(c for c in calls if c[1] == "search")
    assert "-a" in search_cmd, (
        "search must request the alignment backtrace (-a) so convertalis can "
        "emit alntmscore/qcov/tcov without a 'cigar is missing' failure"
    )
