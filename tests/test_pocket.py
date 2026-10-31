"""Tests for the pocket catalytic-template evidence channel.

The channel depends on P2Rank, foldseek and a prebuilt M-CSA template
database, none of which are needed here. PDB/P2Rank parsing, alignment
mapping and template evaluation are tested on synthetic inputs, and the
online orchestration is exercised with a monkeypatched subprocess that emits
deterministic P2Rank and foldseek output.
"""

import json
from pathlib import Path

import pytest

from pufprobe.externals import ToolInfo
from pufprobe.models import EvidenceChannel, Family
from pufprobe.search import pocket


# ---------------------------------------------------------------------------
# PDB fixture helpers
# ---------------------------------------------------------------------------

def _atom_line(chain, resseq, resname, icode="", serial=1):
    """Build an 80-column PDB ATOM record with exact fixed columns."""
    line = [" "] * 80

    def put(start, text):
        for i, ch in enumerate(text):
            line[start - 1 + i] = ch

    put(1, "ATOM")
    put(7, str(serial).rjust(5))
    put(13, " CA")
    put(17, " ")
    put(18, resname.rjust(3))
    put(22, chain)
    put(23, str(resseq).rjust(4))
    put(27, icode)
    put(31, "0.000".rjust(8))
    put(39, "0.000".rjust(8))
    put(47, "0.000".rjust(8))
    put(55, "1.00".rjust(6))
    put(61, "20.00".rjust(6))
    put(77, resname[:2].rjust(2))
    return "".join(line)


def _chain_pdb(chain, resnames, first_auth=1):
    lines = [
        _atom_line(chain, first_auth + i, name, serial=i + 1)
        for i, name in enumerate(resnames)
    ]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------

def test_parse_pdb_residues_seq_index_and_author_number():
    # Author numbers have a gap; sequence positions stay contiguous.
    text = "\n".join(
        [_atom_line("A", 10, "MET", serial=1),
         _atom_line("A", 12, "GLY", serial=2)]
    ) + "\n"
    residues = pocket.parse_pdb_residues(text)
    assert [r.seq_index for r in residues] == [1, 2]
    assert [r.auth_resid for r in residues] == [10, 12]
    assert [r.resname for r in residues] == ["MET", "GLY"]
    assert all(r.chain == "A" for r in residues)


def test_parse_pdb_residues_independent_chains():
    text = (
        _atom_line("A", 1, "MET", serial=1)
        + "\n"
        + _atom_line("B", 1, "LYS", serial=2)
        + "\n"
    )
    residues = pocket.parse_pdb_residues(text)
    by_chain = {r.chain: r for r in residues}
    assert set(by_chain) == {"A", "B"}
    assert by_chain["A"].seq_index == 1
    assert by_chain["B"].seq_index == 1


def test_parse_pdb_residues_dedupes_multiple_atoms():
    # Two ATOM records of the same residue (CA and CB) count once.
    first = _atom_line("A", 1, "MET", serial=1)
    second = _atom_line("A", 1, "MET", serial=2)
    residues = pocket.parse_pdb_residues(first + "\n" + second + "\n")
    assert len(residues) == 1


def test_parse_pdb_residues_records_insertion_code():
    text = _atom_line("A", 12, "GLY", icode="A", serial=1) + "\n"
    residues = pocket.parse_pdb_residues(text)
    assert residues[0].icode == "A"
    assert residues[0].auth_resid == 12


def test_build_key_map_and_index():
    text = _chain_pdb("A", ["MET", "GLY"])
    residues = pocket.parse_pdb_residues(text)
    key_map = pocket.build_residue_key_map(residues)
    assert key_map[("A", 1, "")] == 1
    assert key_map[("A", 2, "")] == 2
    by_index = pocket.index_by_seq(residues)
    assert by_index[1].resname == "MET"
    assert by_index[2].resname == "GLY"


def test_extract_chain_pdb_keeps_only_target_chain():
    text = (
        _atom_line("A", 1, "MET", serial=1)
        + "\n"
        + _atom_line("B", 1, "LYS", serial=2)
        + "\n"
    )
    extracted = pocket.extract_chain_pdb(text, "A")
    assert "MET" in extracted
    assert "LYS" not in extracted
    assert extracted.startswith("ATOM")
    assert "TER" in extracted


# ---------------------------------------------------------------------------
# P2Rank parsing
# ---------------------------------------------------------------------------

def test_parse_residue_token():
    assert pocket.parse_residue_token("A_123") == ("A", 123, "")
    assert pocket.parse_residue_token("A_123A") == ("A", 123, "A")
    assert pocket.parse_residue_token("bad") is None


def _predictions_csv(path, rows):
    header = (
        "name,rank,score,probability,sas_points,surf_atoms,center_x,"
        "center_y,center_z,residue_ids,surf_atom_ids"
    )
    lines = [header]
    for name, rank, score, prob, ids in rows:
        lines.append(
            f"{name},{rank},{score},{prob},10,20,1.0,2.0,3.0,{ids},1 2 3"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_parse_p2rank_predictions():
    csv_path = Path("p2rank_predictions.csv")
    _predictions_csv(
        csv_path,
        [("prot.pdb", 1, 3.4, 0.9, "A_2 A_5"),
         ("prot.pdb", 2, 2.1, 0.7, "A_8")],
    )
    try:
        pockets = pocket.parse_p2rank_predictions(csv_path)
    finally:
        csv_path.unlink()
    assert [p.rank for p in pockets] == [1, 2]
    assert pockets[0].score == pytest.approx(3.4)
    assert pockets[0].probability == pytest.approx(0.9)
    assert pockets[0].residue_keys == {("A", 2, ""), ("A", 5, "")}
    assert pockets[1].residue_keys == {("A", 8, "")}


def test_parse_p2rank_predictions_missing_file():
    assert pocket.parse_p2rank_predictions(Path("absent.csv")) == []


def test_pocket_seq_sets_truncates_and_drops_unknown_keys():
    pockets = [
        pocket.PocketPrediction(
            "p", rank=i + 1, score=1.0, probability=0.5,
            residue_keys={("A", i + 1, "")},
        )
        for i in range(4)
    ]
    key_map = {("A", i + 1, ""): i + 1 for i in range(4)}
    sets = pocket.pocket_seq_sets(pockets, key_map, max_pockets=2)
    assert len(sets) == 2
    assert sets[0][1] == {1}
    assert sets[1][1] == {2}

    # A key absent from the query map is dropped without error.
    pockets[0].residue_keys.add(("Z", 9, ""))
    sets = pocket.pocket_seq_sets(pockets, key_map, max_pockets=1)
    assert sets[0][1] == {1}


# ---------------------------------------------------------------------------
# Alignment mapping
# ---------------------------------------------------------------------------

def test_map_alignment_no_gaps():
    mapping = pocket.map_alignment("AAAA", "AAAA", 1, 1)
    assert mapping == {1: 1, 2: 2, 3: 3, 4: 4}


def test_map_alignment_with_gaps():
    # Query "AA-AA" (4 residues), target "AAA-A" (4 residues): each side has a
    # gap in a different column, so gap columns contribute nothing.
    mapping = pocket.map_alignment("AA-AA", "AAA-A", 1, 1)
    # col1/2 map t1->q1, t2->q2; col3 is a query gap; col4 is a target gap;
    # col5 maps t4->q4.
    assert mapping == {1: 1, 2: 2, 4: 4}


def test_map_alignment_start_offsets():
    mapping = pocket.map_alignment("AA", "AA", 5, 7)
    assert mapping == {7: 5, 8: 6}


# ---------------------------------------------------------------------------
# Template evaluation
# ---------------------------------------------------------------------------

def _residues_by_index(resnames):
    text = _chain_pdb("A", resnames)
    return pocket.index_by_seq(pocket.parse_pdb_residues(text))


def _pocket_sets(seq_sets):
    out = []
    for rank, seq_set in enumerate(seq_sets, 1):
        pred = pocket.PocketPrediction(
            "p", rank=rank, score=1.0, probability=0.9,
            residue_keys={("A", pos, "") for pos in seq_set},
        )
        out.append((pred, seq_set))
    return out


def test_evaluate_template_full_hit():
    # Query residues 2 and 5 are Ser and His, matching the template.
    by_index = _residues_by_index(
        ["MET", "SER", "GLY", "GLY", "HIS", "GLY"]
    )
    template = {
        "target": "1abc_A",
        "catalytic_residues": [
            {"resid": 2, "code": "Ser"},
            {"resid": 5, "code": "His"},
        ],
    }
    pocket_sets = _pocket_sets([{2, 5}])
    result = pocket.evaluate_template(
        template, {2: 2, 5: 5}, by_index, pocket_sets,
        tmscore=0.8, prob=99.0, evalue=1e-10, qcov=0.9, tcov=0.8,
        min_pocket_consistency=0.5,
    )
    assert result is not None
    assert result.catalytic_total == 2
    assert result.catalytic_mapped == 2
    assert result.catalytic_in_pocket == 2
    assert result.consistency == pytest.approx(1.0)
    # structure 0.8 x 0.9 x 1.0 = 0.72; score equals it at consistency 1.
    assert result.structure_score == pytest.approx(0.72, abs=1e-3)
    assert result.score == pytest.approx(0.72, abs=1e-3)


def test_evaluate_template_drops_resname_mismatch():
    # Template expects Ser at position 2, but the query has Lys.
    by_index = _residues_by_index(["MET", "LYS", "GLY", "GLY", "HIS"])
    template = {
        "catalytic_residues": [
            {"resid": 2, "code": "Ser"},
            {"resid": 5, "code": "His"},
        ],
    }
    pocket_sets = _pocket_sets([{2, 5}])
    result = pocket.evaluate_template(
        template, {2: 2, 5: 5}, by_index, pocket_sets,
        tmscore=0.8, prob=99.0, evalue=1e-10, qcov=0.9, tcov=0.8,
        min_pocket_consistency=0.5,
    )
    # Only His maps; it is in the pocket, so consistency over the one mapped
    # residue is 1.0 and the hit survives with a single mapped residue.
    assert result is not None
    assert result.catalytic_mapped == 1
    assert result.catalytic_in_pocket == 1


def test_evaluate_template_no_pocket_hit_returns_none():
    by_index = _residues_by_index(["MET", "SER", "GLY", "GLY", "HIS"])
    template = {
        "catalytic_residues": [
            {"resid": 2, "code": "Ser"},
            {"resid": 5, "code": "His"},
        ],
    }
    # Pocket covers unrelated residues.
    pocket_sets = _pocket_sets([{3, 4}])
    result = pocket.evaluate_template(
        template, {2: 2, 5: 5}, by_index, pocket_sets,
        tmscore=0.8, prob=99.0, evalue=1e-10, qcov=0.9, tcov=0.8,
        min_pocket_consistency=0.5,
    )
    assert result is None


def test_evaluate_template_consistency_gate():
    by_index = _residues_by_index(
        ["MET", "SER", "GLY", "GLY", "HIS", "GLY"]
    )
    template = {
        "catalytic_residues": [
            {"resid": 2, "code": "Ser"},
            {"resid": 5, "code": "His"},
        ],
    }
    # Only one of the two mapped residues is in the pocket: consistency 0.5.
    pocket_sets = _pocket_sets([{2}])
    result = pocket.evaluate_template(
        template, {2: 2, 5: 5}, by_index, pocket_sets,
        tmscore=0.8, prob=99.0, evalue=1e-10, qcov=0.9, tcov=0.8,
        min_pocket_consistency=0.6,
    )
    assert result is None


def test_evaluate_template_tie_breaks_to_lower_rank():
    by_index = _residues_by_index(["MET", "SER", "GLY", "GLY", "HIS"])
    template = {
        "catalytic_residues": [
            {"resid": 2, "code": "Ser"},
            {"resid": 5, "code": "His"},
        ],
    }
    # Both pockets contain one mapped residue; rank 1 must win.
    pocket_sets = _pocket_sets([{2}, {5}])
    result = pocket.evaluate_template(
        template, {2: 2, 5: 5}, by_index, pocket_sets,
        tmscore=0.8, prob=99.0, evalue=1e-10, qcov=0.9, tcov=0.8,
        min_pocket_consistency=0.5,
    )
    assert result is not None
    assert result.pocket_rank == 1
    assert result.catalytic_in_pocket == 1


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def test_evalue_score():
    assert pocket._evalue_score(0.0) == 1.0
    assert pocket._evalue_score(1e-10) == 1.0
    assert pocket._evalue_score(1e-5) == pytest.approx(0.5)
    assert pocket._evalue_score(100.0) == 0.0


def test_structure_score():
    assert pocket._structure_score(0.8, 0.9, 1e-10) == pytest.approx(
        0.72, abs=1e-3
    )
    assert pocket._structure_score(0.0, 0.9, 1e-10) == 0.0
    assert pocket._structure_score(1.5, 1.5, 0.0) == 1.0


# ---------------------------------------------------------------------------
# Template index loading
# ---------------------------------------------------------------------------

def test_load_templates_list_and_dict(tmp_path):
    list_path = tmp_path / "list.json"
    list_path.write_text(
        json.dumps([{"target": "1abc_A", "ec": "1.1.1.1"}]),
        encoding="utf-8",
    )
    loaded = pocket._load_templates(list_path)
    assert "1abc_A" in loaded
    assert loaded["1abc_A"]["ec"] == "1.1.1.1"

    dict_path = tmp_path / "dict.json"
    dict_path.write_text(
        json.dumps({"1abc_A": {"ec": "1.1.1.1"}}),
        encoding="utf-8",
    )
    loaded = pocket._load_templates(dict_path)
    assert loaded["1abc_A"]["target"] == "1abc_A"


# ---------------------------------------------------------------------------
# Online orchestration
# ---------------------------------------------------------------------------

def _family(family_id="FAM1"):
    return Family(
        family_id=family_id,
        representative_key=f"{family_id}@g1",
        members=[f"{family_id}@g1"],
        representative_seq="M" * 10,
        size=1,
        lineage_width=1,
        fold_status="folded",
        plddt=90.0,
    )


def _setup_dbdir(dbdir):
    dbdir.mkdir(parents=True, exist_ok=True)
    (dbdir / "pocket_db.index").write_text("idx\n", encoding="utf-8")
    templates = [
        {
            "target": "1abc_A",
            "mcsa_id": 7,
            "ec": "3.5.2.6",
            "enzyme_name": "beta-lactamase",
            "catalytic_residues": [
                {"resid": 2, "code": "Ser"},
                {"resid": 5, "code": "His"},
            ],
        }
    ]
    (dbdir / "pocket_templates.json").write_text(
        json.dumps(templates), encoding="utf-8"
    )


def _write_structure(structures_dir, family_id):
    structures_dir.mkdir(parents=True, exist_ok=True)
    # Residues 2 and 5 are Ser and His to match the template.
    text = _chain_pdb(
        "A",
        ["MET", "SER", "GLY", "GLY", "HIS", "GLY",
         "GLY", "GLY", "GLY", "GLY"],
    )
    (structures_dir / f"{family_id}.pdb").write_text(text, encoding="utf-8")


def test_run_pocket_search_end_to_end(tmp_path, monkeypatch):
    dbdir = tmp_path / "dbdir"
    _setup_dbdir(dbdir)
    structures_dir = tmp_path / "fold" / "structures"
    _write_structure(structures_dir, "FAM1")
    out_dir = tmp_path / "search"
    out_dir.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(
        pocket.externals, "p2rank",
        lambda: ToolInfo(name="p2rank", available=True, path="p2rank"),
    )
    monkeypatch.setattr(
        pocket.externals, "foldseek",
        lambda: ToolInfo(name="foldseek", available=True, path="foldseek"),
    )

    def fake_run(cmd, capture_output=True, text=True):
        if cmd[1] == "predict":
            out = Path(cmd[cmd.index("-o") + 1])
            out.mkdir(parents=True, exist_ok=True)
            _predictions_csv(
                out / "FAM1_predictions.csv",
                [("FAM1", 1, 3.0, 0.9, "A_2 A_5")],
            )
        elif cmd[1] == "convertalis":
            raw_tsv = Path(cmd[5])
            raw_tsv.parent.mkdir(parents=True, exist_ok=True)
            # 11 columns; query/target 10 residues aligned without gaps.
            row = "\t".join(
                ["FAM1", "1abc_A", "1", "1", "AAAAAAAAAA", "AAAAAAAAAA",
                 "0.8", "99.0", "1e-10", "0.9", "0.8"]
            )
            raw_tsv.write_text(row + "\n", encoding="utf-8")

        class Completed:
            returncode = 0
        return Completed()

    monkeypatch.setattr(pocket.subprocess, "run", fake_run)

    evidence, raw = pocket.run_pocket_search(
        [_family()], structures_dir, dbdir, out_dir, min_plddt=70.0,
    )
    assert raw is not None and raw.exists()
    rep_key = "FAM1@g1"
    assert rep_key in evidence
    hit = evidence[rep_key][0]
    assert hit.channel == EvidenceChannel.POCKET
    assert hit.target_id == "EC:3.5.2.6"
    assert hit.target_label == "beta-lactamase"
    assert hit.details["db"] == "ec"
    assert hit.details["catalytic_in_pocket"] == 2
    assert hit.details["pocket_consistency"] == pytest.approx(1.0)


def test_run_pocket_search_missing_p2rank(tmp_path, monkeypatch):
    monkeypatch.setattr(
        pocket.externals, "p2rank",
        lambda: ToolInfo(name="p2rank", available=False),
    )
    monkeypatch.setattr(
        pocket.externals, "foldseek",
        lambda: ToolInfo(name="foldseek", available=True, path="foldseek"),
    )
    evidence, raw = pocket.run_pocket_search(
        [], Path("x"), None, tmp_path, min_plddt=70.0,
    )
    assert evidence == {}
    assert raw is None


def test_run_pocket_search_missing_database(tmp_path, monkeypatch):
    monkeypatch.setattr(
        pocket.externals, "p2rank",
        lambda: ToolInfo(name="p2rank", available=True, path="p2rank"),
    )
    monkeypatch.setattr(
        pocket.externals, "foldseek",
        lambda: ToolInfo(name="foldseek", available=True, path="foldseek"),
    )
    dbdir = tmp_path / "dbdir"
    dbdir.mkdir()
    evidence, raw = pocket.run_pocket_search(
        [_family()], tmp_path / "structures", dbdir, tmp_path / "out",
        min_plddt=70.0,
    )
    assert evidence == {}
    assert raw is None


def test_run_pocket_search_external_failure_skips(tmp_path, monkeypatch):
    dbdir = tmp_path / "dbdir"
    _setup_dbdir(dbdir)
    structures_dir = tmp_path / "structures"
    _write_structure(structures_dir, "FAM1")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    monkeypatch.setattr(
        pocket.externals, "p2rank",
        lambda: ToolInfo(name="p2rank", available=True, path="p2rank"),
    )
    monkeypatch.setattr(
        pocket.externals, "foldseek",
        lambda: ToolInfo(name="foldseek", available=True, path="foldseek"),
    )

    def failing_run(cmd, capture_output=True, text=True):
        import subprocess as sp
        raise sp.CalledProcessError(1, cmd, "", "boom")

    monkeypatch.setattr(pocket.subprocess, "run", failing_run)

    evidence, raw = pocket.run_pocket_search(
        [_family()], structures_dir, dbdir, out_dir, min_plddt=70.0,
    )
    assert evidence == {}
    assert raw is None
