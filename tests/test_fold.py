"""fold stage: precomputed PDB pLDDT parsing, quality gate, ModelUnavailable path."""

import json
from pathlib import Path

import pytest

from pufprobe.config import load_config
from pufprobe.fold.esmfold import ModelUnavailable, fold_sequence, parse_plddt_from_pdb
from pufprobe.fold.pipeline import run as fold_run
from pufprobe.io_utils import read_json, write_json
from pufprobe.models import Family
from pufprobe.workdir import WorkDir


def _pdb_with_bfactors(values) -> str:
    lines = []
    for i, bf in enumerate(values, start=1):
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
            f"   1.000   2.000   3.000  1.00{bf:6.2f}\n"
        )
    return "".join(lines)


def test_parse_plddt_auto_detects_scale() -> None:
    # 0-100 scale
    pdb100 = _pdb_with_bfactors([80.0, 90.0, 95.0, 60.0])
    assert parse_plddt_from_pdb(pdb100) == pytest.approx(85.0)
    # 0-1 scale -> auto-scaled to 0-100
    pdb01 = _pdb_with_bfactors([0.80, 0.90, 0.95, 0.60])
    assert parse_plddt_from_pdb(pdb01) == pytest.approx(85.0)


def _write_families(wd: WorkDir) -> None:
    families = [
        Family(family_id="FAM00001", representative_key="sA@g1",
               members=["sA@g1"], representative_seq="M" * 50,
               size=1, lineage_width=1),
        Family(family_id="FAM00002", representative_key="sB@g2",
               members=["sB@g2"], representative_seq="M" * 50,
               size=1, lineage_width=1),
    ]
    write_json(wd.path("cluster", "families.json"),
               [json.loads(json.dumps(f.__dict__)) for f in families])


def test_fold_precomputed_quality_gate(tmp_path: Path) -> None:
    wd = WorkDir(tmp_path / "run")
    config = load_config(None)
    _write_families(wd)

    pre = tmp_path / "pre"
    pre.mkdir()
    # high confidence family
    (pre / "FAM00001.pdb").write_text(_pdb_with_bfactors([80, 90, 95, 88]),
                                      encoding="utf-8")
    # low confidence family
    (pre / "FAM00002.pdb").write_text(_pdb_with_bfactors([30, 40, 45, 35]),
                                      encoding="utf-8")

    fold_run(wd, config, precomputed=pre)

    q = read_json(wd.path("cluster", "families.json"))  # type: ignore[arg-type]
    by_id = {f["family_id"]: f for f in q}
    assert by_id["FAM00001"]["fold_status"] == "folded"
    assert by_id["FAM00002"]["fold_status"] == "low_confidence"
    assert by_id["FAM00001"]["plddt"] == pytest.approx(89.0)
    assert by_id["FAM00001"]["structure_path"].endswith("FAM00001.pdb")

    tsv = (wd.path("fold", "fold_quality.tsv")).read_text(encoding="utf-8")
    assert "FAM00001\tprecomputed\tfolded" in tsv
    assert "FAM00002\tprecomputed\tlow_confidence" in tsv
    assert wd.marker("fold").is_file()


def test_af2_engine_marks_skipped(tmp_path: Path) -> None:
    wd = WorkDir(tmp_path / "run")
    _write_families(wd)
    fold_run(wd, load_config(None), engine="af2")
    q = read_json(wd.path("cluster", "families.json"))  # type: ignore[arg-type]
    assert all(f["fold_status"] == "skipped" for f in q)


def test_esmfold_unavailable_without_torch() -> None:
    # This host has no torch; fold_sequence must raise ModelUnavailable.
    with pytest.raises(ModelUnavailable):
        fold_sequence("M" * 100)
