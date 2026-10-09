"""Tests for the rank stage: rule fusion, channel-diversity cap, calibration."""

import json
from copy import deepcopy
from pathlib import Path

import pytest

from pufprobe.config import DEFAULTS, load_config
from pufprobe.models import CDSRecord, Family, GeneStatus, Sample, to_dict
from pufprobe.rank.calibrate import (
    IdentityCalibrator,
    IsotonicCalibrator,
    TemperatureScaling,
    build_calibrator,
)
from pufprobe.rank.fuse import classify_term
from pufprobe.models import Evidence, EvidenceChannel
from pufprobe.workdir import WorkDir


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _cds(sample: str, locus: str) -> CDSRecord:
    return CDSRecord(
        sample_id=sample,
        locus_tag=locus,
        protein_id=None,
        sequence_aa="M" * 120,
        start=1,
        end=360,
        strand="+",
        product="hypothetical protein",
        status=GeneStatus.HYPOTHETICAL,
    )


def _ev(channel: str, target_id: str, label: str, score: float) -> dict:
    return {
        "channel": channel,
        "target_id": target_id,
        "target_label": label,
        "score": score,
        "evalue": None,
        "coverage": None,
        "details": {},
    }


@pytest.fixture()
def rank_workdir(tmp_path: Path) -> WorkDir:
    wd = WorkDir(tmp_path / "run")
    records = [_cds("sA", "g1"), _cds("sA", "g2"), _cds("sA", "g3"), _cds("sA", "g4"), _cds("sB", "g5")]
    (wd.path("prepare", "samples.json")).write_text(
        json.dumps([Sample(sample_id="sA").__dict__, Sample(sample_id="sB").__dict__]),
        encoding="utf-8",
    )
    (wd.path("prepare", "cds_records.json")).write_text(
        json.dumps([to_dict(r) for r in records]), encoding="utf-8"
    )
    family = Family(
        family_id="fam_0",
        representative_key="sA@g4",
        members=["sA@g4"],
        representative_seq="M" * 40,
        size=1,
        lineage_width=1,
        metrics={"sample_count": 1},
    )
    (wd.path("cluster", "families.json")).write_text(
        json.dumps([to_dict(family)]), encoding="utf-8"
    )
    return wd


def _custom_config() -> dict:
    cfg = deepcopy(load_config(None))
    cfg["rank"]["rule_weights"] = {
        "plm_neighbor": 0.0,
        "foldseek": 0.0,
        "remote_homology": 0.0,
        "pocket": 0.0,
        "context": 0.1,
        "database": 0.9,
    }
    return cfg


# ---------------------------------------------------------------------------
# term classification
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "target_id,label,expected",
    [
        ("GO:0003674", "molecular function", ("GO", "GO:0003674")),
        ("EC:1.1.1.1", "alcohol dehydrogenase", ("EC", "EC:1.1.1.1")),
        ("1.1.1.1", "ec number", ("EC", "EC:1.1.1.1")),
        ("PF00001.12", "7tm receptor", ("Pfam", "PF00001.12")),
        ("COG0513", "his kinase", ("COG", "COG0513")),
        ("PDB:1ABC", "crystal structure", ("label", "PDB:1ABC")),
        ("uniref:X9", "some uniprot", ("label", "uniref:X9")),
    ],
)
def test_classify_term(target_id, label, expected):
    ev = Evidence(EvidenceChannel.DATABASE, target_id, label, 0.5)
    assert classify_term(ev) == expected


def test_classify_foldseek_is_structure():
    ev = Evidence(EvidenceChannel.FOLDSEEK, "afdb:Q9XYZ", "AlphaFold hit", 0.7)
    assert classify_term(ev) == ("structure", "afdb:Q9XYZ")


# ---------------------------------------------------------------------------
# fusion rules
# ---------------------------------------------------------------------------
def test_single_channel_capped_and_two_groups_not(rank_workdir):
    from pufprobe.rank import pipeline as rank_pipeline

    search = {
        "sA@g1": [_ev("database", "UniRef90:A", "single-channel term", 1.0)],
        "sA@g2": [
            _ev("database", "UniRef90:B", "two-group term", 1.0),
            _ev("context", "UniRef90:B", "same target from context", 1.0),
        ],
        "sA@g3": [
            _ev("database", "UniRef90:C", "half score", 0.5),
            _ev("context", "UniRef90:C", "same target from context", 0.5),
        ],
    }
    (rank_workdir.path("search", "search_evidence.json")).write_text(
        json.dumps(search), encoding="utf-8"
    )

    rank_pipeline.run(rank_workdir, _custom_config())
    preds = {f"{p['sample_id']}@{p['locus_tag']}": p for p in json.loads(rank_workdir.path("rank", "predictions.json").read_text("utf-8"))}

    g1_term = preds["sA@g1"]["terms"][0]
    # raw would be 0.90 with one evidence group -> capped to high_conf-0.01
    assert g1_term["raw_score"] == pytest.approx(0.79)
    assert "capped" in g1_term["rationale"]

    g2_term = preds["sA@g2"]["terms"][0]
    # two independent groups, no cap
    assert g2_term["raw_score"] == pytest.approx(1.0)
    assert "capped" not in g2_term["rationale"]

    g3_term = preds["sA@g3"]["terms"][0]
    assert g3_term["raw_score"] == pytest.approx(0.9 * 0.5 + 0.1 * 0.5)


def test_calibration_metadata_identity(rank_workdir):
    from pufprobe.rank import pipeline as rank_pipeline

    rank_pipeline.run(rank_workdir, _custom_config())
    cal = json.loads(rank_workdir.path("rank", "calibration.json").read_text("utf-8"))
    assert cal["method"] == "identity"
    assert cal["ece"] is None
    assert cal["note"] == "no held-out calibration set provided"
    assert "rule weights (no learned model)" in cal["fusion"]


def test_priority_scoring(rank_workdir):
    from pufprobe.rank import pipeline as rank_pipeline

    # g4 (family representative) and g5 carry no evidence -> family unanchored.
    rank_pipeline.run(rank_workdir, _custom_config())
    rows = rank_workdir.path("rank", "priorities.tsv").read_text("utf-8").splitlines()
    assert rows[0] == "family_id\tpriority_score\treasons"
    fid, score, reasons = rows[1].split("\t")
    assert fid == "fam_0"
    # novelty 0.5*0.30 + clonal (1 - 1/2 samples)=0.5*0.20 + small 1.0*0.10
    # = 0.15 + 0.10 + 0.10 = 0.35
    assert float(score) == pytest.approx(0.35)
    assert "small protein" in reasons
    assert "structural_novelty" in reasons


# ---------------------------------------------------------------------------
# calibrators
# ---------------------------------------------------------------------------
def test_temperature_scaling_identity_at_t1():
    cal = TemperatureScaling(temperature=1.0)
    assert cal.transform([0.1, 0.9]) == pytest.approx([0.1, 0.9])


def test_isotonic_pava_monotone():
    iso = IsotonicCalibrator().fit([0.1, 0.4, 0.5, 0.9], [0.0, 1.0, 0.0, 1.0])
    out = iso.transform([0.0, 0.5, 1.0])
    assert out == sorted(out)  # monotone non-decreasing
    assert all(0.0 <= v <= 1.0 for v in out)


def test_build_calibrator_identity():
    cal, meta = build_calibrator(None)
    assert isinstance(cal, IdentityCalibrator)
    assert meta["method"] == "identity"
