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
    # Default config; fusion uses the group noisy-OR path with the default
    # channel reliability priors.
    return deepcopy(load_config(None))


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
def _fuse(evidence: dict, cfg: dict = None, pfam_names: dict = None):
    """Run fuse_proteins on a small set of records and return {key: terms}."""
    from pufprobe.rank.fuse import fuse_proteins

    keys = sorted({k for k in evidence})
    records = []
    for key in keys:
        sample, locus = key.split("@", 1)
        records.append(_cds(sample, locus))
    from pufprobe.models import Evidence as _E, EvidenceChannel as _EC

    ev_map: dict = {}
    for key, rows in evidence.items():
        ev_map[key] = [
            _E(
                channel=_EC(r["channel"]),
                target_id=r["target_id"],
                target_label=r["target_label"],
                score=float(r["score"]),
                evalue=r.get("evalue"),
                coverage=r.get("coverage"),
                details=r.get("details", {}),
            )
            for r in rows
        ]
    preds = fuse_proteins(
        records,
        ev_map,
        cfg if cfg is not None else _custom_config(),
        pfam_names=pfam_names,
    )
    return {p.key: p.terms for p in preds}


def test_single_group_capped_medium_confidence():
    # A single strong database hit is capped at single_group_cap (0.60).
    terms = _fuse({"sA@g1": [_ev("database", "UniRef90:A", "single group", 1.0)]})
    term = terms["sA@g1"][0]
    assert term.raw_score == pytest.approx(0.60)
    assert term.confidence == pytest.approx(0.60)
    assert term.evidence_groups == ["database"]
    assert "capped" in term.rationale


def test_two_direct_groups_noisy_or():
    # Two direct groups (database + remote_homology), same target.
    terms = _fuse(
        {
            "sA@g2": [
                _ev("database", "UniRef90:B", "two-group strong", 1.0),
                _ev("remote_homology", "UniRef90:B", "same target from homology", 1.0),
            ],
            "sA@g3": [
                _ev("database", "UniRef90:C", "two-group weak", 0.5),
                _ev("remote_homology", "UniRef90:C", "same target from homology", 0.5),
            ],
        }
    )
    # strong: noisy-OR of p=0.90 (database) and p=0.95 (remote_homology)
    # = 1 - 0.10*0.05 = 0.995
    t2 = terms["sA@g2"][0]
    assert t2.raw_score == pytest.approx(0.995, abs=1e-6)
    assert t2.raw_score >= 0.80
    assert "capped" not in t2.rationale
    # weak: p=0.45 and p=0.475 -> 1 - 0.55*0.525 = 0.71125
    t3 = terms["sA@g3"][0]
    assert t3.raw_score == pytest.approx(0.71125, abs=1e-6)
    assert t3.raw_score < 0.80


def test_context_does_not_unlock_high_confidence():
    # database + context: context is indirect and does not count as a group.
    terms = _fuse(
        {
            "sA@g2": [
                _ev("database", "UniRef90:B", "strong", 1.0),
                _ev("context", "UniRef90:B", "same target from context", 1.0),
            ],
        }
    )
    t2 = terms["sA@g2"][0]
    # noisy-OR would be 0.96 but only one direct group -> capped 0.60
    assert t2.raw_score == pytest.approx(0.60)
    assert "capped" in t2.rationale
    assert "indirect" in t2.rationale


def test_structure_plus_sequence_ontology_reaches_high_confidence():
    # Foldseek (structure) and remote homology (sequence_profile) both map to
    # the same Pfam family -> two orthogonal groups -> high confidence.
    terms = _fuse(
        {
            "sA@g1": [
                _ev("foldseek", "PF00589", "phage integrase core", 1.0),
                _ev("remote_homology", "PF00589", "phage integrase core", 1.0),
            ],
        }
    )
    term = terms["sA@g1"][0]
    assert term.term_type == "Pfam"
    assert term.term_id == "PF00589"
    assert sorted(term.evidence_groups) == ["sequence_profile", "structure"]
    # p=0.85 and p=0.95 -> 1 - 0.15*0.05 = 0.9925
    assert term.raw_score == pytest.approx(0.9925, abs=1e-6)
    assert term.raw_score >= 0.80


def test_strong_sequence_weak_structure_is_capped():
    # Strong sequence profile (group strength 0.95) plus a near-noise structure
    # hit (0.85*0.15=0.1275, below the 0.40 group gate). The weak structure
    # group still appears in the evidence, but it cannot be a second
    # independent witness, so the term is capped at medium confidence.
    terms = _fuse(
        {
            "sA@g1": [
                _ev("remote_homology", "PF00001", "strong sequence", 1.0),
                _ev("foldseek", "PF00001", "weak structure", 0.15),
            ],
        }
    )
    term = terms["sA@g1"][0]
    assert sorted(term.evidence_groups) == ["sequence_profile", "structure"]
    assert term.confidence == pytest.approx(0.60)
    assert "capped" in term.rationale
    assert "structure" in term.rationale


def test_strong_structure_weak_sequence_is_capped():
    # Mirror case: strong structure (0.85) but a weak sequence profile
    # (0.95*0.3=0.285, below the gate) is only one strong group -> capped.
    terms = _fuse(
        {
            "sA@g1": [
                _ev("foldseek", "PF00001", "strong structure", 1.0),
                _ev("remote_homology", "PF00001", "weak sequence", 0.3),
            ],
        }
    )
    term = terms["sA@g1"][0]
    assert sorted(term.evidence_groups) == ["sequence_profile", "structure"]
    assert term.confidence == pytest.approx(0.60)
    assert "capped" in term.rationale


def test_two_groups_above_gate_reach_high_confidence():
    # Sequence 0.95 and structure 0.85*0.7=0.595 both clear the 0.40 gate:
    # 1 - 0.05*0.405 = 0.97975, high confidence.
    terms = _fuse(
        {
            "sA@g1": [
                _ev("remote_homology", "PF00001", "strong sequence", 1.0),
                _ev("foldseek", "PF00001", "solid structure", 0.7),
            ],
        }
    )
    term = terms["sA@g1"][0]
    assert sorted(term.evidence_groups) == ["sequence_profile", "structure"]
    assert term.confidence == pytest.approx(0.97975, abs=1e-4)
    assert "capped" not in term.rationale


def test_pdb_template_label_replaced_for_pfam_term():
    # foldseek maps the PDB-chain template 6jbu_B to PF13855 and scores higher
    # than the sequence hit; the term must display the functional name rather
    # than the bare PDB id.
    terms = _fuse(
        {
            "sA@g1": [
                _ev("foldseek", "PF13855", "6jbu_B", 0.8),
                _ev("remote_homology", "PF13855", "Leucine-rich repeat protein", 0.5),
            ],
        }
    )
    term = terms["sA@g1"][0]
    assert term.term_id == "PF13855"
    assert term.term_name == "Leucine-rich repeat protein"


def test_pfam_with_only_pdb_label_falls_back_to_accession():
    # A foldseek-only Pfam hit with no functional label anywhere shows the
    # Pfam accession, which is at least the correct family identifier.
    terms = _fuse(
        {
            "sA@g1": [
                _ev("foldseek", "PF13855", "6jbu_B", 0.9),
            ],
        }
    )
    term = terms["sA@g1"][0]
    assert term.term_type == "Pfam"
    assert term.term_name == "PF13855"


def test_pfam_official_description_overrides_homolog_and_pdb():
    # With the Pfam name index the official family description names the
    # term, regardless of the PDB template or a distant eukaryotic homolog.
    pfam_names = {
        "PF13855": {"id": "LRR_8", "name": "Leucine rich repeat", "type": "Repeat"},
    }
    terms = _fuse(
        {
            "sA@g1": [
                _ev("foldseek", "PF13855", "6jbu_B", 0.9),
                _ev("remote_homology", "PF13855", "Protein lap4 OS=Drosophila", 0.6),
            ],
        },
        pfam_names=pfam_names,
    )
    term = terms["sA@g1"][0]
    assert term.term_id == "PF13855"
    assert term.term_name == "Leucine rich repeat"


def test_pfam_description_falls_back_to_family_id():
    # A record without a description still resolves to its short identifier.
    pfam_names = {"PF00005": {"id": "ABC_tran", "type": "Domain"}}
    terms = _fuse(
        {
            "sA@g1": [
                _ev("foldseek", "PF00005", "9iqf_B", 0.9),
            ],
        },
        pfam_names=pfam_names,
    )
    assert terms["sA@g1"][0].term_name == "ABC_tran"


def test_unmapped_structure_term_keeps_pdb_label():
    # A foldseek hit without a SIFTS Pfam mapping is a structure term; the PDB
    # template is the result itself and its label must be preserved.
    terms = _fuse(
        {
            "sA@g1": [
                _ev("foldseek", "4fqe", "4fqe", 0.9),
            ],
        }
    )
    term = terms["sA@g1"][0]
    assert term.term_type == "PDB"
    assert term.term_name == "4fqe"


def test_within_group_channels_not_double_counted():
    # foldseek and pocket both belong to the structure group; they are not
    # independent and must not be combined as a second group.
    terms = _fuse(
        {
            "sA@g1": [
                _ev("foldseek", "PF00067", "some fold", 1.0),
                _ev("pocket", "PF00067", "pocket in same fold", 1.0),
            ],
        }
    )
    term = terms["sA@g1"][0]
    assert term.evidence_groups == ["structure"]
    # group strength is max(0.85, 0.60) = 0.85; single group -> capped 0.60
    assert term.raw_score == pytest.approx(0.60)
    assert "capped" in term.rationale


def test_weighted_sum_legacy_method_still_available():
    cfg = _custom_config()
    cfg["rank"]["fusion_method"] = "weighted_sum"
    terms = _fuse(
        {"sA@g1": [_ev("database", "UniRef90:A", "single channel", 1.0)]},
        cfg=cfg,
    )
    term = terms["sA@g1"][0]
    # raw weighted sum = 0.12; single group cap to high_conf - 0.01 = 0.79
    # leaves 0.12 unchanged.
    assert term.raw_score == pytest.approx(0.12)


def test_fusion_via_pipeline(rank_workdir):
    from pufprobe.rank import pipeline as rank_pipeline

    search = {
        "sA@g1": [_ev("database", "UniRef90:A", "single group", 1.0)],
        "sA@g2": [
            _ev("database", "UniRef90:B", "two direct groups", 1.0),
            _ev("remote_homology", "UniRef90:B", "same target from homology", 1.0),
        ],
    }
    (rank_workdir.path("search", "search_evidence.json")).write_text(
        json.dumps(search), encoding="utf-8"
    )

    rank_pipeline.run(rank_workdir, _custom_config())
    preds = {f"{p['sample_id']}@{p['locus_tag']}": p for p in json.loads(rank_workdir.path("rank", "predictions.json").read_text("utf-8"))}

    assert preds["sA@g1"]["terms"][0]["raw_score"] == pytest.approx(0.60)
    assert preds["sA@g2"]["terms"][0]["raw_score"] == pytest.approx(0.995, abs=1e-6)


def test_calibration_metadata_identity(rank_workdir):
    from pufprobe.rank import pipeline as rank_pipeline

    rank_pipeline.run(rank_workdir, _custom_config())
    cal = json.loads(rank_workdir.path("rank", "calibration.json").read_text("utf-8"))
    assert cal["method"] == "identity"
    assert cal["ece"] is None
    assert cal["note"] == "no held-out calibration set provided"
    assert "no learned model" in cal["fusion"]
    assert cal["fusion_method"] == "group_noisy_or"
    assert cal["single_group_cap"] == 0.60
    assert cal["go_true_path_propagation"] == "disabled"
    assert cal["n_go_propagated"] == 0


_GO_OBO = """\
format-version: 1.2

[Term]
id: GO:0000001
name: root one
namespace: biological_process

[Term]
id: GO:0000002
name: child two
namespace: biological_process
is_a: GO:0000001 ! root one
"""


def test_go_true_path_propagation_via_pipeline(rank_workdir, tmp_path: Path):
    from pufprobe.rank import pipeline as rank_pipeline

    dbdir = tmp_path / "dbdir"
    dbdir.mkdir()
    (dbdir / "go-basic.obo").write_text(_GO_OBO, encoding="utf-8")

    search = {
        "sA@g1": [
            {
                "channel": "eggnog",
                "target_id": "GO:0000002",
                "target_label": "child two",
                "score": 1.0,
                "evalue": None,
                "coverage": None,
                "details": {"db": "go"},
            }
        ],
    }
    (rank_workdir.path("search", "search_evidence.json")).write_text(
        json.dumps(search), encoding="utf-8"
    )

    rank_pipeline.run(rank_workdir, _custom_config(), dbdir=dbdir)
    preds = {f"{p['sample_id']}@{p['locus_tag']}": p
             for p in json.loads(rank_workdir.path("rank", "predictions.json").read_text("utf-8"))}
    term_ids = {t["term_id"] for t in preds["sA@g1"]["terms"]}
    # Direct child term and its propagated ancestor are both present.
    assert {"GO:0000002", "GO:0000001"} <= term_ids

    cal = json.loads(rank_workdir.path("rank", "calibration.json").read_text("utf-8"))
    assert cal["go_true_path_propagation"] == "applied"
    assert cal["n_go_propagated"] == 1


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
