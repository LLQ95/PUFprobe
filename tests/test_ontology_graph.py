"""Tests for the GO ontology graph and true-path propagation."""

from pathlib import Path

from pufprobe.models import Evidence, EvidenceChannel
from pufprobe.ontology_graph import load_go_obo, parse_go_obo

# A minimal OBO fragment covering: a three-level is_a chain, a part_of edge,
# a regulates edge that must NOT propagate, an alt_id and an obsolete term.
_MINI_OBO = """\
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

[Term]
id: GO:0000003
name: grandchild three
namespace: biological_process
is_a: GO:0000002 ! child two
relationship: part_of GO:0000010 ! part parent

[Term]
id: GO:0000010
name: part parent
namespace: biological_process

[Term]
id: GO:0000004
name: regulator of grandchild
namespace: biological_process
relationship: regulates GO:0000003 ! grandchild three

[Term]
id: GO:0000005
name: alt holder
namespace: molecular_function
alt_id: GO:0000099

[Term]
id: GO:0000006
name: old term
namespace: biological_process
is_obsolete: true

[Typedef]
id: part_of
name: part of
"""


def _ontology():
    return parse_go_obo(_MINI_OBO.splitlines(keepends=True))


def _go(channel: str, target: str, score: float) -> Evidence:
    return Evidence(
        channel=EvidenceChannel(channel),
        target_id=target,
        target_label=target,
        score=score,
        details={"db": "go"},
    )


# ---------------------------------------------------------------------------
# parsing and graph structure
# ---------------------------------------------------------------------------
def test_parse_counts_and_names():
    onto = _ontology()
    # Six live terms; GO:0000006 is obsolete and excluded.
    assert onto.n_terms == 6
    assert onto.name("GO:0000003") == "grandchild three"
    assert "GO:0000006" in onto.obsolete


def test_isa_and_part_of_closure():
    onto = _ontology()
    closure = set(onto.ancestors("GO:0000003"))
    # is_a chain 3 -> 2 -> 1 and the part_of edge 3 -> 10.
    assert closure == {"GO:0000003", "GO:0000002", "GO:0000001", "GO:0000010"}


def test_regulates_is_not_propagated():
    onto = _ontology()
    # GO:0000004 regulates GO:0000003; this is not a true-path edge, so the
    # regulator has no parents and the regulated term is not an ancestor.
    assert set(onto.ancestors("GO:0000004")) == {"GO:0000004"}
    assert "GO:0000004" not in set(onto.ancestors("GO:0000003"))


def test_alt_id_resolves():
    onto = _ontology()
    assert onto.resolve("GO:0000099") == "GO:0000005"
    assert set(onto.ancestors("GO:0000099")) == {"GO:0000005"}


def test_obsolete_and_unknown_terms():
    onto = _ontology()
    assert onto.ancestors("GO:0000006") == tuple()
    # Unknown GO id is treated as an isolated leaf.
    assert set(onto.ancestors("GO:9999999")) == {"GO:9999999"}


# ---------------------------------------------------------------------------
# evidence propagation
# ---------------------------------------------------------------------------
def test_propagation_adds_only_ancestors():
    onto = _ontology()
    ev_map = {"sA@g1": [_go("eggnog", "GO:0000003", 1.0)]}
    out, added = onto.propagate_evidence(ev_map)
    assert added == 3
    targets = {e.target_id for e in out["sA@g1"]}
    assert targets == {"GO:0000003", "GO:0000002", "GO:0000001", "GO:0000010"}
    for e in out["sA@g1"]:
        if e.target_id != "GO:0000003":
            assert e.details["propagated"] is True
            assert e.details["propagated_from"] == "GO:0000003"
            assert e.channel == EvidenceChannel.EGGNOG
            assert e.target_label == onto.name(e.target_id)


def test_propagated_evidence_keeps_channel_and_score():
    onto = _ontology()
    source = _go("eggnog", "GO:0000003", 0.7)
    out, _ = onto.propagate_evidence({"sA@g1": [source]})
    parent = next(e for e in out["sA@g1"] if e.target_id == "GO:0000001")
    assert parent.score == 0.7
    assert parent.channel == EvidenceChannel.EGGNOG


def test_propagation_deduplicates_same_channel_keeps_strongest():
    onto = _ontology()
    # One strong deep hit and one weak shallow hit in the same channel; the
    # shared root must receive only one propagated record, the strongest.
    ev_map = {
        "sA@g1": [
            _go("eggnog", "GO:0000003", 1.0),
            _go("eggnog", "GO:0000002", 0.5),
        ]
    }
    out, added = onto.propagate_evidence(ev_map)
    root = [e for e in out["sA@g1"] if e.target_id == "GO:0000001"]
    assert len(root) == 1
    assert root[0].score == 1.0
    assert root[0].details["propagated_from"] == "GO:0000003"


def test_direct_annotation_preferred_when_at_least_as_strong():
    onto = _ontology()
    # Direct 0.9 to the parent beats a weak 0.7 propagated from the child.
    ev_map = {
        "sA@g1": [
            _go("eggnog", "GO:0000003", 0.7),
            _go("eggnog", "GO:0000002", 0.9),
        ]
    }
    out, added = onto.propagate_evidence(ev_map)
    middle = [e for e in out["sA@g1"] if e.target_id == "GO:0000002"]
    assert len(middle) == 1
    assert middle[0].score == 0.9
    assert not middle[0].details.get("propagated", False)


def test_propagation_lets_channels_converge_on_ancestor():
    onto = _ontology()
    # eggNOG supports the deep child; a profile hit directly supports the
    # middle term. After propagation the shared ancestor GO:0000002 carries
    # two independent channels (database + sequence_profile).
    ev_map = {
        "sA@g1": [
            _go("eggnog", "GO:0000003", 1.0),
            _go("remote_homology", "GO:0000002", 0.9),
        ]
    }
    out, _ = onto.propagate_evidence(ev_map)
    middle = [e for e in out["sA@g1"] if e.target_id == "GO:0000002"]
    channels = {e.channel for e in middle}
    assert channels == {EvidenceChannel.EGGNOG, EvidenceChannel.REMOTE_HOMOLOGY}


def test_non_go_evidence_is_untouched():
    onto = _ontology()
    pfam = Evidence(
        channel=EvidenceChannel.FOLDSEEK,
        target_id="PF00001",
        target_label="7tm1",
        score=0.8,
        details={"db": "pfam"},
    )
    out, added = onto.propagate_evidence({"sA@g1": [pfam]})
    assert added == 0
    assert out["sA@g1"] == [pfam]


def test_alt_id_evidence_propagates_through_primary():
    onto = _ontology()
    # Hitting the alt id should resolve to the primary term (no crash; it is a
    # leaf here so nothing is added, but the id is recognised).
    ev_map = {"sA@g1": [_go("eggnog", "GO:0000099", 1.0)]}
    out, added = onto.propagate_evidence(ev_map)
    assert added == 0
    assert any(e.target_id == "GO:0000099" for e in out["sA@g1"])


# ---------------------------------------------------------------------------
# loader
# ---------------------------------------------------------------------------
def test_load_go_obo_missing_returns_none(tmp_path: Path):
    assert load_go_obo(tmp_path / "absent.obo") is None
    assert load_go_obo(None) is None


def test_load_go_obo_from_file(tmp_path: Path):
    obo = tmp_path / "go-basic.obo"
    obo.write_text(_MINI_OBO, encoding="utf-8")
    onto = load_go_obo(obo)
    assert onto is not None
    assert onto.n_terms == 6
    assert set(onto.ancestors("GO:0000003")) == {
        "GO:0000003", "GO:0000002", "GO:0000001", "GO:0000010",
    }
