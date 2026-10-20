"""Ontology Dbxref parsing."""

from pufprobe.ontology import parse_dbxrefs


def test_parse_dbxref_tokens() -> None:
    out = parse_dbxrefs(
        [
            "GO:0003735",
            "Pfam:PF00042",
            "EC:3.4.24.77",
            "KEGG:K01977",
            "RefSeq:WP_000001.1",
        ]
    )
    assert out["Pfam"] == ["PF00042"]
    assert out["EC"] == ["EC:3.4.24.77"]
    assert out["GO"] == ["GO:0003735"]


def test_parse_bare_and_incomplete_identifiers() -> None:
    out = parse_dbxrefs(
        [
            "PF00959",
            "EC:3.4.24.-",
            "GO:0005840",
            "COG:COG0001",
        ]
    )
    assert out["Pfam"] == ["PF00959"]
    assert out["EC"] == ["EC:3.4.24.-"]
    assert out["GO"] == ["GO:0005840"]


def test_empty_and_unrelated() -> None:
    out = parse_dbxrefs(["", "SO:0000755", "garbage token"])
    assert out == {"Pfam": [], "EC": [], "GO": []}
