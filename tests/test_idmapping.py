"""search/idmapping: Swiss-Prot cross-reference parsing and evidence expansion."""

import io

from pufprobe.models import Evidence, EvidenceChannel
from pufprobe.search import idmapping

DAT = """ID   DBHA_ECOLI             STANDARD;      PRT; 90 AA.
AC   P0ACF0; Q9LEU2;
DT   01-JAN-1990, integrated into UniProtKB/Swiss-Prot.
DE   RecName: Full=DNA-binding protein HU-alpha;
DE            EC=3.4.24.89;
GN   Name=hupA;
DR   Pfam; PF00216; Bacterial DNA-binding protein; 1.
DR   GO; GO:0003677; F:DNA binding; IEA:UniProtKB.
//
ID   NOX_HUMAN              STANDARD;      PRT; 200 AA.
AC   Q0Q012;
DT   01-JAN-1990, integrated into UniProtKB/Swiss-Prot.
DE   RecName: Full=Some protein;
DR   Pfam; PF00001; Some Domain; 1.
//
"""


def test_parse_sprot_dat() -> None:
    xref = idmapping.parse_sprot_dat(io.StringIO(DAT))
    # Primary and secondary accession map to the same entry.
    assert xref["P0ACF0"] == xref["Q9LEU2"]
    assert xref["P0ACF0"]["Pfam"] == ["PF00216"]
    assert xref["P0ACF0"]["EC"] == ["3.4.24.89"]
    assert xref["P0ACF0"]["GO"] == ["GO:0003677"]
    # Entry without EC/GO still keeps Pfam and omits the empty keys.
    assert xref["Q0Q012"]["Pfam"] == ["PF00001"]
    assert "EC" not in xref["Q0Q012"]
    assert "GO" not in xref["Q0Q012"]


def test_accessions_from_target() -> None:
    assert idmapping.accessions_from_target("sp|P0ACF0|DBHA_ECOLI") == ["P0ACF0"]
    assert idmapping.accessions_from_target("P0ACF0") == ["P0ACF0"]
    assert idmapping.accessions_from_target("  ") == []


def test_expand_evidence() -> None:
    xref = idmapping.parse_sprot_dat(io.StringIO(DAT))
    hit = Evidence(
        channel=EvidenceChannel.REMOTE_HOMOLOGY,
        target_id="sp|P0ACF0|DBHA_ECOLI",
        target_label="DNA-binding protein HU-alpha",
        score=0.8,
        evalue=1e-20,
        coverage=0.9,
        details={},
    )
    expanded, added = idmapping.expand_evidence({"k": [hit]}, xref)
    ids = sorted(e.target_id for e in expanded["k"])
    # Original label plus the three ontology terms.
    assert added == 3
    assert "sp|P0ACF0|DBHA_ECOLI" in ids
    assert "PF00216" in ids
    assert "EC:3.4.24.89" in ids
    assert "GO:0003677" in ids
    # Channel stays remote_homology; the ontology hint is in details.
    ontology = [e for e in expanded["k"] if e.target_id == "PF00216"][0]
    assert ontology.channel == EvidenceChannel.REMOTE_HOMOLOGY
    assert ontology.details["db"] == "pfam"
    assert ontology.details["source_accession"] == "P0ACF0"


def test_expand_evidence_no_xref() -> None:
    hit = Evidence(
        channel=EvidenceChannel.REMOTE_HOMOLOGY,
        target_id="sp|P0ACF0|DBHA_ECOLI",
        target_label="x",
        score=0.8,
        evalue=1e-20,
        coverage=0.9,
        details={},
    )
    expanded, added = idmapping.expand_evidence({"k": [hit]}, {})
    assert added == 0
    assert len(expanded["k"]) == 1
