"""ontology_names: Pfam family name parsing and term-name resolution."""

import io
import json

from pufprobe import ontology_names

DAT = """# STOCKHOLM 1.0
#=GF ID   ABC_tran
#=GF AC   PF00005.34
#=GF DE   ABC transporter
#=GF TP   Domain
//
# STOCKHOLM 1.0
#=GF ID   Some_Fam
#=GF AC   PF99999.1
#=GF TP   Family
//
"""


def test_strip_version() -> None:
    assert ontology_names._strip_version("PF00005.34") == "PF00005"
    assert ontology_names._strip_version("PF00005") == "PF00005"


def test_parse_pfam_dat() -> None:
    names = ontology_names.parse_pfam_dat(io.StringIO(DAT))
    # Accession is stored without its version number.
    assert names["PF00005"] == {
        "id": "ABC_tran",
        "name": "ABC transporter",
        "type": "Domain",
    }
    # A record without a description keeps its id and type.
    assert names["PF99999"] == {"id": "Some_Fam", "type": "Family"}


def test_pfam_term_name() -> None:
    names = ontology_names.parse_pfam_dat(io.StringIO(DAT))
    # Official description is preferred.
    assert ontology_names.pfam_term_name(names, "PF00005") == "ABC transporter"
    # Versioned accessions resolve to the same record.
    assert ontology_names.pfam_term_name(names, "PF00005.34") == "ABC transporter"
    # Without a description the short identifier is used.
    assert ontology_names.pfam_term_name(names, "PF99999") == "Some_Fam"
    # Unknown accessions return None so the caller can fall back.
    assert ontology_names.pfam_term_name(names, "PF00000") is None
    assert ontology_names.pfam_term_name({}, "PF00005") is None


def test_build_and_load_pfam_names(tmp_path) -> None:
    dat_path = tmp_path / "Pfam-A.hmm.dat"
    out_path = tmp_path / "nested" / "pfam_names.json"
    dat_path.write_text(DAT, encoding="utf-8")

    names = ontology_names.build_pfam_names(dat_path, out_path)
    assert out_path.exists()
    assert set(names) == {"PF00005", "PF99999"}

    loaded = ontology_names.load_pfam_names(out_path)
    assert loaded == names
    # Missing and malformed inputs degrade to an empty dict.
    assert ontology_names.load_pfam_names(None) == {}
    assert ontology_names.load_pfam_names(tmp_path / "missing.json") == {}
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(["not", "a", "dict"]), encoding="utf-8")
    assert ontology_names.load_pfam_names(bad) == {}
