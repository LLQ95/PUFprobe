"""Tests for report HTML escaping and annotate backfill behaviour."""

import hashlib
import json
from pathlib import Path

import pytest

from pufprobe.config import load_config
from pufprobe.models import CDSRecord, GeneStatus, to_dict
from pufprobe.workdir import WorkDir


GFF3 = """\
##gff-version 3
sA\tprodigal\tregion\t1\t5000\t.\t+\t.\tID=region1
sA\tprodigal\tCDS\t100\t400\t.\t+\t0\tID=sA_0001;locus_tag=sA_0001;product=hypothetical protein
sA\tprodigal\tCDS\t500\t800\t.\t-\t0\tID=sA_0002;locus_tag=sA_0002;product=old enzyme
"""

FAA = """\
>sA_0001 hypothetical protein
MSEQUENCEAA
>sA_0002 old enzyme
MKAAAAAA
"""

PREDICTIONS = [
    {
        "sample_id": "sA",
        "locus_tag": "sA_0001",
        "terms": [
            {
                "term_type": "Pfam",
                "term_id": "PF00010",
                "term_name": "DUF1000 <putative> & co",
                "raw_score": 0.95,
                "confidence": 0.95,
                "evidence": [
                    {
                        "channel": "database",
                        "target_id": "UniRef90:X9",
                        "target_label": "DUF1000 <putative> & co",
                        "score": 0.95,
                        "evalue": 1e-5,
                        "coverage": 0.8,
                        "details": {},
                    }
                ],
                "rationale": "score: test",
            }
        ],
    },
    {"sample_id": "sA", "locus_tag": "sA_0002", "terms": []},
]


@pytest.fixture()
def workdir(tmp_path: Path) -> WorkDir:
    src = tmp_path / "sources"
    src.mkdir()
    gff = src / "sampleA.gff3"
    faa = src / "sampleA.faa"
    gff.write_text(GFF3, encoding="utf-8", newline="\n")
    faa.write_text(FAA, encoding="utf-8", newline="\n")

    wd = WorkDir(tmp_path / "run")
    recs = [
        CDSRecord(
            sample_id="sA",
            locus_tag=locus,
            protein_id=None,
            sequence_aa="M" * 100,
            start=100,
            end=400,
            strand="+",
            product=prod,
            status=GeneStatus.HYPOTHETICAL,
        )
        for locus, prod in [("sA_0001", "hypothetical protein"), ("sA_0002", "old enzyme")]
    ]
    (wd.path("prepare", "cds_records.json")).write_text(
        json.dumps([to_dict(r) for r in recs]), encoding="utf-8"
    )
    (wd.path("prepare", "samples.json")).write_text(
        json.dumps(
            [
                {
                    "sample_id": "sA",
                    "gff_path": str(gff),
                    "gbk_path": None,
                    "faa_path": str(faa),
                    "fna_path": None,
                    "source": "generic",
                }
            ]
        ),
        encoding="utf-8",
    )
    (wd.path("rank", "predictions.json")).write_text(
        json.dumps(PREDICTIONS), encoding="utf-8"
    )
    wd._gff = gff
    return wd


def test_report_html_escaping_and_sections(workdir):
    from pufprobe.report.pipeline import run

    run(workdir, load_config(None), ["json", "html", "md"])

    html = workdir.path("report", "report.html").read_text("utf-8")
    assert "PUFprobe evidence report" in html
    assert "Overview" in html
    # special characters must be escaped
    assert "&lt;putative&gt; &amp; co" in html
    assert "DUF1000 <putative>" not in html
    # auditable sections
    assert "Per-protein evidence" in html
    assert "Channel" in html and "evalue" in html

    evidence = json.loads(workdir.path("report", "evidence.json").read_text("utf-8"))
    assert len(evidence["predictions"]) == 2
    assert len(evidence["index"]) == 2
    assert evidence["overview"]["proteins_with_terms"] == 1

    md = workdir.path("report", "report.md").read_text("utf-8")
    assert "PUFprobe evidence report" in md


def test_annotate_backfill_only_matched_rows(workdir):
    from pufprobe.annotate.pipeline import run

    original = workdir._gff.read_bytes()
    original_hash = hashlib.sha256(original).hexdigest()

    run(workdir, load_config(None))

    out_gff = workdir.path("annotate", "sA.pufprobe.gff3").read_text("utf-8").splitlines()
    # header + region line untouched
    assert out_gff[0] == "##gff-version 3"
    assert out_gff[1].endswith("ID=region1")
    # unmatched CDS untouched
    assert out_gff[3].endswith("product=old enzyme")

    matched = out_gff[2]
    assert "product=DUF1000 <putative> & co" in matched
    assert "Dbxref=Pfam:PF00010" in matched
    assert "PUFprobe:confidence=0.950" in matched

    # original file untouched
    assert hashlib.sha256(workdir._gff.read_bytes()).hexdigest() == original_hash

    changes = workdir.path("annotate", "changes.tsv").read_text("utf-8").splitlines()
    assert changes[0] == "sample\tlocus\told_product\tnew_product\tconfidence\tsupporting terms"
    assert changes[1].startswith("sA\tsA_0001\thypothetical protein\t")
    assert "Pfam:PF00010" in changes[1]
    assert len(changes) == 2  # header + exactly one edit

    faa = workdir.path("annotate", "sA.pufprobe.faa").read_text("utf-8")
    assert ">sA_0001 DUF1000 <putative> & co" in faa
    assert ">sA_0002 old enzyme" in faa


# ---------------------------------------------------------------------------
# Product-name granularity: shallow GO terms must not be written back.
# ---------------------------------------------------------------------------
def _term(ttype, tid, name, conf, go_depth=None):
    from pufprobe.models import TermPrediction

    return TermPrediction(
        term_type=ttype, term_id=tid, term_name=name,
        raw_score=conf, confidence=conf, go_depth=go_depth,
    )


def _prediction(sample, locus, terms):
    from pufprobe.models import ProteinPrediction

    return ProteinPrediction(sample_id=sample, locus_tag=locus, terms=terms)


def test_select_terms_filters_shallow_go_terms():
    from pufprobe.annotate.backfill import select_terms

    pred = _prediction("sA", "g1", [
        _term("GO", "GO:0003674", "molecular_function", 0.9, go_depth=0),
        _term("GO", "GO:0003824", "catalytic activity", 0.9, go_depth=1),
        _term("GO", "GO:0016829", "specific transferase", 0.9, go_depth=3),
    ])
    accepted = select_terms([pred], 0.8, 2)
    assert list(accepted) == ["sA@g1"]
    edit = accepted["sA@g1"]
    assert edit["new_product"] == "specific transferase"
    assert edit["term_ids"] == ["GO:GO:0016829"]


def test_select_terms_skips_protein_with_only_generic_terms():
    from pufprobe.annotate.backfill import select_terms

    pred = _prediction("sA", "g1", [
        _term("GO", "GO:0005575", "cellular_component", 0.9, go_depth=0),
        _term("GO", "GO:0016020", "membrane", 0.9, go_depth=1),
    ])
    assert select_terms([pred], 0.8, 2) == {}


def test_select_terms_keeps_unknown_depth_and_non_go():
    from pufprobe.annotate.backfill import select_terms

    pred = _prediction("sA", "g1", [
        _term("GO", "GO:0012345", "unknown depth term", 0.9, go_depth=None),
        _term("Pfam", "PF00001", "some family", 0.9),
        _term("EC", "EC:1.2.3.4", "some enzyme", 0.9),
    ])
    edit = select_terms([pred], 0.8, 2)["sA@g1"]
    assert edit["term_ids"] == [
        "EC:EC:1.2.3.4", "GO:GO:0012345", "Pfam:PF00001",
    ]


def test_select_terms_deduplicates_names():
    from pufprobe.annotate.backfill import select_terms

    pred = _prediction("sA", "g1", [
        _term("GO", "GO:0000001", "shared name", 0.9, go_depth=2),
        _term("GO", "GO:0000002", "shared name", 0.9, go_depth=3),
    ])
    edit = select_terms([pred], 0.8, 2)["sA@g1"]
    assert edit["new_product"] == "shared name"
