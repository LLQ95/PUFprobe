"""Tests for HPDB SQLite build and query support."""

import json
from pathlib import Path

import pytest

from pufprobe.config import load_config
from pufprobe.models import CDSRecord, GeneStatus, to_dict
from pufprobe.workdir import WorkDir


SAMPLE_ID = "Escherichia_coli__GCF_000005845.2"


def _cds(locus_tag, product, status, family_id=None):
    return CDSRecord(
        sample_id=SAMPLE_ID,
        locus_tag=locus_tag,
        protein_id=None,
        sequence_aa="M" * 40,
        start=100,
        end=220,
        strand="+",
        product=product,
        status=status,
        contig="NC_000913.3",
        family_id=family_id,
    )


def _evidence_payload():
    return {
        "generated": "2026-10-19T10:00:00",
        "overview": {},
        "predictions": [
            {
                "sample_id": SAMPLE_ID,
                "locus_tag": "b0001",
                "terms": [
                    {
                        "term_type": "Pfam",
                        "term_id": "PF00010",
                        "term_name": "Helix-turn-helix",
                        "raw_score": 0.9,
                        "confidence": 0.9,
                        "go_depth": None,
                        "rationale": "test",
                        "evidence": [
                            {
                                "channel": "foldseek",
                                "target_id": "PDB:1ABC",
                                "target_label": "structure hit",
                                "score": 0.85,
                                "evalue": 1e-4,
                                "coverage": 0.9,
                                "details": {"tmscore": 0.62},
                            }
                        ],
                    }
                ],
                "priority_score": 0.5,
                "priority_reasons": [],
            }
        ],
        "index": [],
    }


def _manifest():
    return [
        {
            "species": "Escherichia coli",
            "taxid": "562",
            "gram": "negative",
            "phylum": "Pseudomonadota",
            "accession": "GCF_000005845.2",
            "sample_dir": str(Path("/tmp/genomes") / SAMPLE_ID),
        }
    ]


@pytest.fixture()
def workdir(tmp_path: Path) -> WorkDir:
    wd = WorkDir(tmp_path / "run")
    records = [
        _cds("b0001", "hypothetical protein", GeneStatus.HYPOTHETICAL, "FAM_0001"),
        _cds("b0002", "old enzyme", GeneStatus.ANNOTATED, None),
    ]
    wd.path("prepare", "cds_records.json").write_text(
        json.dumps([to_dict(r) for r in records]), encoding="utf-8"
    )
    family = {
        "family_id": "FAM_0001",
        "representative_key": f"{SAMPLE_ID}@b0001",
        "members": [f"{SAMPLE_ID}@b0001"],
        "representative_seq": "M" * 40,
        "size": 1,
        "lineage_width": 1,
        "metrics": {"engine": "esmfold"},
        "structure_path": None,
        "plddt": 82.0,
        "fold_status": "folded",
        "novel_score": 0.3,
        "priority_score": 0.5,
    }
    wd.path("cluster", "families.json").write_text(
        json.dumps([family]), encoding="utf-8"
    )
    wd.path("report", "evidence.json").write_text(
        json.dumps(_evidence_payload()), encoding="utf-8"
    )
    neighbours = (
        "sample_id\tlocus_tag\tfamily_id\tneighbour_locus\tneighbour_product\t"
        "neighbour_status\tposition\tintergenic_bp\tstrand_match\tcoupled\n"
        f"{SAMPLE_ID}\tb0001\tFAM_0001\tb0002\told enzyme\tannotated\t"
        "downstream\t12\tTrue\t1\n"
    )
    wd.path("context", "neighbours.tsv").write_text(neighbours, encoding="utf-8")
    (wd.root / "metadata.json").write_text(
        json.dumps(
            {
                "pufprobe_version": "0.10.1",
                "created": "2026-10-19T09:00:00",
                "python": "3.11.9",
                "platform": "Linux",
                "cpu_mode": False,
                "threads": 8,
            }
        ),
        encoding="utf-8",
    )
    return wd


def test_build_populates_all_tables(workdir, tmp_path):
    from pufprobe.db.build import build_database

    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(_manifest()), encoding="utf-8")
    db_path = tmp_path / "hpdb.sqlite"

    counts = build_database(
        workdir.root, db_path, manifest_path=manifest, config=load_config(None)
    )

    assert counts["proteins"] == 2
    assert counts["families"] == 1
    assert counts["structures"] == 1
    assert counts["functions"] == 1
    assert counts["evidence"] == 1
    assert counts["genomic_context"] == 1
    assert counts["runs"] == 1


def test_proteins_get_hp_ids_and_taxonomy(workdir, tmp_path):
    from pufprobe.db.build import build_database
    from pufprobe.db.query import get_protein

    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(_manifest()), encoding="utf-8")
    db_path = tmp_path / "hpdb.sqlite"
    build_database(workdir.root, db_path, manifest_path=manifest)

    protein = get_protein(db_path, "HP-000001")
    assert protein is not None
    assert protein["hp_id"] == "HP-000001"
    assert protein["species"] == "Escherichia coli"
    assert protein["phylum"] == "Pseudomonadota"
    assert protein["gram"] == "negative"
    assert protein["length"] == 40


def test_functions_and_evidence_queries(workdir, tmp_path):
    from pufprobe.db.build import build_database
    from pufprobe.db.query import get_evidence, get_functions

    db_path = tmp_path / "hpdb.sqlite"
    build_database(workdir.root, db_path)

    global_key = f"{SAMPLE_ID}@b0001"
    functions = get_functions(db_path, global_key)
    assert len(functions) == 1
    assert functions[0]["term_id"] == "PF00010"

    evidence = get_evidence(db_path, global_key)
    assert len(evidence) == 1
    assert evidence[0]["channel"] == "foldseek"
    assert evidence[0]["channel_group"] == "structure"
    assert json.loads(evidence[0]["details_json"])["tmscore"] == 0.62


def test_text_search_and_summary(workdir, tmp_path):
    from pufprobe.db.build import build_database
    from pufprobe.db.query import search_text, summary

    db_path = tmp_path / "hpdb.sqlite"
    build_database(workdir.root, db_path)

    hits = search_text(db_path, "hypothetical")
    assert len(hits) == 1
    assert hits[0]["locus_tag"] == "b0001"

    result = summary(db_path)
    assert result["tables"]["proteins"] == 2
    assert result["by_status"]["annotated"] == 1


def test_build_is_deterministic(workdir, tmp_path):
    from pufprobe.db.build import build_database

    first = tmp_path / "a.sqlite"
    second = tmp_path / "b.sqlite"
    build_database(workdir.root, first)
    build_database(workdir.root, second)

    assert first.read_bytes() == second.read_bytes()


def test_run_records_config_and_environment(workdir, tmp_path):
    import sqlite3

    from pufprobe.db.build import build_database

    db_path = tmp_path / "hpdb.sqlite"
    build_database(
        workdir.root, db_path, config=load_config(None)
    )

    connection = sqlite3.connect(db_path)
    row = connection.execute(
        "SELECT config_json, environment_json FROM runs"
    ).fetchone()
    connection.close()

    config_json, environment_json = row
    assert json.loads(config_json)["fold"]["engine"] == "esmfold"
    environment = json.loads(environment_json)
    assert "tools" in environment
    assert "python_modules" in environment


def test_build_without_manifest_leaves_taxonomy_empty(workdir, tmp_path):
    from pufprobe.db.build import build_database
    from pufprobe.db.query import get_protein

    db_path = tmp_path / "hpdb.sqlite"
    build_database(workdir.root, db_path)

    protein = get_protein(db_path, "HP-000001")
    assert protein["species"] is None
    assert protein["gram"] is None


def test_missing_cds_records_raises(tmp_path):
    from pufprobe.db.build import build_database

    empty = WorkDir(tmp_path / "run")
    with pytest.raises(FileNotFoundError):
        build_database(empty.root, tmp_path / "out.sqlite")


def test_partial_workdir_with_only_cds(workdir, tmp_path):
    from pufprobe.db.build import build_database

    # Remove the optional stage outputs to confirm the build still works.
    workdir.path("cluster", "families.json").unlink()
    workdir.path("report", "evidence.json").unlink()
    workdir.path("context", "neighbours.tsv").unlink()

    counts = build_database(workdir.root, tmp_path / "out.sqlite")
    assert counts["proteins"] == 2
    assert counts["families"] == 0
    assert counts["functions"] == 0


def test_dedupe_protein_records_keeps_longest_with_sequence():
    from pufprobe.db.build import _dedupe_protein_records

    repeated = [
        {
            "sample_id": SAMPLE_ID,
            "locus_tag": "b0009",
            "protein_id": "NP_SHORT",
            "sequence_aa": "M" * 20,
            "start": 300,
            "end": 360,
        },
        {
            "sample_id": SAMPLE_ID,
            "locus_tag": "b0009",
            "protein_id": "NP_LONG",
            "sequence_aa": "M" * 80,
            "start": 100,
            "end": 340,
        },
    ]
    other = {
        "sample_id": SAMPLE_ID,
        "locus_tag": "b0010",
        "protein_id": None,
        "sequence_aa": "M" * 10,
        "start": 1,
        "end": 30,
    }
    result = _dedupe_protein_records(repeated + [other])
    assert len(result) == 2
    chosen = next(r for r in result if r["locus_tag"] == "b0009")
    assert chosen["protein_id"] == "NP_LONG"
    assert len(chosen["sequence_aa"]) == 80


def test_dedupe_empty_pseudogene_keeps_earliest_segment():
    from pufprobe.db.build import _dedupe_protein_records

    records = [
        {
            "sample_id": SAMPLE_ID,
            "locus_tag": "b0020",
            "protein_id": None,
            "sequence_aa": "",
            "start": 500,
            "end": 560,
        },
        {
            "sample_id": SAMPLE_ID,
            "locus_tag": "b0020",
            "protein_id": None,
            "sequence_aa": "",
            "start": 100,
            "end": 160,
        },
    ]
    result = _dedupe_protein_records(records)
    assert len(result) == 1
    assert result[0]["start"] == 100


def test_dedupe_predictions_keeps_one_per_global_key():
    from pufprobe.db.build import _dedupe_predictions

    predictions = [
        {"sample_id": SAMPLE_ID, "locus_tag": "b0001", "terms": []},
        {"sample_id": SAMPLE_ID, "locus_tag": "b0001", "terms": []},
        {"sample_id": SAMPLE_ID, "locus_tag": "b0002", "terms": []},
    ]
    result = _dedupe_predictions(predictions)
    assert len(result) == 2
    assert [p["locus_tag"] for p in result] == ["b0001", "b0002"]


def test_build_deduplicates_repeated_global_keys(tmp_path):
    from pufprobe.db.build import build_database

    wd = WorkDir(tmp_path / "run")
    records = [
        {
            "sample_id": SAMPLE_ID,
            "locus_tag": "b0001",
            "protein_id": "NP_1",
            "sequence_aa": "M" * 40,
            "start": 100,
            "end": 220,
            "strand": "+",
            "product": "enzyme",
            "status": "annotated",
            "contig": "NC_1",
        },
        {
            "sample_id": SAMPLE_ID,
            "locus_tag": "b0002",
            "protein_id": "NP_2",
            "sequence_aa": "M" * 40,
            "start": 300,
            "end": 420,
            "strand": "+",
            "product": "hypothetical protein",
            "status": "hypothetical",
            "contig": "NC_1",
        },
        {
            "sample_id": SAMPLE_ID,
            "locus_tag": "b0002",
            "protein_id": "NP_2",
            "sequence_aa": "M" * 40,
            "start": 421,
            "end": 500,
            "strand": "+",
            "product": "hypothetical protein",
            "status": "hypothetical",
            "contig": "NC_1",
        },
    ]
    wd.path("prepare", "cds_records.json").write_text(
        json.dumps(records), encoding="utf-8"
    )

    term = {
        "term_type": "label",
        "term_id": "b0002",
        "term_name": "hypothetical protein",
        "raw_score": 1.0,
        "confidence": 0.6,
        "go_depth": None,
        "rationale": "",
        "evidence": [
            {
                "channel": "context",
                "target_id": "b0002",
                "target_label": "hypothetical protein",
                "score": 1.0,
                "evalue": None,
                "coverage": None,
                "details": {"kind": "operon_coupling", "intergenic_bp": 1},
            }
        ],
    }
    payload = {
        "predictions": [
            {"sample_id": SAMPLE_ID, "locus_tag": "b0002", "terms": [dict(term)]},
            {"sample_id": SAMPLE_ID, "locus_tag": "b0002", "terms": [dict(term)]},
        ]
    }
    wd.path("report", "evidence.json").write_text(
        json.dumps(payload), encoding="utf-8"
    )

    neighbours = (
        "sample_id\tlocus_tag\tfamily_id\tneighbour_locus\tneighbour_product\t"
        "neighbour_status\tposition\tintergenic_bp\tstrand_match\tcoupled\n"
        f"{SAMPLE_ID}\tb0002\t\tb0001\tenzyme\tannotated\tdownstream\t12\tTrue\t1\n"
        f"{SAMPLE_ID}\tb0002\t\tb0001\tenzyme\tannotated\tdownstream\t200\tTrue\t0\n"
    )
    wd.path("context", "neighbours.tsv").write_text(
        neighbours, encoding="utf-8"
    )

    counts = build_database(wd.root, tmp_path / "hpdb.sqlite")
    assert counts["proteins"] == 2
    assert counts["functions"] == 1
    assert counts["evidence"] == 1
    assert counts["genomic_context"] == 1


def test_get_family_lists_members(workdir, tmp_path):
    from pufprobe.db.build import build_database
    from pufprobe.db.query import get_family

    db_path = tmp_path / "hpdb.sqlite"
    build_database(workdir.root, db_path)

    result = get_family(db_path, "FAM_0001")
    assert result is not None
    assert result["family"]["family_id"] == "FAM_0001"
    assert len(result["members"]) == 1
    assert result["members"][0]["locus_tag"] == "b0001"
    assert get_family(db_path, "MISSING") is None


def test_search_terms_finds_protein(workdir, tmp_path):
    from pufprobe.db.build import build_database
    from pufprobe.db.query import search_terms

    db_path = tmp_path / "hpdb.sqlite"
    build_database(workdir.root, db_path)

    by_id = search_terms(db_path, "PF00010")
    assert len(by_id) == 1
    assert by_id[0]["locus_tag"] == "b0001"
    by_name = search_terms(db_path, "helix-turn-helix")
    assert len(by_name) == 1
    assert by_name[0]["term_id"] == "PF00010"


def test_get_neighbours_lists_relationships(workdir, tmp_path):
    from pufprobe.db.build import build_database
    from pufprobe.db.query import get_neighbours

    db_path = tmp_path / "hpdb.sqlite"
    build_database(workdir.root, db_path)

    global_key = f"{SAMPLE_ID}@b0001"
    neighbours = get_neighbours(db_path, global_key)
    assert len(neighbours) == 1
    assert neighbours[0]["neighbour_locus"] == "b0002"
    assert neighbours[0]["neighbour_product"] == "old enzyme"


def test_db_query_cli_lookups(workdir, tmp_path):
    from click.testing import CliRunner

    from pufprobe.cli import cli
    from pufprobe.db.build import build_database

    db_path = tmp_path / "hpdb.sqlite"
    build_database(workdir.root, db_path)
    runner = CliRunner()
    global_key = f"{SAMPLE_ID}@b0001"

    cases = [
        ["db", "query", "--db", str(db_path), "--summary"],
        ["db", "query", "--db", str(db_path), "--id", "HP-000001"],
        ["db", "query", "--db", str(db_path), "--text", "hypothetical"],
        ["db", "query", "--db", str(db_path), "--family", "FAM_0001"],
        ["db", "query", "--db", str(db_path), "--term", "PF00010"],
        ["db", "query", "--db", str(db_path), "--neighbours", global_key],
    ]
    for args in cases:
        result = runner.invoke(cli, args)
        assert result.exit_code == 0, f"{args}: {result.output}"

    missing = runner.invoke(
        cli, ["db", "query", "--db", str(db_path), "--family", "NOPE"]
    )
    assert missing.exit_code != 0

    no_lookup = runner.invoke(cli, ["db", "query", "--db", str(db_path)])
    assert no_lookup.exit_code != 0


def test_export_fasta_writes_sequences(workdir, tmp_path):
    from pufprobe.db.build import build_database
    from pufprobe.db.query import export_fasta

    db_path = tmp_path / "hpdb.sqlite"
    build_database(workdir.root, db_path)
    fasta_path = tmp_path / "proteins.fasta"

    count = export_fasta(db_path, fasta_path)
    assert count == 2
    text = fasta_path.read_text(encoding="utf-8")
    assert ">HP-000001" in text
    assert f"{SAMPLE_ID}@b0001" in text
    # Sequences are wrapped at 60 residues; the 40-residue test fits one line.
    assert "M" * 40 + "\n" in text
