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
