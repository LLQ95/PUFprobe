"""Build a HPDB SQLite file from a completed PUFprobe working directory.

The builder reads the frozen stage outputs rather than recomputing anything.
Each stage file is optional except the prepare CDS records, so a partial run
still populates the tables whose inputs are available. Rows are inserted in a
single transaction and every ordering is deterministic, so rebuilding the same
working directory twice yields the same database content.
"""

from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional

from pufprobe.db.schema import connect
from pufprobe.io_utils import read_json
from pufprobe.logging_utils import get_logger
from pufprobe.models import CHANNEL_GROUPS, gkey
from pufprobe.workdir import STAGE_DIRS

log = get_logger("db.build")


def _stage_file(root: Path, stage: str, name: str) -> Path:
    return root / STAGE_DIRS[stage] / name


def _load_optional_json(path: Path) -> Optional[Any]:
    if not path.exists():
        return None
    return read_json(path)


def load_taxonomy_map(manifest_path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    """Map the PUFprobe sample id to taxonomy fields from a sampling manifest."""
    if manifest_path is None or not manifest_path.exists():
        return {}
    rows = read_json(manifest_path)
    mapping: Dict[str, Dict[str, str]] = {}
    for row in rows:  # type: ignore[union-attr]
        sample_dir = row.get("sample_dir")
        sample_id = Path(sample_dir).name if sample_dir else None
        if not sample_id:
            continue
        mapping[sample_id] = {
            "species": row.get("species", ""),
            "phylum": row.get("phylum", ""),
            "gram": row.get("gram", ""),
            "taxid": str(row.get("taxid", "")),
        }
    return mapping


def _capture_environment() -> Optional[Dict[str, Any]]:
    """Snapshot external tool versions; never let detection break the build."""
    try:
        from pufprobe.externals import summary

        return summary()
    except Exception as exc:  # pragma: no cover - detection is best effort
        log.warning("External environment snapshot skipped: %s", exc)
        return None


def _insert_runs(connection: sqlite3.Connection, root: Path, config: Optional[Dict[str, Any]]) -> None:
    metadata = _load_optional_json(root / "metadata.json")
    if not metadata:
        return
    run_id = metadata.get("created") or "run"
    environment = _capture_environment()
    connection.execute(
        """
        INSERT OR REPLACE INTO runs
            (run_id, pufprobe_version, python_version, platform, cpu_mode,
             threads, created, config_json, environment_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            run_id,
            metadata.get("pufprobe_version"),
            metadata.get("python"),
            metadata.get("platform"),
            int(bool(metadata.get("cpu_mode"))),
            metadata.get("threads"),
            metadata.get("created"),
            json.dumps(config, sort_keys=True) if config is not None else None,
            json.dumps(environment, sort_keys=True) if environment is not None else None,
        ),
    )


def _insert_proteins(
    connection: sqlite3.Connection,
    records: List[dict],
    taxonomy: Dict[str, Dict[str, str]],
) -> Dict[str, str]:
    ordered = sorted(records, key=lambda r: gkey(r["sample_id"], r["locus_tag"]))
    protein_rows = []
    key_to_hp: Dict[str, str] = {}
    for index, rec in enumerate(ordered, start=1):
        global_key = gkey(rec["sample_id"], rec["locus_tag"])
        hp_id = f"HP-{index:06d}"
        key_to_hp[global_key] = hp_id
        tax = taxonomy.get(rec["sample_id"], {})
        sequence = rec.get("sequence_aa") or ""
        protein_rows.append(
            (
                hp_id,
                global_key,
                rec["sample_id"],
                rec["locus_tag"],
                rec.get("protein_id"),
                sequence,
                len(sequence),
                rec.get("product"),
                rec.get("status"),
                rec.get("verdict"),
                rec.get("qc_score"),
                rec.get("family_id"),
                tax.get("species"),
                tax.get("phylum"),
                tax.get("gram"),
                rec.get("contig"),
                rec.get("start"),
                rec.get("end"),
                rec.get("strand"),
            )
        )
    connection.executemany(
        """
        INSERT INTO proteins (
            hp_id, global_key, sample_id, locus_tag, protein_id, sequence, length,
            product, status, verdict, qc_score, family_id, species, phylum, gram,
            contig, start, end, strand)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        protein_rows,
    )
    return key_to_hp


def _insert_families(connection: sqlite3.Connection, families: List[dict]) -> None:
    family_rows = []
    structure_rows = []
    for fam in sorted(families, key=lambda f: f["family_id"]):
        family_rows.append(
            (
                fam["family_id"],
                fam.get("representative_key"),
                fam.get("size"),
                fam.get("lineage_width"),
                fam.get("fold_status"),
                fam.get("plddt"),
                fam.get("structure_path"),
                fam.get("novel_score"),
                fam.get("priority_score"),
            )
        )
        structure_rows.append(
            (
                fam["family_id"],
                fam.get("metrics", {}).get("engine") if fam.get("metrics") else None,
                fam.get("fold_status"),
                fam.get("plddt"),
                fam.get("structure_path"),
            )
        )
    connection.executemany(
        """
        INSERT INTO families (
            family_id, representative_key, size, lineage_width, fold_status,
            plddt, structure_path, novel_score, priority_score)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        family_rows,
    )
    connection.executemany(
        """
        INSERT OR REPLACE INTO structures (
            family_id, engine, fold_status, plddt, structure_path)
        VALUES (?, ?, ?, ?, ?)
        """,
        structure_rows,
    )


def _insert_predictions(connection: sqlite3.Connection, payload: Dict[str, Any]) -> None:
    function_rows = []
    evidence_rows = []
    seen_functions = set()
    for prediction in payload.get("predictions", []):
        global_key = gkey(prediction["sample_id"], prediction["locus_tag"])
        for term in prediction.get("terms", []):
            func_key = (global_key, term["term_type"], term["term_id"])
            if func_key not in seen_functions:
                seen_functions.add(func_key)
                function_rows.append(
                    (
                        global_key,
                        term["term_type"],
                        term["term_id"],
                        term.get("term_name"),
                        term.get("raw_score"),
                        term.get("confidence"),
                        term.get("go_depth"),
                        term.get("rationale"),
                    )
                )
            for ev in term.get("evidence", []):
                channel = ev.get("channel", "")
                evidence_rows.append(
                    (
                        global_key,
                        term["term_type"],
                        term["term_id"],
                        channel,
                        CHANNEL_GROUPS.get(channel),
                        ev.get("target_id"),
                        ev.get("target_label"),
                        ev.get("score"),
                        ev.get("evalue"),
                        ev.get("coverage"),
                        json.dumps(ev.get("details", {}), sort_keys=True),
                    )
                )
    connection.executemany(
        """
        INSERT INTO functions (
            global_key, term_type, term_id, term_name, raw_score, confidence,
            go_depth, rationale)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        function_rows,
    )
    connection.executemany(
        """
        INSERT INTO evidence (
            global_key, term_type, term_id, channel, channel_group, target_id,
            target_label, score, evalue, coverage, details_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        evidence_rows,
    )


def _insert_genomic_context(connection: sqlite3.Connection, path: Path) -> None:
    if not path.exists():
        return
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            global_key = gkey(row["sample_id"], row["locus_tag"])
            intergenic = row.get("intergenic_bp", "")
            rows.append(
                (
                    global_key,
                    row.get("sample_id"),
                    row.get("locus_tag"),
                    row.get("neighbour_locus"),
                    row.get("neighbour_product"),
                    row.get("neighbour_status"),
                    row.get("position"),
                    int(intergenic) if intergenic.lstrip("-").isdigit() else None,
                    row.get("strand_match"),
                    int(row.get("coupled", 0)) if str(row.get("coupled")).isdigit() else 0,
                )
            )
    connection.executemany(
        """
        INSERT INTO genomic_context (
            global_key, sample_id, locus_tag, neighbour_locus, neighbour_product,
            neighbour_status, position, intergenic_bp, strand_match, coupled)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )


def build_database(
    workdir: Path,
    db_path: Path,
    manifest_path: Optional[Path] = None,
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, int]:
    """Build the SQLite database and return per-table row counts."""
    root = Path(workdir)
    cds_path = _stage_file(root, "prepare", "cds_records.json")
    records = _load_optional_json(cds_path)
    if records is None:
        raise FileNotFoundError(
            f"Missing {cds_path}; run the prepare stage before building the database."
        )

    connection = connect(db_path)
    # The build is a controlled bulk import and supports partial runs, so a
    # protein can reference a family whose own stage output was not produced.
    # Foreign key DDL is still declared for querying clients; enforcement is
    # disabled only during the import.
    connection.execute("PRAGMA foreign_keys = OFF")
    try:
        with connection:
            taxonomy = load_taxonomy_map(manifest_path)
            _insert_runs(connection, root, config)

            # Parent tables first: families precede proteins, which precede the
            # term, evidence and context tables that reference them.
            families = _load_optional_json(_stage_file(root, "cluster", "families.json"))
            if families:
                _insert_families(connection, families)

            _insert_proteins(connection, records, taxonomy)

            evidence_payload = _load_optional_json(
                _stage_file(root, "report", "evidence.json")
            )
            if evidence_payload:
                _insert_predictions(connection, evidence_payload)

            _insert_genomic_context(
                connection, _stage_file(root, "context", "neighbours.tsv")
            )
    finally:
        connection.close()

    counts = table_counts(db_path)
    log.info("HPDB built at %s: %s", db_path, counts)
    return counts


def table_counts(db_path: Path) -> Dict[str, int]:
    """Return the row count for every table in the built database."""
    tables = [
        "runs",
        "proteins",
        "families",
        "functions",
        "evidence",
        "structures",
        "genomic_context",
    ]
    connection = sqlite3.connect(str(db_path))
    try:
        return {
            table: connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            for table in tables
        }
    finally:
        connection.close()
