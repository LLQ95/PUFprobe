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


def _protein_priority(rec: dict) -> tuple:
    """Deterministic ordering for choosing one row that shares a global key.

    A record with a translated protein sequence is preferred, then the longest
    sequence, which keeps the full-length isoform rather than an alternative
    translation start, then the earliest start and widest span. The protein id
    is a final tie breaker so the chosen row is identical on every run.
    """
    sequence = rec.get("sequence_aa") or ""
    return (
        0 if sequence else 1,
        -len(sequence),
        int(rec.get("start") or 0),
        -int(rec.get("end") or 0),
        str(rec.get("protein_id") or ""),
    )


def _dedupe_protein_records(records: List[dict]) -> List[dict]:
    """Collapse rows sharing a global key to one deterministic record.

    RefSeq GFF files can list several CDS rows for one locus tag. Rows with the
    same GFF id represent a multi-segment or interrupted pseudogene, while rows
    with distinct ids represent alternative translation starts. They describe
    the same protein locus, so the proteins table keeps one row per global key.
    """
    best: Dict[str, dict] = {}
    for rec in records:
        key = gkey(rec["sample_id"], rec["locus_tag"])
        current = best.get(key)
        if current is None or _protein_priority(rec) < _protein_priority(current):
            best[key] = rec
    return list(best.values())


def _dedupe_predictions(predictions: List[dict]) -> List[dict]:
    """Keep one prediction per global key.

    Duplicate CDS rows produce identical term predictions, so the extra copies
    would only add redundant function and evidence rows.
    """
    seen: set = set()
    unique: List[dict] = []
    for prediction in predictions:
        key = gkey(prediction["sample_id"], prediction["locus_tag"])
        if key in seen:
            continue
        seen.add(key)
        unique.append(prediction)
    return unique


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
        engine = fam.get("metrics", {}).get("engine") if fam.get("metrics") else None
        # The fold stage does not persist the engine in families.json. The GPU
        # pipeline only runs ESMFold, so a produced structure is attributed to
        # it here; a skipped family keeps a null engine.
        if engine is None and fam.get("fold_status") in ("folded", "low_confidence"):
            engine = "esmfold"
        structure_rows.append(
            (
                fam["family_id"],
                engine,
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
    seen_evidence = set()
    for prediction in _dedupe_predictions(payload.get("predictions", [])):
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
                evidence_tuple = (
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
                # Multi-segment CDS rows can repeat the same evidence object.
                if evidence_tuple in seen_evidence:
                    continue
                seen_evidence.add(evidence_tuple)
                evidence_rows.append(evidence_tuple)
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


def _intergenic_rank(raw: Any) -> int:
    """Numeric intergenic distance for ordering; a missing value sorts last."""
    text = str(raw)
    if text.lstrip("-").isdigit():
        return int(text)
    return 10 ** 12


def _insert_genomic_context(connection: sqlite3.Connection, path: Path) -> None:
    if not path.exists():
        return
    # A multi-segment CDS yields one neighbour row per segment with a different
    # intergenic distance. Keep one row per neighbour relationship, choosing the
    # nearest segment boundary, which is the geometry relevant to operon calls.
    chosen: Dict[tuple, tuple] = {}
    ordering: Dict[tuple, tuple] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            global_key = gkey(row["sample_id"], row["locus_tag"])
            intergenic_raw = row.get("intergenic_bp", "")
            dedup_key = (global_key, row.get("neighbour_locus"), row.get("position"))
            record = (
                global_key,
                row.get("sample_id"),
                row.get("locus_tag"),
                row.get("neighbour_locus"),
                row.get("neighbour_product"),
                row.get("neighbour_status"),
                row.get("position"),
                int(intergenic_raw)
                if str(intergenic_raw).lstrip("-").isdigit()
                else None,
                row.get("strand_match"),
                int(row.get("coupled", 0))
                if str(row.get("coupled")).isdigit()
                else 0,
            )
            rank = (
                _intergenic_rank(intergenic_raw),
                0 if str(row.get("coupled")).lower() == "true" else 1,
                str(row.get("neighbour_product") or ""),
            )
            current_rank = ordering.get(dedup_key)
            if current_rank is None or rank < current_rank:
                ordering[dedup_key] = rank
                chosen[dedup_key] = record
    connection.executemany(
        """
        INSERT INTO genomic_context (
            global_key, sample_id, locus_tag, neighbour_locus, neighbour_product,
            neighbour_status, position, intergenic_bp, strand_match, coupled)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        list(chosen.values()),
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
    records = _dedupe_protein_records(records)

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
