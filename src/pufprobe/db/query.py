"""Read-only query helpers for a built HPDB SQLite file.

The functions return plain dictionaries so they can be printed as JSON or used
by the command line. They cover the common lookups; arbitrary SQL is available
through the standard sqlite3 module and through Datasette.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional


def _open(db_path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(str(Path(db_path)))
    connection.row_factory = sqlite3.Row
    return connection


def get_protein(db_path: Path, identifier: str) -> Optional[Dict[str, Any]]:
    """Fetch one protein by HP id, global key, locus tag or protein id."""
    connection = _open(db_path)
    try:
        row = connection.execute(
            """
            SELECT * FROM proteins
            WHERE hp_id = ? OR global_key = ? OR locus_tag = ? OR protein_id = ?
            LIMIT 1
            """,
            (identifier, identifier, identifier, identifier),
        ).fetchone()
        return dict(row) if row else None
    finally:
        connection.close()


def get_functions(db_path: Path, global_key: str) -> List[Dict[str, Any]]:
    """Return the accepted functional terms for one protein."""
    connection = _open(db_path)
    try:
        rows = connection.execute(
            """
            SELECT term_type, term_id, term_name, raw_score, confidence,
                   go_depth, rationale
            FROM functions WHERE global_key = ?
            ORDER BY term_type, term_id
            """,
            (global_key,),
        ).fetchall()
        return [dict(row) for row in rows]
    finally:
        connection.close()


def get_evidence(db_path: Path, global_key: str) -> List[Dict[str, Any]]:
    """Return the raw evidence rows supporting one protein."""
    connection = _open(db_path)
    try:
        rows = connection.execute(
            """
            SELECT channel, channel_group, target_id, target_label, score,
                   evalue, coverage, details_json
            FROM evidence WHERE global_key = ?
            ORDER BY channel, target_id
            """,
            (global_key,),
        ).fetchall()
        return [dict(row) for row in rows]
    finally:
        connection.close()


def search_text(db_path: Path, text: str, limit: int = 100) -> List[Dict[str, Any]]:
    """Search product and term labels with a case-insensitive substring."""
    connection = _open(db_path)
    try:
        pattern = f"%{text}%"
        rows = connection.execute(
            """
            SELECT hp_id, global_key, sample_id, locus_tag, product, status,
                   family_id, species
            FROM proteins
            WHERE product LIKE ? OR global_key LIKE ? OR species LIKE ?
            ORDER BY global_key
            LIMIT ?
            """,
            (pattern, pattern, pattern, limit),
        ).fetchall()
        return [dict(row) for row in rows]
    finally:
        connection.close()


def summary(db_path: Path) -> Dict[str, Any]:
    """Return table counts and a breakdown by gene status and gram stain."""
    connection = _open(db_path)
    try:
        table_names = [
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table' ORDER BY name"
            ).fetchall()
        ]
        table_counts = {
            name: connection.execute(f"SELECT COUNT(*) FROM {name}").fetchone()[0]
            for name in table_names
        }
        by_status = {
            row[0]: row[1]
            for row in connection.execute(
                "SELECT status, COUNT(*) FROM proteins GROUP BY status ORDER BY status"
            ).fetchall()
        }
        by_gram = {
            row[0]: row[1]
            for row in connection.execute(
                "SELECT gram, COUNT(*) FROM proteins GROUP BY gram ORDER BY gram"
            ).fetchall()
        }
        return {
            "tables": table_counts,
            "by_status": by_status,
            "by_gram": by_gram,
        }
    finally:
        connection.close()
