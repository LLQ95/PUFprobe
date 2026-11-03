"""SQLite schema for the hypothetical protein database (HPDB).

The schema is expressed as ordered CREATE statements so a fresh database can be
created with the standard library sqlite3 module. Foreign keys are declared for
documentation and for clients that enable PRAGMA foreign_keys; the build itself
inserts parent rows before child rows and does not rely on cascade behaviour.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import List

CREATE_TABLES: List[str] = [
    """
    CREATE TABLE IF NOT EXISTS runs (
        run_id          TEXT PRIMARY KEY,
        pufprobe_version TEXT,
        python_version  TEXT,
        platform        TEXT,
        cpu_mode        INTEGER,
        threads         INTEGER,
        created         TEXT,
        config_json     TEXT,
        environment_json TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS proteins (
        hp_id       TEXT PRIMARY KEY,
        global_key  TEXT NOT NULL UNIQUE,
        sample_id   TEXT NOT NULL,
        locus_tag   TEXT NOT NULL,
        protein_id  TEXT,
        sequence    TEXT,
        length      INTEGER,
        product     TEXT,
        status      TEXT,
        verdict     TEXT,
        qc_score    REAL,
        family_id   TEXT,
        species     TEXT,
        phylum      TEXT,
        gram        TEXT,
        contig      TEXT,
        start       INTEGER,
        end         INTEGER,
        strand      TEXT,
        FOREIGN KEY (family_id) REFERENCES families (family_id)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS families (
        family_id         TEXT PRIMARY KEY,
        representative_key TEXT,
        size              INTEGER,
        lineage_width     INTEGER,
        fold_status       TEXT,
        plddt             REAL,
        structure_path    TEXT,
        novel_score       REAL,
        priority_score    REAL
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS functions (
        id          INTEGER PRIMARY KEY AUTOINCREMENT,
        global_key  TEXT NOT NULL,
        term_type   TEXT NOT NULL,
        term_id     TEXT NOT NULL,
        term_name   TEXT,
        raw_score   REAL,
        confidence  REAL,
        go_depth    INTEGER,
        rationale   TEXT,
        UNIQUE (global_key, term_type, term_id),
        FOREIGN KEY (global_key) REFERENCES proteins (global_key)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS evidence (
        id             INTEGER PRIMARY KEY AUTOINCREMENT,
        global_key     TEXT NOT NULL,
        term_type      TEXT,
        term_id        TEXT,
        channel        TEXT NOT NULL,
        channel_group  TEXT,
        target_id      TEXT,
        target_label   TEXT,
        score          REAL,
        evalue         REAL,
        coverage       REAL,
        details_json   TEXT,
        FOREIGN KEY (global_key) REFERENCES proteins (global_key)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS structures (
        family_id      TEXT PRIMARY KEY,
        engine         TEXT,
        fold_status    TEXT,
        plddt          REAL,
        structure_path TEXT,
        FOREIGN KEY (family_id) REFERENCES families (family_id)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS genomic_context (
        id                INTEGER PRIMARY KEY AUTOINCREMENT,
        global_key        TEXT NOT NULL,
        sample_id         TEXT,
        locus_tag         TEXT,
        neighbour_locus   TEXT,
        neighbour_product TEXT,
        neighbour_status  TEXT,
        position          TEXT,
        intergenic_bp     INTEGER,
        strand_match      TEXT,
        coupled           INTEGER,
        FOREIGN KEY (global_key) REFERENCES proteins (global_key)
    )
    """,
]

CREATE_INDEXES: List[str] = [
    "CREATE INDEX IF NOT EXISTS idx_proteins_sample ON proteins (sample_id)",
    "CREATE INDEX IF NOT EXISTS idx_proteins_family ON proteins (family_id)",
    "CREATE INDEX IF NOT EXISTS idx_proteins_status ON proteins (status)",
    "CREATE INDEX IF NOT EXISTS idx_proteins_species ON proteins (species)",
    "CREATE INDEX IF NOT EXISTS idx_functions_term ON functions (term_type, term_id)",
    "CREATE INDEX IF NOT EXISTS idx_functions_key ON functions (global_key)",
    "CREATE INDEX IF NOT EXISTS idx_evidence_key ON evidence (global_key)",
    "CREATE INDEX IF NOT EXISTS idx_evidence_channel ON evidence (channel)",
    "CREATE INDEX IF NOT EXISTS idx_context_key ON genomic_context (global_key)",
]


def create_schema(connection: sqlite3.Connection) -> None:
    """Create every table and index on an open connection."""
    for statement in CREATE_TABLES:
        connection.execute(statement)
    for statement in CREATE_INDEXES:
        connection.execute(statement)
    connection.commit()


def connect(path: Path) -> sqlite3.Connection:
    """Open a database at path and ensure the schema exists."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(str(path))
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    create_schema(connection)
    return connection
