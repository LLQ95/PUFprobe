"""MMseqs2 sequence search against a reference sequence database.

Runs ``createdb`` on the representative FASTA, ``search`` against the MMseqs2
database under ``dbdir`` (convention ``dbdir/uniref`` or ``dbdir/mmseqs_db``)
and ``convertalis`` to a TSV. Parsed hits become
``Evidence(channel=REMOTE_HOMOLOGY)`` with a 0-1 score derived from the
e-value. A missing binary or database disables this channel (logged) rather
than raising.
"""

from __future__ import annotations

import math
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from pufprobe import externals
from pufprobe.logging_utils import get_logger
from pufprobe.models import Evidence, EvidenceChannel

logger = get_logger(__name__)

INSTALL_HINT = "See docs/server_setup.md for downloading the UniRef MMseqs2 database."
_DB_CANDIDATES = ("uniref90", "uniref50", "uniref100", "uniref", "swissprot", "mmseqs_db")

COLS = "query,target,pident,evalue,qstart,qend,tstart,tend,qcov,tcov,theader"


def run_mmseqs_search(
    reps_faa: Path,
    dbdir: Optional[Path],
    out_dir: Path,
    threads: int = 1,
    report_evalue: float = 1e-3,
) -> Tuple[Dict[str, List[Evidence]], Optional[Path]]:
    """Search representatives; return ``{query_key: [Evidence]}`` and raw TSV path."""
    info = externals.mmseqs()
    if not info.available:
        logger.warning("mmseqs not found; skipping sequence-search channel.")
        return {}, None

    db_path = _resolve_db(dbdir)
    if db_path is None:
        logger.warning("No sequence database under dbdir (looked for %s); skipping mmseqs channel. %s",
                       ", ".join(_DB_CANDIDATES), INSTALL_HINT)
        return {}, None

    exe = info.path or "mmseqs"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="mmseqs_search_", dir=str(out_dir)))
    qdb = tmp / "qdb"
    res = tmp / "res"
    raw_tsv = out_dir / "mmseqs_raw.tsv"

    try:
        _run([exe, "createdb", str(reps_faa), str(qdb)])
        _run([exe, "search", str(qdb), str(db_path), str(res), str(tmp / "stmp"),
              "--threads", str(threads), "-e", str(report_evalue)])
        _run([exe, "convertalis", str(qdb), str(db_path), str(res), str(raw_tsv),
              "--format-output", COLS])
    except (subprocess.SubprocessError, OSError) as exc:
        logger.warning("MMseqs2 search failed (%s); channel skipped.", exc)
        return {}, None

    evidence = parse_convertalis(raw_tsv, report_evalue)
    return evidence, raw_tsv


def _resolve_db(dbdir: Optional[Path]) -> Optional[Path]:
    if dbdir is None:
        return None
    for name in _DB_CANDIDATES:
        candidate = dbdir / name
        # mmseqs DBs are stored as a set of files sharing a prefix.
        if candidate.with_suffix(".index").exists() or Path(str(candidate) + ".index").exists():
            return candidate
    return None


def parse_convertalis(tsv_path: Path, report_evalue: float) -> Dict[str, List[Evidence]]:
    evidence: Dict[str, List[Evidence]] = {}
    if not tsv_path.exists():
        return evidence
    with tsv_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) < 11:
                continue
            query, target, pident, evalue_s, qs, qe, ts, te, qcov, tcov, theader = row[:11]
            try:
                evalue = float(evalue_s)
                qcov_f = float(qcov)
                pident_f = float(pident)
            except ValueError:
                continue
            if evalue > report_evalue:
                continue
            score = _evalue_score(evalue)
            ev = Evidence(
                channel=EvidenceChannel.REMOTE_HOMOLOGY,
                target_id=target,
                target_label=_description(target, theader),
                score=score,
                evalue=evalue,
                coverage=qcov_f,
                details={
                    "pident": pident_f, "qcov": qcov_f, "tcov": float(tcov),
                    "qstart": qs, "qend": qe, "tstart": ts, "tend": te,
                    "database": "mmseqs",
                },
            )
            evidence.setdefault(query, []).append(ev)
    return evidence


def _evalue_score(evalue: float) -> float:
    if evalue <= 0.0:
        return 1.0
    return max(0.0, min(1.0, -math.log10(evalue) / 10.0))


def _description(target_id: str, header: str) -> str:
    """Return the free-text description part of a full header line.

    The ``theader`` column contains the identifier token followed by the
    description, e.g. ``sp|P0AG14|SOHB_ECOLI Probable protease SohB ...``.
    The first whitespace-delimited token is removed; the full header is
    returned when it contains no description.
    """
    header = header.strip()
    if not header:
        return target_id
    parts = header.split(None, 1)
    return parts[1].strip() if len(parts) == 2 else header


def _run(cmd: list) -> None:
    completed = subprocess.run(cmd, capture_output=True, text=True)
    if completed.returncode != 0:
        # Surface the program's own stderr so the cause is visible in the log;
        # CalledProcessError would only print the command and the exit code.
        detail = (completed.stderr or completed.stdout or "").strip()
        raise subprocess.SubprocessError(
            f"{cmd[0]} {cmd[1]} failed (exit {completed.returncode}): {detail}"
        )
