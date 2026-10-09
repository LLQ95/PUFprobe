"""HHblits profile search against an HH-suite3 database.

Runs ``hhblits`` per representative sequence into ``dbdir/hhblits_db`` and
parses the ``.hhr`` result tables into ``Evidence(REMOTE_HOMOLOGY)`` with
score = Prob/100. A missing binary or database disables the channel (logged).
The optional ``tmvec`` channel is detected but skipped when absent.
"""

from __future__ import annotations

import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from pufprobe import externals
from pufprobe.io_utils import read_fasta, write_fasta
from pufprobe.logging_utils import get_logger
from pufprobe.models import Evidence, EvidenceChannel

logger = get_logger(__name__)

INSTALL_HINT = "See docs/server_setup.md for building an HH-suite3 database (e.g. BFD/MGnify)."
_DB_NAME = "hhblits_db"

# hhr result-table row: leading rank, target name, Prob%, E-value, P-value, Score...
_HIT_ROW = re.compile(
    r"^\s*(?P<rank>\d+)\s+(?P<name>\S.*?)\s{2,}(?P<prob>\d+(?:\.\d+)?)\s+"
    r"(?P<evalue>\S+)\s+(?P<pvalue>\S+)\s+(?P<score>[\d.]+)"
)


def run_hhblits_search(
    reps_faa: Path,
    dbdir: Optional[Path],
    out_dir: Path,
    iterations: int = 2,
    threads: int = 1,
) -> Tuple[Dict[str, List[Evidence]], Optional[Path]]:
    info = externals.hhblits()
    if not info.available:
        logger.warning("hhblits not found; skipping hhblits channel.")
        return {}, None

    db_path = dbdir / _DB_NAME if dbdir else None
    if db_path is None or not Path(str(db_path) + ".hhr").exists() \
            and not Path(str(db_path) + "_hhm.ffdata").exists():
        logger.warning("No hhblits database at dbdir/%s; skipping hhblits channel. %s",
                       _DB_NAME, INSTALL_HINT)
        return {}, None

    exe = info.path or "hhblits"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="hhblits_", dir=str(out_dir)))
    raw_tsv = out_dir / "hhblits_raw.tsv"

    evidence: Dict[str, List[Evidence]] = {}
    rows: List[str] = ["query\ttarget\tprob\tevalue\tscore"]
    for header, seq in read_fasta(reps_faa):
        query_key = header.split()[0]
        q_fasta = tmp / f"{query_key.replace('@', '_')}.faa"
        out_hhr = tmp / f"{query_key.replace('@', '_')}.hhr"
        write_fasta(q_fasta, [(query_key, seq)])
        cmd = [exe, "-i", str(q_fasta), "-d", str(db_path),
               "-n", str(iterations), "-o", str(out_hhr),
               "-cpu", str(threads)]
        try:
            completed = subprocess.run(cmd, capture_output=True, text=True)
            if completed.returncode != 0 or not out_hhr.exists():
                logger.warning("hhblits failed for %s (%s); skipping query.",
                               query_key, completed.stderr.strip()[:200])
                continue
        except (subprocess.SubprocessError, OSError) as exc:
            logger.warning("hhblits error for %s (%s); skipping query.", query_key, exc)
            continue

        for name, prob, evalue, score in parse_hhr(out_hhr):
            rows.append(f"{query_key}\t{name}\t{prob:.2f}\t{evalue:.2e}\t{score:.2f}")
            evidence.setdefault(query_key, []).append(Evidence(
                channel=EvidenceChannel.REMOTE_HOMOLOGY,
                target_id=name,
                target_label=name,
                score=max(0.0, min(1.0, prob / 100.0)),
                evalue=evalue,
                details={"prob": prob, "score": score, "database": "hhblits"},
            ))

    raw_tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return evidence, raw_tsv


def parse_hhr(hhr_path: Path) -> List[Tuple[str, float, float, float]]:
    """Return [(target_name, prob_pct, evalue, score)] from an .hhr file."""
    hits: List[Tuple[str, float, float, float]] = []
    in_table = False
    with hhr_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("No Hit"):
                in_table = True
                continue
            if in_table and line.startswith(">"):
                break
            if not in_table:
                continue
            match = _HIT_ROW.match(line)
            if not match:
                continue
            try:
                hits.append((
                    match.group("name").strip(),
                    float(match.group("prob")),
                    float(match.group("evalue")),
                    float(match.group("score")),
                ))
            except ValueError:
                continue
    return hits
