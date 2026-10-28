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
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from pufprobe import externals
from pufprobe.io_utils import read_fasta, write_fasta
from pufprobe.logging_utils import get_logger
from pufprobe.models import Evidence, EvidenceChannel

logger = get_logger(__name__)

INSTALL_HINT = "See docs/server_setup.md for building an HH-suite3 database (e.g. BFD/MGnify)."
_DB_NAME = "hhblits_db"

# hhr result-table row. The "No Hit" table is fixed-column, but the target
# name is space-padded and often runs long (UniRef30), leaving a single space
# before Prob. The numeric fields are anchored from the fixed tail (SS, Cols,
# query/template spans, template length) so a long name cannot shift
# prob/evalue/score. The table header carries one leading space; parse_hhr
# strips it before matching.
_HIT_ROW = re.compile(
    r"^\s*(?P<rank>\d+)\s+(?P<name>.*?)\s+"
    r"(?P<prob>\d{1,3}\.\d+)\s+"
    r"(?P<evalue>\S+)\s+(?P<pvalue>\S+)\s+"
    r"(?P<score>[\d.]+)\s+(?P<ss>[\d.]+)\s+(?P<cols>\d+)\s+"
    r"(?P<qrange>\d+-\d+)\s+(?P<trange>\d+-\d+)\s+\((?P<tlen>\d+)\)\s*$"
)


def run_hhblits_search(
    reps_faa: Path,
    dbdir: Optional[Path],
    out_dir: Path,
    iterations: int = 2,
    threads: int = 1,
    jobs: int = 1,
    query_timeout: int = 300,
) -> Tuple[Dict[str, List[Evidence]], Optional[Path]]:
    info = externals.hhblits()
    if not info.available:
        logger.warning("hhblits not found; skipping hhblits channel.")
        return {}, None

    db_path = dbdir / _DB_NAME if dbdir else None
    db_markers = ("_hhm.ffdata", "_a3m.ffdata", "_msa.ffdata")
    if db_path is None or not any(
        Path(str(db_path) + marker).exists() for marker in db_markers
    ):
        logger.warning(
            "No hh-suite database at dbdir/%s (looked for files %s); skipping hhblits channel. %s",
            _DB_NAME,
            ", ".join(f"hhblits_db{m}" for m in db_markers),
            INSTALL_HINT,
        )
        return {}, None

    exe = info.path or "hhblits"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="hhblits_", dir=str(out_dir)))
    raw_tsv = out_dir / "hhblits_raw.tsv"

    queries = [(header.split()[0], seq) for header, seq in read_fasta(reps_faa)]
    workers = max(1, min(int(jobs), len(queries))) if queries else 1
    # Split the CPU budget across the concurrent searches. A single query
    # scales sublinearly with -cpu, so running several independent searches
    # uses the cores more efficiently than one -cpu=$threads invocation.
    cpu_per_worker = max(1, int(threads) // workers)

    def _run_one(item: Tuple[str, str]) -> Tuple[str, List[Tuple[str, float, float, float]]]:
        query_key, seq = item
        safe = query_key.replace("@", "_")
        q_fasta = tmp / f"{safe}.faa"
        out_hhr = tmp / f"{safe}.hhr"
        write_fasta(q_fasta, [(query_key, seq)])
        cmd = [exe, "-i", str(q_fasta), "-d", str(db_path),
               "-n", str(iterations), "-o", str(out_hhr),
               "-cpu", str(cpu_per_worker)]
        try:
            completed = subprocess.run(
                cmd, capture_output=True, text=True, timeout=query_timeout
            )
            if completed.returncode != 0 or not out_hhr.exists():
                logger.warning("hhblits failed for %s (%s); skipping query.",
                               query_key, completed.stderr.strip()[:200])
                return query_key, []
        except subprocess.TimeoutExpired:
            logger.warning("hhblits timed out for %s after %ss; skipping query.",
                           query_key, query_timeout)
            if out_hhr.exists():
                out_hhr.unlink()
            return query_key, []
        except (subprocess.SubprocessError, OSError) as exc:
            logger.warning("hhblits error for %s (%s); skipping query.", query_key, exc)
            return query_key, []
        return query_key, parse_hhr(out_hhr)

    if workers > 1:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            # pool.map preserves input order, so the raw TSV stays deterministic.
            results = list(pool.map(_run_one, queries))
    else:
        results = [_run_one(item) for item in queries]

    evidence: Dict[str, List[Evidence]] = {}
    rows: List[str] = ["query\ttarget\tprob\tevalue\tscore"]
    for query_key, hits in results:
        for name, prob, evalue, score in hits:
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
            stripped = line.lstrip()
            if stripped.startswith("No Hit"):
                in_table = True
                continue
            if in_table and stripped.startswith(">"):
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
