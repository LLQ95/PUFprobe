"""Foldseek structure search against a structural database (AFDB/PDB).

Only families with a usable model (``fold_status == folded``, or
``low_confidence`` but pLDDT at/above the gate) contribute structures. When
the foldseek binary, a usable structure or the structural database is missing
the channel is skipped (logged) rather than raised. Hits become
``Evidence(channel=FOLDSEEK)``.
"""

from __future__ import annotations

import math
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from pufprobe import externals
from pufprobe.logging_utils import get_logger
from pufprobe.models import Family, Evidence, EvidenceChannel

logger = get_logger(__name__)

INSTALL_HINT = "See docs/server_setup.md for setting up an AlphaFoldDB/PDB foldseek database."
_DB_NAME = "foldseek_db"

COLS = "query,target,alntmscore,prob,evalue,qcov,tcov"


def usable_structures(
    families: List[Family], min_plddt: float
) -> List[Family]:
    """Families whose model may be used for structural search."""
    out: List[Family] = []
    for fam in families:
        if fam.fold_status == "folded":
            out.append(fam)
        elif fam.fold_status == "low_confidence" and fam.plddt is not None \
                and fam.plddt >= min_plddt:
            out.append(fam)
    return out


def run_foldseek_search(
    families: List[Family],
    structures_dir: Path,
    dbdir: Optional[Path],
    out_dir: Path,
    min_plddt: float,
    sensitivity: float = 9.5,
    threads: int = 1,
) -> Tuple[Dict[str, List[Evidence]], Optional[Path]]:
    info = externals.foldseek()
    if not info.available:
        logger.warning("foldseek not found; skipping structure-search channel.")
        return {}, None

    usable = usable_structures(families, min_plddt)
    present = [f for f in usable if f.structure_path and Path(f.structure_path).exists()]
    if not present:
        logger.warning("No usable folded structures; skipping foldseek channel.")
        return {}, None

    db_path = dbdir / _DB_NAME if dbdir else None
    if db_path is None or not Path(str(db_path) + ".index").exists():
        logger.warning("No foldseek database at dbdir/%s; skipping foldseek channel. %s",
                       _DB_NAME, INSTALL_HINT)
        return {}, None

    # write a structure-only query dir (foldseek createdb expects a dir of PDBs)
    query_dir = out_dir / "_foldseek_query"
    query_dir.mkdir(parents=True, exist_ok=True)
    fam_to_rep = {}
    for fam in present:
        src = structures_dir / f"{fam.family_id}.pdb"
        if src.exists():
            (query_dir / f"{fam.family_id}.pdb").write_text(
                src.read_text(encoding="utf-8", errors="replace"), encoding="utf-8"
            )
            fam_to_rep[fam.family_id] = fam.representative_key

    exe = info.path or "foldseek"
    tmp = Path(tempfile.mkdtemp(prefix="foldseek_", dir=str(out_dir)))
    qdb = tmp / "qdb"
    res = tmp / "res"
    raw_tsv = out_dir / "foldseek_raw.tsv"

    try:
        _run([exe, "createdb", str(query_dir), str(qdb)])
        _run([exe, "search", str(qdb), str(db_path), str(res), str(tmp / "stmp"),
              "-s", str(sensitivity), "--threads", str(threads)])
        _run([exe, "convertalis", str(qdb), str(db_path), str(res), str(raw_tsv),
              "--format-output", COLS])
    except (subprocess.SubprocessError, OSError) as exc:
        logger.warning("Foldseek search failed (%s); channel skipped.", exc)
        return {}, None

    evidence = parse_foldseek_tsv(raw_tsv, fam_to_rep)
    return evidence, raw_tsv


def parse_foldseek_tsv(
    tsv_path: Path, fam_to_rep: Dict[str, str]
) -> Dict[str, List[Evidence]]:
    evidence: Dict[str, List[Evidence]] = {}
    if not tsv_path.exists():
        return evidence
    with tsv_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) < 7:
                continue
            query_fam, target, tmscore, prob, evalue_s, qcov, tcov = row[:7]
            rep_key = fam_to_rep.get(query_fam, query_fam)
            try:
                evalue = float(evalue_s)
                prob_f = float(prob)
                tm_f = float(tmscore)
            except ValueError:
                continue
            score = max(0.0, min(1.0, prob_f / 100.0)) if prob_f > 1.0 else prob_f
            ev = Evidence(
                channel=EvidenceChannel.FOLDSEEK,
                target_id=target,
                target_label=target,
                score=score,
                evalue=evalue,
                coverage=float(qcov),
                details={"tmscore": tm_f, "qcov": float(qcov), "tcov": float(tcov),
                         "family_id": query_fam},
            )
            evidence.setdefault(rep_key, []).append(ev)
    return evidence


def _run(cmd: list) -> None:
    completed = subprocess.run(cmd, capture_output=True, text=True)
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(
            completed.returncode, cmd, completed.stdout, completed.stderr
        )
