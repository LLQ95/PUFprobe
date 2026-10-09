"""MMseqs2 sequence clustering (``createdb`` / ``cluster`` / ``createtsv``).

When the ``mmseqs`` binary is detected on PATH this runs greedy cascaded
clustering in a throw-away temp directory and returns a unified
``{global_key: cluster_index}`` mapping. When the binary is missing the
function returns ``None`` (after logging install guidance) so the caller can
drop to the pure-numpy :mod:`pufprobe.clustering.pycluster` path.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Optional

from pufprobe import externals
from pufprobe.logging_utils import get_logger

logger = get_logger(__name__)

INSTALL_HINT = "Install MMseqs2 with `conda install -c conda-forge mmseqs2` (see https://mmseqs.com)."


def run_mmseqs_cluster(
    fasta_path: Path,
    work_dir: Path,
    identity: float,
    coverage: float,
    threads: int = 1,
) -> Optional[Dict[str, int]]:
    """Cluster sequences in ``fasta_path``; return ``{key: cluster_index}``.

    Returns ``None`` when the mmseqs binary is unavailable so the caller can
    fall back to pycluster. FASTA headers must be exactly the global keys.
    """
    info = externals.mmseqs()
    if not info.available:
        logger.warning(
            "mmseqs binary not found; skipping MMseqs2 clustering and falling "
            "back to pure-numpy clustering. %s", INSTALL_HINT,
        )
        return None

    exe = info.path or "mmseqs"
    work_dir.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="mmseqs_cluster_", dir=str(work_dir)))
    qdb = tmp / "qdb"
    cdb = tmp / "clu"
    tsv = tmp / "clu.tsv"

    cmd_base = [exe]
    try:
        _run(cmd_base + ["createdb", str(fasta_path), str(qdb)])
        _run(
            cmd_base
            + [
                "cluster",
                str(qdb),
                str(cdb),
                str(tmp / "mm_tmp"),
                "--min-seq-id",
                str(identity),
                "-c",
                str(coverage),
                "--threads",
                str(threads),
            ]
        )
        _run(cmd_base + ["createtsv", str(qdb), str(qdb), str(cdb), str(tsv)])
    except (subprocess.SubprocessError, OSError) as exc:
        logger.warning("MMseqs2 clustering failed (%s); falling back to pycluster.", exc)
        return None

    return parse_clusters_tsv(tsv)


def parse_clusters_tsv(tsv_path: Path) -> Dict[str, int]:
    """Parse ``createtsv`` output (rep, member per line) into key->index."""
    mapping: Dict[str, int] = {}
    cluster_order: Dict[str, int] = {}
    next_index = 0
    if not tsv_path.exists():
        return mapping
    with tsv_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            if rep not in cluster_order:
                cluster_order[rep] = next_index
                next_index += 1
            mapping[member] = cluster_order[rep]
    return mapping


def _run(cmd: list) -> None:
    completed = subprocess.run(cmd, capture_output=True, text=True, check=True)
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(
            completed.returncode, cmd, completed.stdout, completed.stderr
        )
