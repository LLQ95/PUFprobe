"""Working-directory layout, checkpointing and run metadata.

Each stage writes under a numbered directory and records a ``_DONE`` marker on
success, so re-running a command resumes rather than repeating completed work.
"""

import json
import platform
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from pufprobe import __version__
from pufprobe.io_utils import read_json, write_json
from pufprobe.logging_utils import setup_logging

STAGE_DIRS: Dict[str, str] = {
    "prepare": "00_prepare",
    "qc": "01_qc",
    "embed": "02_embed",
    "cluster": "03_cluster",
    "fold": "04_fold",
    "search": "05_search",
    "context": "06_context",
    "rank": "07_rank",
    "report": "08_report",
    "annotate": "09_annotate",
}


class WorkDir:
    """Filesystem state for one PUFprobe run."""

    def __init__(self, root: Path, cpu: bool = False, threads: int = 1, log_level: str = "INFO"):
        self.root = Path(root).resolve()
        self.cpu = cpu
        self.threads = threads
        self.log_level = log_level
        self.logs_dir = self.root / "logs"
        for dirname in STAGE_DIRS.values():
            (self.root / dirname).mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        setup_logging(self.logs_dir / "pufprobe.log", level=log_level)
        self.metadata_path = self.root / "metadata.json"

    # ------------------------------------------------------------------
    # Paths
    # ------------------------------------------------------------------
    def stage_dir(self, stage: str) -> Path:
        if stage not in STAGE_DIRS:
            raise KeyError(f"Unknown stage: {stage}")
        path = self.root / STAGE_DIRS[stage]
        path.mkdir(parents=True, exist_ok=True)
        return path

    def path(self, stage: str, name: str) -> Path:
        return self.stage_dir(stage) / name

    # ------------------------------------------------------------------
    # Checkpointing
    # ------------------------------------------------------------------
    def marker(self, stage: str) -> Path:
        return self.stage_dir(stage) / "_DONE"

    def is_complete(self, stage: str) -> bool:
        return self.marker(stage).is_file()

    def mark_complete(self, stage: str, payload: Optional[dict] = None) -> None:
        marker = self.marker(stage)
        marker.write_text(
            json.dumps(
                {"stage": stage, "finished": datetime.now().isoformat(timespec="seconds"), **(payload or {})},
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )

    # ------------------------------------------------------------------
    # Metadata
    # ------------------------------------------------------------------
    def write_metadata(self, extra: Optional[Dict[str, Any]] = None) -> None:
        payload: Dict[str, Any] = {
            "pufprobe_version": __version__,
            "created": datetime.now().isoformat(timespec="seconds"),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "cpu_mode": self.cpu,
            "threads": self.threads,
        }
        if extra:
            payload.update(extra)
        write_json(self.metadata_path, payload)

    def update_metadata(self, values: Dict[str, Any]) -> None:
        current: Dict[str, Any] = {}
        if self.metadata_path.exists():
            current = read_json(self.metadata_path)  # type: ignore[assignment]
        current.update(values)
        write_json(self.metadata_path, current)
