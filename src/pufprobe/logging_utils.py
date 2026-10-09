"""Logging configuration shared by every subcommand."""

import logging
from pathlib import Path

_LOG_FORMAT = "%(asctime)s %(levelname)-7s %(name)s | %(message)s"
_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logging(log_file: Path | None = None, level: str = "INFO") -> None:
    """Configure root logging once; optionally mirror records to a file."""
    root = logging.getLogger("pufprobe")
    root.setLevel(getattr(logging, level.upper(), logging.INFO))
    root.handlers.clear()
    formatter = logging.Formatter(_LOG_FORMAT, datefmt=_DATE_FORMAT)
    stream = logging.StreamHandler()
    stream.setFormatter(formatter)
    root.addHandler(stream)
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setFormatter(formatter)
        root.addHandler(file_handler)
    root.propagate = False


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(f"pufprobe.{name}")
