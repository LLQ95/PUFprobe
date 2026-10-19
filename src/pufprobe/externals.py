"""Availability detection for external binaries and optional Python modules.

Every heavy dependency is imported or located lazily through this module. The
pipeline follows a detect-warn-skip policy: a missing tool disables one
evidence channel but never crashes the run.
"""

import importlib.util
import os
import shutil
import subprocess
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional


@dataclass
class ToolInfo:
    name: str
    available: bool
    path: Optional[str] = None
    version: Optional[str] = None


def _extra_search_dirs() -> list[str]:
    """Common bin directories that hold conda/mamba-installed tools.

    These cover the case where pufprobe runs inside a conda environment but the
    binary was installed into an environment that is not currently activated
    on PATH (operators typically invoke pufprobe with ``conda run -n``).
    """
    dirs: list[str] = []
    for env_var in ("CONDA_PREFIX", "MAMBA_ROOT_PREFIX"):
        prefix = os.environ.get(env_var)
        if prefix:
            dirs.append(str(Path(prefix) / "bin"))
    return dirs


def _which(name: str) -> Optional[str]:
    path = shutil.which(name)
    if path is not None:
        return path
    for directory in _extra_search_dirs():
        candidate = Path(directory) / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return None


def detect_tool(name: str, version_args: tuple[str, ...] = ("--version",)) -> ToolInfo:
    path = _which(name)
    if path is None:
        return ToolInfo(name=name, available=False)
    version: Optional[str] = None
    try:
        completed = subprocess.run(
            [path, *version_args],
            capture_output=True,
            text=True,
            timeout=30,
        )
        raw = (completed.stdout or completed.stderr).strip().splitlines()
        if raw:
            version = raw[0].strip()
    except (subprocess.SubprocessError, OSError):
        version = None
    return ToolInfo(name=name, available=True, path=path, version=version)


@lru_cache(maxsize=None)
def mmseqs() -> ToolInfo:
    return detect_tool("mmseqs")


@lru_cache(maxsize=None)
def foldseek() -> ToolInfo:
    return detect_tool("foldseek")


@lru_cache(maxsize=None)
def hhblits() -> ToolInfo:
    return detect_tool("hhblits")


@lru_cache(maxsize=None)
def hhsearch() -> ToolInfo:
    return detect_tool("hhsearch")


@lru_cache(maxsize=None)
def p2rank() -> ToolInfo:
    # P2Rank is not distributed on bioconda; the setup script installs the
    # portable GitHub release under thirdparty/p2rank (launcher: distro/prank).
    # Resolution order: P2RANK_HOME, repo-local thirdparty, then PATH.
    candidates = []
    env_home = os.environ.get("P2RANK_HOME")
    if env_home:
        candidates.append(os.path.join(env_home, "distro", "prank"))
    repo_root = Path(__file__).resolve().parents[2]
    candidates.append(str(repo_root / "thirdparty" / "p2rank" / "distro" / "prank"))
    for candidate in candidates:
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return ToolInfo(name="p2rank", available=True, path=candidate)
    path = shutil.which("p2rank") or shutil.which("prank")
    if path is None:
        return ToolInfo(name="p2rank", available=False)
    return ToolInfo(name="p2rank", available=True, path=path)


@lru_cache(maxsize=None)
def tmvec() -> ToolInfo:
    return detect_tool("tmvec", version_args=("--help",))


def has_module(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


def torch_info() -> dict:
    """Return torch/CUDA availability without importing torch at module load."""
    if not has_module("torch"):
        return {"torch": False, "cuda": False, "cuda_devices": 0, "version": None}
    import torch  # type: ignore

    return {
        "torch": True,
        "cuda": bool(torch.cuda.is_available()),
        "cuda_devices": torch.cuda.device_count() if torch.cuda.is_available() else 0,
        "version": torch.__version__,
    }


def require_tool(info: ToolInfo, install_hint: str) -> str:
    """Return the executable path or raise an informative RuntimeError."""
    if not info.available:
        raise RuntimeError(
            f"Required tool '{info.name}' was not found on PATH. {install_hint}"
        )
    return info.path or info.name


def summary() -> dict:
    """Plain JSON snapshot used by run metadata and ``pufprobe doctor``."""
    tools = {
        "mmseqs": mmseqs().__dict__,
        "foldseek": foldseek().__dict__,
        "hhblits": hhblits().__dict__,
        "hhsearch": hhsearch().__dict__,
        "p2rank": p2rank().__dict__,
    }
    return {
        "tools": tools,
        "python_modules": {
            "torch": has_module("torch"),
            "sklearn": has_module("sklearn"),
            "igraph": has_module("igraph"),
            "leidenalg": has_module("leidenalg"),
        },
        "torch": torch_info(),
    }
