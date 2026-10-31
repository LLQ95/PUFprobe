"""Keep the package version string aligned with the build metadata."""
import re
from pathlib import Path

import pufprobe


def test_version_matches_pyproject():
    root = Path(__file__).resolve().parents[1]
    text = (root / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    assert match is not None
    assert pufprobe.__version__ == match.group(1)
