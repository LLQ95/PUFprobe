"""Tests for external binary detection, focusing on P2Rank layout handling."""

import stat
import sys

import pytest

from pufprobe import externals


def _make_launcher(path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture(autouse=True)
def _clear_p2rank_cache():
    externals.p2rank.cache_clear()
    yield
    externals.p2rank.cache_clear()


def test_p2rank_top_level_layout(tmp_path, monkeypatch):
    base = tmp_path / "p2rank"
    _make_launcher(base / "prank")
    (base / "bin").mkdir(parents=True)
    (base / "bin" / "p2rank.jar").write_text("jar", encoding="utf-8")
    monkeypatch.setenv("P2RANK_HOME", str(base))

    info = externals.p2rank()

    assert info.available
    assert info.path == str(base / "prank")


def test_p2rank_distro_layout(tmp_path, monkeypatch):
    base = tmp_path / "p2rank"
    _make_launcher(base / "distro" / "prank")
    (base / "distro" / "bin").mkdir(parents=True)
    (base / "distro" / "bin" / "p2rank.jar").write_text("jar", encoding="utf-8")
    monkeypatch.setenv("P2RANK_HOME", str(base))

    info = externals.p2rank()

    assert info.available
    assert info.path == str(base / "distro" / "prank")


def test_p2rank_requires_jar(tmp_path, monkeypatch):
    base = tmp_path / "p2rank"
    _make_launcher(base / "prank")  # launcher without the P2Rank jar
    monkeypatch.setenv("P2RANK_HOME", str(base))
    monkeypatch.setenv("PATH", str(tmp_path / "nobin"))

    info = externals.p2rank()

    assert not info.available


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="an extension-less launcher is not on Windows PATHEXT",
)
def test_p2rank_on_path_by_name(tmp_path, monkeypatch):
    bindir = tmp_path / "bin"
    _make_launcher(bindir / "p2rank")
    monkeypatch.delenv("P2RANK_HOME", raising=False)
    monkeypatch.setenv("PATH", str(bindir))

    info = externals.p2rank()

    assert info.available
    assert info.path == str(bindir / "p2rank")


def test_bioconda_prank_not_misclassified(tmp_path, monkeypatch):
    bindir = tmp_path / "bin"
    _make_launcher(bindir / "prank")  # bioconda alignment tool, no jar nearby
    monkeypatch.delenv("P2RANK_HOME", raising=False)
    monkeypatch.setenv("PATH", str(bindir))

    info = externals.p2rank()

    assert not info.available
