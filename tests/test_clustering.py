"""clustering: pure-numpy pycluster (union-find) and mmseqs absence fallback."""

from pathlib import Path

import numpy as np

from pufprobe.clustering import mmseqs as mmseqs_mod
from pufprobe.clustering import pycluster
from pufprobe.externals import ToolInfo


def test_pycluster_groups_close_vectors() -> None:
    # Group A: three near-identical unit-ish vectors along x.
    # Group B: two near-identical along y.
    vectors = {
        "a1": np.array([1.0, 0.0, 0.0], dtype=np.float32),
        "a2": np.array([1.0, 0.02, 0.0], dtype=np.float32),
        "a3": np.array([1.0, -0.02, 0.0], dtype=np.float32),
        "b1": np.array([0.0, 1.0, 0.0], dtype=np.float32),
        "b2": np.array([0.0, 1.0, 0.02], dtype=np.float32),
    }
    mapping = pycluster.cluster_vectors(vectors, cosine_threshold=0.9, batch_size=2)
    assert mapping["a1"] == mapping["a2"] == mapping["a3"]
    assert mapping["b1"] == mapping["b2"]
    assert mapping["a1"] != mapping["b1"]


def test_pycluster_empty_and_singleton() -> None:
    assert pycluster.cluster_vectors({}, 0.8) == {}
    single = pycluster.cluster_vectors({"only": np.array([1.0, 0.0])}, 0.8)
    assert single == {"only": 0}


def test_mmseqs_absent_returns_none(tmp_path: Path, monkeypatch) -> None:
    from pufprobe import externals
    monkeypatch.setattr(externals, "mmseqs", lambda: ToolInfo(name="mmseqs", available=False))
    (tmp_path / "q.faa").write_text(">k1\nMAAA\n", encoding="utf-8")
    result = mmseqs_mod.run_mmseqs_cluster(
        tmp_path / "q.faa", tmp_path / "work", identity=0.3, coverage=0.8
    )
    assert result is None


def test_representative_selection_longest_tie_key() -> None:
    from pufprobe.embed.pipeline import _pick_representative
    seqs = {"k_short": "AAA", "k_longer": "MMMMM", "k_longest": "KKKKKKK"}
    assert _pick_representative(list(seqs), seqs) == "k_longest"
    # tie length -> lexicographically smaller key
    tied = {"b_long": "MMMMM", "a_long": "MMMMM"}
    assert _pick_representative(list(tied), tied) == "a_long"
