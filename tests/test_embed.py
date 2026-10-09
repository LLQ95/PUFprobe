"""embed stage: fallback vector shape/manifest, candidate filtering, bookkeeping."""

import json
from pathlib import Path

import numpy as np

from pufprobe.config import load_config
from pufprobe.embed.fallback import compute_fallback_embeddings, fallback_vector
from pufprobe.embed.pipeline import run as embed_run
from pufprobe.io_utils import read_json
from pufprobe.models import CDSRecord, GeneStatus, GeneVerdict, to_dict
from pufprobe.workdir import WorkDir


def _rec(sample: str, locus: str, seq: str,
         status: str = GeneStatus.HYPOTHETICAL.value,
         verdict=None) -> CDSRecord:
    return CDSRecord(
        sample_id=sample, locus_tag=locus, protein_id=None, sequence_aa=seq,
        start=1, end=len(seq) * 3, strand="+", product="hypothetical protein",
        status=GeneStatus(status), verdict=verdict,
    )


def _write_records(wd: WorkDir, records) -> None:
    rows = [to_dict(r) for r in records]
    wd.path("prepare", "cds_records.json").write_text(
        json.dumps(rows), encoding="utf-8"
    )


def test_fallback_vector_pads_to_config_dim() -> None:
    config = load_config(None)
    dim = config["embed"]["fallback_feature_dim"]  # 48
    vec, filled = fallback_vector("MKMKMKMKGP", dim)
    assert vec.shape == (dim,)
    assert filled == 30  # 20 composition + 10 physicochemical
    # M is index 10 in AA_ORDER="ACDEFGHIKLMNPQRSTVWY" -> nonzero there
    assert vec[10] > 0
    assert np.all(vec[filled:] == 0)


def test_embed_pipeline_fallback_end_to_end(tmp_path: Path) -> None:
    wd = WorkDir(tmp_path / "run")
    config = load_config(None)
    dim = config["embed"]["fallback_feature_dim"]

    # 4 candidates: two near-identical hydrophobic (cluster together), one
    # distinct, one annotated (excluded), one likely_false_orf (excluded).
    seq_h1 = "M" * 40 + "AKL"
    seq_h2 = "M" * 40 + "AKM"
    seq_diff = "N" * 40 + "D"
    records = [
        _rec("sA", "g1", seq_h1),
        _rec("sA", "g2", seq_h2),
        _rec("sB", "g3", seq_diff),
        _rec("sB", "g4", seq_h1, status=GeneStatus.ANNOTATED.value),
        _rec("sB", "g5", seq_h1, verdict=GeneVerdict.LIKELY_FALSE_ORF),
    ]
    _write_records(wd, records)

    embed_run(wd, config, fallback_only=True)

    # manifest
    manifest = read_json(wd.path("embed", "manifest.json"))  # type: ignore[arg-type]
    assert manifest["source"] == "fallback"
    assert manifest["dim"] == dim
    assert manifest["n_records"] == 3  # annotated + false-orf excluded
    assert manifest["padded"] is True  # 30 -> 48

    # embeddings npz round-trip
    npz = np.load(wd.path("embed", "embeddings.npz"), allow_pickle=True)
    assert set(npz.files) == {"sA@g1", "sA@g2", "sB@g3"}
    assert npz["sA@g1"].shape == (dim,)

    # families: g1/g2 cosine-similar, g3 separate
    families = read_json(wd.path("cluster", "families.json"))  # type: ignore[arg-type]
    assert len(families) == 2
    by_id = {f["family_id"]: f for f in families}
    assert by_id["FAM00001"]["size"] == 2 or by_id["FAM00002"]["size"] == 2
    pair_fam = by_id["FAM00001"] if by_id["FAM00001"]["size"] == 2 else by_id["FAM00002"]
    assert set(pair_fam["members"]) == {"sA@g1", "sA@g2"}
    # lineage width = distinct samples within cluster
    assert pair_fam["lineage_width"] == 1
    # representative = longest, ties by key
    assert pair_fam["representative_key"] in {"sA@g1", "sA@g2"}

    # cds_records.json back-wrote only family_id / is_representative
    rows = read_json(wd.path("prepare", "cds_records.json"))  # type: ignore[arg-type]
    by_key = {f"{r['sample_id']}@{r['locus_tag']}": r for r in rows}
    assert by_key["sA@g1"]["family_id"] == pair_fam["family_id"]
    assert by_key["sA@g2"]["family_id"] == pair_fam["family_id"]
    assert by_key["sB@g3"]["family_id"] != pair_fam["family_id"]
    # annotated + false-orf records got no family
    assert by_key["sB@g4"]["family_id"] is None
    assert by_key["sB@g5"]["family_id"] is None
    # exactly one representative per family
    reps = [k for k, r in by_key.items() if r["is_representative"]]
    assert len(reps) == 2

    # completion markers
    assert wd.marker("embed").is_file()
    assert wd.marker("cluster").is_file()


def test_esm_unavailable_triggers_fallback(tmp_path: Path, monkeypatch) -> None:
    wd = WorkDir(tmp_path / "run")
    config = load_config(None)
    _write_records(wd, [_rec("sA", "g1", "M" * 30)])

    # Force the ESM path to raise ModelUnavailable; pipeline must fall back.
    from pufprobe.embed import pipeline as embed_pkg
    from pufprobe.embed.esm import ModelUnavailable

    def _boom(*args, **kwargs):
        raise ModelUnavailable("no torch here")

    monkeypatch.setattr(embed_pkg, "compute_embeddings", _boom)
    embed_run(wd, config)  # no fallback_only, but ESM blows up -> fallback
    manifest = read_json(wd.path("embed", "manifest.json"))  # type: ignore[arg-type]
    assert manifest["source"] == "fallback"
    assert "esm_error" in manifest
