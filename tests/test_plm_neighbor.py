"""Tests for the pLM embedding nearest-neighbour evidence channel.

A small synthetic reference library (no torch / network) is built on disk so
the loading, cosine voting and end-to-end evidence generation can be tested
deterministically.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from pufprobe.io_utils import write_fasta
from pufprobe.models import EvidenceChannel
from pufprobe.search import plm_neighbor as plm

_DIM = 8
_MODEL = "esm2_test_dummy"


def _unit(values):
    vector = np.asarray(values, dtype=np.float32)
    return vector / np.linalg.norm(vector)


# Reference proteins expressed in the first three coordinates. Query q1 is the
# e1 axis, so cosine with each reference is its e1 component.
_REF_PROTEINS = [
    ("R1", [0.99, 0.141, 0.0, 0.0, 0, 0, 0, 0], {"GO": ["GO:0000001"]}),
    ("R2", [0.80, 0.60, 0.0, 0.0, 0, 0, 0, 0], {"GO": ["GO:0000001"]}),
    ("R3", [0.00, 1.00, 0.0, 0.0, 0, 0, 0, 0], {"GO": ["GO:0000002"]}),
    ("R4", [0.60, 0.00, 0.80, 0.0, 0, 0, 0, 0], {"Pfam": ["PF00001"]}),
    ("R5", [0.70, 0.00, 0.714, 0.0, 0, 0, 0, 0], {"EC": ["1.1.1.1"]}),
]


def _write_reference(ref_dir: Path, model: str = _MODEL):
    ref_dir.mkdir(parents=True, exist_ok=True)
    ids = [pid for pid, _, _ in _REF_PROTEINS]
    vectors = np.stack([_unit(vec) for _, vec, _ in _REF_PROTEINS]).astype(np.float16)
    terms = {pid: term for pid, _, term in _REF_PROTEINS}
    np.save(ref_dir / plm.REF_VECTORS, vectors)
    (ref_dir / plm.REF_IDS).write_text(json.dumps(ids), encoding="utf-8")
    (ref_dir / plm.REF_TERMS).write_text(json.dumps(terms), encoding="utf-8")
    (ref_dir / plm.REF_MANIFEST).write_text(
        json.dumps({"model": model, "dim": _DIM, "n": len(ids), "normalized": True}),
        encoding="utf-8",
    )


def _write_query(run_dir: Path, model: str = _MODEL, vector=None):
    q = _unit(vector if vector is not None else [1, 0, 0, 0, 0, 0, 0, 0])
    reps = run_dir / "representatives.faa"
    write_fasta(reps, [("q1", "M" * 40)])
    emb_dir = run_dir / "embed"
    emb_dir.mkdir(parents=True, exist_ok=True)
    np.savez(emb_dir / "embeddings.npz", q1=q)
    (emb_dir / "manifest.json").write_text(
        json.dumps({"model": model, "source": "esm", "dim": _DIM}), encoding="utf-8"
    )
    return reps, emb_dir / "embeddings.npz", emb_dir / "manifest.json"


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def test_load_reference_roundtrip_and_normalises(tmp_path: Path):
    ref_dir = tmp_path / "plm_ref"
    _write_reference(ref_dir)
    library = plm.load_reference(ref_dir)
    assert library is not None
    assert library.vectors.shape == (5, _DIM)
    assert library.model == _MODEL
    # Stored float16 vectors are re-normalised on load.
    norms = np.linalg.norm(library.vectors, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-3)


def test_load_reference_missing_returns_none(tmp_path: Path):
    assert plm.load_reference(None) is None
    assert plm.load_reference(tmp_path / "absent") is None


# ---------------------------------------------------------------------------
# voting
# ---------------------------------------------------------------------------
def _neighbors(library):
    q = _unit([1, 0, 0, 0, 0, 0, 0, 0])
    hits = plm._batch_topk(library.vectors, q[None, :], ["q1"], k=5)
    return hits["q1"]


def test_vote_single_strong_neighbor(tmp_path: Path):
    library = plm.load_reference(_write_reference_and_return(tmp_path))
    idx, cos = _neighbors(library)
    # R1 alone: cos 0.99 -> weight 0.98, require that one voter.
    cands = plm.vote_terms(
        idx[:1], cos[:1], library.ids, library.terms,
        cos_floor=0.5, cos_min=0.5, min_term_score=0.5, min_voters=1,
    )
    by_term = {c["term_id"]: c for c in cands}
    assert "GO:0000001" in by_term
    assert by_term["GO:0000001"]["score"] == pytest.approx(0.98, abs=0.02)
    assert by_term["GO:0000001"]["n_voters"] == 1


def test_vote_agreement_accumulates(tmp_path: Path):
    library = plm.load_reference(_write_reference_and_return(tmp_path))
    idx, cos = _neighbors(library)
    single = plm.vote_terms(
        idx[:1], cos[:1], library.ids, library.terms, min_term_score=0.0,
    )
    both = plm.vote_terms(
        idx[:2], cos[:2], library.ids, library.terms, min_term_score=0.0,
    )
    single_score = {c["term_id"]: c["score"] for c in single}["GO:0000001"]
    both_score = {c["term_id"]: c["score"] for c in both}["GO:0000001"]
    assert both_score > single_score
    assert both_score == pytest.approx(0.992, abs=0.01)


def test_vote_filters_orthogonal_and_weak(tmp_path: Path):
    library = plm.load_reference(_write_reference_and_return(tmp_path))
    idx, cos = _neighbors(library)
    cands = plm.vote_terms(
        idx, cos, library.ids, library.terms,
        cos_floor=0.5, cos_min=0.5, min_term_score=0.5, min_voters=1,
    )
    terms = {c["term_id"] for c in cands}
    # R3 is orthogonal (cos 0) -> GO:0000002 absent; R4/R5 weight too low.
    assert "GO:0000002" not in terms
    assert "PF00001" not in terms
    assert "1.1.1.1" not in terms
    assert "EC:1.1.1.1" not in terms


def test_vote_min_voters_gate(tmp_path: Path):
    library = plm.load_reference(_write_reference_and_return(tmp_path))
    idx, cos = _neighbors(library)
    cands = plm.vote_terms(
        idx, cos, library.ids, library.terms,
        min_term_score=0.0, min_voters=2,
    )
    # GO:0000001 has two voters (R1,R2) and survives; every other term has one.
    assert {c["term_id"] for c in cands} == {"GO:0000001"}


def test_vote_normalises_ec_and_lower_threshold(tmp_path: Path):
    library = plm.load_reference(_write_reference_and_return(tmp_path))
    idx, cos = _neighbors(library)
    cands = plm.vote_terms(
        idx, cos, library.ids, library.terms,
        min_term_score=0.3, min_voters=1,
    )
    by_term = {c["term_id"]: c for c in cands}
    # R5 carries EC at cos 0.7 -> weight 0.4; emitted as EC:1.1.1.1.
    assert "EC:1.1.1.1" in by_term
    assert by_term["EC:1.1.1.1"]["kind"] == "ec"


# ---------------------------------------------------------------------------
# end-to-end run
# ---------------------------------------------------------------------------
def test_run_end_to_end(tmp_path: Path):
    ref_dir = tmp_path / "plm_ref"
    _write_reference(ref_dir)
    run_dir = tmp_path / "run"
    reps, emb, manifest = _write_query(run_dir)
    out_dir = run_dir / "search"
    out_dir.mkdir(parents=True, exist_ok=True)

    evidence, raw = plm.run_plm_neighbor_search(
        reps, emb, None, out_dir, ref_dir=ref_dir,
        query_manifest_path=manifest,
        top_k=5, cos_floor=0.5, cos_min=0.5, min_term_score=0.5,
    )
    assert raw is not None and raw.exists()
    term_ids = [e.target_id for e in evidence["q1"]]
    assert "GO:0000001" in term_ids
    go_ev = next(e for e in evidence["q1"] if e.target_id == "GO:0000001")
    assert go_ev.channel == EvidenceChannel.PLM_NEIGHBOR
    assert go_ev.details["db"] == "go"
    assert go_ev.details["n_voters"] == 2


def test_run_model_mismatch_disables_channel(tmp_path: Path):
    ref_dir = tmp_path / "plm_ref"
    _write_reference(ref_dir, model="esm2_t33_650M_UR50D")
    run_dir = tmp_path / "run"
    reps, emb, manifest = _write_query(run_dir, model="esm2_t12_35M_UR50D")
    out_dir = run_dir / "search"
    out_dir.mkdir(parents=True, exist_ok=True)

    evidence, raw = plm.run_plm_neighbor_search(
        reps, emb, None, out_dir, ref_dir=ref_dir, query_manifest_path=manifest,
    )
    assert evidence == {}
    assert raw is None


def test_run_finds_default_library_under_dbdir(tmp_path: Path):
    dbdir = tmp_path / "dbdir"
    _write_reference(dbdir / "plm_ref")
    run_dir = tmp_path / "run"
    reps, emb, manifest = _write_query(run_dir)
    out_dir = run_dir / "search"
    out_dir.mkdir(parents=True, exist_ok=True)

    evidence, _ = plm.run_plm_neighbor_search(
        reps, emb, dbdir, out_dir, ref_dir=None, query_manifest_path=manifest,
    )
    assert "q1" in evidence


def _write_reference_and_return(tmp_path: Path) -> Path:
    ref_dir = tmp_path / "plm_ref"
    _write_reference(ref_dir)
    return ref_dir
