"""Tests for the learning-to-rank fusion framework."""

from __future__ import annotations

import numpy as np
import pytest

from pufprobe.models import (
    CDSRecord,
    Evidence,
    EvidenceChannel,
    GeneStatus,
)
from pufprobe.rank import ltr
from pufprobe.rank.fuse import fuse_proteins


def _evidence(channel: str, tid: str, score: float,
              evalue: float | None = None, coverage: float | None = None) -> Evidence:
    return Evidence(
        channel=EvidenceChannel(channel),
        target_id=tid,
        target_label=tid,
        score=score,
        evalue=evalue,
        coverage=coverage,
        details={"db": "pfam"},
    )


def _synthetic_problem(n_proteins: int = 12):
    """True term PF00001 backed by two strong groups; decoy PF00002 weak.

    Returns evidence_map, labels and the ordered candidate term ids.
    """
    evidence_map = {}
    labels = {}
    for i in range(n_proteins):
        key = f"sp@g{i}"
        evidence_map[key] = [
            _evidence("remote_homology", "PF00001", 0.95, 1e-60, 0.92),
            _evidence("eggnog", "PF00001", 0.95),
            _evidence("remote_homology", "PF00002", 0.15),
        ]
        labels[key] = {"Pfam": ["PF00001"]}
    return evidence_map, labels


def test_feature_vector_layout():
    gevs = [
        _evidence("remote_homology", "PF00001", 0.9, 1e-50, 0.8),
        _evidence("eggnog", "PF00001", 0.8),
    ]
    x = ltr.candidate_feature_vector(gevs, ltr.DEFAULT_RELIABILITY, "PF00001")
    assert x.shape == (len(ltr.FEATURE_NAMES),)
    # Raw channel maxima are exposed in the first channel block.
    i_remote = ltr.FEATURE_NAMES.index("remote_homology")
    assert x[i_remote] == pytest.approx(0.9)
    i_groups = ltr.FEATURE_NAMES.index("n_groups")
    # remote_homology (sequence_profile) and eggnog (database) are two groups.
    assert x[i_groups] == 2.0
    i_strong = ltr.FEATURE_NAMES.index("n_strong_direct_groups")
    assert x[i_strong] == 2.0
    i_ev = ltr.FEATURE_NAMES.index("evalue_signal")
    assert x[i_ev] > 0
    i_cov = ltr.FEATURE_NAMES.index("mean_coverage")
    assert x[i_cov] == pytest.approx(0.8)


def test_build_dataset_labels_and_groups():
    evidence_map, labels = _synthetic_problem(6)
    rows, groups = ltr.build_dataset(evidence_map, labels)
    # Two candidate terms per protein.
    assert groups == [2] * 6
    assert len(rows) == 12
    by_key = {}
    for r in rows:
        by_key.setdefault(r.protein_key, []).append(r)
    for key, candidates in by_key.items():
        truth = {r.term_id: r.y for r in candidates}
        assert truth["PF00001"] == 1
        assert truth["PF00002"] == 0


def test_temporal_split_is_disjoint_and_date_ordered():
    keys = [f"p{i}" for i in range(20)]
    dates = {f"p{i}": f"2020-{1 + i % 12:02d}-01" for i in range(20)}
    train, val, test = ltr.temporal_group_split(keys, dates=dates,
                                                test_frac=0.2, val_frac=0.1)
    assert not (set(train) & set(test))
    assert not (set(train) & set(val))
    assert not (set(val) & set(test))
    # The most recent proteins must be in the temporal test set.
    assert train and test
    assert max(dates[k] for k in test) >= max(dates[k] for k in train)
    # A seeded split without dates is still group-disjoint.
    t2, v2, e2 = ltr.temporal_group_split(keys, dates=None, seed=1)
    assert not (set(t2) & set(e2))


def test_logreg_ranks_true_term_above_decoy():
    evidence_map, labels = _synthetic_problem(12)
    rows, groups = ltr.build_dataset(evidence_map, labels)
    X, y = ltr.rows_to_matrix(rows)
    ranker = ltr.fit_ranker(X, y, groups, backend="logreg")
    assert ranker.backend == "logreg"
    scores = ranker.score(X)
    for protein in sorted(set(r.protein_key for r in rows)):
        idx = [k for k, r in enumerate(rows) if r.protein_key == protein]
        pair = {rows[k].term_id: scores[k] for k in idx}
        assert pair["PF00001"] > pair["PF00002"]


def test_fmax_aupr_smin_perfect_separation():
    evidence_map, labels = _synthetic_problem(12)
    rows, groups = ltr.build_dataset(evidence_map, labels)
    X, y = ltr.rows_to_matrix(rows)
    ranker = ltr.fit_ranker(X, y, groups, backend="logreg")
    scores = ranker.score(X)
    f_best, _, precision, recall = ltr.fmax(rows, scores)
    assert f_best >= 0.99
    assert precision >= 0.99 and recall >= 0.99
    assert ltr.micro_aupr(rows, scores) >= 0.99
    ic = {"PF00001": 6.0, "PF00002": 5.0}
    s_best, _ = ltr.smin(rows, scores, ic)
    assert s_best == pytest.approx(0.0, abs=1e-6)


def test_apply_term_stats_fills_prior_columns():
    evidence_map, labels = _synthetic_problem(3)
    rows, _ = ltr.build_dataset(evidence_map, labels)
    i_family = ltr.FEATURE_NAMES.index("family_size_norm")
    assert all(r.x[i_family] == 0.0 for r in rows)
    stats = {"PF00001": {"family_size": 3.0, "lineage_breadth": 2.0}}
    ltr.apply_term_stats(rows, stats)
    for r in rows:
        if r.term_id == "PF00001":
            assert r.x[i_family] > 0
        else:
            assert r.x[i_family] == 0.0


def test_save_load_round_trip(tmp_path):
    evidence_map, labels = _synthetic_problem(6)
    rows, groups = ltr.build_dataset(evidence_map, labels)
    X, y = ltr.rows_to_matrix(rows)
    ranker = ltr.fit_ranker(X, y, groups, backend="logreg")
    ranker.term_stats = {"PF00001": {"family_size": 6.0, "lineage_breadth": 6.0}}
    path = tmp_path / "ltr.pkl"
    ranker.save(str(path))
    reloaded = ltr.LTRRanker.load(str(path))
    assert reloaded.backend == "logreg"
    assert reloaded.term_stats["PF00001"]["family_size"] == 6.0
    np.testing.assert_allclose(reloaded.score(X), ranker.score(X), rtol=1e-9)


def _cds(locus: str) -> CDSRecord:
    return CDSRecord(
        sample_id="sp",
        locus_tag=locus,
        protein_id=None,
        sequence_aa="M" * 120,
        start=1,
        end=360,
        strand="+",
        product="hypothetical protein",
        status=GeneStatus.HYPOTHETICAL,
    )


def _load_trainset_module():
    import importlib.util
    from pathlib import Path

    script = Path(__file__).resolve().parents[1] / "server" / "build_ltr_trainset.py"
    spec = importlib.util.spec_from_file_location("build_ltr_trainset", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_trainset_parse_accession():
    mod = _load_trainset_module()
    assert mod.parse_accession("sp|P12345|NAME some desc") == "P12345"
    assert mod.parse_accession("tr|Q00000|NAME") == "Q00000"
    assert mod.parse_accession("P99999 a plain header") == "P99999"


def test_trainset_remove_self_hits():
    mod = _load_trainset_module()
    evidence = {
        "holdout@P12345": [
            {"target_id": "sp|P12345|NAME", "score": 1.0},
            {"target_id": "sp|P99999|OTHER", "score": 0.5},
            {"target_id": "PF00001", "score": 0.8},
        ]
    }
    cleaned, removed = mod.remove_self_hits(evidence)
    assert removed == 1
    targets = [r["target_id"] for r in cleaned["holdout@P12345"]]
    assert targets == ["sp|P99999|OTHER", "PF00001"]


def test_trainset_build_labels_and_dates():
    mod = _load_trainset_module()
    proteins = [("P12345", "M" * 10)]
    xref = {"P12345": {"GO": ["GO:1"], "EC": ["1.1.1.1"],
                       "Pfam": ["PF00001"], "COG": ["COG0001"]}}
    labels = mod.build_labels(proteins, xref)
    assert labels["holdout@P12345"]["Pfam"] == ["PF00001"]
    dates = mod.align_dates(proteins, {"P12345": "2020-01-01"})
    assert dates["holdout@P12345"] == "2020-01-01"


def test_fuse_proteins_uses_ltr_ranker():
    evidence_map, labels = _synthetic_problem(4)
    rows, groups = ltr.build_dataset(evidence_map, labels)
    X, y = ltr.rows_to_matrix(rows)
    ranker = ltr.fit_ranker(X, y, groups, backend="logreg")
    records = [_cds(f"g{i}") for i in range(4)]
    config = {"rank": {"channel_reliability": ltr.DEFAULT_RELIABILITY,
                       "require_independent_groups": 2,
                       "fusion_method": "group_noisy_or",
                       "single_group_cap": 0.6,
                       "group_min_strength": 0.4,
                       "high_confidence": 0.8,
                       "rule_weights": {}}}
    predictions = fuse_proteins(records, evidence_map, config, ltr_model=ranker)
    for pred in predictions:
        assert pred.terms, "each protein should have candidates"
        assert pred.terms[0].term_id == "PF00001"
        assert "LTR" in pred.terms[0].rationale
        assert pred.terms[0].raw_score > pred.terms[-1].raw_score
