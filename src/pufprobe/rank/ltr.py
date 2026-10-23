"""Learning-to-rank fusion of multi-channel evidence.

The rule-based fusion in :mod:`pufprobe.rank.fuse` (group noisy-OR) is a fixed
prior over how channels combine. This module learns that combination from data,
in the style of the CAFA challenges: a set of already-annotated proteins is
treated as if its function were unknown, the evidence channels are run against
the reference databases, and a ranker is trained to put the protein's true
annotations above the rest of its candidate terms.

Each training/test row is one ``(protein, candidate_term)`` pair. Rows that
belong to the same protein form one ranking group, so a LambdaRank objective
optimises the within-protein ordering directly. The feature vector combines:

* the best raw score from every evidence channel for that term,
* the reliability-discounted strength of every independent evidence group,
* cross-channel agreement and coverage statistics,
* a handful of meta features: mean hit coverage, a -log10 E-value signal, the
  term's family size and taxonomic breadth (priors learned on the training
  set and stored inside the model), and the protein's pLDDT.

Three backends are supported and selected automatically. ``lightgbm`` uses a
true pairwise LambdaRank model and is the recommended backend; ``sklearn``
falls back to a pointwise gradient-boosted classifier; ``logreg`` is a
dependency-free numpy L2 logistic ranker kept so the framework trains and
runs on a minimal install. The model (including the term priors) is stored
with pickle, so a LightGBM model needs lightgbm at inference time; the
numpy fallback ranker needs nothing beyond numpy.

CAFA-style evaluation returns protein-centric Fmax, micro-average AUPR and,
when a term information-content table is supplied, semantic distance Smin.
"""

from __future__ import annotations

import math
import pickle
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from pufprobe.logging_utils import get_logger
from pufprobe.models import CHANNEL_GROUPS, Evidence, EvidenceChannel
from pufprobe.rank.fuse import INDIRECT_GROUPS, classify_term

logger = get_logger("rank.ltr")

# Keep the same default reliabilities as the rule config so group-strength
# features are on the same scale regardless of where the matrix is built.
DEFAULT_RELIABILITY: Dict[str, float] = {
    EvidenceChannel.PLM_NEIGHBOR.value: 0.75,
    EvidenceChannel.FOLDSEEK.value: 0.85,
    EvidenceChannel.REMOTE_HOMOLOGY.value: 0.95,
    EvidenceChannel.POCKET.value: 0.60,
    EvidenceChannel.CONTEXT.value: 0.60,
    EvidenceChannel.DATABASE.value: 0.90,
    EvidenceChannel.EGGNOG.value: 0.90,
}

# Fixed feature ordering. The seven raw channels, the five group strengths,
# then the meta features. The model stores this list and rejects matrices that
# do not match.
CHANNEL_FEATURES: List[str] = [ch.value for ch in EvidenceChannel]
GROUP_FEATURES: List[str] = [
    "sequence", "structure", "sequence_profile", "context", "database",
]
META_FEATURES: List[str] = [
    "n_channels",
    "n_groups",
    "n_strong_direct_groups",
    "n_evidence",
    "neighbor_consistency",
    "evalue_signal",
    "mean_coverage",
    "family_size_norm",
    "lineage_breadth_norm",
    "protein_plddt",
]
FEATURE_NAMES: List[str] = CHANNEL_FEATURES + GROUP_FEATURES + META_FEATURES

_STRONG_GATE = 0.40
_EVALUE_SAT = 300.0
_FAMILY_SAT = 10000.0
_LINEAGE_SAT = 1000.0


@dataclass
class CandidateRow:
    """One feature row for a single candidate term of one protein."""

    protein_key: str
    term_type: str
    term_id: str
    x: np.ndarray
    y: int


# ---------------------------------------------------------------------------
# Feature construction
# ---------------------------------------------------------------------------

def _sigmoid(z: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(z, -60.0, 60.0)))


def _term_prior(term_stats: Dict[str, Dict[str, float]], term_id: str,
                key: str, sat: float) -> float:
    entry = term_stats.get(term_id)
    if not entry:
        return 0.0
    value = float(entry.get(key, 0.0))
    return min(math.log1p(max(value, 0.0)) / math.log1p(sat), 1.0)


def candidate_feature_vector(
    gevs: List[Evidence],
    reliability: Dict[str, float],
    term_id: str,
    protein_plddt: float = 0.0,
    term_stats: Optional[Dict[str, Dict[str, float]]] = None,
    strong_gate: float = _STRONG_GATE,
) -> np.ndarray:
    """Build the fixed-length feature vector for one candidate term."""
    best_per_channel: Dict[str, float] = {}
    coverages: List[float] = []
    evalue_signal = 0.0
    for ev in gevs:
        ch = ev.channel.value
        best_per_channel[ch] = max(best_per_channel.get(ch, ev.score), ev.score)
        if ev.coverage is not None:
            coverages.append(float(ev.coverage))
        if ev.evalue is not None:
            try:
                ev_f = float(ev.evalue)
                if ev_f > 0.0:
                    evalue_signal = max(evalue_signal,
                                        min(-math.log10(ev_f) / _EVALUE_SAT, 1.0))
            except (TypeError, ValueError):
                pass

    # Reliability-discounted group strengths, mirroring rule fusion.
    group_strength: Dict[str, float] = {}
    for ch, score in best_per_channel.items():
        grp = CHANNEL_GROUPS[ch]
        discounted = float(reliability.get(ch, 0.5)) * float(score)
        group_strength[grp] = max(group_strength.get(grp, 0.0), discounted)

    present_groups = {CHANNEL_GROUPS[ev.channel.value] for ev in gevs}
    direct_groups = present_groups - INDIRECT_GROUPS
    strong_direct = [g for g in direct_groups
                     if group_strength.get(g, 0.0) >= strong_gate]

    # Cross-channel agreement: geometric mean of group strengths rewards terms
    # that are supported consistently across independent witnesses.
    if group_strength:
        logs = [math.log(max(s, 1e-6)) for s in group_strength.values()]
        neighbor_consistency = math.exp(sum(logs) / len(logs))
    else:
        neighbor_consistency = 0.0

    term_stats = term_stats or {}
    x = np.zeros(len(FEATURE_NAMES), dtype=np.float64)
    for i, ch in enumerate(CHANNEL_FEATURES):
        x[i] = float(best_per_channel.get(ch, 0.0))
    base = len(CHANNEL_FEATURES)
    for j, grp in enumerate(GROUP_FEATURES):
        x[base + j] = float(group_strength.get(grp, 0.0))
    m = base + len(GROUP_FEATURES)
    meta = {
        "n_channels": float(len(best_per_channel)),
        "n_groups": float(len(present_groups)),
        "n_strong_direct_groups": float(len(strong_direct)),
        "n_evidence": float(len(gevs)),
        "neighbor_consistency": float(neighbor_consistency),
        "evalue_signal": float(evalue_signal),
        "mean_coverage": float(np.mean(coverages)) if coverages else 0.0,
        "family_size_norm": _term_prior(term_stats, term_id,
                                        "family_size", _FAMILY_SAT),
        "lineage_breadth_norm": _term_prior(term_stats, term_id,
                                            "lineage_breadth", _LINEAGE_SAT),
        "protein_plddt": float(protein_plddt) / 100.0
        if protein_plddt and protein_plddt > 1.5 else float(protein_plddt or 0.0),
    }
    for k, name in enumerate(META_FEATURES):
        x[m + k] = meta[name]
    return x


def build_dataset(
    evidence_map: Dict[str, List[Evidence]],
    labels: Dict[str, Dict[str, Sequence[str]]],
    term_stats: Optional[Dict[str, Dict[str, float]]] = None,
    protein_meta: Optional[Dict[str, Dict[str, float]]] = None,
    reliability: Optional[Dict[str, float]] = None,
) -> Tuple[List[CandidateRow], List[int]]:
    """Build ranking rows and per-protein group sizes.

    ``labels`` maps a protein key to a mapping of term type -> iterable of true
    term ids. GO labels must already be extended to true-path ancestors if the
    inference ranker propagates GO evidence. Only proteins present in both
    ``labels`` and ``evidence_map`` contribute rows.
    """
    reliability = reliability or DEFAULT_RELIABILITY
    protein_meta = protein_meta or {}
    rows: List[CandidateRow] = []
    groups: List[int] = []

    for key in sorted(labels):
        evs = evidence_map.get(key)
        if not evs:
            continue
        truth_by_type: Dict[str, set] = {}
        for ttype, ids in labels[key].items():
            truth_by_type.setdefault(str(ttype), set()).update(str(t) for t in ids)

        grouped: Dict[Tuple[str, str], List[Evidence]] = {}
        for ev in evs:
            ttype, tid = classify_term(ev)
            grouped.setdefault((ttype, tid), []).append(ev)

        meta = protein_meta.get(key, {})
        plddt = float(meta.get("plddt", 0.0) or 0.0)
        n_rows = 0
        for (ttype, tid), gevs in sorted(grouped.items()):
            x = candidate_feature_vector(
                gevs, reliability, tid,
                protein_plddt=plddt, term_stats=term_stats,
            )
            y = int(tid in truth_by_type.get(ttype, set()))
            rows.append(CandidateRow(key, ttype, tid, x, y))
            n_rows += 1
        if n_rows:
            groups.append(n_rows)
    return rows, groups


def rows_to_matrix(rows: Sequence[CandidateRow]) -> Tuple[np.ndarray, np.ndarray]:
    X = np.stack([r.x for r in rows]).astype(np.float64)
    y = np.asarray([r.y for r in rows], dtype=np.int32)
    return X, y


def apply_term_stats(rows: Sequence[CandidateRow],
                     term_stats: Dict[str, Dict[str, float]]) -> None:
    """Fill the family-size and lineage-breadth feature columns in place.

    Used after the priors are estimated on the training split, so every split
    (train/validation/test) sees the same training-derived priors and no test
    label frequencies leak into the features.
    """
    i_family = FEATURE_NAMES.index("family_size_norm")
    i_lineage = FEATURE_NAMES.index("lineage_breadth_norm")
    for r in rows:
        r.x[i_family] = _term_prior(term_stats, r.term_id,
                                    "family_size", _FAMILY_SAT)
        r.x[i_lineage] = _term_prior(term_stats, r.term_id,
                                     "lineage_breadth", _LINEAGE_SAT)


# ---------------------------------------------------------------------------
# Term priors (family size / lineage breadth), learned on the training split
# ---------------------------------------------------------------------------

def compute_term_stats(
    rows: Sequence[CandidateRow],
    term_taxa: Optional[Dict[str, set]] = None,
) -> Dict[str, Dict[str, float]]:
    """Count how often each true term occurs (family size).

    Family size is the number of proteins that carry the term; lineage breadth
    is the number of distinct taxa it spans. ``term_taxa`` optionally maps
    ``protein_key -> taxon`` so breadth can be computed; without it breadth is
    left equal to family size. Only positive (true) rows are counted.
    """
    proteins_by_term: Dict[str, set] = {}
    taxa_by_term: Dict[str, set] = {}
    for r in rows:
        if r.y != 1:
            continue
        proteins_by_term.setdefault(r.term_id, set()).add(r.protein_key)
        if term_taxa and r.protein_key in term_taxa:
            taxa_by_term.setdefault(r.term_id, set()).add(
                term_taxa[r.protein_key])
    stats: Dict[str, Dict[str, float]] = {}
    for term, proteins in proteins_by_term.items():
        stats[term] = {
            "family_size": float(len(proteins)),
            "lineage_breadth": float(len(taxa_by_term.get(term, proteins))),
        }
    return stats


# ---------------------------------------------------------------------------
# CAFA-style temporal / group split
# ---------------------------------------------------------------------------

def temporal_group_split(
    keys: Sequence[str],
    dates: Optional[Dict[str, str]] = None,
    test_frac: float = 0.2,
    val_frac: float = 0.1,
    seed: int = 42,
) -> Tuple[List[str], List[str], List[str]]:
    """Split proteins into train/val/test without leaking a protein across sets.

    When ``dates`` (protein key -> ISO date) are given the split is a CAFA
    temporal split: the most recent proteins are held out for testing, the
    preceding band is the validation set and the oldest proteins are used for
    training. Without dates the proteins are shuffled by group with a fixed
    seed.
    """
    keys = sorted(set(keys))
    n = len(keys)
    if n == 0:
        return [], [], []
    if dates:
        ordered = sorted(keys, key=lambda k: (dates.get(k, ""), k))
        n_test = max(1, int(round(n * test_frac))) if n > 2 else 1
        n_val = max(1, int(round(n * val_frac))) if n > 4 else 0
        test = ordered[-n_test:]
        val = ordered[-(n_test + n_val):-n_test] if n_val else []
        train = ordered[: n - n_test - n_val]
        return train, val, test
    rng = np.random.default_rng(seed)
    shuffled = list(keys)
    rng.shuffle(shuffled)
    n_test = max(1, int(round(n * test_frac))) if n > 2 else 1
    n_val = max(1, int(round(n * val_frac))) if n > 4 else 0
    test = shuffled[:n_test]
    val = shuffled[n_test:n_test + n_val]
    train = shuffled[n_test + n_val:]
    return train, val, test


def select_rows(rows: Sequence[CandidateRow], keys: Sequence[str]
                ) -> List[CandidateRow]:
    keep = set(keys)
    return [r for r in rows if r.protein_key in keep]


# ---------------------------------------------------------------------------
# Numpy dependency-free logistic ranker
# ---------------------------------------------------------------------------

class _LogisticRanker:
    """L2 logistic regression with class-balancing (pointwise ranking)."""

    def __init__(self, lr: float = 0.2, n_iter: int = 600,
                 l2: float = 1e-3, seed: int = 42):
        self.lr = lr
        self.n_iter = n_iter
        self.l2 = l2
        self.seed = seed
        self.w: Optional[np.ndarray] = None
        self.mean: Optional[np.ndarray] = None
        self.std: Optional[np.ndarray] = None

    def _standardise_fit(self, X: np.ndarray) -> np.ndarray:
        self.mean = X.mean(axis=0)
        self.std = X.std(axis=0)
        self.std[self.std < 1e-8] = 1.0
        return (X - self.mean) / self.std

    def standardise(self, X: np.ndarray) -> np.ndarray:
        if self.mean is None or self.std is None:
            return X
        return (X - self.mean) / self.std

    def fit(self, X: np.ndarray, y: np.ndarray) -> "_LogisticRanker":
        Xs = self._standardise_fit(X)
        Xa = np.hstack([Xs, np.ones((Xs.shape[0], 1))])
        n, d = Xa.shape
        pos = max(int(y.sum()), 1)
        neg = max(n - pos, 1)
        w_pos = n / (2.0 * pos)
        w_neg = n / (2.0 * neg)
        sw = np.where(y == 1, w_pos, w_neg)
        rng = np.random.default_rng(self.seed)
        self.w = rng.normal(0.0, 0.01, size=d)
        for _ in range(self.n_iter):
            p = _sigmoid(Xa @ self.w)
            grad = Xa.T @ (sw * (p - y)) / n
            grad[:-1] += self.l2 * self.w[:-1]
            self.w -= self.lr * grad
        return self

    def decision(self, X: np.ndarray) -> np.ndarray:
        Xa = np.hstack([self.standardise(X), np.ones((X.shape[0], 1))])
        return Xa @ self.w

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        return _sigmoid(self.decision(X))


# ---------------------------------------------------------------------------
# Public ranker wrapper
# ---------------------------------------------------------------------------

@dataclass
class LTRRanker:
    """A trained fusion ranker with its feature spec and term priors."""

    backend: str = "logreg"
    feature_names: List[str] = field(default_factory=lambda: list(FEATURE_NAMES))
    term_stats: Dict[str, Dict[str, float]] = field(default_factory=dict)
    model: Any = None
    logreg: Optional[_LogisticRanker] = None

    def score(self, X: np.ndarray) -> np.ndarray:
        """Return ranking scores (larger is more relevant)."""
        if list(self.feature_names) != list(FEATURE_NAMES):
            raise ValueError("feature order mismatch: model was trained on a "
                             "different feature set than the runtime builder")
        if self.backend == "lightgbm":
            raw = np.asarray(self.model.predict(X), dtype=np.float64)
            return _sigmoid(raw)
        if self.backend == "sklearn":
            return np.asarray(self.model.predict_proba(X)[:, 1], dtype=np.float64)
        if self.logreg is not None:
            return self.logreg.predict_proba(X)
        raise RuntimeError("LTR ranker has no fitted backend")

    def save(self, path: str) -> None:
        with open(path, "wb") as handle:
            pickle.dump(self, handle)

    @staticmethod
    def load(path: str) -> "LTRRanker":
        with open(path, "rb") as handle:
            obj = pickle.load(handle)
        if not isinstance(obj, LTRRanker):
            raise ValueError(f"{path} does not contain an LTRRanker")
        return obj


def _available_backend(preferred: str) -> str:
    order = {
        "auto": ["lightgbm", "sklearn", "logreg"],
        "lightgbm": ["lightgbm"],
        "sklearn": ["sklearn"],
        "logreg": ["logreg"],
    }[preferred]
    for backend in order:
        if backend == "lightgbm":
            try:
                import lightgbm  # noqa: F401
                return "lightgbm"
            except Exception:
                continue
        if backend == "sklearn":
            try:
                import sklearn  # noqa: F401
                return "sklearn"
            except Exception:
                continue
        return "logreg"
    return "logreg"


def fit_ranker(
    X: np.ndarray,
    y: np.ndarray,
    groups: Sequence[int],
    X_val: Optional[np.ndarray] = None,
    y_val: Optional[np.ndarray] = None,
    groups_val: Optional[Sequence[int]] = None,
    backend: str = "auto",
    seed: int = 42,
) -> LTRRanker:
    """Fit a ranker to one training ranking set."""
    chosen = _available_backend(backend)
    ranker = LTRRanker(backend=chosen)

    if chosen == "lightgbm":
        from lightgbm import LGBMRanker
        model = LGBMRanker(
            objective="lambdarank",
            metric="ndcg",
            n_estimators=400,
            learning_rate=0.05,
            num_leaves=31,
            min_child_samples=20,
            subsample=0.9,
            colsample_bytree=0.9,
            random_state=seed,
            n_jobs=-1,
        )
        fit_kwargs: Dict[str, Any] = dict(
            X=X, y=y, group=np.asarray(groups, dtype=np.int32))
        if X_val is not None and y_val is not None and groups_val:
            fit_kwargs["eval_set"] = [(X_val, y_val)]
            fit_kwargs["eval_group"] = [np.asarray(groups_val, dtype=np.int32)]
            fit_kwargs["eval_metric"] = "ndcg"
        model.fit(**fit_kwargs)
        ranker.model = model
        logger.info("Fitted LightGBM LambdaRank on %d rows / %d groups.",
                    len(y), len(groups))
    elif chosen == "sklearn":
        from sklearn.ensemble import GradientBoostingClassifier
        model = GradientBoostingClassifier(
            n_estimators=200, max_depth=3, learning_rate=0.05,
            subsample=0.9, random_state=seed)
        model.fit(X, y)
        ranker.model = model
        logger.info("Fitted sklearn gradient-boosted classifier on %d rows.", len(y))
    else:
        ranker.logreg = _LogisticRanker(seed=seed).fit(X, y)
        logger.info("Fitted numpy logistic ranker on %d rows (no ML backend).", len(y))
    return ranker


# ---------------------------------------------------------------------------
# CAFA-style metrics
# ---------------------------------------------------------------------------

def _score_predictions(
    rows: Sequence[CandidateRow], scores: np.ndarray,
) -> Tuple[Dict[str, List[Tuple[str, float]]], Dict[str, set]]:
    pred: Dict[str, List[Tuple[str, float]]] = {}
    truth: Dict[str, set] = {}
    for r, s in zip(rows, scores):
        pred.setdefault(r.protein_key, []).append((r.term_id, float(s)))
        if r.y == 1:
            truth.setdefault(r.protein_key, set()).add(r.term_id)
    for key in truth:
        pred.setdefault(key, [])
    return pred, truth


def fmax(
    rows: Sequence[CandidateRow], scores: np.ndarray,
    thresholds: Optional[np.ndarray] = None,
) -> Tuple[float, float, float, int]:
    """Protein-centric Fmax; returns (fmax, best_threshold, precision, recall).

    Recall (denominator) covers every protein that has a true annotation;
    proteins without a prediction contribute zero recall at that threshold.
    The count returned is the number of proteins with true labels.
    """
    pred, truth = _score_predictions(rows, scores)
    if not truth:
        return 0.0, 0.0, 0.0, 0
    if thresholds is None:
        thresholds = np.linspace(0.0, 1.0, 101)
    best = (0.0, 0.0, 0.0, 0.0)
    for thr in thresholds:
        tp = fp = 0
        n_predicted = 0
        for key, true_terms in truth.items():
            hits = {tid for tid, s in pred.get(key, []) if s >= thr}
            if hits:
                n_predicted += 1
            tp += len(hits & true_terms)
            fp += len(hits - true_terms)
        total_truth = sum(len(t) for t in truth.values())
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / total_truth if total_truth else 0.0
        f1 = (2 * precision * recall / (precision + recall)
              if precision + recall else 0.0)
        if f1 > best[0]:
            best = (f1, float(thr), precision, recall)
    return best[0], best[1], best[2], best[3]


def micro_aupr(rows: Sequence[CandidateRow], scores: np.ndarray) -> float:
    """Micro-average precision-recall AUC over every candidate row."""
    y = np.asarray([r.y for r in rows], dtype=np.int32)
    if y.sum() == 0 or y.sum() == len(y):
        return 0.0
    order = np.argsort(-scores, kind="mergesort")
    y_sorted = y[order]
    tp = np.cumsum(y_sorted)
    fp = np.cumsum(1 - y_sorted)
    precision = tp / np.maximum(tp + fp, 1)
    recall = tp / float(y.sum())
    # Trapezoidal integration over recall.
    aupr = float(np.sum(np.diff(recall, prepend=0.0) * precision))
    return aupr


def smin(
    rows: Sequence[CandidateRow], scores: np.ndarray,
    information_content: Dict[str, float],
    thresholds: Optional[np.ndarray] = None,
) -> Tuple[float, float]:
    """Minimum semantic distance (Smin) using information content.

    ``information_content`` maps a term id to its information content. The
    protein-aggregated remaining uncertainty (RU) and missing information (MI)
    follow the standard CAFA definitions; Smin is their Euclidean minimum over
    thresholds. Returns (smin, threshold).
    """
    pred, truth = _score_predictions(rows, scores)
    if not truth or not information_content:
        return float("nan"), float("nan")
    if thresholds is None:
        thresholds = np.linspace(0.0, 1.0, 51)

    def ic_sum(terms: Iterable[str]) -> float:
        return sum(float(information_content.get(t, 0.0)) for t in terms)

    best = (float("inf"), float("nan"))
    n = len(truth)
    for thr in thresholds:
        ru = mi = 0.0
        for key, true_terms in truth.items():
            hits = {tid for tid, s in pred.get(key, []) if s >= thr}
            ru += ic_sum(true_terms - hits)
            mi += ic_sum(hits - true_terms)
        ru /= n
        mi /= n
        s = math.sqrt(ru * ru + mi * mi)
        if s < best[0]:
            best = (s, float(thr))
    return best[0], best[1]
