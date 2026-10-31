"""Post-hoc calibration of fused raw scores.

Three calibrators are provided:

* ``IdentityCalibrator``  — passthrough, used when no held-out calibration set
  is available.
* ``TemperatureScaling``  — a simple logistic map in logit space
  (``p_cal = sigmoid(logit(p_raw) / T)``).  T is fitted on a labelled set by
  minimising the binary negative log-likelihood;  T = 1 reproduces the input.
* ``IsotonicCalibrator``  — isotonic (pool-adjacent-violators) regression.
  scikit-learn is used when it is importable; otherwise a small dependency-free
  PAVA implementation is used so the whole stage still runs in a minimal
  click+numpy+pytest environment.

``build_calibrator`` returns ``(calibrator, metadata)``; the metadata dict is
serialised into ``07_rank/calibration.json``.
"""

from __future__ import annotations

import math
import pickle
from typing import Any, List, Optional, Sequence

from pufprobe.logging_utils import get_logger

log = get_logger("rank.calibrate")

_EPS = 1e-6


def _sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def _logit(p: float) -> float:
    p = min(max(p, _EPS), 1.0 - _EPS)
    return math.log(p / (1.0 - p))


class IdentityCalibrator:
    """Passthrough calibration."""

    def transform(self, scores: Sequence[float]) -> List[float]:
        return [float(s) for s in scores]


def expected_calibration_error(scores: Sequence[float], labels: Sequence[int], n_bins: int = 10) -> float:
    """Standard 10-bin expected calibration error."""
    if not scores:
        return float("nan")
    bin_edges = [i / n_bins for i in range(n_bins + 1)]
    n = len(scores)
    ece = 0.0
    for b in range(n_bins):
        lo, hi = bin_edges[b], bin_edges[b + 1]
        idx = [i for i, s in enumerate(scores) if lo <= s <= hi] if b == n_bins - 1 else [
            i for i, s in enumerate(scores) if lo <= s < hi
        ]
        if not idx:
            continue
        avg_conf = sum(scores[i] for i in idx) / len(idx)
        avg_acc = sum(labels[i] for i in idx) / len(idx)
        ece += (len(idx) / n) * abs(avg_conf - avg_acc)
    return float(ece)


class TemperatureScaling:
    """Logistic temperature map fitted by NLL grid search."""

    def __init__(self, temperature: float = 1.0):
        self.temperature = float(temperature)

    def fit(self, scores: Sequence[float], labels: Sequence[int]) -> "TemperatureScaling":
        if not scores:
            return self
        best_t, best_nll = self.temperature, float("inf")
        for t in [round(0.2 + 0.05 * i, 3) for i in range(0, 161)]:  # 0.2 .. 8.0
            nll = 0.0
            for s, y in zip(scores, labels):
                p = _sigmoid(_logit(s) / t)
                p = min(max(p, _EPS), 1.0 - _EPS)
                nll -= y * math.log(p) + (1 - y) * math.log(1 - p)
            if nll < best_nll:
                best_nll, best_t = nll, t
        self.temperature = best_t
        return self

    def transform(self, scores: Sequence[float]) -> List[float]:
        return [_sigmoid(_logit(s) / self.temperature) for s in scores]


class IsotonicCalibrator:
    """Isotonic regression; sklearn when available, PAVA fallback otherwise."""

    def __init__(self) -> None:
        self._backend = "identity"
        self._model = IdentityCalibrator()
        self._xs: List[float] = []
        self._ys: List[float] = []

    def fit(self, scores: Sequence[float], labels: Sequence[int]) -> "IsotonicCalibrator":
        if not scores:
            return self
        try:
            from sklearn.isotonic import IsotonicRegression  # type: ignore

            self._model = IsotonicRegression(out_of_bounds="clip")
            self._model.fit(list(scores), list(labels))
            self._backend = "sklearn.isotonic"
            return self
        except Exception:  # pragma: no cover - exercised only without sklearn
            log.info("sklearn unavailable; using built-in PAVA isotonic regression.")
        # Pool Adjacent Violators Algorithm (equal weights).
        order = sorted(range(len(scores)), key=lambda i: scores[i])
        xs = [float(scores[i]) for i in order]
        ys = [float(labels[i]) for i in order]
        # Block start indices and running sums.
        block_start: List[int] = [0]
        block_sum: List[float] = [ys[0]]
        block_len: List[int] = [1]
        for i in range(1, len(xs)):
            block_start.append(len(block_sum))
            block_sum.append(ys[i])
            block_len.append(1)
            while len(block_sum) > 1 and block_sum[-2] / block_len[-2] > block_sum[-1] / block_len[-1]:
                block_sum[-2] += block_sum[-1]
                block_len[-2] += block_len[-1]
                block_sum.pop()
                block_len.pop()
                block_start.pop()
        self._xs = [xs[s] for s in block_start]
        self._ys = [block_sum[i] / block_len[i] for i in range(len(block_sum))]
        self._backend = "pava"
        return self

    def transform(self, scores: Sequence[float]) -> List[float]:
        if self._backend.startswith("sklearn"):
            return [float(v) for v in self._model.transform(list(scores))]  # type: ignore[union-attr]
        if not self._xs:
            return [float(s) for s in scores]
        out: List[float] = []
        for s in scores:
            if s <= self._xs[0]:
                out.append(self._ys[0])
            elif s >= self._xs[-1]:
                out.append(self._ys[-1])
            else:
                lo, hi = 0, len(self._xs) - 1
                while hi - lo > 1:
                    mid = (lo + hi) // 2
                    if self._xs[mid] <= s:
                        lo = mid
                    else:
                        hi = mid
                out.append(self._ys[lo])
        return out


def build_calibrator(calibration_set: Optional[Sequence[dict]] = None):
    """Return ``(calibrator, metadata)``.

    ``calibration_set`` is an optional iterable of ``{"score": float,
    "label": 0|1}`` records.  Without it the identity calibration is used and
    ``calibration.json`` records ``method="identity"``.
    """
    if not calibration_set:
        return IdentityCalibrator(), {
            "method": "identity",
            "ece": None,
            "note": "no held-out calibration set provided",
        }
    scores = [float(r["score"]) for r in calibration_set]
    labels = [int(r["label"]) for r in calibration_set]
    iso = IsotonicCalibrator().fit(scores, labels)
    cal_scores = iso.transform(scores)
    meta = {
        "method": "isotonic",
        "backend": iso._backend,
        "ece": expected_calibration_error(cal_scores, labels),
        "note": f"isotonic regression fitted on {len(scores)} held-out records",
    }
    return iso, meta


def save_calibrator(path: str, calibrator: Any) -> None:
    """Persist a fitted calibrator (temperature, isotonic or identity)."""
    with open(path, "wb") as handle:
        pickle.dump(calibrator, handle)


def load_calibrator(path: str) -> Any:
    """Load a calibrator previously stored with :func:`save_calibrator`."""
    with open(path, "rb") as handle:
        return pickle.load(handle)
