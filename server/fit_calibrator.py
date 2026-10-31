#!/usr/bin/env python
"""Fit a post-hoc confidence calibrator for a trained LTR fusion ranker.

The ranker produced by ``train_ltr.py`` applies a sigmoid to its LightGBM
output so that scores lie in [0, 1], but the sigmoid is only a fixed map. The
values are optimised for within-protein ranking, not for absolute probability,
so a confidence of 0.80 does not yet imply an 80 percent chance of being
correct. This script rebuilds the same CAFA-style temporal split, scores the
validation and test proteins with the saved ranker, and fits temperature and
isotonic calibrators.

Calibration is fitted on the validation split only. The validation proteins
are not used to fit ranker parameters, and the test split is never used for
fitting, so it stays an independent check of expected calibration error,
Brier score, negative log-likelihood and per-bin reliability.

Examples
--------
::

    python -u fit_calibrator.py \
        --evidence ltr_dataset/search_evidence_full4_white.json \
        --labels   ltr_dataset/labels.json \
        --dates    ltr_dataset/dates.json \
        --go-obo   dbs/dbdir/go-basic.obo \
        --model    dbs/dbdir/ltr_fusion_white.pkl \
        --calibrator-out dbs/dbdir/calibrator_white.pkl \
        --metrics-out ltr_dataset/metrics_calibration_white.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pufprobe.logging_utils import setup_logging, get_logger  # noqa: E402
from pufprobe.rank import ltr  # noqa: E402
from pufprobe.rank.calibrate import (  # noqa: E402
    TemperatureScaling,
    IsotonicCalibrator,
    expected_calibration_error,
    save_calibrator,
)
from train_ltr import (  # noqa: E402
    load_json,
    normalize_labels,
    align_to_keys,
    propagate_go_labels,
    load_evidence,
)

logger = get_logger("fit_calibrator")

_EPS = 1e-6


def brier_score(probs: Sequence[float], labels: Sequence[int]) -> float:
    """Mean squared difference between predicted probability and outcome."""
    p = np.asarray(probs, dtype=np.float64)
    y = np.asarray(labels, dtype=np.float64)
    return float(np.mean((p - y) ** 2))


def negative_log_likelihood(probs: Sequence[float], labels: Sequence[int]) -> float:
    """Mean binary negative log-likelihood."""
    total = 0.0
    for p, y in zip(probs, labels):
        p = min(max(float(p), _EPS), 1.0 - _EPS)
        total -= y * np.log(p) + (1 - y) * np.log(1.0 - p)
    return float(total / len(labels))


def reliability_table(
    probs: Sequence[float], labels: Sequence[int], n_bins: int = 10,
) -> List[Dict[str, Any]]:
    """Per-bin confidence, accuracy and count, for a reliability diagram."""
    p = np.asarray(probs, dtype=np.float64)
    y = np.asarray(labels, dtype=np.int32)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    table: List[Dict[str, Any]] = []
    for b in range(n_bins):
        lo, hi = edges[b], edges[b + 1]
        if b == n_bins - 1:
            idx = np.where((p >= lo) & (p <= hi))[0]
        else:
            idx = np.where((p >= lo) & (p < hi))[0]
        if idx.size == 0:
            continue
        table.append({
            "bin": b,
            "range": [round(float(lo), 2), round(float(hi), 2)],
            "count": int(idx.size),
            "mean_confidence": round(float(p[idx].mean()), 4),
            "accuracy": round(float(y[idx].mean()), 4),
        })
    return table


def score_split(
    ranker: ltr.LTRRanker,
    evidence_map: Dict[str, Any],
    labels: Dict[str, Dict[str, List[str]]],
    keys: Sequence[str],
    protein_meta: Optional[Dict[str, Any]],
) -> Tuple[np.ndarray, np.ndarray]:
    """Return raw ranker scores and binary labels for one protein split."""
    split_labels = {k: labels[k] for k in keys}
    rows, _ = ltr.build_dataset(
        evidence_map, split_labels,
        term_stats=getattr(ranker, "term_stats", None),
        protein_meta=protein_meta,
    )
    X, y = ltr.rows_to_matrix(rows)
    scores = ranker.score(X)
    return np.asarray(scores, dtype=np.float64), y


def method_metrics(
    probs: Sequence[float], labels: Sequence[int],
) -> Dict[str, float]:
    """Calibration summaries for one set of probabilities."""
    return {
        "nll": round(negative_log_likelihood(probs, labels), 5),
        "brier": round(brier_score(probs, labels), 5),
        "ece": round(expected_calibration_error(probs, labels), 5),
    }


def main() -> int:
    setup_logging()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, nargs="+", required=True,
                        help="One or more search_evidence.json files.")
    parser.add_argument("--labels", type=Path, required=True,
                        help="Protein key -> true GO/EC/Pfam/COG term lists.")
    parser.add_argument("--dates", type=Path, default=None,
                        help="Protein key -> ISO date, for a temporal split.")
    parser.add_argument("--plddt", type=Path, default=None,
                        help="Protein key -> mean pLDDT (meta feature).")
    parser.add_argument("--go-obo", type=Path, default=None,
                        help="go-basic.obo, used to close GO labels (true path).")
    parser.add_argument("--model", type=Path, required=True,
                        help="Trained LTR ranker saved by train_ltr.py.")
    parser.add_argument("--calibrator-out", type=Path, required=True,
                        help="Output path for the selected fitted calibrator.")
    parser.add_argument("--metrics-out", type=Path, default=None,
                        help="Optional JSON file for calibration metrics.")
    parser.add_argument("--method", choices=["auto", "temperature", "isotonic"],
                        default="auto",
                        help="auto selects by validation NLL.")
    parser.add_argument("--test-frac", type=float, default=0.2)
    parser.add_argument("--val-frac", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    logger.info("Loading evidence from %d file(s) ...", len(args.evidence))
    evidence_map = load_evidence(args.evidence)
    evidence_keys = sorted(evidence_map)
    logger.info("Evidence loaded for %d proteins.", len(evidence_keys))

    labels = normalize_labels(load_json(args.labels))
    labels = align_to_keys(labels, evidence_keys)
    labels = {k: v for k, v in labels.items() if v}
    n_go_added = propagate_go_labels(labels, args.go_obo)
    if n_go_added:
        logger.info("GO label true-path closure added %d ancestor terms.", n_go_added)

    dates = None
    if args.dates:
        dates = align_to_keys(load_json(args.dates), evidence_keys, default="")
    protein_meta = None
    if args.plddt:
        raw_plddt = align_to_keys(load_json(args.plddt), evidence_keys, default=0.0)
        protein_meta = {k: {"plddt": float(v)} for k, v in raw_plddt.items()}

    usable = sorted(set(labels) & set(evidence_keys))
    if not usable:
        raise SystemExit("No protein has both evidence and labels; cannot fit.")
    labels = {k: labels[k] for k in usable}

    train_keys, val_keys, test_keys = ltr.temporal_group_split(
        usable, dates=dates, test_frac=args.test_frac,
        val_frac=args.val_frac, seed=args.seed,
    )
    logger.info("Split: %d train / %d val / %d test proteins.",
                len(train_keys), len(val_keys), len(test_keys))

    ranker = ltr.LTRRanker.load(str(args.model))
    logger.info("Loaded LTR ranker from %s (backend=%s).", args.model, ranker.backend)

    val_raw, val_y = score_split(
        ranker, evidence_map, labels, val_keys, protein_meta)
    test_raw, test_y = score_split(
        ranker, evidence_map, labels, test_keys, protein_meta)
    logger.info("Validation rows: %d, test rows: %d.", len(val_y), len(test_y))

    # Fit both calibrators on validation only.
    temperature = TemperatureScaling().fit(val_raw, val_y)
    isotonic = IsotonicCalibrator().fit(val_raw, val_y)
    val_temperature = temperature.transform(val_raw)
    val_isotonic = isotonic.transform(val_raw)

    val_summary: Dict[str, Any] = {
        "identity": method_metrics(val_raw, val_y),
        "temperature": method_metrics(val_temperature, val_y),
        "isotonic": method_metrics(val_isotonic, val_y),
    }
    logger.info("Validation NLL: identity=%.4f temperature=%.4f isotonic=%.4f",
                val_summary["identity"]["nll"],
                val_summary["temperature"]["nll"],
                val_summary["isotonic"]["nll"])

    if args.method == "temperature":
        chosen_name, chosen = "temperature", temperature
    elif args.method == "isotonic":
        chosen_name, chosen = "isotonic", isotonic
    else:
        candidates = {
            "temperature": val_summary["temperature"]["nll"],
            "isotonic": val_summary["isotonic"]["nll"],
        }
        chosen_name = min(candidates, key=candidates.get)
        chosen = temperature if chosen_name == "temperature" else isotonic

    args.calibrator_out.parent.mkdir(parents=True, exist_ok=True)
    save_calibrator(str(args.calibrator_out), chosen)
    logger.info("Saved %s calibrator to %s.", chosen_name, args.calibrator_out)

    # Evaluate every method on the untouched test split.
    test_temperature = temperature.transform(test_raw)
    test_isotonic = isotonic.transform(test_raw)
    test_summary: Dict[str, Any] = {
        "identity": method_metrics(test_raw, test_y),
        "temperature": method_metrics(test_temperature, test_y),
        "isotonic": method_metrics(test_isotonic, test_y),
    }
    chosen_test_probs = (
        test_temperature if chosen_name == "temperature" else test_isotonic
    )

    metrics: Dict[str, Any] = {
        "model": str(args.model),
        "ranker_backend": ranker.backend,
        "split_kind": "temporal" if dates else "seeded_group",
        "n_train_proteins": len(train_keys),
        "n_val_proteins": len(val_keys),
        "n_test_proteins": len(test_keys),
        "n_val_rows": int(len(val_y)),
        "n_test_rows": int(len(test_y)),
        "selected_method": chosen_name,
        "validation": val_summary,
        "test": test_summary,
        "temperature": round(float(temperature.temperature), 4),
        "test_reliability_identity": reliability_table(test_raw, test_y),
        "test_reliability_selected": reliability_table(
            chosen_test_probs, test_y),
    }

    if args.metrics_out:
        args.metrics_out.parent.mkdir(parents=True, exist_ok=True)
        args.metrics_out.write_text(
            json.dumps(metrics, indent=2), encoding="utf-8")

    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
