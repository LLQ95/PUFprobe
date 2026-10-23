#!/usr/bin/env python
"""Train the learning-to-rank fusion model on an annotated hold-out set.

The training set is produced by running the PUFprobe search channels against a
collection of proteins whose annotations are already known (the reviewed
Swiss-Prot hold-out used by ``baseline_holdout.py`` is one such source). Each
query's own annotations are hidden during the search, and the evidence found
for every candidate term is recorded in ``search_evidence.json``. This script
rebuilds the per-(protein, term) feature matrix, labels each candidate against
the hidden truth, splits proteins by a CAFA-style temporal cut (or a seeded
group split when no dates are available), fits a ranker and reports Fmax,
AUPR and, optionally, semantic distance Smin.

Examples
--------
Temporal split with GO true-path labels and a LightGBM LambdaRank model::

    python -u train_ltr.py \
        --evidence workdir/05_search/search_evidence.json \
        --labels   holdout_labels.json \
        --dates    holdout_dates.json \
        --go-obo   dbdir/go-basic.obo \
        --plddt    holdout_plddt.json \
        --model-out dbdir/ltr_fusion.pkl \
        --backend auto

A minimal dependency-free run (numpy logistic ranker) is identical without
``--backend`` and without lightgbm installed.

The label file maps a protein key (or a bare UniProt accession, matched to the
``sample@locus`` tail) to true term lists, e.g.::

    {"P12345": {"GO": ["GO:0003674"], "EC": ["1.2.3.4"],
                "Pfam": ["PF00001"], "COG": []}}
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pufprobe.logging_utils import setup_logging, get_logger  # noqa: E402
from pufprobe.ontology_graph import load_go_obo  # noqa: E402
from pufprobe.rank.fuse import build_evidence_map  # noqa: E402
from pufprobe.rank import ltr  # noqa: E402

logger = get_logger("train_ltr")


def load_json(path: Path) -> Any:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _normalise_ec(term: str) -> str:
    t = str(term).strip()
    if t and not t.startswith("EC:") and "." in t and t[0].isdigit():
        return f"EC:{t}"
    return t


def normalize_labels(raw: Dict[str, Any]) -> Dict[str, Dict[str, List[str]]]:
    """Accept either the nested label format or a Swiss-Prot xref mapping."""
    out: Dict[str, Dict[str, List[str]]] = {}
    for key, value in raw.items():
        if not isinstance(value, dict):
            continue
        by_type: Dict[str, List[str]] = defaultdict(list)
        for ttype, terms in value.items():
            ttype = str(ttype)
            if ttype not in ("GO", "EC", "Pfam", "COG"):
                continue
            for t in terms or []:
                t = _normalise_ec(t) if ttype == "EC" else str(t)
                if t and t not in by_type[ttype]:
                    by_type[ttype].append(t)
        out[str(key)] = dict(by_type)
    return out


def align_to_keys(
    mapping: Dict[str, Any], evidence_keys: Sequence[str], default: Any = None,
) -> Dict[str, Any]:
    """Map a keyed file onto the actual evidence keys, tolerating accessions."""
    out: Dict[str, Any] = {}
    for ek in evidence_keys:
        if ek in mapping:
            out[ek] = mapping[ek]
            continue
        accession = ek.split("@")[-1]
        if accession in mapping:
            out[ek] = mapping[accession]
        elif default is not None:
            out[ek] = default
    return out


def propagate_go_labels(
    labels: Dict[str, Dict[str, List[str]]], obo_path: Optional[Path],
) -> int:
    """Extend every GO label to its true-path closure. Returns added count."""
    if obo_path is None:
        return 0
    obo = load_go_obo(obo_path)
    if obo is None:
        logger.warning("Could not load GO ontology from %s; labels kept as-is.", obo_path)
        return 0
    added = 0
    for by_type in labels.values():
        direct = set(by_type.get("GO", []))
        closed: set = set()
        for term in direct:
            closed.update(obo.ancestors(term))
        by_type["GO"] = sorted(closed)
        added += len(closed - direct)
    return added


def load_evidence(paths: Sequence[Path]) -> Dict[str, Any]:
    merged: Dict[str, Any] = {}
    for path in paths:
        part = build_evidence_map(Path(path), None)
        for key, evs in part.items():
            merged.setdefault(key, []).extend(evs)
    return merged


def groups_for_rows(rows: Sequence[ltr.CandidateRow], keys: Sequence[str]
                    ) -> List[int]:
    keep = set(keys)
    counts: Dict[str, int] = defaultdict(int)
    for r in rows:
        if r.protein_key in keep:
            counts[r.protein_key] += 1
    return [counts[k] for k in keys if counts.get(k, 0) > 0]


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
    parser.add_argument("--taxa", type=Path, default=None,
                        help="Protein key -> taxon id, for lineage breadth.")
    parser.add_argument("--go-obo", type=Path, default=None,
                        help="go-basic.obo, used to close GO labels (true path).")
    parser.add_argument("--ic", type=Path, default=None,
                        help="Term id -> information content, for Smin.")
    parser.add_argument("--model-out", type=Path, required=True)
    parser.add_argument("--metrics-out", type=Path, default=None)
    parser.add_argument("--backend", choices=["auto", "lightgbm", "sklearn",
                                              "logreg"], default="auto")
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
    term_taxa = None
    if args.taxa:
        term_taxa = align_to_keys(load_json(args.taxa), evidence_keys, default="")

    # Keep only proteins that carry both evidence and truth.
    usable = sorted(set(labels) & set(evidence_keys))
    if not usable:
        raise SystemExit("No protein has both evidence and labels; cannot train.")
    labels = {k: labels[k] for k in usable}
    logger.info("Training proteins with evidence + labels: %d.", len(usable))

    train_keys, val_keys, test_keys = ltr.temporal_group_split(
        usable, dates=dates, test_frac=args.test_frac,
        val_frac=args.val_frac, seed=args.seed,
    )
    logger.info("Split: %d train / %d val / %d test proteins (%s).",
                len(train_keys), len(val_keys), len(test_keys),
                "temporal" if dates else "seeded group")

    # Term priors are learned on the training split only, then applied to every
    # split so test frequencies never leak into the features.
    train_labels = {k: labels[k] for k in train_keys}
    train_rows_pre, _ = ltr.build_dataset(
        evidence_map, train_labels, term_stats=None,
        protein_meta=protein_meta,
    )
    term_stats = ltr.compute_term_stats(train_rows_pre, term_taxa=term_taxa)
    logger.info("Term priors estimated for %d terms from the train split.",
                len(term_stats))

    def _rows_for(keys: Sequence[str]) -> List[ltr.CandidateRow]:
        rows, _ = ltr.build_dataset(
            evidence_map, {k: labels[k] for k in keys},
            term_stats=term_stats, protein_meta=protein_meta,
        )
        ltr.apply_term_stats(rows, term_stats)
        return rows

    train_rows = _rows_for(train_keys)
    val_rows = _rows_for(val_keys)
    test_rows = _rows_for(test_keys)
    if not test_rows:
        raise SystemExit("Test split produced no candidate rows; lower test_frac.")

    X_train, y_train = ltr.rows_to_matrix(train_rows)
    X_val, y_val = (ltr.rows_to_matrix(val_rows) if val_rows else (None, None))
    X_test, y_test = ltr.rows_to_matrix(test_rows)
    g_train = groups_for_rows(train_rows, train_keys)
    g_val = groups_for_rows(val_rows, val_keys) if val_rows else None
    pos_rate = float(y_train.mean()) if len(y_train) else 0.0
    logger.info("Feature matrix: %d train rows (positive rate %.3f), %d test rows.",
                len(y_train), pos_rate, len(y_test))

    ranker = ltr.fit_ranker(
        X_train, y_train, g_train,
        X_val=X_val, y_val=y_val, groups_val=g_val,
        backend=args.backend, seed=args.seed,
    )
    ranker.term_stats = term_stats

    test_scores = ranker.score(X_test)
    f_best, thr, precision, recall = ltr.fmax(test_rows, test_scores)
    aupr = ltr.micro_aupr(test_rows, test_scores)
    metrics: Dict[str, Any] = {
        "backend": ranker.backend,
        "feature_names": ltr.FEATURE_NAMES,
        "n_train_proteins": len(train_keys),
        "n_val_proteins": len(val_keys),
        "n_test_proteins": len(test_keys),
        "n_train_rows": int(len(y_train)),
        "n_test_rows": int(len(y_test)),
        "train_positive_rate": round(pos_rate, 5),
        "test_fmax": round(float(f_best), 4),
        "fmax_threshold": round(float(thr), 4),
        "test_precision_at_fmax": round(float(precision), 4),
        "test_recall_at_fmax": round(float(recall), 4),
        "test_micro_aupr": round(float(aupr), 4),
        "split_kind": "temporal" if dates else "seeded_group",
        "n_term_priors": len(term_stats),
    }

    if args.ic:
        ic = load_json(args.ic)
        ic = {str(k): float(v) for k, v in ic.items()}
        s_best, s_thr = ltr.smin(test_rows, test_scores, ic)
        if np.isfinite(s_best):
            metrics["test_smin"] = round(float(s_best), 4)
            metrics["smin_threshold"] = round(float(s_thr), 4)

    if ranker.backend == "lightgbm" and hasattr(ranker.model, "feature_importances_"):
        importance = sorted(
            zip(ltr.FEATURE_NAMES,
                (int(v) for v in ranker.model.feature_importances_)),
            key=lambda kv: -kv[1],
        )
        metrics["feature_importance"] = [
            {"feature": name, "importance": value} for name, value in importance
        ]

    args.model_out.parent.mkdir(parents=True, exist_ok=True)
    ranker.save(str(args.model_out))
    logger.info("Model saved to %s (backend=%s).", args.model_out, ranker.backend)

    if args.metrics_out:
        args.metrics_out.parent.mkdir(parents=True, exist_ok=True)
        args.metrics_out.write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
