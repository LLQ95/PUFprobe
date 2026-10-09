"""Evidence fusion: turn per-channel evidence into candidate function terms.

For every CDS record, all incoming Evidence objects are grouped into candidate
terms.  A term's raw score is the channel-weighted sum of the best evidence per
channel::

    raw = sum_c  weight_c * max(score_c for this term)

Channels with no evidence for the term contribute nothing.  A diversity rule
then caps terms that are backed by too few independent evidence groups below
``rank.high_confidence``.

Term type is inferred from the evidence target:

* ``GO:xxxxxxx`` ids                      -> GO
* EC numbers (optionally ``EC:`` prefixed) -> EC
* ``PFxxxxx`` accessions                   -> Pfam
* ``COGxxxx`` accessions                   -> COG
* foldseek/structure-derived targets       -> PDB / structure
* anything else                            -> generic ``label`` (target id/name)
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from pufprobe.io_utils import read_json
from pufprobe.logging_utils import get_logger
from pufprobe.models import (
    CHANNEL_GROUPS,
    CDSRecord,
    Evidence,
    EvidenceChannel,
    GeneStatus,
    GeneVerdict,
    ProteinPrediction,
    TermPrediction,
)

log = get_logger("rank.fuse")

CHANNEL_ORDER = [
    "plm_neighbor",
    "foldseek",
    "remote_homology",
    "pocket",
    "context",
    "database",
]

_GO_RE = re.compile(r"^GO:\d{7}$")
_EC_RE = re.compile(r"^(EC:)?\d+\.\d+\.\d+\.[\dR]+(?:\.[A-Z])?$")
_PFAM_RE = re.compile(r"^PF\d{5}")
_COG_RE = re.compile(r"^COG\d{4}")
_PDB_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}(?:[._-][A-Za-z0-9]+)?$")

# AMR / virulence neighbourhood keywords shared with priority scoring.
AMR_VIRULENCE_KEYWORDS = (
    "card", "vfdb", "aroe", "amr", "antibiotic", "resistance", "virulence",
    "toxin", "tox", "beta-lactam", "bl2", "bla", "tet", "sul", "dfr", "qnr",
    "erm", "aph", "aac", "ant", "strA", "strB",
)


def classify_term(ev: Evidence) -> tuple[str, str]:
    """Return ``(term_type, normalized_term_id)`` for one evidence hit."""
    tid = (ev.target_id or "").strip()
    details = ev.details or {}
    db = str(details.get("db", "")).lower()

    if tid.startswith("GO:") or db == "go":
        return "GO", tid
    if _EC_RE.match(tid) or db == "ec":
        norm = tid if tid.startswith("EC:") else f"EC:{tid}"
        return "EC", norm
    if _PFAM_RE.match(tid) or db == "pfam":
        return "Pfam", tid
    if _COG_RE.match(tid) or db == "cog":
        return "COG", tid
    if (
        ev.channel == EvidenceChannel.FOLDSEEK
        or details.get("hit_type") == "structure"
        or _PDB_RE.match(tid)
    ):
        return ("PDB", tid) if _PDB_RE.match(tid) else ("structure", tid)
    return "label", tid


def build_evidence_map(search_path: Optional[Path], context_path: Optional[Path]) -> Dict[str, List[Evidence]]:
    """Merge ``search_evidence.json`` and ``context_evidence.json`` per global key."""
    merged: Dict[str, List[Evidence]] = {}

    def _load(path: Optional[Path]) -> Dict[str, List[Evidence]]:
        if path is None or not path.exists():
            return {}
        raw = read_json(path)
        out: Dict[str, List[Evidence]] = {}
        assert isinstance(raw, dict)
        for key, rows in raw.items():
            out[str(key)] = [_evidence_from_row(r) for r in rows]
        return out

    for source in (_load(search_path), _load(context_path)):
        for key, evs in source.items():
            merged.setdefault(key, []).extend(evs)
    return merged


def _evidence_from_row(row: dict) -> Evidence:
    row = dict(row)
    row["channel"] = EvidenceChannel(row["channel"])
    return Evidence(**row)


def term_from_row(row: dict) -> TermPrediction:
    """Reconstruct a TermPrediction without models.from_dict (typing-compat)."""
    return TermPrediction(
        term_type=row["term_type"],
        term_id=row["term_id"],
        term_name=row.get("term_name", row["term_id"]),
        raw_score=float(row["raw_score"]),
        confidence=float(row.get("confidence", row["raw_score"])),
        evidence=[_evidence_from_row(e) for e in row.get("evidence", [])],
        rationale=row.get("rationale", ""),
    )


def cds_record_from_row(row: dict) -> CDSRecord:
    """Reconstruct a CDSRecord without models.from_dict (typing-compat)."""
    row = dict(row)
    row["status"] = GeneStatus(row["status"])
    verdict = row.get("verdict")
    if verdict is not None:
        row["verdict"] = GeneVerdict(verdict)
    return CDSRecord(**row)


def prediction_from_row(row: dict) -> ProteinPrediction:
    return ProteinPrediction(
        sample_id=row["sample_id"],
        locus_tag=row["locus_tag"],
        terms=[term_from_row(t) for t in row.get("terms", [])],
        priority_score=row.get("priority_score"),
        priority_reasons=list(row.get("priority_reasons", [])),
    )


def _term_label(evs: List[Evidence]) -> str:
    """Pick the most descriptive label from the highest-scoring evidence."""
    best = max(evs, key=lambda e: e.score)
    return (best.target_label or best.target_id).strip()


def fuse_proteins(
    records: List[CDSRecord],
    evidence_map: Dict[str, List[Evidence]],
    config: Dict[str, Any],
    weights_model: Any = None,
) -> List[ProteinPrediction]:
    """Fuse evidence for every CDS record into ranked terms."""
    rank_cfg = config["rank"]
    rule_weights: Dict[str, float] = rank_cfg["rule_weights"]
    required_groups: int = int(rank_cfg["require_independent_groups"])
    high_conf: float = float(rank_cfg["high_confidence"])

    predictions: List[ProteinPrediction] = []
    seen_go = False

    for rec in records:
        evs = evidence_map.get(rec.key, [])
        groups: Dict[tuple[str, str], List[Evidence]] = {}
        for ev in evs:
            ttype, tid = classify_term(ev)
            if ttype == "GO" and not seen_go:
                log.info(
                    "GO true-path propagation not applied: no ontology library loaded; "
                    "GO candidates are kept as annotated."
                )
                seen_go = True
            groups.setdefault((ttype, tid), []).append(ev)

        terms: List[TermPrediction] = []
        for (ttype, tid), gevs in sorted(groups.items()):
            # best score per channel for this term
            best_per_channel: Dict[str, float] = {}
            for ev in gevs:
                ch = ev.channel.value
                best_per_channel[ch] = max(best_per_channel.get(ch, ev.score), ev.score)

            if weights_model is not None:
                try:
                    raw = _learned_score(weights_model, best_per_channel, gevs)
                    rationale = "score: learned fusion model"
                except Exception as exc:  # pragma: no cover - defensive
                    log.warning("learned fusion failed (%s); falling back to rule weights.", exc)
                    raw = 0.0
                    rationale = ""
                    weights_model = None
                else:
                    raw = min(max(raw, 0.0), 1.0)
                    terms.append(
                        TermPrediction(
                            term_type=ttype,
                            term_id=tid,
                            term_name=_term_label(gevs),
                            raw_score=raw,
                            confidence=raw,
                            evidence=list(gevs),
                            rationale=rationale,
                        )
                    )
                    continue
            raw = sum(rule_weights.get(ch, 0.0) * score for ch, score in best_per_channel.items())
            raw = min(max(raw, 0.0), 1.0)
            rationale = "score: weighted sum of channel maxima (rule weights)"

            group_set = sorted({CHANNEL_GROUPS[ev.channel.value] for ev in gevs})
            if len(group_set) < required_groups:
                capped = min(raw, high_conf - 0.01)
                rationale += (
                    f"; capped at {capped:.3f}: only {len(group_set)} independent "
                    f"evidence group(s) ({', '.join(group_set)}), need {required_groups}"
                )
                raw = capped

            terms.append(
                TermPrediction(
                    term_type=ttype,
                    term_id=tid,
                    term_name=_term_label(gevs),
                    raw_score=raw,
                    confidence=raw,
                    evidence=list(gevs),
                    rationale=rationale,
                )
            )

        terms.sort(key=lambda t: t.raw_score, reverse=True)
        predictions.append(
            ProteinPrediction(
                sample_id=rec.sample_id,
                locus_tag=rec.locus_tag,
                terms=terms,
            )
        )
    return predictions


def _learned_score(model: Any, best_per_channel: Dict[str, float], gevs: List[Evidence]) -> float:
    """Build a fixed feature vector and ask the loaded model for a probability."""
    import numpy as np  # local import; numpy is a hard dependency

    vec = [best_per_channel.get(ch, 0.0) for ch in CHANNEL_ORDER]
    vec.append(float(len({CHANNEL_GROUPS[ev.channel.value] for ev in gevs})))
    vec.append(float(len(gevs)))
    x = np.asarray([vec], dtype=float)
    if hasattr(model, "predict_proba"):
        proba = model.predict_proba(x)
        return float(proba[0][-1])
    return float(model.predict(x)[0])


def load_weights_model(weights_path: Optional[Path]) -> tuple[Any, str]:
    """Lazily load a learned fusion model.

    Returns ``(model, note)``.  Any failure (missing joblib/sklearn, bad file)
    yields ``(None, "fusion: rule weights (no learned model)")``.
    """
    if weights_path is None:
        return None, "fusion: rule weights (no learned model)"
    try:
        try:
            import joblib  # type: ignore
        except ImportError:
            log.warning("--weights given but joblib is not installed; using rule weights.")
            return None, "fusion: rule weights (no learned model)"
        model = joblib.load(str(weights_path))
    except Exception as exc:
        log.warning("could not load fusion weights %s (%s); using rule weights.", weights_path, exc)
        return None, "fusion: rule weights (no learned model)"
    return model, f"fusion: learned model loaded from {weights_path}"


def load_predictions(path: Path) -> List[ProteinPrediction]:
    """Read predictions.json written by an earlier rank run."""
    rows = read_json(path)
    assert isinstance(rows, list)
    return [prediction_from_row(row) for row in rows]
