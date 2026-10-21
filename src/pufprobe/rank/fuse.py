"""Evidence fusion: turn per-channel evidence into candidate function terms.

For every CDS record, all incoming Evidence objects are grouped into candidate
terms. Two rule-based fusion methods are supported and selected with
``rank.fusion_method``.

The default ``group_noisy_or`` method treats each independent evidence group
(see ``CHANNEL_GROUPS``) as an approximately independent witness. Within a
group the strongest reliability-discounted hit is taken; channels in the same
group are not counted twice (foldseek and pocket both rely on structure, for
example). Group strengths are then combined with a noisy-OR::

    p_g    = max over channels c in group g of (reliability_c * score_c)
    fused  = 1 - prod_g (1 - p_g)

A term backed by fewer than ``rank.require_independent_groups`` groups is
capped at ``rank.single_group_cap`` (medium confidence, for manual review),
because a single line of evidence, however strong, is not enough for automatic
high-confidence transfer in this second-stage predictor. Two strong,
orthogonal groups reach high confidence.

The legacy ``weighted_sum`` method retains ``raw = sum_c weight_c *
max(score_c for this term)`` and caps terms supported by too few independent
groups below ``rank.high_confidence``; it is kept for comparison.

A learned fusion model (``--weights``) overrides either rule method.

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
from pufprobe.ontology_names import pfam_term_name
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

# Evidence groups that provide indirect support (guilt by association). They
# can raise a score through noisy-OR but are not counted toward
# require_independent_groups, so they cannot by themselves unlock high
# confidence. High-confidence transfer needs two direct evidence groups.
INDIRECT_GROUPS = frozenset({"context"})

_GO_RE = re.compile(r"^GO:\d{7}$")
_EC_RE = re.compile(r"^(EC:)?\d+\.\d+\.\d+\.[\dR]+(?:\.[A-Z])?$")
_PFAM_RE = re.compile(r"^PF\d{5}")
_COG_RE = re.compile(r"^COG\d{4}")
_PDB_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}(?:[._-][A-Za-z0-9]+)?$")
# A bare PDB-chain template label (4fqe, 6jbu_B, 9iqf_B_1, ...); chain and
# domain suffixes may repeat. Used only to tell structure-template labels
# apart from human-readable functional names.
_PDB_LABEL_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}(?:[._-][A-Za-z0-9]+)*$")

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


def _term_label(
    evs: List[Evidence],
    ttype: str,
    term_id: str,
    pfam_names: Optional[Dict[str, Dict[str, str]]] = None,
) -> str:
    """Pick the most descriptive label for a term.

    For a Pfam term the official family description from the Pfam name index
    is the term's own name. A Foldseek hit carries the PDB template, and a
    sequence hit may carry a distant homolog, sometimes from a eukaryote, so
    neither should name the family. The official Pfam description is used
    whenever it is available. When the index is absent we fall back to the
    best non-PDB human-readable label and then the accession. An unmapped
    structure term (``ttype`` of ``structure`` or ``PDB``) keeps the PDB
    label, because the template is itself the result.
    """
    if ttype == "Pfam":
        official = pfam_term_name(pfam_names or {}, term_id)
        if official:
            return official
    best = max(evs, key=lambda e: e.score)
    default = (best.target_label or best.target_id or term_id).strip()
    if ttype in ("structure", "PDB") or not _PDB_LABEL_RE.match(default):
        return default
    for ev in sorted(evs, key=lambda e: e.score, reverse=True):
        label = (ev.target_label or "").strip()
        if label and not _PDB_LABEL_RE.match(label):
            return label
    return term_id.strip()


def _group_strengths(
    best_per_channel: Dict[str, float],
    reliability: Dict[str, float],
) -> Dict[str, float]:
    """Strongest reliability-discounted strength per independent evidence group.

    Channels inside one group are not independent (foldseek and pocket both
    rest on the predicted structure), so the group is represented by its
    strongest discounted hit rather than a sum.
    """
    strength: Dict[str, float] = {}
    for ch, score in best_per_channel.items():
        group = CHANNEL_GROUPS[ch]
        discounted = float(reliability.get(ch, 0.5)) * float(score)
        strength[group] = max(strength.get(group, 0.0), discounted)
    return strength


def _noisy_or(strengths: List[float]) -> float:
    """Probability that at least one independent witness supports the term."""
    miss = 1.0
    for s in strengths:
        miss *= max(0.0, 1.0 - min(max(s, 0.0), 1.0))
    return min(max(1.0 - miss, 0.0), 1.0)


def _rule_fusion(
    best_per_channel: Dict[str, float],
    gevs: List[Evidence],
    rank_cfg: Dict[str, Any],
) -> tuple[float, str]:
    """Fuse one term's channel maxima under the configured rule method."""
    required_groups = int(rank_cfg["require_independent_groups"])
    method = str(rank_cfg.get("fusion_method", "group_noisy_or"))

    all_groups = {CHANNEL_GROUPS[ev.channel.value] for ev in gevs}
    direct_groups = sorted(all_groups - INDIRECT_GROUPS)
    indirect_present = sorted(all_groups & INDIRECT_GROUPS)

    reliability: Dict[str, float] = rank_cfg.get("channel_reliability", {})
    single_group_cap = float(rank_cfg.get("single_group_cap", 0.60))
    group_min_strength = float(rank_cfg.get("group_min_strength", 0.40))
    group_strength = _group_strengths(best_per_channel, reliability)
    # Only direct groups whose strength clears the gate count as independent
    # witnesses. Weaker groups still raise the noisy-OR score, but a near-noise
    # group cannot unlock high confidence on its own.
    strong_direct = sorted(
        g for g in direct_groups if group_strength.get(g, 0.0) >= group_min_strength
    )
    weak_direct = sorted(
        g for g in direct_groups if group_strength.get(g, 0.0) < group_min_strength
    )

    def _cap_note(capped: float) -> str:
        named = ", ".join(strong_direct) if strong_direct else "none"
        note = (
            f"; capped at {capped:.3f}: only {len(strong_direct)} strong direct "
            f"group(s) ({named}) at gate {group_min_strength:.2f}, need "
            f"{required_groups} for high confidence"
        )
        if weak_direct:
            note += (
                f"; weaker group(s) ({', '.join(weak_direct)}) add to the score "
                f"but do not count as independent witnesses"
            )
        if indirect_present:
            note += f"; indirect group(s) {', '.join(indirect_present)} do not count"
        return note

    if method == "weighted_sum":
        rule_weights: Dict[str, float] = rank_cfg["rule_weights"]
        high_conf = float(rank_cfg["high_confidence"])
        raw = sum(rule_weights.get(ch, 0.0) * score for ch, score in best_per_channel.items())
        raw = min(max(raw, 0.0), 1.0)
        rationale = "score: weighted sum of channel maxima (rule weights)"
        if len(strong_direct) < required_groups:
            capped = min(raw, high_conf - 0.01)
            rationale += _cap_note(capped)
            raw = capped
        return raw, rationale

    # Default: group-level noisy-OR.
    ordered = sorted(group_strength)
    fused = _noisy_or([group_strength[g] for g in ordered])
    detail = ", ".join(f"{g}={group_strength[g]:.3f}" for g in ordered)
    rationale = (
        f"score: noisy-OR over {len(ordered)} evidence group(s) "
        f"[{detail}] = {fused:.3f}"
    )
    if len(strong_direct) < required_groups:
        capped = min(fused, single_group_cap)
        rationale += _cap_note(capped)
        fused = capped
    return fused, rationale


def fuse_proteins(
    records: List[CDSRecord],
    evidence_map: Dict[str, List[Evidence]],
    config: Dict[str, Any],
    weights_model: Any = None,
    pfam_names: Optional[Dict[str, Dict[str, str]]] = None,
) -> List[ProteinPrediction]:
    """Fuse evidence for every CDS record into ranked terms."""
    rank_cfg = config["rank"]

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
                            term_name=_term_label(gevs, ttype, tid, pfam_names),
                            raw_score=raw,
                            confidence=raw,
                            evidence=list(gevs),
                            rationale=rationale,
                        )
                    )
                    continue
            raw, rationale = _rule_fusion(best_per_channel, gevs, rank_cfg)

            terms.append(
                TermPrediction(
                    term_type=ttype,
                    term_id=tid,
                    term_name=_term_label(gevs, ttype, tid, pfam_names),
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
