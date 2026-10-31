"""rank stage pipeline.

Contract implemented here:

1. Load CDS records, families, search_evidence.json and context_evidence.json.
2. Collect candidate terms per protein from evidence. Fuse channel scores:
   - when --weights points to a trained fusion model it is loaded lazily
     (sklearn joblib/pickle) and used for learning-to-rank scoring;
   - otherwise config rank.rule_weights give transparent rule scores, and the
     output notes that default rule weights rather than learned weights were
     used.
3. Apply the independent-channel agreement policy: a term needs evidence from
   at least rank.require_independent_groups distinct evidence groups to reach
   high confidence; single-channel terms are capped below the threshold.
4. Calibrate raw scores through pufprobe.rank.calibrate (temperature scaling
   and isotonic regression; sklearn used when present, pool-adjacent-violators
   fallback otherwise). Without a held-out calibration set the calibration is
   identity and is labelled as such.
5. For families with no anchored term, pufprobe.rank.priority computes a
   wet-lab priority score from structural novelty, predicted pockets, clonal
   specificity, AMR/virulence neighbourhood and small-protein status.

Outputs:
- 07_rank/predictions.json  list of ProteinPrediction dicts
- 07_rank/calibration.json  method, ECE placeholder and per-ontology notes
- 07_rank/priorities.tsv    unanchored family ranking
Completion marker written for rank.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional

from pufprobe.io_utils import read_json, write_json
from pufprobe.logging_utils import get_logger
from pufprobe.ontology_graph import load_go_obo
from pufprobe.ontology_names import load_pfam_names
from pufprobe.models import (
    Family,
    Sample,
    families_from_dicts,
    to_dict,
)
from pufprobe.rank.calibrate import (
    IdentityCalibrator,
    IsotonicCalibrator,
    TemperatureScaling,
    build_calibrator,
    load_calibrator,
)
from pufprobe.rank.fuse import build_evidence_map, cds_record_from_row, fuse_proteins, load_weights_model
from pufprobe.rank.priority import compute_family_priority
from pufprobe.workdir import WorkDir

log = get_logger("rank.pipeline")


def run(
    workdir: WorkDir,
    config: Dict[str, Any],
    weights: Optional[Path] = None,
    dbdir: Optional[Path] = None,
) -> None:
    records_path = workdir.path("prepare", "cds_records.json")
    samples_path = workdir.path("prepare", "samples.json")
    families_path = workdir.path("cluster", "families.json")
    search_ev_path = workdir.path("search", "search_evidence.json")
    context_ev_path = workdir.path("context", "context_evidence.json")

    records = [cds_record_from_row(row) for row in read_json(records_path)]  # type: ignore[arg-type]
    samples: List[Sample] = (
        [Sample(**row) for row in read_json(samples_path)]  # type: ignore[arg-type]
        if samples_path.exists()
        else []
    )
    families: List[Family] = (
        families_from_dicts(read_json(families_path))  # type: ignore[arg-type]
        if families_path.exists()
        else []
    )

    evidence_map = build_evidence_map(search_ev_path, context_ev_path)
    log.info(
        "evidence loaded: %d proteins with evidence "
        "(search=%s, context=%s)",
        sum(1 for evs in evidence_map.values() if evs),
        search_ev_path.name if search_ev_path.exists() else "missing",
        context_ev_path.name if context_ev_path.exists() else "missing",
    )

    pfam_names_path = dbdir / "pfam_names.json" if dbdir is not None else None
    pfam_names = load_pfam_names(pfam_names_path)
    if pfam_names:
        log.info(
            "Pfam family name index loaded: %d families (%s)",
            len(pfam_names), pfam_names_path,
        )
    else:
        log.info(
            "Pfam family name index not found (%s); Pfam terms use evidence labels.",
            pfam_names_path,
        )

    # Gene Ontology true-path propagation. The configured path wins; otherwise
    # go-basic.obo (or a gzip copy) is sought in the database directory.
    go_ontology = load_go_obo(_resolve_go_obo_path(config, dbdir))
    go_propagated = go_ontology is not None
    n_go_propagated = 0
    if go_ontology is not None:
        evidence_map, n_go_propagated = go_ontology.propagate_evidence(evidence_map)
        log.info(
            "GO true-path propagation applied: %d ancestor-term evidence entries added.",
            n_go_propagated,
        )

    # Learning-to-rank fusion (optional). Configured by rank.ltr_model_path;
    # a relative path is resolved against the database directory. Any failure
    # falls back to rule fusion without interrupting the run.
    ltr_model = None
    ltr_path_cfg = str(config["rank"].get("ltr_model_path", "") or "")
    if ltr_path_cfg:
        ltr_path = Path(ltr_path_cfg)
        if not ltr_path.is_absolute() and dbdir is not None:
            ltr_path = dbdir / ltr_path
        if ltr_path.exists():
            try:
                from pufprobe.rank.ltr import LTRRanker
                ltr_model = LTRRanker.load(str(ltr_path))
                log.info("LTR fusion model loaded: %s (backend=%s)",
                         ltr_path, ltr_model.backend)
            except Exception as exc:  # pragma: no cover - defensive
                log.warning("could not load LTR model %s (%s); using rule fusion.",
                            ltr_path, exc)
        else:
            log.warning("configured LTR model not found: %s; using rule fusion.",
                        ltr_path)

    # Post-hoc confidence calibration (optional). Configured by
    # rank.calibration_model_path with the same relative-path and fallback
    # behaviour as the LTR model. Without a fitted calibrator the identity
    # calibrator leaves raw scores unchanged.
    fitted_calibrator = None
    cal_path_cfg = str(config["rank"].get("calibration_model_path", "") or "")
    if cal_path_cfg:
        cal_path = Path(cal_path_cfg)
        if not cal_path.is_absolute() and dbdir is not None:
            cal_path = dbdir / cal_path
        if cal_path.exists():
            try:
                fitted_calibrator = load_calibrator(str(cal_path))
                log.info("Confidence calibrator loaded: %s (%s).",
                         cal_path, type(fitted_calibrator).__name__)
            except Exception as exc:  # pragma: no cover - defensive
                log.warning("could not load calibrator %s (%s); using identity.",
                            cal_path, exc)
        else:
            log.warning("configured calibrator not found: %s; using identity.",
                        cal_path)

    model, model_note = load_weights_model(weights)
    if ltr_model is not None:
        fusion_note = f"fusion: learned LTR ranker ({ltr_model.backend}) from {ltr_path_cfg}"
    elif model is None:
        method = config["rank"].get("fusion_method", "group_noisy_or")
        fusion_note = f"fusion: {method} over independent evidence groups (rule-based, no learned model)"
    else:
        fusion_note = model_note
    predictions = fuse_proteins(
        records, evidence_map, config,
        weights_model=model, pfam_names=pfam_names,
        go_propagated=go_propagated,
        ltr_model=ltr_model,
    )

    # Attach GO depth so the annotate stage can skip shallow, generic terms
    # when building product names. Non-GO terms and runs without an ontology
    # keep go_depth None, which the backfill treats as "do not filter".
    if go_ontology is not None:
        for pred in predictions:
            for term in pred.terms:
                if term.term_type == "GO":
                    term.go_depth = go_ontology.depth(term.term_id)

    # Use the fitted calibrator when one loaded; otherwise the identity
    # calibrator keeps the raw scores. No ground truth is available at
    # inference, so ECE comes from the held-out fit rather than this stage.
    if fitted_calibrator is None:
        calibrator, cal_meta = build_calibrator(None)
    else:
        calibrator = fitted_calibrator
        if isinstance(calibrator, TemperatureScaling):
            kind = "temperature"
        elif isinstance(calibrator, IsotonicCalibrator):
            kind = "isotonic"
        else:
            kind = "identity"
        cal_meta = {
            "method": kind,
            "ece": None,
            "note": f"loaded fitted calibrator from {cal_path_cfg}",
        }
    for pred in predictions:
        raw_scores = [t.raw_score for t in pred.terms]
        calibrated = calibrator.transform(raw_scores)
        for term, conf in zip(pred.terms, calibrated):
            term.confidence = float(conf)

    # Families with no anchored term get a wet-lab priority score.
    anchored_keys = {pred.key for pred in predictions if pred.terms}
    total_samples = max(len(samples), len({r.sample_id for r in records}), 1)
    priority_rows = []
    for family in families:
        if any(member in anchored_keys for member in family.members):
            continue
        rep_evs = evidence_map.get(family.representative_key, [])
        score, reasons = compute_family_priority(family, rep_evs, total_samples, config)
        family.priority_score = score
        family.priority_reasons = reasons
        priority_rows.append((family.family_id, score, reasons))
    priority_rows.sort(key=lambda row: row[1], reverse=True)

    write_json(workdir.path("rank", "predictions.json"), [to_dict(p) for p in predictions])

    calibration_payload = {
        **cal_meta,
        "fusion": fusion_note,
        "fusion_method": config["rank"].get("fusion_method", "group_noisy_or"),
        "ltr_backend": getattr(ltr_model, "backend", None),
        "high_confidence": config["rank"]["high_confidence"],
        "single_group_cap": config["rank"].get("single_group_cap"),
        "require_independent_groups": config["rank"]["require_independent_groups"],
        "n_predictions": len(predictions),
        "n_terms": sum(len(p.terms) for p in predictions),
        "n_unanchored_families": len(priority_rows),
        "go_true_path_propagation": "applied" if go_propagated else "disabled",
        "n_go_propagated": n_go_propagated,
    }
    write_json(workdir.path("rank", "calibration.json"), calibration_payload)

    priorities_path = workdir.path("rank", "priorities.tsv")
    with priorities_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("family_id\tpriority_score\treasons\n")
        for fid, score, reasons in priority_rows:
            handle.write(f"{fid}\t{score:.4f}\t{'; '.join(reasons)}\n")

    workdir.mark_complete("rank")
    log.info(
        "rank done: %d predictions, %d terms, %d unanchored families prioritised.",
        len(predictions),
        sum(len(p.terms) for p in predictions),
        len(priority_rows),
    )


def _resolve_go_obo_path(config: Dict[str, Any], dbdir: Optional[Path]) -> Optional[Path]:
    """Locate go-basic.obo.

    An explicit ``rank.go_obo_path`` takes precedence. Otherwise the plain
    and gzip-compressed names are probed in the database directory. Returns
    None when no ontology is configured or present.
    """
    configured = str(config.get("rank", {}).get("go_obo_path", "") or "").strip()
    if configured:
        return Path(configured)
    if dbdir is None:
        return None
    for name in ("go-basic.obo", "go-basic.obo.gz"):
        candidate = dbdir / name
        if candidate.exists():
            return candidate
    return None
