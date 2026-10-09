"""fold stage pipeline.

Contract implemented here:

1. Load families from 03_cluster/families.json. Only representative sequences
   are folded; members inherit results through family assignment.
2. Structures come from --precomputed (a directory of <family_id>.pdb files)
   or from ESMFold through pufprobe.fold.esmfold when torch and weights are
   available. The af2 engine is an accepted option but is not bundled; when
   selected without a local AlphaFold installation the family is marked
   skipped with the reason recorded.
3. Apply the pLDDT quality gate (config fold.min_plddt): families below the
   threshold keep fold_status low_confidence and their structure is not used
   by structure-search channels.

Outputs:
- 04_fold/structures/<family_id>.pdb
- 04_fold/fold_quality.tsv  family_id, engine, status, median_plddt, mean_pae
Family fields structure_path, plddt and fold_status are updated and
families.json is rewritten. Completion marker written for fold.

Every unavailable backend is logged and skipped, never raised.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from pufprobe.fold.esmfold import ModelUnavailable, fold_sequence, parse_plddt_from_pdb
from pufprobe.io_utils import read_json, write_json
from pufprobe.logging_utils import get_logger
from pufprobe.models import Family, families_from_dicts, to_dict
from pufprobe.workdir import WorkDir

logger = get_logger(__name__)


def run(
    workdir: WorkDir,
    config: Dict[str, Any],
    engine: Optional[str] = None,
    precomputed: Optional[Path] = None,
) -> None:
    families_path = workdir.path("cluster", "families.json")
    if not families_path.exists():
        raise FileNotFoundError(f"Missing {families_path}; run `pufprobe embed` first.")

    families: List[Family] = families_from_dicts(read_json(families_path))  # type: ignore[arg-type]
    fold_cfg = config["fold"]
    engine = engine or fold_cfg["engine"]
    min_plddt = float(fold_cfg["min_plddt"])
    max_len = int(fold_cfg["max_length_esmfold"])

    fold_dir = workdir.stage_dir("fold")
    struct_dir = fold_dir / "structures"
    struct_dir.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, Any]] = []
    for fam in families:
        plddt: Optional[float] = None
        mean_pae: Optional[float] = None
        pdb_text: Optional[str] = None
        status = "skipped"
        engine_used = engine

        if precomputed is not None:
            pdb_path = precomputed / f"{fam.family_id}.pdb"
            if pdb_path.exists():
                pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
                plddt = parse_plddt_from_pdb(pdb_text)
                engine_used = "precomputed"
            else:
                logger.warning("Precomputed PDB missing for %s; marked skipped.", fam.family_id)
        elif engine == "af2":
            logger.warning(
                "AlphaFold2 engine is not bundled; %s marked skipped. Install a "
                "local AlphaFold and place <family_id>.pdb under --precomputed.",
                fam.family_id,
            )
        elif engine == "esmfold":
            if len(fam.representative_seq) > max_len:
                logger.warning(
                    "%s representative length %d exceeds max_length_esmfold=%d; skipped.",
                    fam.family_id, len(fam.representative_seq), max_len,
                )
            else:
                try:
                    pdb_text, plddt, mean_pae = fold_sequence(fam.representative_seq)
                except ModelUnavailable as exc:
                    logger.warning("ESMFold unavailable for %s (%s); marked failed.",
                                   fam.family_id, exc)
                    status = "failed"
                except Exception as exc:  # noqa: BLE001 - per-family isolation
                    logger.warning("Folding failed for %s (%s); marked failed.", fam.family_id, exc)
                    status = "failed"
        else:
            logger.warning("Unknown fold engine '%s'; %s marked skipped.", engine, fam.family_id)

        # quality gate
        if pdb_text is not None and plddt is not None:
            status = "folded" if plddt >= min_plddt else "low_confidence"

        structure_rel: Optional[str] = None
        if pdb_text is not None:
            out = struct_dir / f"{fam.family_id}.pdb"
            out.write_text(pdb_text, encoding="utf-8", newline="\n")
            structure_rel = str(out.relative_to(workdir.root))

        fam.structure_path = structure_rel
        fam.plddt = plddt
        fam.fold_status = status

        rows.append({
            "family_id": fam.family_id,
            "engine": engine_used,
            "status": status,
            "median_plddt": plddt,
            "mean_pae": mean_pae,
        })

    _write_quality_tsv(fold_dir / "fold_quality.tsv", rows)
    write_json(families_path, [to_dict(f) for f in families])
    workdir.mark_complete("fold", payload={"n_families": len(families), "engine": engine})
    logger.info("fold done: %d families (engine=%s).", len(families), engine)


def _write_quality_tsv(path: Path, rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("family_id\tengine\tstatus\tmedian_plddt\tmean_pae\n")
        for row in rows:
            handle.write(
                f"{row['family_id']}\t{row['engine']}\t{row['status']}\t"
                f"{_fmt(row['median_plddt'])}\t{_fmt(row['mean_pae'])}\n"
            )


def _fmt(value: Optional[float]) -> str:
    return "" if value is None else f"{value:.2f}"
