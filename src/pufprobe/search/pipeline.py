"""search stage pipeline.

Contract implemented here:

1. Channels default to config search.channels and may be overridden with
   --channels. Supported channels: mmseqs, foldseek, hhblits, tmvec.
2. Each channel first checks pufprobe.externals for its binary (and, for
   foldseek, that 04_fold contains usable structures). A missing binary or
   database disables that channel with a logged message and install guidance;
   other channels still run.
3. Raw outputs are parsed into pufprobe.models.Evidence objects:
   - mmseqs and hhblits hits -> EvidenceChannel.REMOTE_HOMOLOGY
   - foldseek hits -> EvidenceChannel.FOLDSEEK
   - tmvec hits -> EvidenceChannel.FOLDSEEK (structure-derived)
   Target labels carry the database description/name when available.

Outputs:
- 05_search/<channel>_raw.tsv
- 05_search/search_evidence.json   global key -> list of Evidence dicts
Completion marker written for search.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from pufprobe import externals
from pufprobe.io_utils import read_json, write_json
from pufprobe.logging_utils import get_logger
from pufprobe.models import Evidence, families_from_dicts, to_dict
from pufprobe.search import foldseek as foldseek_channel
from pufprobe.search import hhblits as hhblits_channel
from pufprobe.search import idmapping
from pufprobe.search import mmseqs_search as mmseqs_channel
from pufprobe.workdir import WorkDir

logger = get_logger(__name__)


def run(
    workdir: WorkDir,
    config: Dict[str, Any],
    channels: Optional[List[str]] = None,
    dbdir: Optional[Path] = None,
) -> None:
    families_path = workdir.path("cluster", "families.json")
    reps_faa = workdir.path("cluster", "representatives.faa")
    if not families_path.exists() or not reps_faa.exists():
        raise FileNotFoundError("Run `pufprobe embed` before `pufprobe search`.")

    families = families_from_dicts(read_json(families_path))  # type: ignore[arg-type]
    search_cfg = config["search"]
    fold_cfg = config["fold"]
    channels = channels or list(search_cfg["channels"])
    out_dir = workdir.stage_dir("search")
    structures_dir = workdir.stage_dir("fold") / "structures"

    merged: Dict[str, List[Evidence]] = {}

    for channel in channels:
        channel = channel.strip().lower()
        if channel == "mmseqs":
            evidence, _ = mmseqs_channel.run_mmseqs_search(
                reps_faa, dbdir, out_dir,
                threads=int(workdir.threads),
                report_evalue=float(search_cfg["report_evalue"]),
            )
        elif channel == "foldseek":
            sifts_path = dbdir / "pdb_chain_pfam.tsv.gz" if dbdir is not None else None
            evidence, _ = foldseek_channel.run_foldseek_search(
                families, structures_dir, dbdir, out_dir,
                min_plddt=float(fold_cfg["min_plddt"]),
                sensitivity=float(search_cfg["foldseek_sensitivity"]),
                threads=int(workdir.threads),
                sifts_path=sifts_path,
            )
        elif channel == "hhblits":
            evidence, _ = hhblits_channel.run_hhblits_search(
                reps_faa, dbdir, out_dir,
                iterations=int(search_cfg["hhblits_iterations"]),
                threads=int(workdir.threads),
            )
        elif channel == "tmvec":
            info = externals.tmvec()
            if not info.available:
                logger.info("tmvec not installed; skipping optional tmvec channel.")
                continue
            logger.info("tmvec detected but not wired in this build; channel skipped.")
            evidence = {}
        else:
            logger.warning("Unknown search channel '%s'; ignored.", channel)
            continue

        for key, ev_list in evidence.items():
            merged.setdefault(key, []).extend(ev_list)
        logger.info("Channel '%s' contributed %d query(s) / %d hits.",
                    channel, len(evidence), sum(len(v) for v in evidence.values()))

    # Optional term normalization: map Swiss-Prot accessions onto Pfam/EC/GO so
    # that independent evidence groups can converge on the same functional term.
    xref_path = dbdir / "swissprot_xref.json" if dbdir is not None else None
    xref = idmapping.load_xref(xref_path)
    if xref:
        merged, n_added = idmapping.expand_evidence(merged, xref)
        logger.info(
            "Term normalization added %d ontology evidence entries (Pfam/EC/GO).",
            n_added,
        )

    payload = {key: [to_dict(ev) for ev in evs] for key, evs in merged.items()}
    write_json(out_dir / "search_evidence.json", payload)
    workdir.mark_complete("search", payload={
        "channels": channels,
        "n_queries": len(merged),
        "n_hits": sum(len(v) for v in merged.values()),
    })
    logger.info("search done: %d queries / %d evidence entries.",
                len(merged), sum(len(v) for v in merged.values()))
