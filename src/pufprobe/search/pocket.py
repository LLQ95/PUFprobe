"""Pocket catalytic-template evidence channel (P2Rank + M-CSA).

The channel asks whether a predicted binding pocket of the query structure
matches a catalytic site curated in the Mechanism and Catalytic Site Atlas
(M-CSA). An offline build step (see ``server/build_pocket_templates.py``)
turns every M-CSA reference structure into a single-chain PDB and a foldseek
database, and records each template's catalytic residues. At search time the
query structure is scanned with P2Rank and aligned against the template
database with foldseek. For every template hit, the foldseek alignment maps
the template catalytic residues onto the query sequence, and the mapped
residues are checked against the P2Rank pockets. A hit is admitted only when
at least one catalytic residue lands inside a predicted pocket.

Each M-CSA entry is associated with a single EC number and a single reference
structure, so an accepted hit becomes ``Evidence(channel=POCKET)`` that votes
on the EC term. The channel belongs to the ``structure`` evidence group, the
same group as foldseek, because both rest on the predicted structure; it is
not counted as an independent witness in fusion. When P2Rank, foldseek, the
template database or a usable query structure is missing the channel is
skipped (logged) rather than raised.
"""

from __future__ import annotations

import csv
import math
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from pufprobe import externals
from pufprobe.logging_utils import get_logger
from pufprobe.models import Evidence, EvidenceChannel, Family

logger = get_logger(__name__)

INSTALL_HINT = "See docs/server_setup.md for building the M-CSA pocket template database."
_TEMPLATE_NAME = "pocket_templates.json"
_POCKET_DB_NAME = "pocket_db"

# foldseek convertalis columns. qaln/taln carry the gapped alignment and need
# the search backtrace, so the search runs with -a. qstart/tstart are the
# 1-based positions of the first aligned residue; qaln/taln are the gapped
# sequences and share one column per aligned position.
COLS = "query,target,qstart,tstart,qaln,taln,alntmscore,prob,evalue,qcov,tcov"
_NCOL = 11

# Structural quality gates, aligned with the foldseek channel: a TM-score
# around 0.5 is the same-fold boundary, and the E-value states significance.
DEFAULT_MIN_TMSCORE = 0.5
DEFAULT_MAX_EVALUE = 1e-2
DEFAULT_SENSITIVITY = 9.5
# A template needs at least this fraction of its aligned catalytic residues
# inside a single predicted pocket, and only the top-ranked pockets are
# considered to limit false positives.
DEFAULT_MIN_POCKET_CONSISTENCY = 0.5
DEFAULT_MAX_POCKETS = 5

_RESIDUE_TOKEN_RE = re.compile(r"^([A-Za-z0-9])_(\d+)([A-Za-z]?)$")


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class PdbResidue:
    """One residue in a PDB file, with chain-local sequence position."""

    chain: str
    auth_resid: int
    icode: str
    resname: str
    seq_index: int


@dataclass
class PocketPrediction:
    """One P2Rank pocket with the author-numbered residues it contains."""

    name: str
    rank: int
    score: float
    probability: float
    residue_keys: Set[Tuple[str, int, str]] = field(default_factory=set)


@dataclass
class TemplateMatch:
    """Evaluated match between a query structure and one pocket template."""

    target: str
    mcsa_id: int
    ec: str
    enzyme_name: str
    pocket_rank: int
    tmscore: float
    prob: float
    evalue: float
    qcov: float
    tcov: float
    catalytic_total: int
    catalytic_mapped: int
    catalytic_in_pocket: int
    consistency: float
    structure_score: float
    score: float


# ---------------------------------------------------------------------------
# PDB parsing (shared with the offline template build)
# ---------------------------------------------------------------------------

def parse_pdb_residues(pdb_text: str) -> List[PdbResidue]:
    """Parse ATOM records into per-chain residues with 1-based sequence positions.

    Residues are indexed in file order independently for each chain, which is
    the position foldseek assigns when it reads a PDB. The PDB author number
    (``auth_resid``) may contain gaps and is kept separately.
    """
    residues: List[PdbResidue] = []
    chain_count: Dict[str, int] = {}
    seen: Set[Tuple[str, int, str]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        # PDB fixed columns: resname 18-20, chain 22, resid 23-26, icode 27
        # (0-based slices 17:20, 21, 22:26, 26).
        resname = line[17:20].strip()
        chain = line[21].strip()
        icode = line[26].strip()
        try:
            auth_resid = int(line[22:26])
        except ValueError:
            continue
        key = (chain, auth_resid, icode)
        if key in seen:
            continue  # multiple atoms of the same residue
        seen.add(key)
        seq_index = chain_count.get(chain, 0) + 1
        chain_count[chain] = seq_index
        residues.append(
            PdbResidue(
                chain=chain,
                auth_resid=auth_resid,
                icode=icode,
                resname=resname,
                seq_index=seq_index,
            )
        )
    return residues


def build_residue_key_map(
    residues: List[PdbResidue],
) -> Dict[Tuple[str, int, str], int]:
    """Map author-numbered residue keys to chain-local sequence positions."""
    return {
        (r.chain, r.auth_resid, r.icode): r.seq_index for r in residues
    }


def index_by_seq(residues: List[PdbResidue]) -> Dict[int, PdbResidue]:
    """Map chain-local sequence positions to residues.

    Fold outputs are single-chain structures, so the sequence position is
    unique across the parsed residues.
    """
    return {r.seq_index: r for r in residues}


def extract_chain_pdb(pdb_text: str, chain: str) -> str:
    """Extract one chain's ATOM and TER records into a single-chain PDB."""
    kept: List[str] = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "TER")) and line[21].strip() == chain:
            kept.append(line)
    if kept and not kept[-1].startswith("TER"):
        kept.append(_ter_record(chain))
    return "\n".join(kept) + "\n"


def _ter_record(chain: str) -> str:
    return "TER   {:>5}      {:>3} {:>1}".format("", "", chain)


# ---------------------------------------------------------------------------
# P2Rank output parsing
# ---------------------------------------------------------------------------

def parse_residue_token(token: str) -> Optional[Tuple[str, int, str]]:
    """Parse a P2Rank residue id (``A_123`` or ``A_123A``).

    P2Rank prints the chain, an underscore and the PDB author residue number,
    followed by the insertion code when one exists.
    """
    match = _RESIDUE_TOKEN_RE.match(token.strip())
    if match is None:
        return None
    chain, resid_s, icode = match.groups()
    return chain, int(resid_s), icode


def parse_p2rank_predictions(csv_path: Path) -> List[PocketPrediction]:
    """Parse a P2Rank ``predictions.csv`` into pocket records."""
    pockets: List[PocketPrediction] = []
    if not csv_path.exists():
        return pockets
    with csv_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                rank = int(row["rank"])
                score = float(row["score"])
                probability = float(row["probability"])
            except (KeyError, ValueError):
                continue
            residue_keys: Set[Tuple[str, int, str]] = set()
            for token in (row.get("residue_ids") or "").split():
                parsed = parse_residue_token(token)
                if parsed is not None:
                    residue_keys.add(parsed)
            pockets.append(
                PocketPrediction(
                    name=(row.get("name") or "").strip(),
                    rank=rank,
                    score=score,
                    probability=probability,
                    residue_keys=residue_keys,
                )
            )
    pockets.sort(key=lambda p: p.rank)
    return pockets


def pocket_seq_sets(
    pockets: List[PocketPrediction],
    key_map: Dict[Tuple[str, int, str], int],
    max_pockets: int,
) -> List[Tuple[PocketPrediction, Set[int]]]:
    """Translate each pocket's author-numbered keys to sequence positions.

    Only the top-ranked pockets are returned. Residue keys that cannot be
    located in the query structure (for example a chain P2Rank did not load)
    are dropped rather than treated as misses.
    """
    out: List[Tuple[PocketPrediction, Set[int]]] = []
    for pocket in pockets[:max_pockets]:
        seq_set = {
            key_map[key] for key in pocket.residue_keys if key in key_map
        }
        out.append((pocket, seq_set))
    return out


# ---------------------------------------------------------------------------
# foldseek alignment mapping and scoring
# ---------------------------------------------------------------------------

def map_alignment(
    qaln: str, taln: str, qstart: int, tstart: int
) -> Dict[int, int]:
    """Map target to query sequence positions from a gapped alignment.

    Both returned positions are 1-based. A column with a gap on either side
    contributes no mapping, so residues aligned only against a gap are left
    unmapped.
    """
    mapping: Dict[int, int] = {}
    qpos = int(qstart)
    tpos = int(tstart)
    for q_char, t_char in zip(qaln, taln):
        q_cur: Optional[int] = None
        t_cur: Optional[int] = None
        if q_char != "-":
            q_cur = qpos
            qpos += 1
        if t_char != "-":
            t_cur = tpos
            tpos += 1
        if q_cur is not None and t_cur is not None:
            mapping[t_cur] = q_cur
    return mapping


def evaluate_template(
    template: dict,
    mapping: Dict[int, int],
    query_by_index: Dict[int, PdbResidue],
    pocket_sets: List[Tuple[PocketPrediction, Set[int]]],
    tmscore: float,
    prob: float,
    evalue: float,
    qcov: float,
    tcov: float,
    min_pocket_consistency: float,
) -> Optional[TemplateMatch]:
    """Evaluate one aligned template against the query pockets.

    Catalytic residues are kept only when the alignment maps them to a query
    residue of the same amino acid. The pocket with the most mapped catalytic
    residues is selected; ties resolve to the lower rank. Returns ``None`` when
    no catalytic residue is aligned, lands in a pocket, or the consistency
    gate is not met.
    """
    catalytic = template.get("catalytic_residues", [])
    total = len(catalytic)
    mapped_positions: List[int] = []
    for residue in catalytic:
        try:
            target_resid = int(residue["resid"])
            code = str(residue["code"])
        except (KeyError, ValueError, TypeError):
            continue
        query_resid = mapping.get(target_resid)
        if query_resid is None:
            continue
        query_res = query_by_index.get(query_resid)
        if query_res is None:
            continue
        if query_res.resname.upper() != code.upper():
            continue
        mapped_positions.append(query_resid)

    mapped = len(mapped_positions)
    if total == 0 or mapped == 0:
        return None

    best_pocket: Optional[PocketPrediction] = None
    best_hits = 0
    for pocket, seq_set in pocket_sets:
        hits = sum(1 for pos in mapped_positions if pos in seq_set)
        if hits > best_hits or (
            hits == best_hits
            and hits > 0
            and (best_pocket is None or pocket.rank < best_pocket.rank)
        ):
            best_hits = hits
            best_pocket = pocket

    if best_hits == 0 or best_pocket is None:
        return None

    consistency = best_hits / mapped
    if consistency < min_pocket_consistency:
        return None

    structure_score = _structure_score(tmscore, qcov, evalue)
    score = min(1.0, max(0.0, structure_score * consistency))
    return TemplateMatch(
        target=str(template.get("target", "")),
        mcsa_id=int(template.get("mcsa_id", 0)),
        ec=str(template.get("ec", "")),
        enzyme_name=str(template.get("enzyme_name", "")),
        pocket_rank=best_pocket.rank,
        tmscore=tmscore,
        prob=prob,
        evalue=evalue,
        qcov=qcov,
        tcov=tcov,
        catalytic_total=total,
        catalytic_mapped=mapped,
        catalytic_in_pocket=best_hits,
        consistency=consistency,
        structure_score=structure_score,
        score=score,
    )


def _evalue_score(evalue: float) -> float:
    if evalue <= 0.0:
        return 1.0
    return max(0.0, min(1.0, -math.log10(evalue) / 10.0))


def _structure_score(tmscore: float, qcov: float, evalue: float) -> float:
    """Function-transfer strength of the structural alignment.

    Same formula as the foldseek channel, combining TM-score, query coverage
    and E-value significance so the structure group carries graded strength.
    """
    fold = max(0.0, min(1.0, tmscore))
    coverage = max(0.0, min(1.0, qcov))
    significance = _evalue_score(evalue)
    return min(1.0, max(0.0, fold * coverage * significance))


# ---------------------------------------------------------------------------
# Online search orchestration
# ---------------------------------------------------------------------------

def run_pocket_search(
    families: List[Family],
    structures_dir: Path,
    dbdir: Optional[Path],
    out_dir: Path,
    min_plddt: float,
    sensitivity: float = DEFAULT_SENSITIVITY,
    threads: int = 1,
    min_tmscore: float = DEFAULT_MIN_TMSCORE,
    max_evalue: float = DEFAULT_MAX_EVALUE,
    min_pocket_consistency: float = DEFAULT_MIN_POCKET_CONSISTENCY,
    max_pockets: int = DEFAULT_MAX_POCKETS,
) -> Tuple[Dict[str, List[Evidence]], Optional[Path]]:
    """Run the pocket channel for every usable folded family."""
    p2rank_info = externals.p2rank()
    if not p2rank_info.available:
        logger.warning("P2Rank not found; skipping pocket channel.")
        return {}, None
    foldseek_info = externals.foldseek()
    if not foldseek_info.available:
        logger.warning("foldseek not found; skipping pocket channel.")
        return {}, None

    templates_path = dbdir / _TEMPLATE_NAME if dbdir else None
    pocket_db = dbdir / _POCKET_DB_NAME if dbdir else None
    if templates_path is None or not templates_path.exists():
        logger.warning("No %s at dbdir; skipping pocket channel. %s",
                       _TEMPLATE_NAME, INSTALL_HINT)
        return {}, None
    if pocket_db is None or not Path(str(pocket_db) + ".index").exists():
        logger.warning("No foldseek pocket database at dbdir/%s; skipping. %s",
                       _POCKET_DB_NAME, INSTALL_HINT)
        return {}, None

    templates = _load_templates(templates_path)
    if not templates:
        logger.warning("Empty pocket template index; skipping pocket channel.")
        return {}, None

    present = _select_query_structures(families, structures_dir, min_plddt)
    if not present:
        logger.warning("No usable folded structures; skipping pocket channel.")
        return {}, None

    work = Path(tempfile.mkdtemp(prefix="pocket_", dir=str(out_dir)))
    raw_tsv = out_dir / "pocket_raw.tsv"
    try:
        p2rank_dir = work / "p2rank"
        p2rank_dir.mkdir(parents=True, exist_ok=True)
        query_dir = work / "query"
        query_dir.mkdir(parents=True, exist_ok=True)

        query_pockets: Dict[str, List[PocketPrediction]] = {}
        query_residues: Dict[str, List[PdbResidue]] = {}
        for fam in present:
            src = structures_dir / f"{fam.family_id}.pdb"
            pdb_text = src.read_text(encoding="utf-8", errors="replace")
            (query_dir / f"{fam.family_id}.pdb").write_text(
                pdb_text, encoding="utf-8", newline="\n"
            )
            residues = parse_pdb_residues(pdb_text)
            query_residues[fam.family_id] = residues
            predictions_csv = _run_p2rank(
                p2rank_info.path or "p2rank",
                query_dir / f"{fam.family_id}.pdb",
                p2rank_dir / fam.family_id,
                threads,
            )
            if predictions_csv is not None:
                query_pockets[fam.family_id] = parse_p2rank_predictions(
                    predictions_csv
                )

        foldseek_tsv = _run_foldseek(
            foldseek_info.path or "foldseek",
            query_dir,
            pocket_db,
            work,
            raw_tsv,
            sensitivity,
            threads,
        )
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or "").strip()
        step = exc.cmd[1] if len(exc.cmd) > 1 else "?"
        logger.warning(
            "Pocket channel %s failed (exit %d)%s; channel skipped.",
            step, exc.returncode,
            (": " + detail[-1500:]) if detail else "",
        )
        return {}, None
    except OSError as exc:
        logger.warning("Pocket channel could not launch a dependency (%s); skipped.", exc)
        return {}, None

    evidence = parse_pocket_tsv(
        foldseek_tsv,
        templates,
        query_residues,
        query_pockets,
        family_representative={f.family_id: f.representative_key for f in present},
        min_tmscore=min_tmscore,
        max_evalue=max_evalue,
        min_pocket_consistency=min_pocket_consistency,
        max_pockets=max_pockets,
    )
    shutil.rmtree(work, ignore_errors=True)
    return evidence, raw_tsv


def _select_query_structures(
    families: List[Family], structures_dir: Path, min_plddt: float
) -> List[Family]:
    """Families with a present, quality-gated structure for pocket search."""
    out: List[Family] = []
    for fam in families:
        if not (structures_dir / f"{fam.family_id}.pdb").exists():
            continue
        if fam.fold_status == "folded":
            out.append(fam)
        elif fam.fold_status == "low_confidence" and fam.plddt is not None \
                and fam.plddt >= min_plddt:
            out.append(fam)
    return out


def _load_templates(path: Path) -> Dict[str, dict]:
    """Load the template index keyed by the foldseek target name."""
    from pufprobe.io_utils import read_json

    payload = read_json(path)
    out: Dict[str, dict] = {}
    if isinstance(payload, list):
        for entry in payload:
            if isinstance(entry, dict) and "target" in entry:
                out[str(entry["target"])] = entry
    elif isinstance(payload, dict):
        for key, entry in payload.items():
            if isinstance(entry, dict):
                entry = dict(entry)
                entry.setdefault("target", key)
                out[str(key)] = entry
    return out


def _run_p2rank(
    exe: str, pdb_path: Path, out_dir: Path, threads: int
) -> Optional[Path]:
    """Run P2Rank on one structure and locate its predictions CSV.

    Visualizations are disabled because the channel only consumes the tabular
    output and the server has no PyMOL display. The output label is not
    assumed, so the produced ``*_predictions.csv`` is located by glob.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        exe,
        "predict",
        "-f",
        str(pdb_path),
        "-o",
        str(out_dir),
        "-visualizations",
        "0",
        "-threads",
        str(threads),
    ]
    _run(cmd)
    matches = sorted(out_dir.glob("*_predictions.csv"))
    return matches[0] if matches else None


def _run_foldseek(
    exe: str,
    query_dir: Path,
    pocket_db: Path,
    work: Path,
    raw_tsv: Path,
    sensitivity: float,
    threads: int,
) -> Path:
    """Build the query database, search the pocket templates and convert output."""
    qdb = work / "qdb"
    res = work / "res"
    stmp = work / "stmp"
    _run([exe, "createdb", str(query_dir), str(qdb)])
    # -a stores the backtrace so convertalis can emit qaln/taln.
    _run([
        exe, "search", str(qdb), str(pocket_db), str(res), str(stmp),
        "-a", "-s", str(sensitivity), "--threads", str(threads),
    ])
    _run([
        exe, "convertalis", str(qdb), str(pocket_db), str(res), str(raw_tsv),
        "--format-output", COLS,
    ])
    return raw_tsv


def parse_pocket_tsv(
    tsv_path: Path,
    templates: Dict[str, dict],
    query_residues: Dict[str, List[PdbResidue]],
    query_pockets: Dict[str, List[PocketPrediction]],
    family_representative: Dict[str, str],
    min_tmscore: float,
    max_evalue: float,
    min_pocket_consistency: float,
    max_pockets: int,
) -> Dict[str, List[Evidence]]:
    """Parse foldseek hits, evaluate each template and build EC evidence."""
    evidence: Dict[str, List[Evidence]] = {}
    if not tsv_path.exists():
        return evidence

    context: Dict[str, Tuple[Dict[int, PdbResidue], list]] = {}
    for family_id, residues in query_residues.items():
        key_map = build_residue_key_map(residues)
        pocket_sets = pocket_seq_sets(
            query_pockets.get(family_id, []), key_map, max_pockets
        )
        context[family_id] = (index_by_seq(residues), pocket_sets)

    with tsv_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) < _NCOL:
                continue
            (query_fam, target, qstart_s, tstart_s, qaln, taln,
             tmscore_s, prob_s, evalue_s, qcov_s, tcov_s) = row[:_NCOL]
            template = templates.get(target)
            if template is None or query_fam not in context:
                continue
            try:
                qstart = int(qstart_s)
                tstart = int(tstart_s)
                tmscore = float(tmscore_s)
                prob = float(prob_s)
                evalue = float(evalue_s)
                qcov = float(qcov_s)
                tcov = float(tcov_s)
            except ValueError:
                continue
            if tmscore < min_tmscore or evalue > max_evalue:
                continue

            mapping = map_alignment(qaln, taln, qstart, tstart)
            query_by_index, pocket_sets = context[query_fam]
            match = evaluate_template(
                template,
                mapping,
                query_by_index,
                pocket_sets,
                tmscore,
                prob,
                evalue,
                qcov,
                tcov,
                min_pocket_consistency,
            )
            if match is None:
                continue
            rep_key = family_representative.get(query_fam, query_fam)
            evidence.setdefault(rep_key, []).append(
                _evidence_from_match(match, target, query_fam)
            )
    return evidence


def _evidence_from_match(
    match: TemplateMatch, target: str, family_id: str
) -> Evidence:
    details = {
        "db": "ec",
        "mcsa_id": match.mcsa_id,
        "ec": match.ec,
        "enzyme_name": match.enzyme_name,
        "pdb_target": target,
        "pocket_rank": match.pocket_rank,
        "tmscore": match.tmscore,
        "prob": match.prob,
        "qcov": match.qcov,
        "tcov": match.tcov,
        "catalytic_total": match.catalytic_total,
        "catalytic_mapped": match.catalytic_mapped,
        "catalytic_in_pocket": match.catalytic_in_pocket,
        "pocket_consistency": match.consistency,
        "family_id": family_id,
    }
    return Evidence(
        channel=EvidenceChannel.POCKET,
        target_id=f"EC:{match.ec}",
        target_label=match.enzyme_name or f"EC:{match.ec}",
        score=match.score,
        evalue=match.evalue,
        coverage=match.qcov,
        details=details,
    )


def _run(cmd: list) -> None:
    completed = subprocess.run(cmd, capture_output=True, text=True)
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(
            completed.returncode, cmd, completed.stdout, completed.stderr
        )
