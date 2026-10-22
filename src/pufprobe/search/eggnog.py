"""eggNOG functional annotation via eggNOG-mapper.

Runs ``emapper`` (eggNOG 7) on the representative sequences and parses the
``.emapper.annotations`` table into ``Evidence(EGGNOG)``. The seed-ortholog
E-value supplies the score, and the orthologous group, description, GO, EC and
Pfam fields supply candidate terms. A missing binary or database disables the
channel (logged).

Unlike the per-query profile searches, emapper annotates the whole fasta in a
single invocation and parallelises internally.

eggNOG is treated as a curated database witness (see CHANNEL_GROUPS),
independent of the representative profile searches the pipeline runs itself,
so an eggNOG annotation and a profile hit to the same function count as two
orthogonal evidence groups.
"""

from __future__ import annotations

import math
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from pufprobe import externals
from pufprobe.io_utils import read_fasta, write_fasta
from pufprobe.logging_utils import get_logger
from pufprobe.models import Evidence, EvidenceChannel

logger = get_logger(__name__)

INSTALL_HINT = (
    "Install eggNOG-mapper from bioconda "
    "(`conda install -c bioconda eggnog-mapper`) and download the eggNOG "
    "database (`download_eggnog_data.py`)."
)

# Column names in the .emapper.annotations header for emapper 3.0 / eggNOG 7.
# emapper 2.x uses a different 21-column header (max_annot_lvl, Description,
# no tax_ceiling/farthest_donor_lineage/annotation_confidence); the parser is
# driven by the header it reads, so both layouts are accepted.
ANNOT_FIELDS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
    "tax_ceiling", "farthest_donor_lineage", "COG_category",
    "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway",
    "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
    "annotation_confidence",
]

# One orthologous-group token, e.g. "COG0460@1|root",
# "1MW3H@1224|Proteobacteria" or "DUF874@131567|A-1". A trailing "*" on the
# token marks the primary OG in emapper 3.0.
_OG_TOKEN_RE = re.compile(
    r"^(?P<og>[^@|]+)@(?P<level>[^|]+?)(?:\|(?P<lvl_name>.*?))?(?P<star>\*)?$"
)

# One Pfam hit as written by emapper 3.0: "<short_id>_<start>_<end>",
# e.g. "DUF874_144_191". emapper 2.x writes the plain short id instead.
_PFAM_HIT_RE = re.compile(r"^(?P<id>.+)_(?P<start>\d+)_(?P<end>\d+)$")

# Values emapper uses for missing annotations.
_NA = {"-", "", "NA"}

# Saturation point for the E-value -> score mapping: E <= 1e-200 scores 1.0.
_EVALUE_SAT_LOG10 = 200.0


def evalue_to_score(evalue_text: str) -> float:
    """Map a seed-ortholog E-value to a bounded score.

    -log10(E) is normalised so that E <= 1e-200 reaches 1.0. An E-value of 0
    (an effectively exact match) scores 1.0.
    """
    try:
        evalue = float(evalue_text)
    except (TypeError, ValueError):
        return 0.0
    if evalue <= 0.0:
        return 1.0
    return max(0.0, min(1.0, -math.log10(evalue) / _EVALUE_SAT_LOG10))


def _split_list(text: str) -> List[str]:
    """Split a comma-separated emapper column, dropping missing values."""
    out: List[str] = []
    for token in (text or "").split(","):
        token = token.strip()
        if token and token not in _NA:
            out.append(token)
    return out


def parse_ogs(ogs_text: str) -> List[Dict[str, str]]:
    """Parse the eggNOG_OGs column into ordered OG tokens."""
    out: List[Dict[str, str]] = []
    for token in (ogs_text or "").split(","):
        match = _OG_TOKEN_RE.match(token.strip())
        if not match:
            continue
        out.append({
            "og": match.group("og").strip(),
            "level": match.group("level").strip(),
            "level_name": (match.group("lvl_name") or "").strip(),
            "primary": "1" if match.group("star") else "0",
        })
    return out


def parse_pfam_hits(pfams_text: str) -> List[str]:
    """Parse the PFAMs column into unique Pfam short ids.

    emapper 3.0 writes each hit as "<id>_<start>_<end>" (one domain region),
    so the coordinate suffix is stripped and repeated ids are collapsed;
    emapper 2.x writes the plain short id.
    """
    ids: List[str] = []
    for token in _split_list(pfams_text):
        match = _PFAM_HIT_RE.match(token)
        short = match.group("id") if match else token
        if short and short not in ids:
            ids.append(short)
    return ids


def _best_og(ogs: List[Dict[str, str]]) -> Optional[Dict[str, str]]:
    """Pick the term-naming OG.

    The primary OG marked by emapper 3.0 ("*") is preferred, then a COG id
    (the broad, stable identifier that lets homologs converge across
    species); otherwise the most specific (last) OG is used.
    """
    if not ogs:
        return None
    for og in ogs:
        if og.get("primary") == "1":
            return og
    for og in ogs:
        if og["og"].startswith("COG"):
            return og
    return ogs[-1]


def build_pfam_id_index(
    pfam_names: Dict[str, Dict[str, str]],
) -> Dict[str, str]:
    """Build a Pfam short-name -> accession index from the name table.

    ``pfam_names`` maps accession -> {id, name, type}. The short identifier
    (``id``, e.g. "DUF874" or "ABC_tran") is the token used in the eggNOG
    PFAMs column, so reversing the index names a Pfam term by accession and
    lets it converge with the sequence channels.
    """
    index: Dict[str, str] = {}
    for accession, record in pfam_names.items():
        short = (record.get("id") or "").strip()
        if short:
            index[short] = accession
    return index


def _evidence_payloads(
    rec: Dict[str, str],
    pfam_id_index: Optional[Dict[str, str]] = None,
) -> List[Dict[str, Any]]:
    """Build raw evidence payloads for one parsed annotation row."""
    score = evalue_to_score(rec.get("evalue", ""))
    ogs = parse_ogs(rec.get("eggNOG_OGs", ""))
    # Description exists only in emapper 2.x; 3.0 drops the column.
    desc = (rec.get("Description") or "").strip()
    pref = (rec.get("Preferred_name") or "").strip()
    desc = "" if desc in _NA else desc
    pref = "" if pref in _NA else pref

    payloads: List[Dict[str, Any]] = []

    best = _best_og(ogs)
    if best is not None:
        og_id = best["og"]
        payloads.append({
            "target_id": og_id,
            "target_label": desc or pref or og_id,
            "score": score,
            "kind": "og",
            "og_level": best["level"],
            "has_description": bool(desc),
        })

    for go in _split_list(rec.get("GOs", "")):
        if go.startswith("GO:"):
            payloads.append({
                "target_id": go, "target_label": go,
                "score": score, "kind": "go",
            })

    for ec in _split_list(rec.get("EC", "")):
        norm = ec if ec.startswith("EC:") else f"EC:{ec}"
        payloads.append({
            "target_id": norm, "target_label": ec,
            "score": score, "kind": "ec",
        })

    for short in parse_pfam_hits(rec.get("PFAMs", "")):
        accession = (pfam_id_index or {}).get(short)
        if accession is not None:
            payloads.append({
                "target_id": accession, "target_label": accession,
                "score": score, "kind": "pfam",
            })
        else:
            # Keep the short name as a generic label when no accession exists.
            payloads.append({
                "target_id": short, "target_label": short,
                "score": score, "kind": "pfam_name",
            })

    return payloads


def parse_emapper(
    annot_path: Path,
    pfam_id_index: Optional[Dict[str, str]] = None,
) -> Dict[str, List[Dict[str, Any]]]:
    """Parse an .emapper.annotations file into per-query evidence payloads."""
    rows: Dict[str, List[Dict[str, Any]]] = {}
    header: Optional[List[str]] = None
    with annot_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line[1:].lstrip().split("\t")
                continue
            if header is None:
                continue
            values = line.split("\t")
            if len(values) < len(header):
                values += [""] * (len(header) - len(values))
            rec = dict(zip(header, values))
            qid = rec.get("query", "").strip()
            if qid:
                rows[qid] = _evidence_payloads(rec, pfam_id_index)
    return rows


def run_eggnog_search(
    reps_faa: Path,
    dbdir: Optional[Path],
    out_dir: Path,
    threads: int = 1,
    data_dir: Optional[str] = None,
    search_mode: str = "diamond",
    pfam_names: Optional[Dict[str, Dict[str, str]]] = None,
    emapper_path: Optional[str] = None,
) -> Tuple[Dict[str, List[Evidence]], Optional[Path]]:
    """Run emapper on the representative fasta and return Evidence."""
    # Resolve the interpreter and environment. An absolute emapper path often
    # points into a dedicated conda env (emapper.py is not on the active env's
    # PATH). Invoke it with that env's python and prepend its bin directory to
    # PATH so the bundled diamond/MMseqs2 are found even when PUFprobe itself
    # runs from a different conda env.
    run_env = os.environ.copy()
    if emapper_path:
        emapper_exe = Path(emapper_path)
        if not emapper_exe.exists():
            logger.warning(
                "configured emapper path does not exist: %s; skipping eggnog.",
                emapper_path,
            )
            return {}, None
        env_bin = emapper_exe.parent
        env_python = env_bin / "python"
        if env_python.exists():
            cmd_head = [str(env_python), str(emapper_exe)]
        else:
            cmd_head = [str(emapper_exe)]
        run_env["PATH"] = str(env_bin) + os.pathsep + run_env.get("PATH", "")
    else:
        info = externals.emapper()
        if not info.available:
            logger.warning(
                "emapper not found; skipping eggnog channel. %s", INSTALL_HINT)
            return {}, None
        cmd_head = [info.path or "emapper"]

    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="eggnog_", dir=str(out_dir)))
    raw_tsv = out_dir / "eggnog_raw.tsv"

    queries = [(header.split()[0], seq) for header, seq in read_fasta(reps_faa)]
    safe_map: Dict[str, str] = {}
    safe_records: List[Tuple[str, str]] = []
    for global_key, seq in queries:
        safe = global_key.replace("@", "_")
        safe_map[safe] = global_key
        safe_records.append((safe, seq))
    in_fasta = tmp / "reps.faa"
    write_fasta(in_fasta, safe_records)

    prefix = "eggnog"
    cmd = cmd_head + [
        "-i", str(in_fasta),
        "--output", prefix,
        "--output_dir", str(tmp),
        "--cpu", str(max(1, int(threads))),
        "-m", search_mode,
        "--override",
    ]
    if data_dir:
        cmd += ["--data_dir", data_dir]

    try:
        completed = subprocess.run(
            cmd, capture_output=True, text=True, env=run_env)
        annot = tmp / f"{prefix}.emapper.annotations"
        if completed.returncode != 0 or not annot.exists():
            logger.warning(
                "emapper failed (returncode %s); skipping eggnog channel.\n"
                "STDOUT: %s\nSTDERR: %s\n%s",
                completed.returncode,
                (completed.stdout or "").strip()[-800:],
                (completed.stderr or "").strip()[-800:],
                INSTALL_HINT,
            )
            return {}, None
    except (subprocess.SubprocessError, OSError) as exc:
        logger.warning("emapper error (%s); skipping eggnog channel.", exc)
        return {}, None

    pfam_id_index = build_pfam_id_index(pfam_names or {})
    parsed = parse_emapper(annot, pfam_id_index if pfam_names else None)

    evidence: Dict[str, List[Evidence]] = {}
    rows = ["query\tkind\ttarget\tscore"]
    for safe, payloads in parsed.items():
        global_key = safe_map.get(safe, safe)
        for p in payloads:
            rows.append(
                f"{global_key}\t{p['kind']}\t{p['target_id']}\t{p['score']:.2f}"
            )
            details: Dict[str, Any] = {"database": "eggnog", "kind": p["kind"]}
            if "og_level" in p:
                details["og_level"] = p["og_level"]
            evidence.setdefault(global_key, []).append(Evidence(
                channel=EvidenceChannel.EGGNOG,
                target_id=p["target_id"],
                target_label=p["target_label"],
                score=p["score"],
                details=details,
            ))

    raw_tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")
    logger.info(
        "eggnog annotated %d / %d representatives (%d evidence entries).",
        len(parsed), len(queries),
        sum(len(v) for v in evidence.values()),
    )
    return evidence, raw_tsv
