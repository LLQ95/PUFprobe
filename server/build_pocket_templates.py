#!/usr/bin/env python
"""Build the M-CSA pocket template database for the pocket evidence channel.

The M-CSA reference structures are turned into single-chain PDB files and a
foldseek database, together with a JSON index that records each template's
catalytic residues. At search time foldseek aligns query structures against
the database and the mapped catalytic residues are checked against P2Rank
pockets.

Usage:
    python build_pocket_templates.py <mcsa_dir> <work_dir> <dbdir> \
        [--foldseek PATH] [--threads N] [--rebuild-db]

Inputs in <mcsa_dir> (from 12_download_evidence_dbs.sh mcsa):
    mcsa_residues.json   flat list of catalytic residues
    mcsa_entries.json    first API page (the API caps page size at 100, so
                         the full entry set is paged here and cached)

Intermediate caches in <work_dir>:
    raw_pdb/             downloaded reference PDB files
    chain_pdb/           single-chain PDB files fed to foldseek
    entries_all.json     merged full entry set

Final artifacts installed in <dbdir>:
    pocket_db*           foldseek database
    pocket_templates.json

The build is idempotent and resumable: downloaded PDB files, extracted
chains and the foldseek database are reused when present, so an interrupted
run can be restarted with the same arguments.
"""

from __future__ import annotations

import argparse
import json
import shutil
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from pufprobe.search.pocket import (
    extract_chain_pdb,
    parse_pdb_residues,
)
from pufprobe.search.pocket_templates import (
    collect_reference_sites,
    entries_by_id,
    make_template_record,
)

_RCSB_URL = "https://files.rcsb.org/download/{pdb}.pdb"
_EBI_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb}.pdb"
_PDB_SOURCES = (_RCSB_URL, _EBI_URL)
_MCSA_ENTRIES_NAME = "mcsa_entries.json"
_MCSA_RESIDUES_NAME = "mcsa_residues.json"
_ENTRIES_ALL_NAME = "entries_all.json"
_TEMPLATE_NAME = "pocket_templates.json"
_DB_NAME = "pocket_db"
# Retry budget for transient network failures.
_MAX_RETRIES = 4
_RETRY_WAIT = 5


def _fetch(url: str, timeout: int = 90) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "pufprobe-build"})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read()


def _fetch_with_retry(url: str) -> bytes:
    last_error: Optional[Exception] = None
    for attempt in range(_MAX_RETRIES):
        try:
            return _fetch(url)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = exc
            if attempt + 1 < _MAX_RETRIES:
                time.sleep(_RETRY_WAIT * (attempt + 1))
    raise RuntimeError(f"failed to download {url}: {last_error}")


def ensure_full_entries(mcsa_dir: Path, work_dir: Path) -> Dict[int, dict]:
    """Load the full M-CSA entry set, paging the API and caching the result."""
    cache_path = work_dir / _ENTRIES_ALL_NAME
    if cache_path.exists():
        with cache_path.open("r", encoding="utf-8") as handle:
            cached = json.load(handle)
        if isinstance(cached, list) and cached:
            return entries_by_id(cached)

    first_page_path = mcsa_dir / _MCSA_ENTRIES_NAME
    if not first_page_path.exists():
        raise SystemExit(f"missing {first_page_path}; run the mcsa download step")
    with first_page_path.open("r", encoding="utf-8") as handle:
        first_page = json.load(handle)
    if isinstance(first_page, list):
        all_entries = first_page
    elif isinstance(first_page, dict):
        all_entries = list(first_page.get("results", []))
        next_url = first_page.get("next")
        page = 1
        while next_url:
            page += 1
            payload = json.loads(_fetch_with_retry(next_url).decode("utf-8"))
            all_entries.extend(payload.get("results", []))
            next_url = payload.get("next")
    else:
        raise SystemExit(f"unrecognized entries payload in {first_page_path}")

    work_dir.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(all_entries, handle, ensure_ascii=False, indent=2)
    return entries_by_id(all_entries)


def download_pdb(pdb_id: str, raw_dir: Path) -> Path:
    """Download a reference PDB once, trying RCSB then the EBI mirror."""
    target = raw_dir / f"{pdb_id}.pdb"
    if target.exists() and target.stat().st_size > 0:
        return target
    raw_dir.mkdir(parents=True, exist_ok=True)
    last_error: Optional[Exception] = None
    for template in _PDB_SOURCES:
        url = template.format(pdb=pdb_id)
        try:
            data = _fetch_with_retry(url)
        except RuntimeError as exc:
            last_error = exc
            continue
        if data:
            target.write_bytes(data)
            return target
    raise RuntimeError(f"could not download {pdb_id} from any source: {last_error}")


def build_chain_pdb(
    pdb_path: Path,
    chain: str,
    target: str,
    chain_dir: Path,
    expected: Dict[int, str],
) -> List[dict]:
    """Extract one chain, write it and verify its catalytic residues.

    Returns the residues verified against the structure at their chain-local
    sequence position. Residues that are missing, out of range or whose
    three-letter code does not match M-CSA are dropped and reported.
    """
    out_path = chain_dir / f"{target}.pdb"
    if not out_path.exists():
        pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
        chain_text = extract_chain_pdb(pdb_text, chain)
        chain_dir.mkdir(parents=True, exist_ok=True)
        out_path.write_text(chain_text, encoding="utf-8", newline="\n")

    residues = parse_pdb_residues(out_path.read_text(encoding="utf-8"))
    by_index = {r.seq_index: r for r in residues}
    verified: List[dict] = []
    for resid, code in sorted(expected.items()):
        residue = by_index.get(resid)
        if residue is None:
            print(f"  [warn] {target}: catalytic resid {resid} out of range")
            continue
        if residue.resname.upper() != code.upper():
            print(
                f"  [warn] {target}: resid {resid} is {residue.resname}, "
                f"expected {code}"
            )
            continue
        verified.append({"resid": resid, "code": code})
    return verified


def read_foldseek_targets(db_prefix: Path) -> List[str]:
    """Read the entry names foldseek assigned, in database order."""
    header_path = Path(str(db_prefix) + "_h")
    if not header_path.exists():
        return []
    names: List[str] = []
    text = header_path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        line = line.strip()
        if line:
            names.append(line.split(" ", 1)[0])
    return names


def run_foldseek_createdb(
    foldseek: str, chain_dir: Path, db_prefix: Path, threads: int, rebuild: bool
) -> None:
    """Build the foldseek database unless a valid one already exists."""
    index_path = Path(str(db_prefix) + ".index")
    if index_path.exists() and not rebuild:
        print(f"foldseek database exists at {db_prefix}; reusing")
        return
    import subprocess

    cmd = [
        foldseek,
        "createdb",
        str(chain_dir),
        str(db_prefix),
        "--threads",
        str(threads),
    ]
    print("running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mcsa_dir", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("dbdir", type=Path)
    parser.add_argument("--foldseek", default=None)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--rebuild-db", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    foldseek = args.foldseek or shutil.which("foldseek")
    if not foldseek:
        raise SystemExit(
            "foldseek not found; pass --foldseek PATH or put it on PATH"
        )

    mcsa_dir: Path = args.mcsa_dir
    work_dir: Path = args.work_dir
    dbdir: Path = args.dbdir
    raw_dir = work_dir / "raw_pdb"
    chain_dir = work_dir / "chain_pdb"

    residues_path = mcsa_dir / _MCSA_RESIDUES_NAME
    if not residues_path.exists():
        raise SystemExit(f"missing {residues_path}; run the mcsa download step")
    with residues_path.open("r", encoding="utf-8") as handle:
        residues = json.load(handle)

    print("loading M-CSA entries...")
    entries = ensure_full_entries(mcsa_dir, work_dir)
    sites = collect_reference_sites(residues)
    print(f"reference sites: {len(sites)}; entries: {len(entries)}")

    # Build and verify every single-chain template. Sites whose structure or
    # catalytic residues cannot be verified are skipped but do not abort.
    prepared: List[Tuple[str, dict, List[dict]]] = []
    for index, ((pdb_id, chain), site) in enumerate(sorted(sites.items()), 1):
        target = f"{pdb_id}_{chain}"
        print(f"[{index}/{len(sites)}] {target}")
        try:
            pdb_path = download_pdb(pdb_id, raw_dir)
            verified = build_chain_pdb(
                pdb_path, chain, target, chain_dir, site["residues"]
            )
        except RuntimeError as exc:
            print(f"  [warn] skipping {target}: {exc}")
            continue
        if not verified:
            print(f"  [warn] skipping {target}: no verified catalytic residue")
            continue
        prepared.append((target, site, verified))

    db_prefix = dbdir / _DB_NAME
    dbdir.mkdir(parents=True, exist_ok=True)
    run_foldseek_createdb(
        foldseek, chain_dir, db_prefix, args.threads, args.rebuild_db
    )
    actual_targets = set(read_foldseek_targets(db_prefix))

    templates: List[dict] = []
    for target, site, verified in prepared:
        if target not in actual_targets:
            print(f"  [warn] {target} absent from foldseek headers; skipped")
            continue
        entry = entries.get(int(site["mcsa_id"]))
        templates.append(make_template_record(target, site, entry, verified))

    templates_path = dbdir / _TEMPLATE_NAME
    with templates_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(templates, handle, ensure_ascii=False, indent=2)

    with_ec = sum(1 for record in templates if record["ec"])
    print(
        f"wrote {len(templates)} templates ({with_ec} with EC) to "
        f"{templates_path}"
    )


if __name__ == "__main__":
    main()
