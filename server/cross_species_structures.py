#!/usr/bin/env python3
"""Fetch precomputed structures for family representatives and search them.

Cross-species family representatives are plain UniProt accessions. For each
representative the script first tries to download an AlphaFold DB model (no GPU
cost, already computed for most UniProt entries). Accessions absent from
AlphaFold DB can optionally be folded with ESMFold on the local GPU. Every
representative with a model is then searched against the foldseek structural
database, and the hits are mapped to Pfam families through SIFTS.

The download, the foldseek database build and the search are all idempotent:
completed models and the raw search table are reused on re-run, so the step can
be interrupted and resumed without re-downloading or re-folding.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

# pufprobe is importable when installed editable (the conda env on the server).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pufprobe.fold.esmfold import (  # noqa: E402
    ModelUnavailable,
    fold_with_model,
    load_esmfold,
    parse_plddt_from_pdb,
)
from pufprobe.search.sifts import load_pdb_chain_pfam, pdb_chain_key  # noqa: E402

AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
USER_AGENT = "PUFprobe-cross-species-structures/1.0"
FOLDSEEK_COLS = "query,target,alntmscore,prob,evalue,qcov,tcov"


def _http_get(url: str, retries: int = 3, accept: str = "*/*"):
    request = urllib.request.Request(
        url, headers={"User-Agent": USER_AGENT, "Accept": accept})
    last_error: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (400, 404):
                raise  # deterministic client error: accession absent or invalid
            last_error = exc
        except Exception as exc:  # noqa: BLE001 - transient network errors
            last_error = exc
        if attempt < retries:
            time.sleep(2 * attempt)
    raise RuntimeError(f"download failed for {url}: {last_error}")


def download_one(acc: str, dest: Path, retries: int = 3
                 ) -> tuple[str, Optional[str], Optional[float]]:
    """Download one AlphaFold model via the prediction API.

    The API returns the current PDB URL and mean pLDDT, so the script stays
    valid as the model version advances (v4, v6, ...). Returns
    (status, version, mean_plddt). A 404 from the API means the accession is
    not modelled in AlphaFold DB.
    """
    target = dest / f"{acc}.pdb"
    if target.exists() and target.stat().st_size > 0:
        text = target.read_text(encoding="utf-8", errors="replace")
        return "cached", "cached", parse_plddt_from_pdb(text)

    try:
        meta_raw = _http_get(AF_API.format(acc=acc), retries=retries,
                             accept="application/json")
        records = json.loads(meta_raw.decode("utf-8"))
        record = records[0] if isinstance(records, list) else records
    except urllib.error.HTTPError as exc:
        if exc.code in (400, 404):
            return "missing", None, None
        return "error", None, None
    except Exception:  # noqa: BLE001
        return "error", None, None

    pdb_url = record.get("pdbUrl")
    if not pdb_url:
        return "missing", None, None
    try:
        raw = _http_get(pdb_url, retries=retries)
        if raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)
        text = raw.decode("utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        return "error", None, None

    target.write_text(text, encoding="utf-8")
    # Median per-residue pLDDT from the B-factor column, matching the ESMFold
    # path so structure confidence is comparable across both sources.
    return "alphafold", str(record.get("latestVersion", "unknown")), \
        parse_plddt_from_pdb(text)


def fetch_alphafold(accessions: List[str], struct_dir: Path, workers: int
                    ) -> Dict[str, dict]:
    struct_dir.mkdir(parents=True, exist_ok=True)
    manifest: Dict[str, dict] = {}
    missing: List[str] = []
    errors: List[str] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download_one, acc, struct_dir): acc
                   for acc in accessions}
        done = 0
        for future in concurrent.futures.as_completed(futures):
            acc = futures[future]
            status, version, plddt = future.result()
            manifest[acc] = {"source": status, "version": version, "plddt": plddt}
            if status == "missing":
                missing.append(acc)
            elif status == "error":
                errors.append(acc)
            done += 1
            if done % 200 == 0 or done == len(accessions):
                print(f"    AlphaFold: {done}/{len(accessions)} processed "
                      f"({len(missing)} missing, {len(errors)} errors)")
    (struct_dir / "missing_alphafold.txt").write_text(
        "\n".join(sorted(missing)) + ("\n" if missing else ""), encoding="utf-8")
    print(f"  AlphaFold models: {len(accessions) - len(missing) - len(errors)} "
          f"available, {len(missing)} absent, {len(errors)} transient errors.")
    return manifest


def read_fasta_headers(path: Path) -> List[str]:
    accessions = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                accessions.append(line[1:].strip().split()[0])
    return accessions


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    key = None
    chunks: List[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if key is not None:
                records[key] = "".join(chunks)
            key = line[1:].strip().split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if key is not None:
        records[key] = "".join(chunks)
    return records


def fold_missing(accessions: List[str], sequences: Dict[str, str],
                 struct_dir: Path, manifest: Dict[str, dict]) -> None:
    """Fold representatives absent from AlphaFold DB with ESMFold."""
    to_fold = [acc for acc in accessions
               if manifest.get(acc, {}).get("source") == "missing"]
    if not to_fold:
        print("  No representatives require ESMFold fallback.")
        return
    try:
        loaded = load_esmfold()
    except ModelUnavailable as exc:
        print(f"  ESMFold unavailable ({exc}); leaving {len(to_fold)} "
              "representatives without models.", file=sys.stderr)
        return
    folded = 0
    for acc in to_fold:
        target = struct_dir / f"{acc}.pdb"
        if target.exists() and target.stat().st_size > 0:
            manifest[acc] = {"source": "esmfold", "version": "esmfold_v1",
                             "plddt": parse_plddt_from_pdb(
                                 target.read_text(encoding="utf-8", errors="replace"))}
            folded += 1
            continue
        sequence = sequences.get(acc, "")
        if not sequence:
            continue
        try:
            pdb_str, plddt, _ = fold_with_model(loaded, sequence)
            target.write_text(pdb_str, encoding="utf-8")
            manifest[acc] = {"source": "esmfold", "version": "esmfold_v1",
                             "plddt": plddt}
            folded += 1
            if folded % 100 == 0:
                print(f"    ESMFold: {folded}/{len(to_fold)} folded")
        except ModelUnavailable as exc:
            print(f"  ESMFold stopped ({exc}); {len(to_fold) - folded} remain.",
                  file=sys.stderr)
            break
    if loaded.device == "cuda":
        loaded.torch.cuda.empty_cache()
    print(f"  ESMFold folded {folded} additional representatives.")


def run(cmd: List[str]) -> None:
    completed = subprocess.run(cmd, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(cmd)}\n"
            f"{completed.stdout}\n{completed.stderr}"
        )


def foldseek_search(struct_dir: Path, dbdir: Path, db_name: str,
                    out_dir: Path, threads: int, sensitivity: float) -> Path:
    foldseek = shutil.which("foldseek")
    if foldseek is None:
        raise RuntimeError("foldseek binary not found on PATH; run 06_update_software.sh.")
    db_prefix = dbdir / db_name
    if not Path(str(db_prefix) + ".index").exists():
        raise RuntimeError(
            f"foldseek database {db_prefix}.index not found; run 01_download_dbs.sh."
        )
    raw_tsv = out_dir / "cross_species_foldseek_raw.tsv"
    if raw_tsv.exists() and raw_tsv.stat().st_size > 0:
        return raw_tsv
    import tempfile
    tmp = Path(tempfile.mkdtemp(prefix="foldseek_cross_", dir=str(out_dir)))
    qdb = tmp / "qdb"
    res = tmp / "res"
    run([foldseek, "createdb", str(struct_dir), str(qdb)])
    run([foldseek, "search", str(qdb), str(db_prefix), str(res), str(tmp / "stmp"),
         "-a", "-s", str(sensitivity), "--threads", str(threads)])
    run([foldseek, "convertalis", str(qdb), str(db_prefix), str(res), str(raw_tsv),
         "--format-output", FOLDSEEK_COLS])
    shutil.rmtree(tmp, ignore_errors=True)
    return raw_tsv


def summarize_hits(raw_tsv: Path, sifts_path: Optional[Path],
                   min_tm: float, max_evalue: float) -> List[dict]:
    sifts = None
    if sifts_path is not None and sifts_path.exists():
        sifts = load_pdb_chain_pfam(sifts_path)
    best: Dict[str, dict] = {}
    if raw_tsv.exists():
        with raw_tsv.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                row = line.rstrip("\n").split("\t")
                if len(row) < 7:
                    continue
                query, target, tm_s, prob_s, evalue_s, qcov_s, tcov_s = row[:7]
                acc = Path(query).stem
                try:
                    tm, qcov, evalue = float(tm_s), float(qcov_s), float(evalue_s)
                except ValueError:
                    continue
                if tm < min_tm or evalue > max_evalue:
                    continue
                pfams = []
                if sifts:
                    pfams = sifts.get(pdb_chain_key(target), [])
                score = tm * qcov
                current = best.get(acc)
                if current is None or score > current["score"]:
                    best[acc] = {
                        "accession": acc, "best_target": target, "tmscore": tm,
                        "qcov": qcov, "evalue": evalue, "score": round(score, 4),
                        "pfams": ";".join(pfams) if pfams else "",
                    }
    return list(best.values())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=Path, required=True,
                        help="representatives.faa from stage 08")
    parser.add_argument("--struct-dir", type=Path, required=True,
                        help="directory for downloaded/predicted PDB models")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="directory for structure search outputs")
    parser.add_argument("--dbdir", type=Path, required=True)
    parser.add_argument("--foldseek-db", default="foldseek_db")
    parser.add_argument("--sifts", type=Path, default=None)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--sensitivity", type=float, default=9.5)
    parser.add_argument("--min-tm", type=float, default=0.5)
    parser.add_argument("--max-evalue", type=float, default=10.0)
    parser.add_argument("--esmfold-fallback", action="store_true")
    parser.add_argument("--max-reps", type=int, default=0,
                        help=">0 limits representatives (smoke test); 0 = all")
    args = parser.parse_args()

    args.struct_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    accessions = read_fasta_headers(args.reps)
    if args.max_reps > 0:
        accessions = accessions[:args.max_reps]
    print(f"Family representatives: {len(accessions)}")

    manifest = fetch_alphafold(accessions, args.struct_dir, args.workers)

    if args.esmfold_fallback:
        sequences = read_fasta(args.reps)
        fold_missing(accessions, sequences, args.struct_dir, manifest)

    (args.out_dir / "structure_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")

    print("Running foldseek structural search ...")
    try:
        raw_tsv = foldseek_search(args.struct_dir, args.dbdir, args.foldseek_db,
                                  args.out_dir, args.threads, args.sensitivity)
    except RuntimeError as exc:
        print(f"Foldseek search skipped: {exc}", file=sys.stderr)
        return 1

    hits = summarize_hits(raw_tsv, args.sifts, args.min_tm, args.max_evalue)
    hits_path = args.out_dir / "structure_hits.tsv"
    with hits_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("accession\tbest_target\ttmscore\tqcov\tevalue\tscore\tpfams\n")
        for hit in sorted(hits, key=lambda h: -h["score"]):
            handle.write(f"{hit['accession']}\t{hit['best_target']}\t"
                         f"{hit['tmscore']:.3f}\t{hit['qcov']:.3f}\t"
                         f"{hit['evalue']:.2e}\t{hit['score']}\t{hit['pfams']}\n")

    n_with_pfam = sum(1 for h in hits if h["pfams"])
    print(f"Structure search done: {len(hits)} representatives with a confident "
          f"fold-level hit; {n_with_pfam} mapped to Pfam families.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
