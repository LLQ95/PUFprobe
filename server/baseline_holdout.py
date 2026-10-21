#!/usr/bin/env python3
"""Swiss-Prot hold-out benchmark for sequence, structure and fused prediction.

A set of reviewed bacterial proteins with curated Pfam annotations is treated
as if their function were unknown. Their own annotations are hidden (the
self-hit is removed from sequence results, and each prediction is scored
against the curated Pfam set). Three channels are scored independently:

  * sequence homology  - best non-self MMseqs2 hit against Swiss-Prot,
  * structure          - Foldseek hit over the AlphaFold model, mapped to Pfam
                         through SIFTS,
  * fusion             - the union of both channels (recall) and the subset
                         supported by both (consensus precision).

Metrics are broken down by phylum, by sequence identity to the nearest
non-self homologue, and by structural TM-score. The identity breakdown is the
key result: the structure channel is expected to recover annotations in the
low-identity range where sequence homology stops being informative.

External annotators (Bakta, Prokka, InterProScan, eggNOG, DeepFRI) are not
run here; the per-protein table exposes the same ground truth so their
predictions can be scored by the same routine when those tools are available.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional

# Reuse the validated UniProt paging and phylum parser.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch_reference_proteomes import request_page, phylum_of, is_hypothetical  # noqa: E402

import numpy as np  # noqa: E402

PAGE_SIZE = 500
FIELDS = ("accession,protein_name,sequence,length,organism_id,organism_name,"
          "lineage,xref_pfam,annotation_score")


def phylum_from_ranked_lineage(lineage: str) -> str:
    """Extract the phylum from a UniProtKB lineage string with rank suffixes.

    Unlike the proteomes endpoint, UniProtKB returns entries such as
    "Bacillota (phylum)", so the explicit rank tag is the most reliable
    marker. Falls back to the unranked phylum set when no tag is present.
    """
    if not lineage:
        return ""
    for token in lineage.split(","):
        token = token.strip()
        if token.endswith("(phylum)"):
            return token[: -len("(phylum)")].strip()
    cleaned = ", ".join(
        part.split(" (")[0].strip() for part in lineage.split(",") if part.strip())
    return phylum_of(cleaned)


def sample_holdout(target: int, seed: int, max_pages: int,
                   min_length: int, max_length: int) -> List[dict]:
    """Sample reviewed bacterial proteins with curated Pfam annotations."""
    query = (f"(reviewed:true) AND (taxonomy_id:2) AND "
             f"(length:[{min_length} TO {max_length}])")
    chosen: Dict[str, dict] = {}
    by_phylum: Counter = Counter()
    cursor: Optional[str] = None
    pages = 0
    while len(chosen) < target and pages < max_pages:
        params = {"query": query, "fields": FIELDS, "format": "tsv",
                  "size": str(PAGE_SIZE), "compressed": "true"}
        if cursor:
            params["cursor"] = cursor
        rows, cursor = request_page("uniprotkb", params)
        pages += 1
        if not rows:
            break
        for row in rows:
            acc = row.get("Entry", "").strip()
            pfam_raw = row.get("Pfam", "").strip()
            pfams = sorted({p.strip() for p in pfam_raw.split(";") if p.strip()})
            if not acc or not pfams:
                continue
            name = row.get("Protein names", "")
            if is_hypothetical(name):
                continue
            phylum = phylum_from_ranked_lineage(row.get("Taxonomic lineage", ""))
            if not phylum:
                continue
            record = {
                "accession": acc,
                "protein_name": name,
                "sequence": row.get("Sequence", "").replace(" ", ""),
                "length": row.get("Length", ""),
                "phylum": phylum,
                "organism": row.get("Organism", ""),
                "pfams": pfams,
            }
            # Keep a phylum quota so the set is not dominated by proteobacteria.
            if by_phylum[phylum] >= max(20, target // 6) and acc not in chosen:
                continue
            chosen[acc] = record
            by_phylum[phylum] += 1
            if len(chosen) >= target:
                break
        if not cursor:
            break
    records = list(chosen.values())
    rng = np.random.default_rng(seed)
    rng.shuffle(records)
    print(f"Hold-out set: {len(records)} reviewed proteins across "
          f"{len(by_phylum)} phyla ({pages} pages).")
    return records


def write_fasta(records: List[dict], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for rec in records:
            handle.write(f">{rec['accession']}\n{rec['sequence']}\n")


def run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(cmd)}\n{result.stderr[-2000:]}")


def mmseqs_sequence_hits(query_fasta: Path, swissprot_db: Path,
                         out_dir: Path, threads: int, sensitivity: float
                         ) -> Dict[str, dict]:
    """Best non-self MMseqs2 hit per query, mapped to Pfam via xref."""
    work = out_dir / "mmseqs_work"
    work.mkdir(parents=True, exist_ok=True)
    qdb = work / "qdb"
    if not (work / "qdb.index").exists():
        run(["mmseqs", "createdb", str(query_fasta), str(qdb)])
    aln = work / "aln"
    if not (work / "aln.index").exists():
        run(["mmseqs", "search", str(qdb), str(swissprot_db), str(aln),
             str(work / "tmp"), "--threads", str(threads), "-s", str(sensitivity)])
    raw = out_dir / "sequence_hits.tsv"
    if not raw.exists():
        run(["mmseqs", "convertalis", str(qdb), str(swissprot_db), str(aln),
             str(raw), "--format-output",
             "query,target,pident,evalue,qcov,tcov"])
    from pufprobe.search.idmapping import load_xref, accessions_from_target
    xref = load_xref(swissprot_db.parent / "swissprot_xref.json")
    best: Dict[str, dict] = {}
    with raw.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) < 6:
                continue
            query, target, pident, evalue, qcov, tcov = row[:6]
            qacc = query.split()[0].split("|")[-1].strip()
            target_accs = accessions_from_target(target)
            if qacc in target_accs:
                continue  # remove the self-hit
            try:
                pident_f, evalue_f, qcov_f = float(pident), float(evalue), float(qcov)
            except ValueError:
                continue
            current = best.get(qacc)
            if current is None or evalue_f < current["evalue"]:
                pfams = set()
                for tacc in target_accs:
                    pfams.update(xref.get(tacc, {}).get("Pfam", []))
                best[qacc] = {
                    "target": target, "pident": pident_f, "evalue": evalue_f,
                    "qcov": qcov_f, "pfams": sorted(pfams),
                }
    return best


def foldseek_structure_hits(records: List[dict], struct_dir: Path,
                            dbdir: Path, foldseek_db: str,
                            out_dir: Path, sifts: Path, threads: int,
                            workers: int, sensitivity: float,
                            min_tm: float) -> Dict[str, dict]:
    """Download AlphaFold models and map the best Foldseek hit to Pfam."""
    import cross_species_structures as cs

    accessions = [r["accession"] for r in records]
    manifest = cs.fetch_alphafold(accessions, struct_dir, workers)
    raw = cs.foldseek_search(struct_dir, dbdir, foldseek_db, out_dir,
                             threads, sensitivity)
    rows = cs.summarize_hits(raw, sifts, min_tm=min_tm, max_evalue=10.0)
    hits = {r["accession"]: r for r in rows}
    for acc in accessions:
        if acc in hits:
            hits[acc]["plddt"] = manifest.get(acc, {}).get("plddt")
    return hits


def identity_bin(pident: Optional[float]) -> str:
    if pident is None:
        return "no hit"
    if pident < 20:
        return "<20%"
    if pident < 30:
        return "20-30%"
    if pident < 50:
        return "30-50%"
    return ">=50%"


def tm_bin(tm: Optional[float]) -> str:
    if tm is None:
        return "no hit"
    if tm < 0.5:
        return "<0.5"
    if tm < 0.7:
        return "0.5-0.7"
    if tm < 0.9:
        return "0.7-0.9"
    return ">=0.9"


def score_channel(pred: List[str], truth: List[str]) -> tuple:
    truth_set, pred_set = set(truth), set(pred)
    if not pred_set:
        return 0, 0.0, 0.0, False
    overlap = pred_set & truth_set
    precision = len(overlap) / len(pred_set)
    recall = len(overlap) / len(truth_set) if truth_set else 0.0
    return len(overlap), precision, recall, bool(overlap)


def summarize(results: List[dict], key: str) -> dict:
    scored = [r for r in results if r.get("seq_pfams") is not None or
              r.get("struct_pfams") is not None]
    has_pred = [r for r in scored if r.get(f"{key}_pfams")]
    correct = [r for r in has_pred if r.get(f"{key}_correct")]
    precisions = [r[f"{key}_precision"] for r in has_pred]
    recalls = [r[f"{key}_recall"] for r in has_pred]
    return {
        "channel": key,
        "n": len(scored),
        "coverage": round(len(has_pred) / len(scored), 4) if scored else 0,
        "any_correct_rate": round(len(correct) / len(has_pred), 4) if has_pred else 0,
        "mean_precision": round(float(np.mean(precisions)), 4) if precisions else 0,
        "mean_recall": round(float(np.mean(recalls)), 4) if recalls else 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dbdir", type=Path, required=True)
    parser.add_argument("--swissprot-db", default="swissprot")
    parser.add_argument("--foldseek-db", default="foldseek_db")
    parser.add_argument("--sifts", default=None)
    parser.add_argument("--sample", type=int, default=20000)
    parser.add_argument("--max-pages", type=int, default=120)
    parser.add_argument("--min-length", type=int, default=50)
    parser.add_argument("--max-length", type=int, default=1500)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seq-sensitivity", type=float, default=7.5)
    parser.add_argument("--struct-sensitivity", type=float, default=9.5)
    parser.add_argument("--min-tm", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--no-structures", action="store_true")
    parser.add_argument("--cache", default=None,
                        help="hold-out sample cache TSV (reused if present)")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    struct_dir = args.output / "structures"
    struct_dir.mkdir(parents=True, exist_ok=True)
    sifts = Path(args.sifts) if args.sifts else args.dbdir / "pdb_chain_pfam.tsv.gz"
    swissprot_db = (
        args.dbdir / args.swissprot_db if not Path(args.swissprot_db).is_absolute()
        else Path(args.swissprot_db))

    cache_path = Path(args.cache) if args.cache else args.output / "holdout_sample.tsv"
    if cache_path.exists():
        records = []
        with cache_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row["pfams"] = row["pfams"].split(";") if row["pfams"] else []
                records.append(row)
        print(f"Loaded cached hold-out sample: {len(records)} proteins.")
    else:
        records = sample_holdout(args.sample, args.seed, args.max_pages,
                                 args.min_length, args.max_length)
        with cache_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["accession", "protein_name", "sequence", "length",
                             "phylum", "organism", "pfams"])
            for rec in records:
                writer.writerow([rec["accession"], rec["protein_name"],
                                 rec["sequence"], rec["length"], rec["phylum"],
                                 rec["organism"], ";".join(rec["pfams"])])
    if args.sample and len(records) > args.sample:
        records = records[:args.sample]

    query_fasta = args.output / "holdout.faa"
    write_fasta(records, query_fasta)

    print("Sequence channel: MMseqs2 against Swiss-Prot ...")
    seq_hits = mmseqs_sequence_hits(query_fasta, swissprot_db, args.output,
                                    args.threads, args.seq_sensitivity)
    print(f"  non-self sequence hits: {len(seq_hits)}")

    struct_hits: Dict[str, dict] = {}
    if not args.no_structures:
        print("Structure channel: AlphaFold + Foldseek ...")
        struct_hits = foldseek_structure_hits(
            records, struct_dir, args.dbdir, args.foldseek_db,
            args.output / "structure", sifts, args.threads, args.workers,
            args.struct_sensitivity, args.min_tm)
        print(f"  structural hits: {len(struct_hits)}")

    # ---- per-protein scoring -------------------------------------------------
    results = []
    for rec in records:
        acc = rec["accession"]
        truth = rec["pfams"]
        seq = seq_hits.get(acc)
        seq_pfams = seq["pfams"] if seq else []
        struct = struct_hits.get(acc)
        struct_pfams = sorted(set((struct or {}).get("pfams", "").split(";")) - {""})
        union_pfams = sorted(set(seq_pfams) | set(struct_pfams))
        consensus_pfams = sorted(set(seq_pfams) & set(struct_pfams))
        row = {
            "accession": acc, "phylum": rec["phylum"], "length": rec["length"],
            "truth_pfams": ";".join(truth),
            "seq_pident": seq["pident"] if seq else None,
            "seq_pfams": ";".join(seq_pfams),
            "struct_tm": struct.get("tmscore") if struct else None,
            "struct_plddt": struct.get("plddt") if struct else None,
            "struct_pfams": ";".join(struct_pfams),
            "union_pfams": ";".join(union_pfams),
            "consensus_pfams": ";".join(consensus_pfams),
        }
        for key, pred in (("seq", seq_pfams), ("struct", struct_pfams),
                          ("union", union_pfams), ("consensus", consensus_pfams)):
            _, prec, reca, correct = score_channel(pred, truth)
            row[f"{key}_precision"] = round(prec, 3)
            row[f"{key}_recall"] = round(reca, 3)
            row[f"{key}_correct"] = correct
        row["identity_bin"] = identity_bin(row["seq_pident"])
        row["tm_bin"] = tm_bin(row["struct_tm"])
        results.append(row)

    fields = ["accession", "phylum", "length", "truth_pfams",
              "seq_pident", "seq_pfams", "struct_tm", "struct_plddt",
              "struct_pfams", "union_pfams", "consensus_pfams",
              "identity_bin", "tm_bin"] + \
             [f"{k}_{m}" for k in ("seq", "struct", "union", "consensus")
              for m in ("precision", "recall", "correct")]
    with (args.output / "holdout_results.tsv").open(
            "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    summary = [summarize(results, k) for k in
               ("seq", "struct", "union", "consensus")]
    with (args.output / "benchmark_summary.tsv").open(
            "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle,
                                fieldnames=list(summary[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    # Breakdown by identity: the gain of the structure channel in the remote
    # homology regime.
    by_identity: Dict[str, list] = defaultdict(list)
    for row in results:
        by_identity[row["identity_bin"]].append(row)
    with (args.output / "metrics_by_identity.tsv").open(
            "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["identity_bin", "n", "seq_any_correct",
                         "struct_any_correct", "union_any_correct"])
        for bin_name in ("<20%", "20-30%", "30-50%", ">=50%", "no hit"):
            rows = by_identity.get(bin_name, [])
            if not rows:
                continue
            n = len(rows)
            writer.writerow([
                bin_name, n,
                round(sum(r["seq_correct"] for r in rows) / n, 3),
                round(sum(r["struct_correct"] for r in rows) / n, 3),
                round(sum(r["union_correct"] for r in rows) / n, 3),
            ])

    by_phylum: Dict[str, list] = defaultdict(list)
    for row in results:
        by_phylum[row["phylum"]].append(row)
    with (args.output / "metrics_by_phylum.tsv").open(
            "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["phylum", "n", "seq_correct", "struct_correct",
                         "union_correct"])
        for phylum, rows in sorted(by_phylum.items(),
                                   key=lambda kv: -len(kv[1])):
            n = len(rows)
            writer.writerow([
                phylum, n,
                round(sum(r["seq_correct"] for r in rows) / n, 3),
                round(sum(r["struct_correct"] for r in rows) / n, 3),
                round(sum(r["union_correct"] for r in rows) / n, 3),
            ])

    (args.output / "benchmark_summary.json").write_text(
        json.dumps({"n_proteins": len(results), "channels": summary},
                   indent=2), encoding="utf-8")
    print("Benchmark summary:")
    for item in summary:
        print(f"  {item['channel']:9s} coverage={item['coverage']:.3f} "
              f"any_correct={item['any_correct_rate']:.3f} "
              f"precision={item['mean_precision']:.3f} "
              f"recall={item['mean_recall']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
