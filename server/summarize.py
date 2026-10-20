#!/usr/bin/env python3
"""Summarize one or more PUFprobe working directories.

Reads:
  07_rank/predictions.json, 07_rank/priorities.tsv, 07_rank/calibration.json
  05_search/search_evidence.json, 06_context/context_evidence.json
  03_cluster/families.json

Writes per-run CSV tables and a combined Markdown overview:
  per_sample_summary.csv, channel_hits.csv, family_overview.csv,
  wetlab_priorities_topN.csv, summary.md
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List


def load_json(path: Path, default: Any) -> Any:
    if not path.exists():
        return default
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def iter_workdirs(paths: List[str]) -> List[Path]:
    workdirs: List[Path] = []
    for item in paths:
        path = Path(item)
        if path.is_dir() and (path / "07_rank" / "predictions.json").exists():
            workdirs.append(path)
        elif path.is_dir():
            workdirs.extend(
                sorted(path.parent for path in path.glob("*/07_rank/predictions.json"))
            )
    return workdirs


def count_channels(evidence_map: Dict[str, List[dict]]) -> Counter:
    counts: Counter = Counter()
    for evidences in evidence_map.values():
        for evidence in evidences:
            counts[evidence.get("channel", "unknown")] += 1
    return counts


CHANNEL_GROUPS = {
    "plm_neighbor": "sequence",
    "foldseek": "structure",
    "remote_homology": "sequence_profile",
    "pocket": "structure",
    "context": "context",
    "database": "database",
}


def groups_of(term: dict) -> str:
    groups = sorted(
        {
            CHANNEL_GROUPS.get(item.get("channel", ""), item.get("channel", ""))
            for item in term.get("evidence", [])
        }
    )
    return ",".join(group for group in groups if group)


def summarize_run(workdir: Path, high_confidence: float) -> Dict[str, Any]:
    predictions = load_json(workdir / "07_rank" / "predictions.json", [])
    search_evidence = load_json(workdir / "05_search" / "search_evidence.json", {})
    context_evidence = load_json(workdir / "06_context" / "context_evidence.json", {})
    families = load_json(workdir / "03_cluster" / "families.json", [])
    calibration = load_json(workdir / "07_rank" / "calibration.json", {})
    threshold = calibration.get("high_confidence", high_confidence)

    per_sample: Dict[str, Dict[str, int]] = defaultdict(
        lambda: {"n_proteins": 0, "n_terms": 0, "n_high_conf": 0}
    )
    high_conf_rows: List[dict] = []
    for protein in predictions:
        sample = protein.get("sample_id", "?")
        terms = protein.get("terms", [])
        if not terms:
            continue
        row = per_sample[sample]
        row["n_proteins"] += 1
        row["n_terms"] += len(terms)
        for term in terms:
            if term.get("confidence", 0.0) >= threshold:
                row["n_high_conf"] += 1
                high_conf_rows.append(
                    {
                        "run": workdir.name,
                        "sample_id": sample,
                        "locus_tag": protein.get("locus_tag", ""),
                        "term_type": term.get("term_type", ""),
                        "term_id": term.get("term_id", ""),
                        "term_name": term.get("term_name", ""),
                        "confidence": f"{term.get('confidence', 0.0):.3f}",
                        "raw_score": f"{term.get('raw_score', 0.0):.3f}",
                        "evidence_groups": groups_of(term),
                    }
                )

    channel_counts = count_channels(search_evidence)
    channel_counts.update(count_channels(context_evidence))

    family_rows = []
    for family in families:
        family_rows.append(
            {
                "run": workdir.name,
                "family_id": family.get("family_id", ""),
                "size": family.get("size", 0),
                "lineage_width": family.get("lineage_width", 0),
                "fold_status": family.get("fold_status", ""),
                "plddt": (
                    f"{family['plddt']:.1f}" if family.get("plddt") is not None else ""
                ),
                "priority_score": (
                    f"{family['priority_score']:.4f}"
                    if family.get("priority_score") is not None
                    else ""
                ),
            }
        )

    priority_rows = []
    priorities_path = workdir / "07_rank" / "priorities.tsv"
    if priorities_path.exists():
        with priorities_path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                priority_rows.append(
                    {
                        "run": workdir.name,
                        "family_id": row.get("family_id", ""),
                        "priority_score": row.get("priority_score", ""),
                        "reasons": row.get("reasons", ""),
                    }
                )

    return {
        "workdir": workdir,
        "threshold": threshold,
        "per_sample": per_sample,
        "high_conf_rows": high_conf_rows,
        "channel_counts": channel_counts,
        "family_rows": family_rows,
        "priority_rows": priority_rows,
        "calibration": calibration,
        "n_predictions": len(predictions),
    }


def write_csv(path: Path, fieldnames: List[str], rows: Iterable[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize PUFprobe runs")
    parser.add_argument(
        "--workdir",
        action="append",
        default=None,
        help="PUFprobe working directory, or a directory containing runs; repeatable",
    )
    parser.add_argument("--runs-root", default=None, help="Directory holding all runs")
    parser.add_argument("--out", required=True, help="Output directory for tables")
    parser.add_argument("--top", type=int, default=30, help="Top-N wet-lab priorities")
    parser.add_argument(
        "--high-confidence",
        type=float,
        default=0.80,
        help="Fallback threshold if calibration.json lacks one",
    )
    args = parser.parse_args()

    targets = list(args.workdir or [])
    if args.runs_root:
        targets.append(args.runs_root)
    if not targets:
        parser.error("provide at least one --workdir or --runs-root")

    workdirs = iter_workdirs(targets)
    if not workdirs:
        raise SystemExit(f"No completed runs found in: {targets}")

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_rows: List[dict] = []
    channel_rows: List[dict] = []
    family_rows: List[dict] = []
    priority_rows: List[dict] = []
    high_conf_rows: List[dict] = []
    md_lines = ["# PUFprobe run summary", ""]

    for workdir in workdirs:
        result = summarize_run(workdir, args.high_confidence)
        run_name = workdir.name
        for sample, counts in sorted(result["per_sample"].items()):
            sample_rows.append({"run": run_name, "sample_id": sample, **counts})
        for channel, count in sorted(result["channel_counts"].items()):
            channel_rows.append(
                {"run": run_name, "channel": channel, "n_evidence": count}
            )
        family_rows.extend(result["family_rows"])
        high_conf_rows.extend(result["high_conf_rows"])
        top = result["priority_rows"][: args.top]
        priority_rows.extend(top)

        md_lines.append(f"## {run_name}")
        cal = result["calibration"]
        md_lines.append("")
        md_lines.append(
            f"- Proteins with predictions: {result['n_predictions']}; "
            f"total functional terms: {cal.get('n_terms', 'NA')}; "
            f"anchorless families: {cal.get('n_unanchored_families', 'NA')}"
        )
        md_lines.append(
            f"- High-confidence threshold: {result['threshold']}; "
            f"independent evidence groups required: {cal.get('require_independent_groups', 'NA')}"
        )
        md_lines.append("- Evidence count per channel: " + (
            ", ".join(
                f"{channel} {count}"
                for channel, count in sorted(result["channel_counts"].items())
            )
            or "none"
        ))
        md_lines.append("")
        md_lines.append("| Strain | Proteins with predictions | Functional terms | High-confidence terms |")
        md_lines.append("| --- | --- | --- | --- |")
        for sample, counts in sorted(result["per_sample"].items()):
            md_lines.append(
                f"| {sample} | {counts['n_proteins']} | {counts['n_terms']} "
                f"| {counts['n_high_conf']} |"
            )
        md_lines.append("")
        if top:
            md_lines.append(f"Top {len(top)} wet-lab priority candidates:")
            md_lines.append("")
            md_lines.append("| Family | Size | Priority | Reasons |")
            md_lines.append("| --- | --- | --- | --- |")
            family_size = {row["family_id"]: row["size"] for row in result["family_rows"]}
            for row in top:
                md_lines.append(
                    f"| {row['family_id']} | {family_size.get(row['family_id'], '')} "
                    f"| {row['priority_score']} | {row['reasons']} |"
                )
        md_lines.append("")

    write_csv(
        out_dir / "per_sample_summary.csv",
        ["run", "sample_id", "n_proteins", "n_terms", "n_high_conf"],
        sample_rows,
    )
    write_csv(
        out_dir / "channel_hits.csv",
        ["run", "channel", "n_evidence"],
        channel_rows,
    )
    write_csv(
        out_dir / "family_overview.csv",
        [
            "run",
            "family_id",
            "size",
            "lineage_width",
            "fold_status",
            "plddt",
            "priority_score",
        ],
        family_rows,
    )
    write_csv(
        out_dir / "high_confidence_predictions.csv",
        [
            "run",
            "sample_id",
            "locus_tag",
            "term_type",
            "term_id",
            "term_name",
            "confidence",
            "raw_score",
            "evidence_groups",
        ],
        high_conf_rows,
    )
    write_csv(
        out_dir / f"wetlab_priorities_top{args.top}.csv",
        ["run", "family_id", "priority_score", "reasons"],
        priority_rows,
    )

    (out_dir / "summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"Summary written to {out_dir}")
    for name in (
        "per_sample_summary.csv",
        "channel_hits.csv",
        "family_overview.csv",
        "high_confidence_predictions.csv",
        f"wetlab_priorities_top{args.top}.csv",
        "summary.md",
    ):
        print(f"  - {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
