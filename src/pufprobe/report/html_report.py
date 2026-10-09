"""Self-contained HTML/Markdown report rendering.

The HTML uses only inline CSS: no external stylesheets, no JavaScript, no
frameworks.  Every interpolated string is passed through :func:`html.escape`
so database labels containing ``<``, ``&`` or quotes can never break the page.
"""

from __future__ import annotations

import html
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from pufprobe.models import CDSRecord, ProteinPrediction

_CSS = """
body { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
       margin: 2rem; color: #1f2933; line-height: 1.45; }
h1 { font-size: 1.5rem; border-bottom: 2px solid #2c5f8a; padding-bottom: .3rem; }
h2 { font-size: 1.15rem; margin-top: 2rem; color: #2c5f8a; }
h3 { font-size: 1rem; margin-top: 1.5rem; }
table { border-collapse: collapse; width: 100%; margin: .6rem 0 1.4rem; font-size: .85rem; }
th, td { border: 1px solid #cbd2d9; padding: .35rem .5rem; text-align: left;
         vertical-align: top; }
th { background: #eef3f8; }
tr:nth-child(even) td { background: #fafbfc; }
.num { text-align: right; font-variant-numeric: tabular-nums; }
.muted { color: #7b8794; }
.tag { display: inline-block; background: #e3edf6; border-radius: 3px;
       padding: 0 .4rem; margin-right: .2rem; }
"""


def _e(text: Any) -> str:
    return html.escape(str(text if text is not None else ""), quote=True)


def _fmt(value: Any, nd: int = 3) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.{nd}f}"
    return str(value)


def compute_overview(
    records: Sequence[CDSRecord],
    predictions: Sequence[ProteinPrediction],
    high_confidence: float,
) -> Dict[str, Any]:
    """Aggregate counts used by both the HTML and Markdown reports."""
    verdict_counts: Dict[str, int] = {}
    puf = 0
    for rec in records:
        if rec.status.value in ("hypothetical", "conserved_hypothetical"):
            puf += 1
        if rec.verdict is not None:
            verdict_counts[rec.verdict.value] = verdict_counts.get(rec.verdict.value, 0) + 1
    proteins_with_terms = sum(1 for p in predictions if p.terms)
    high_conf_terms = sum(1 for p in predictions for t in p.terms if t.confidence >= high_confidence)
    return {
        "samples": len({r.sample_id for r in records}),
        "cds_total": len(records),
        "puf_total": puf,
        "verdict_counts": verdict_counts,
        "proteins_with_terms": proteins_with_terms,
        "high_confidence_terms": high_conf_terms,
    }


def _record_map(records: Sequence[CDSRecord]) -> Dict[str, CDSRecord]:
    return {r.key: r for r in records}


def render_html(
    overview: Dict[str, Any],
    predictions: Sequence[ProteinPrediction],
    records: Sequence[CDSRecord],
    priority_rows: List[Tuple[str, float, List[str]]],
    high_confidence: float,
    generated_at: Optional[str] = None,
) -> str:
    """Return a complete self-contained HTML document."""
    generated_at = generated_at or datetime.now().isoformat(timespec="seconds")
    rmap = _record_map(records)
    parts: List[str] = [
        "<!DOCTYPE html>",
        '<html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>PUFprobe evidence report</title><style>{_CSS}</style></head><body>",
        "<h1>PUFprobe evidence report</h1>",
        f'<p class="muted">generated: {_e(generated_at)} &middot; '
        f"high-confidence threshold: {_e(f'{high_confidence:.2f}')}</p>",
        "<h2>Overview</h2><table><tr><th>Metric</th><th class='num'>Value</th></tr>",
    ]

    def _row(label: str, value: Any) -> None:
        parts.append(f"<tr><td>{_e(label)}</td><td class='num'>{_e(value)}</td></tr>")

    _row("Samples", overview["samples"])
    _row("CDS records", overview["cds_total"])
    _row("PUF records (hypothetical / conserved hypothetical)", overview["puf_total"])
    for verdict, count in sorted(overview["verdict_counts"].items()):
        _row(f"Verdict: {verdict}", count)
    _row("Proteins with at least one proposed term", overview["proteins_with_terms"])
    _row(f"High-confidence terms (&gt;= {high_confidence:.2f})", overview["high_confidence_terms"])
    parts.append("</table>")

    anchored = [p for p in predictions if p.terms]
    parts.append("<h2>Per-protein evidence</h2>")
    if not anchored:
        parts.append("<p class='muted'>No protein received a proposed term.</p>")
    for pred in anchored:
        rec = rmap.get(pred.key)
        product = rec.product if rec else ""
        status = rec.status.value if rec else ""
        parts.append(
            f"<h3>{_e(pred.key)} <span class='tag'>{_e(status)}</span></h3>"
        )
        if product:
            parts.append(f"<p class='muted'>original product: {_e(product)}</p>")
        parts.append(
            "<table><tr><th>Type</th><th>ID</th><th>Name</th><th class='num'>raw</th>"
            "<th class='num'>confidence</th><th class='num'>evidence groups</th></tr>"
        )
        for term in pred.terms:
            groups = ", ".join(term.evidence_groups)
            parts.append(
                "<tr>"
                f"<td>{_e(term.term_type)}</td>"
                f"<td>{_e(term.term_id)}</td>"
                f"<td>{_e(term.term_name)}</td>"
                f"<td class='num'>{_fmt(term.raw_score)}</td>"
                f"<td class='num'>{_fmt(term.confidence)}</td>"
                f"<td class='num'>{_e(groups)}</td>"
                "</tr>"
            )
        parts.append("</table>")
        parts.append(
            "<table><tr><th>Channel</th><th>Group</th><th>Target id</th><th>Label</th>"
            "<th class='num'>score</th><th class='num'>evalue</th><th class='num'>coverage</th></tr>"
        )
        seen = set()
        for term in pred.terms:
            for ev in term.evidence:
                key = (ev.channel.value, ev.target_id)
                if key in seen:
                    continue
                seen.add(key)
                parts.append(
                    "<tr>"
                    f"<td>{_e(ev.channel.value)}</td>"
                    f"<td>{_e(ev.group)}</td>"
                    f"<td>{_e(ev.target_id)}</td>"
                    f"<td>{_e(ev.target_label)}</td>"
                    f"<td class='num'>{_fmt(ev.score)}</td>"
                    f"<td class='num'>{_fmt(ev.evalue, 1) if ev.evalue is not None else ''}</td>"
                    f"<td class='num'>{_fmt(ev.coverage)}</td>"
                    "</tr>"
                )
        parts.append("</table>")

    parts.append("<h2>Unanchored families (wet-lab priority)</h2>")
    if priority_rows:
        parts.append(
            "<table><tr><th>Family</th><th class='num'>priority</th><th>Reasons</th></tr>"
        )
        for fid, score, reasons in priority_rows:
            parts.append(
                "<tr>"
                f"<td>{_e(fid)}</td>"
                f"<td class='num'>{_fmt(score, 4)}</td>"
                f"<td>{_e('; '.join(reasons))}</td>"
                "</tr>"
            )
        parts.append("</table>")
    else:
        parts.append("<p class='muted'>All families are anchored by at least one term.</p>")

    parts.append("</body></html>")
    return "\n".join(parts)


def render_markdown(
    overview: Dict[str, Any],
    predictions: Sequence[ProteinPrediction],
    records: Sequence[CDSRecord],
    priority_rows: List[Tuple[str, float, List[str]]],
    high_confidence: float,
    generated_at: Optional[str] = None,
) -> str:
    """Same content as :func:`render_html`, in GitHub-flavoured Markdown."""
    generated_at = generated_at or datetime.now().isoformat(timespec="seconds")
    rmap = _record_map(records)
    lines: List[str] = [
        "# PUFprobe evidence report",
        "",
        f"generated: {generated_at}  ",
        f"high-confidence threshold: {high_confidence:.2f}",
        "",
        "## Overview",
        "",
        "| Metric | Value |",
        "| --- | --- |",
        f"| Samples | {overview['samples']} |",
        f"| CDS records | {overview['cds_total']} |",
        f"| PUF records | {overview['puf_total']} |",
    ]
    for verdict, count in sorted(overview["verdict_counts"].items()):
        lines.append(f"| Verdict: {verdict} | {count} |")
    lines += [
        f"| Proteins with at least one proposed term | {overview['proteins_with_terms']} |",
        f"| High-confidence terms | {overview['high_confidence_terms']} |",
        "",
        "## Per-protein evidence",
        "",
    ]
    anchored = [p for p in predictions if p.terms]
    if not anchored:
        lines.append("_No protein received a proposed term._")
    for pred in anchored:
        rec = rmap.get(pred.key)
        lines.append(f"### {pred.key}")
        if rec is not None:
            lines.append(f"- status: `{rec.status.value}`")
            if rec.product:
                lines.append(f"- original product: {rec.product}")
        lines += [
            "",
            "| Type | ID | Name | raw | confidence | evidence groups |",
            "| --- | --- | --- | ---: | ---: | --- |",
        ]
        for term in pred.terms:
            lines.append(
                f"| {term.term_type} | {term.term_id} | {term.term_name} "
                f"| {term.raw_score:.3f} | {term.confidence:.3f} "
                f"| {', '.join(term.evidence_groups)} |"
            )
        lines += [
            "",
            "| Channel | Group | Target id | Label | score | evalue | coverage |",
            "| --- | --- | --- | ---: | ---: | ---: |",
        ]
        seen = set()
        for term in pred.terms:
            for ev in term.evidence:
                key = (ev.channel.value, ev.target_id)
                if key in seen:
                    continue
                seen.add(key)
                evalue = f"{ev.evalue:.1e}" if ev.evalue is not None else ""
                coverage = f"{ev.coverage:.3f}" if ev.coverage is not None else ""
                lines.append(
                    f"| {ev.channel.value} | {ev.group} | {ev.target_id} | {ev.target_label} "
                    f"| {ev.score:.3f} | {evalue} | {coverage} |"
                )
        lines.append("")

    lines += ["## Unanchored families (wet-lab priority)", ""]
    if priority_rows:
        lines += [
            "| Family | priority | reasons |",
            "| --- | ---: | --- |",
        ]
        for fid, score, reasons in priority_rows:
            lines.append(f"| {fid} | {score:.4f} | {'; '.join(reasons)} |")
    else:
        lines.append("_All families are anchored by at least one term._")
    lines.append("")
    return "\n".join(lines)
