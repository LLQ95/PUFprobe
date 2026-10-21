#!/usr/bin/env python3
"""Fetch cross-species hypothetical proteins from UniProt reference proteomes.

The script lists bacterial reference proteomes, samples them by phylum so the
major bacterial lineages are covered, downloads each sampled proteome, and keeps
the proteins that are still uncharacterized. Every proteome is stored in a
separate file, so the download is resumable and can be re-run safely.

Only protein sequences and metadata are downloaded. No assemblies or reads are
fetched, and predicted structures are not computed here.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import random
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

REST_BASE = "https://rest.uniprot.org"
USER_AGENT = "PUFprobe-reference-fetch/1.0"
PAGE_SIZE = 500

GENERIC_NAMES = {
    "hypothetical protein",
    "uncharacterized protein",
    "putative protein",
    "predicted protein",
    "conserved protein",
    "conserved hypothetical protein",
    "uncharacterized conserved protein",
    "protein of unknown function",
}

DUF_PATTERN = re.compile(r"\b(?:duf|upf)\d+\b", re.IGNORECASE)

# Recognised bacterial phylum names, covering both GTDB (-ota) and legacy
# (-ia/-ae) nomenclature. The phylum is located by matching the lineage against
# this set, because the new taxonomy inserts a subkingdom rank (e.g.
# Pseudomonadati) and occasionally an intermediate "group" rank between
# Bacteria and the phylum.
KNOWN_PHYLA = {
    # GTDB / current names
    "acidobacteriota", "actinobacteriota", "actinomycetota", "aquificota",
    "armatimonadota", "atribacterota", "bdellovibrionota", "bacteroidota",
    "binatota", "caldisericota", "calditrichota", "campylobacterota",
    "chlamydiota", "chlorobiota", "chloroflexota", "coprothermobacterota",
    "cyanobacteriota", "deferribacterota", "deinococcota", "dependentiae",
    "desulfobacterota", "dictyoglomota", "elusimicrobiota", "fibrobacterota",
    "bacillota", "fusobacteriota", "gemmatimonadota", "ignavibacteriota",
    "kiritimatiellota", "krumholzibacteriota", "lentisphaerota",
    "margulisbacteria", "methylomirabilota", "modulibacteria", "myxococcota",
    "nitrospinota", "nitrospirota", "planctomycetota", "poribacteria",
    "pseudomonadota", "rhodothermota", "riflebacteria", "spirochaetota",
    "synergistota", "thermotogota", "verrucomicrobiota", "zixibacteria",
    "saltatorellota", "abditibacteriota", "latescibacterota",
    "coprothermobacterota", "desulfuromonadota", "fervidibacterota",
    "goldbacteria", "krumholzibacteriota", "lindowiibacteriota",
    "muiribacteriota", "neocallimastigomycota", "sumerlaeota",
    "wallbacteria", "worodarchiaota",
    # Legacy names still used by UniProt
    "acidobacteria", "actinobacteria", "aquificae", "bacteroidetes",
    "chlamydiae", "chlorobi", "chloroflexi", "cyanobacteria",
    "deferribacteres", "deinococcus-thermus", "elusimicrobia",
    "fibrobacteres", "firmicutes", "fusobacteria", "gemmatimonadetes",
    "lentisphaerae", "nitrospirae", "planctomycetes", "proteobacteria",
    "spirochaetes", "synergistetes", "thermotogae", "verrucomicrobia",
}

# Subkingdom ranks that appear directly after Bacteria and must be skipped
# when falling back to positional phylum extraction.
SKIP_RANKS = {
    "pseudomonadati", "bacillati", "actinobacteriobionta", "glissandra",
    "rhodobionta", "fusobacteriobionta", "thermotogati", "spirochaetoti",
    "chlamydiobionta", "deinococcobionta", "synergistobionta",
    "chloroflexobionta", "bacteroidobionta", "cyanobacteriobionta",
}


def is_hypothetical(protein_name: str) -> bool:
    """Return True when the protein name describes an uncharacterized product."""
    name = protein_name.lower().strip().rstrip(".")
    if not name:
        return True
    if name in GENERIC_NAMES:
        return True
    if "hypothetical" in name or "uncharacterized" in name:
        return True
    if "unknown function" in name:
        return True
    if DUF_PATTERN.search(name):
        return True
    return False


def request_page(endpoint: str, parameters: dict,
                 retries: int = 4) -> tuple[list[dict], str | None]:
    """Request one page and return the rows together with the next cursor.

    The cursor is read from the Link header (rel="next"). Gzip responses are
    decompressed transparently.
    """
    url = f"{REST_BASE}/{endpoint}/search?" + urllib.parse.urlencode(parameters)
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                raw = response.read()
                # The compressed=true parameter always returns gzip bytes, even
                # when the Content-Encoding header is absent, so detect gzip by
                # its magic number rather than relying on the header.
                if raw[:2] == b"\x1f\x8b":
                    text = gzip.decompress(raw).decode("utf-8", errors="replace")
                else:
                    text = raw.decode("utf-8", errors="replace")
                link = response.headers.get("Link", "")
            break
        except Exception as error:  # noqa: BLE001 - retry any network error
            last_error = error
            if attempt < retries:
                time.sleep(3 * attempt)
            else:
                raise RuntimeError(f"download failed: {url}\n{last_error}") from error

    next_cursor: str | None = None
    match = re.search(r'<([^>]+)>;\s*rel="next"', link)
    if match:
        cursor_match = re.search(r"cursor=([^&>]+)", match.group(1))
        if cursor_match:
            next_cursor = urllib.parse.unquote(cursor_match.group(1))

    handle = io.StringIO(text)
    header_line = handle.readline().rstrip("\n")
    if not header_line:
        return [], next_cursor
    header = header_line.split("\t")
    rows: list[dict] = []
    for line in handle:
        values = line.rstrip("\n").split("\t")
        values += [""] * (len(header) - len(values))
        rows.append(dict(zip(header, values)))
    return rows, next_cursor


def fetch_all(endpoint: str, query: str, fields: str) -> list[dict]:
    """Page through the search endpoint until every row is collected."""
    rows: list[dict] = []
    cursor: str | None = None
    while True:
        parameters = {
            "query": query,
            "fields": fields,
            "format": "tsv",
            "size": str(PAGE_SIZE),
            "compressed": "true",
        }
        if cursor:
            parameters["cursor"] = cursor
        page, next_cursor = request_page(endpoint, parameters)
        rows.extend(page)
        if not next_cursor:
            break
        cursor = next_cursor
    return rows


def phylum_of(lineage: str) -> str:
    """Extract the phylum from a comma-separated lineage string.

    The first rank that matches a known phylum name is returned. If no rank
    matches, the rank immediately after Bacteria (or Archaea) is used, skipping
    the subkingdom rank and any intermediate "group" ranks.
    """
    parts = [part.strip() for part in lineage.split(",") if part.strip()]
    for part in parts:
        if part.lower() in KNOWN_PHYLA:
            return part
    try:
        root_index = next(
            index for index, part in enumerate(parts)
            if part.lower() in {"bacteria", "archaea"}
        )
        for part in parts[root_index + 1:]:
            lowered = part.lower()
            if lowered in SKIP_RANKS or lowered.endswith(" group"):
                continue
            return part
    except StopIteration:
        pass
    return parts[0] if parts else "Unknown"


def sample_proteomes(
    proteomes: list[dict], target: int, min_per_phylum: int, seed: int
) -> list[dict]:
    """Sample proteomes by phylum with a floor for small phyla."""
    by_phylum: dict[str, list[dict]] = {}
    for proteome in proteomes:
        phylum = phylum_of(proteome.get("Taxonomic lineage", ""))
        by_phylum.setdefault(phylum, []).append(proteome)

    rng = random.Random(seed)
    selected: list[dict] = []
    for members in by_phylum.values():
        members = list(members)
        rng.shuffle(members)
        selected.extend(members[:min(min_per_phylum, len(members))])

    chosen_ids = {item.get("Proteome Id") for item in selected}
    remaining = target - len(selected)
    if remaining > 0:
        pool = [item for item in proteomes if item.get("Proteome Id") not in chosen_ids]
        rng.shuffle(pool)
        selected.extend(pool[:remaining])
    return selected[:target]


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(field, "")) for field in fields) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True,
                        help="directory for the cross-species dataset")
    parser.add_argument("--max-proteomes", type=int, default=200,
                        help="number of reference proteomes to sample")
    parser.add_argument("--min-per-phylum", type=int, default=2,
                        help="minimum proteomes sampled for each phylum")
    parser.add_argument("--seed", type=int, default=42,
                        help="random seed for reproducible sampling")
    args = parser.parse_args()

    output = args.output
    per_proteome = output / "per_proteome"
    per_proteome.mkdir(parents=True, exist_ok=True)

    index_path = output / "proteomes_index.tsv"
    sampled_path = output / "sampled_proteomes.tsv"
    index_fields = ["Proteome Id", "Organism", "Organism Id",
                    "Protein count", "Taxonomic lineage"]

    if index_path.exists() and index_path.stat().st_size > 0:
        with index_path.open(encoding="utf-8", newline="") as handle:
            proteomes = list(csv.DictReader(handle, delimiter="\t"))
        print(f"Reusing cached proteome index ({len(proteomes)} entries).")
    else:
        print("Listing bacterial reference proteomes ...")
        query = "(reference:true) AND (superkingdom:Bacteria)"
        proteome_fields = "upid,organism,organism_id,protein_count,lineage"
        proteomes = fetch_all("proteomes", query, proteome_fields)
        if not proteomes:
            print("No reference proteomes returned; check the network connection.",
                  file=sys.stderr)
            return 1
        write_tsv(index_path, proteomes, index_fields)
    print(f"  {len(proteomes)} reference proteomes indexed")

    sampled = sample_proteomes(proteomes, args.max_proteomes,
                               args.min_per_phylum, args.seed)
    write_tsv(sampled_path, sampled, index_fields)
    phylum_counts: dict[str, int] = {}
    for item in sampled:
        phylum = phylum_of(item.get("Taxonomic lineage", ""))
        phylum_counts[phylum] = phylum_counts.get(phylum, 0) + 1
    print(f"  sampled {len(sampled)} proteomes across {len(phylum_counts)} phyla")

    protein_fields = ("accession,protein_name,sequence,length,annotation_score,"
                      "protein_existence,organism_id,organism_name")
    output_fields = ["Entry", "Protein names", "Sequence", "Length",
                     "Annotation", "Protein existence", "Organism (ID)",
                     "Organism", "Proteome", "Phylum"]

    total_hypothetical = 0
    for proteome in sampled:
        upid = proteome.get("Proteome Id")
        phylum = phylum_of(proteome.get("Taxonomic lineage", ""))
        proteome_file = per_proteome / f"{upid}.tsv"
        if proteome_file.exists() and proteome_file.stat().st_size > 0:
            with proteome_file.open(encoding="utf-8") as handle:
                total_hypothetical += sum(1 for _ in handle) - 1
            continue

        proteins = fetch_all("uniprotkb", f"(proteome:{upid})", protein_fields)
        kept: list[dict] = []
        for protein in proteins:
            if is_hypothetical(protein.get("Protein names", "")):
                row = dict(protein)
                row["Proteome"] = upid
                row["Phylum"] = phylum
                kept.append(row)
        write_tsv(proteome_file, kept, output_fields)
        total_hypothetical += len(kept)
        print(f"  {upid} ({proteome.get('Organism', '?')}): "
              f"{len(proteins)} proteins, {len(kept)} hypothetical")

    merged_path = output / "all_hypothetical_proteins.tsv"
    with merged_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(output_fields) + "\n")
        for proteome_file in sorted(per_proteome.glob("*.tsv")):
            with proteome_file.open(encoding="utf-8") as source:
                source.readline()
                for line in source:
                    handle.write(line)

    summary = {
        "reference_proteomes": len(proteomes),
        "sampled_proteomes": len(sampled),
        "phylum_counts": phylum_counts,
        "hypothetical_proteins": total_hypothetical,
        "merged_table": str(merged_path),
    }
    (output / "fetch_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Done. {total_hypothetical} hypothetical proteins written to {merged_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
