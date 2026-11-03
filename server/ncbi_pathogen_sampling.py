#!/usr/bin/env python3
"""Stratified sampling of bacterial pathogen genomes from NCBI.

The script resolves one high-quality RefSeq assembly for each target pathogen
and downloads the genomic sequence, protein sequences and GFF3 annotation.
Targets are chosen so that Gram-positive and Gram-negative organisms are both
represented and so that the major bacterial phyla are covered. The resolved
assembly list is cached, and every assembly is written to its own directory, so
interrupted runs can be resumed and the download can be repeated safely.

Only public RefSeq data are retrieved through the NCBI E-utilities and the NCBI
HTTPS file server. No credentials are embedded in the script. An NCBI API key,
when available, is read from the NCBI_API_KEY environment variable or the
--api-key option.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP_HTTPS_ROOT = "https://ftp.ncbi.nlm.nih.gov"
USER_AGENT = "PUFprobe-ncbi-sampling/1.0"

# RefSeq categories are ranked so that a reference genome is preferred over a
# representative genome, and assemblies without a category are considered last.
CATEGORY_RANK: Dict[str, int] = {
    "reference genome": 0,
    "representative genome": 1,
}
UNRANKED_CATEGORY = 2

# Each entry holds the species name, NCBI taxonomy id, Gram stain and phylum.
# The stain follows the conventional classification used in clinical
# microbiology; acid-fast actinobacteria are grouped with Gram-positive taxa.
TARGETS: List[Dict[str, str]] = [
    # Gram-negative organisms.
    {"species": "Escherichia coli", "taxid": "562",
     "gram": "negative", "phylum": "Pseudomonadota"},
    {"species": "Salmonella enterica", "taxid": "28901",
     "gram": "negative", "phylum": "Pseudomonadota"},
    {"species": "Klebsiella pneumoniae", "taxid": "573",
     "gram": "negative", "phylum": "Pseudomonadota"},
    {"species": "Pseudomonas aeruginosa", "taxid": "287",
     "gram": "negative", "phylum": "Pseudomonadota"},
    {"species": "Acinetobacter baumannii", "taxid": "470",
     "gram": "negative", "phylum": "Pseudomonadota"},
    {"species": "Vibrio cholerae", "taxid": "243277",
     "gram": "negative", "phylum": "Pseudomonadota"},
    {"species": "Neisseria meningitidis", "taxid": "487",
     "gram": "negative", "phylum": "Pseudomonadota"},
    {"species": "Campylobacter jejuni", "taxid": "197",
     "gram": "negative", "phylum": "Campylobacterota"},
    {"species": "Bacteroides fragilis", "taxid": "817",
     "gram": "negative", "phylum": "Bacteroidota"},
    # Gram-positive organisms.
    {"species": "Staphylococcus aureus", "taxid": "1280",
     "gram": "positive", "phylum": "Bacillota"},
    {"species": "Streptococcus pyogenes", "taxid": "1314",
     "gram": "positive", "phylum": "Bacillota"},
    {"species": "Listeria monocytogenes", "taxid": "1639",
     "gram": "positive", "phylum": "Bacillota"},
    {"species": "Enterococcus faecalis", "taxid": "1351",
     "gram": "positive", "phylum": "Bacillota"},
    {"species": "Clostridium perfringens", "taxid": "1502",
     "gram": "positive", "phylum": "Bacillota"},
    {"species": "Mycobacterium tuberculosis", "taxid": "1773",
     "gram": "positive", "phylum": "Actinobacteriota"},
    {"species": "Corynebacterium diphtheriae", "taxid": "1717",
     "gram": "positive", "phylum": "Actinobacteriota"},
]

# Gzip-compressed NCBI files together with the plain output name used locally.
FILE_SPECS: List[Tuple[str, str]] = [
    ("_genomic.fna.gz", "genomic.fna"),
    ("_protein.faa.gz", "protein.faa"),
    ("_genomic.gff.gz", "genomic.gff"),
]


def http_request(url: str, timeout: int = 120,
                 retries: int = 4) -> bytes:
    """Retrieve a URL as bytes with exponential backoff on failure."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Optional[BaseException] = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            last_error = error
            if attempt < retries:
                time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"download failed after {retries} attempts: {url}\n"
                       f"{last_error}")


class EutilsClient:
    """Small NCBI E-utilities client with polite rate limiting."""

    def __init__(self, api_key: Optional[str] = None) -> None:
        self.api_key = api_key
        # With an API key NCBI allows ten requests per second; stay below that.
        self.pause = 0.12 if api_key else 0.35

    def call(self, tool: str, parameters: Dict[str, Any]) -> Dict[str, Any]:
        query = dict(parameters)
        query["retmode"] = "json"
        if self.api_key:
            query["api_key"] = self.api_key
        url = f"{EUTILS_BASE}/{tool}.fcgi?" + urllib.parse.urlencode(query)
        time.sleep(self.pause)
        raw = http_request(url)
        return json.loads(raw.decode("utf-8", errors="replace"))

    def search(self, term: str, retmax: int = 200) -> List[str]:
        result = self.call("esearch", {
            "db": "assembly", "term": term, "retmax": str(retmax)})
        return list(result["esearchresult"].get("idlist", []))

    def summarize(self, uids: Sequence[str]) -> Dict[str, Any]:
        if not uids:
            return {"result": {}}
        result = self.call("esummary", {
            "db": "assembly", "id": ",".join(uids)})
        return result


def category_rank(category: str) -> int:
    return CATEGORY_RANK.get((category or "").strip().lower(), UNRANKED_CATEGORY)


def ftp_to_https(ftp_path: str) -> str:
    """Convert an NCBI FTP URL to the equivalent HTTPS URL."""
    if ftp_path.startswith("ftp://ftp.ncbi.nlm.nih.gov"):
        return FTP_HTTPS_ROOT + ftp_path[len("ftp://ftp.ncbi.nlm.nih.gov"):]
    if ftp_path.startswith("https://"):
        return ftp_path
    return ftp_path


def safe_token(text: str) -> str:
    """Build a filesystem-safe token from a species or assembly name."""
    cleaned = "".join(ch if ch.isalnum() or ch in "._-" else "_"
                      for ch in text)
    while "__" in cleaned:
        cleaned = cleaned.replace("__", "_")
    return cleaned.strip("_")


def resolve_assembly(target: Dict[str, str],
                     client: EutilsClient) -> Dict[str, Any]:
    """Resolve the best RefSeq assembly for one target pathogen.

    Reference genomes are queried first, then representative genomes, and a
    complete RefSeq assembly is used as a final fallback. The default Organism
    field is used rather than the no-expansion form because reference genomes
    are indexed under strain-level taxids below the target species.
    """
    taxid = target["taxid"]
    queries = [
        f'txid{taxid}[Organism] AND "reference genome"[RCAT] '
        f'AND "latest refseq"[filter]',
        f'txid{taxid}[Organism] AND "representative genome"[RCAT] '
        f'AND "latest refseq"[filter]',
        f'txid{taxid}[Organism] AND "latest refseq"[filter] '
        f'AND "complete genome"[filter]',
    ]

    chosen: Optional[Dict[str, Any]] = None
    records: List[Dict[str, Any]] = []
    for query in queries:
        uids = client.search(query, retmax=50)
        if not uids:
            continue
        summary = client.summarize(uids)["result"]
        records = [
            summary[uid] for uid in uids
            if summary.get(uid) and summary[uid].get("assemblyaccession")
        ]
        if records:
            # Deterministic ordering: category rank, then accession, so the
            # earliest established assembly wins within the same category.
            records.sort(key=lambda record: (
                category_rank(str(record.get("refseq_category", ""))),
                str(record.get("assemblyaccession", ""))))
            chosen = records[0]
            break

    if chosen is None:
        raise RuntimeError(f"No RefSeq assembly found for {target['species']} "
                           f"(taxid {taxid})")
    ftp_path = chosen.get("ftppath_refseq") or chosen.get("ftppath_genbank")
    return {
        "species": target["species"],
        "taxid": taxid,
        "gram": target["gram"],
        "phylum": target["phylum"],
        "accession": chosen["assemblyaccession"],
        "assembly_name": chosen.get("assemblyname", ""),
        "refseq_category": chosen.get("refseq_category", ""),
        "assembly_status": chosen.get("assemblystatus", ""),
        "assembly_level": chosen.get("assembly_level", ""),
        "organism": chosen.get("organism", ""),
        "ftp_path": ftp_path,
        "candidate_count": len(records),
    }


def decompress_gzip(gz_path: Path, dest: Path) -> None:
    """Decompress a gzip file to a plain output path."""
    with gzip.open(gz_path, "rb") as source, \
            dest.open("wb") as output:
        shutil.copyfileobj(source, output)


def download_assembly(entry: Dict[str, Any], genomes_root: Path,
                      client: EutilsClient) -> Dict[str, Any]:
    """Download and verify the three files for one resolved assembly."""
    base_url = ftp_to_https(entry["ftp_path"])
    base_name = entry["ftp_path"].rstrip("/").rsplit("/", 1)[-1]
    sample_dir = genomes_root / (
        safe_token(entry["species"]) + "__" + safe_token(entry["accession"]))
    sample_dir.mkdir(parents=True, exist_ok=True)

    files: Dict[str, str] = {}
    for remote_suffix, output_name in FILE_SPECS:
        output_path = sample_dir / output_name
        if output_path.exists() and output_path.stat().st_size > 0:
            files[output_name] = str(output_path)
            continue
        remote_url = f"{base_url}/{base_name}{remote_suffix}"
        gz_path = sample_dir / (output_name + ".gz")
        if not (gz_path.exists() and gz_path.stat().st_size > 0):
            gz_path.write_bytes(http_request(remote_url))
        decompress_gzip(gz_path, output_path)
        gz_path.unlink()
        if output_path.stat().st_size == 0:
            raise RuntimeError(f"empty file downloaded: {remote_url}")
        files[output_name] = str(output_path)

    entry["sample_dir"] = str(sample_dir)
    entry["files"] = files
    return entry


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True,
                        help="directory for the sampled pathogen dataset")
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"),
                        help="NCBI API key (defaults to NCBI_API_KEY)")
    parser.add_argument("--manifest-only", action="store_true",
                        help="resolve assemblies without downloading files")
    args = parser.parse_args()

    output = args.output
    genomes_root = output / "genomes"
    genomes_root.mkdir(parents=True, exist_ok=True)

    manifest_path = output / "pathogen_manifest.json"
    client = EutilsClient(args.api_key)

    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        print(f"Reusing cached manifest ({len(manifest)} assemblies).")
    else:
        manifest = []
        for target in TARGETS:
            print(f"Resolving {target['species']} (taxid {target['taxid']}) ...")
            entry = resolve_assembly(target, client)
            print(f"  {entry['accession']} [{entry['refseq_category'] or 'na'}] "
                  f"from {entry['candidate_count']} candidates")
            manifest.append(entry)
        write_json(manifest_path, manifest)

    if not args.manifest_only:
        for entry in manifest:
            if entry.get("sample_dir") and all(
                    Path(path).exists() for path in entry.get("files", {}).values()):
                print(f"Skip completed {entry['accession']}")
                continue
            print(f"Downloading {entry['species']} {entry['accession']} ...")
            download_assembly(entry, genomes_root, client)
            write_json(manifest_path, manifest)
            print("  done")

    gram_counts: Dict[str, int] = {}
    phylum_counts: Dict[str, int] = {}
    for entry in manifest:
        gram_counts[entry["gram"]] = gram_counts.get(entry["gram"], 0) + 1
        phylum_counts[entry["phylum"]] = (
            phylum_counts.get(entry["phylum"], 0) + 1)
    summary = {
        "pathogens": len(manifest),
        "gram_counts": gram_counts,
        "phylum_counts": phylum_counts,
        "reference_genomes": sum(
            1 for e in manifest
            if e.get("refseq_category") == "reference genome"),
        "representative_genomes": sum(
            1 for e in manifest
            if e.get("refseq_category") == "representative genome"),
        "manifest": str(manifest_path),
    }
    write_json(output / "sampling_summary.json", summary)
    print(f"Done. {len(manifest)} pathogens; gram={gram_counts}; "
          f"phyla={phylum_counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
