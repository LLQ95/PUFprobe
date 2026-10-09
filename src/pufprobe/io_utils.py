"""Small IO helpers: FASTA parsing/writing and JSON (de)serialization."""

import json
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) pairs from a FASTA file."""
    header: str | None = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def write_fasta(path: Path, records: Iterable[Tuple[str, str]], line_width: int = 60) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for i in range(0, len(sequence), line_width):
                handle.write(sequence[i : i + line_width] + "\n")


def read_json(path: Path) -> object:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: object, indent: int = 2) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=indent, ensure_ascii=False)


def read_jsonl(path: Path) -> List[dict]:
    rows: List[dict] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def write_jsonl(path: Path, rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def find_inputs(root: Path) -> Dict[str, List[Path]]:
    """Recursively collect GFF/GBK/FAA/FNA files under a directory."""
    groups: Dict[str, List[Path]] = {"gff": [], "gbk": [], "faa": [], "fna": []}
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        if suffix in (".gff", ".gff3"):
            groups["gff"].append(path)
        elif suffix in (".gbk", ".gbff", ".genbank"):
            groups["gbk"].append(path)
        elif suffix in (".faa", ".fasta", ".fa") and "protein" in path.name.lower():
            groups["faa"].append(path)
        elif suffix == ".faa":
            groups["faa"].append(path)
        elif suffix == ".fna":
            groups["fna"].append(path)
    return groups
