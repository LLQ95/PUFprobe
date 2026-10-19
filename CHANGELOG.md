# Changelog

All notable changes to PUFprobe are documented in this file. The format follows
Keep a Changelog, and the project adheres to semantic versioning.

## [0.2.0] - 2026-10-19

Server-hardening release for the computational deployment phase.

### Added

- server/ directory with a turnkey deployment bundle for the biflab-1 cluster:
  00_setup_env.sh, 01_download_dbs.sh (two tiers), 02_annotate.sh,
  03_run_pufprobe.sh, 04_summary.sh, summarize.py,
  configs/server_config.yaml and the Chinese step-by-step runbook
  README_runbook.md
- All server scripts are LF, `set -euo pipefail`, logged, idempotent and
  resumable; paths, threads, mirrors and database tiers are overridable
  through environment variables
- Automatic GPU detection (nvidia-smi) with separate CUDA and CPU torch paths,
  and China-accessible mirror configuration for conda, pip and database
  downloads
- Real-format regression tests for Prokka .gff (PROKKA_ ids, embedded
  ##FASTA section, percent-encoded products) and Bakta .gff3 outputs, plus
  .gbff discovery coverage

### Changed

- File discovery now recognises the Bakta .gbff extension and *.genbank
  alongside .gff/.gff3/.gbk
- MMseqs2 database resolver accepts the swissprot, uniref90, uniref50 and
  uniref100 prefixes in addition to uniref/mmseqs_db, matching the two-tier
  database layout used by the server scripts
- HHblits database detection recognises _hhm.ffdata, _a3m.ffdata and
  _msa.ffdata ffindex markers instead of relying on a .hhr file
- Annotation source detection reads the GFF column-2 source field
  case-insensitively, so lowercase "bakta" sources classify correctly
- External tool lookup also searches CONDA_PREFIX/bin and the mamba root
  prefix, which covers conda run invocation without an activated shell

### Fixed

- Families that could not be folded because the ESMFold backend is absent are
  recorded as "skipped" (missing tool) rather than "failed" (folding error)

## [0.1.0] - 2026-10-09

First phase release.

### Added

- Standard Python package using a src layout, installable with pip install -e .
- Command line interface with eight subcommands: prepare, embed, fold, search,
  context, rank, report, annotate
- prepare parses Bakta and Prokka GFF3, GenBank GBK and protein FAA files into
  a unified sample/CDS registry, with checkpointing and run metadata
- True-gene quality control that moves likely false ORFs out of the inference
  queue with interpretable, configurable scores and reasons
- embed with optional ESM-2 embeddings, precomputed-embedding support and a
  pure-Python sequence feature fallback
- Family clustering through MMseqs2 when present, with an embedding/graph
  fallback and representative-sequence selection
- fold wrapping optional ESMFold for family representatives, with a pLDDT
  quality gate and precomputed-structure support
- search wrapping MMseqs2, Foldseek and HHblits (TM-Vec optional) into a
  unified evidence format, with graceful per-channel skipping
- context implementing operon/neighbour, phylogenetic profile, cross-strain
  micro-synteny signals and CARD/VFDB/SPI neighbourhood hooks in pure Python
- rank combining evidence through a learning-to-rank/logistic framework with
  rule-based fallback, calibration hooks, independent-channel agreement policy
  and wet-lab priority ranking for unanchored families
- report producing machine-readable JSON and auditable HTML/Markdown reports
- annotate backfilling high-confidence predictions into new GFF3/GBK/FAA files
  without overwriting the originals
- pytest unit tests and a CPU-only end-to-end smoke test on bundled examples
