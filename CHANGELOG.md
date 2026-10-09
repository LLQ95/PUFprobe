# Changelog

All notable changes to PUFprobe are documented in this file. The format follows
Keep a Changelog, and the project adheres to semantic versioning.

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
