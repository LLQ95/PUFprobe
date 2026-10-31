# Changelog

All notable changes to PUFprobe are documented in this file. The format follows
Keep a Changelog, and the project adheres to semantic versioning.

## [Unreleased]

## [0.10.0] - 2026-10-31

### Added

- A pocket catalytic-template evidence channel. P2Rank predicts binding
  pockets on the query structures, and foldseek aligns them against a
  prebuilt database of M-CSA catalytic sites. A match passes the foldseek
  TM-score and E-value gates, requires the mapped catalytic residues to have
  the same amino-acid type as the template, and keeps the hit only when a
  minimum fraction of the catalytic residues fall inside a predicted pocket.
  Each accepted template emits an EC evidence record in the structure group.
- `server/build_pocket_templates.py` builds the template database from M-CSA:
  it pages through the full entry list, downloads each reference PDB structure,
  extracts and validates the reference chain, and creates the foldseek database
  together with a `pocket_templates.json` index. The build is resumable and
  reuses an existing foldseek database unless `--rebuild-db` is given.
- `server/12_download_evidence_dbs.sh` gains a `pocket` target that runs the
  template build, and five `search` config options control the channel
  (`pocket_sensitivity`, `pocket_min_tmscore`, `pocket_max_evalue`,
  `pocket_min_consistency`, `pocket_max_pockets`).

### Notes

- The pocket channel is off by default. It depends on P2Rank and a Java
  runtime and is resource-intensive, so it is enabled by adding `"pocket"` to
  `search.channels`, typically in a server configuration. Missing tools or
  databases disable the channel without stopping the run.

## [0.9.1] - 2026-10-31

### Fixed

- The hard-coded `pufprobe.__version__` string, which had stayed at 0.6.0
  across several releases, is synchronized with the release metadata. A new
  test asserts that `pufprobe.__version__` matches the version in
  `pyproject.toml`, guarding against future drift.
- Product names now exclude true-path ancestor terms supported only by
  propagated evidence. After the 2019 GO restructuring, generic terms such
  as `membrane` and `cell periphery` sit at depth 2 and passed the depth
  gate even though they carry no protein-specific function. The annotate
  stage keeps a GO term only when some evidence is not propagated, applied
  in addition to the depth gate.

## [0.9.0] - 2026-10-31

### Added

- `GoOntology.depth` returns the longest true-path distance from a GO term to
  a root (a root has depth 0), following the `is_a` and `part_of` DAG. Unknown
  and obsolete terms return None. The rank stage attaches this value to every
  GO term as a new optional `go_depth` field, and a `rank.min_term_depth`
  config option (default 2) sets the gate used when writing product names.

### Changed

- Product names no longer include GO roots and shallow generic terms. The
  annotate stage skips GO terms whose known depth is below `min_term_depth`,
  so terms such as `molecular_function`, `catalytic activity` and `membrane`
  are not written back. Non-GO terms and GO terms with an unknown depth are
  kept, duplicate term names are collapsed, and a protein whose accepted
  terms are all generic is left unannotated rather than given a vague name.

## [0.8.1] - 2026-10-31

### Fixed

- Calibrator fitting now accepts numpy arrays. Empty-input checks use
  `len(scores) == 0` rather than `if not scores`, which raised an ambiguous
  truth-value error on the validation and test arrays, and the temperature
  negative log-likelihood grid search is vectorized with numpy instead of
  nested Python loops.
- Isotonic calibration returns an empty list for a prediction with no terms
  instead of passing a zero-length array to sklearn's `IsotonicRegression`,
  which rejected it with a minimum-of-one error. A regression test covers the
  empty term-list case.

## [0.8.0] - 2026-10-31

### Added

- A learning-to-rank fusion framework (`pufprobe.rank.ltr`). Each
  `(protein, candidate term)` pair is described by a fixed feature vector:
  the best raw score from every channel, the reliability-discounted strength
  of every independent evidence group, cross-channel agreement and coverage
  statistics, and meta features (mean coverage, a capped -log10 E-value signal,
  training-derived term family size and taxonomic breadth, and protein pLDDT).
  Rows are grouped by protein and split by a CAFA-style temporal cut (dates
  when provided, otherwise a seeded group split that never leaks a protein
  across sets). A LightGBM LambdaRank model is used when lightgbm is
  installed, with a sklearn gradient-boosting fallback and a dependency-free
  numpy logistic ranker so the framework always runs; the model, its feature
  order and term priors are stored in one pickle. Evaluation reports
  protein-centric Fmax and its threshold, micro-average AUPR, and semantic
  distance Smin when term information content is supplied.
  `server/train_ltr.py` is the training CLI. `server/build_ltr_trainset.py`
  builds the annotated training set: a FASTA of curated proteins is run through
  the normal embed and search stages, each query's self-hit is removed so a
  protein cannot label itself, and the evidence is paired with the curated
  truth from the Swiss-Prot cross-reference index, reproducing the inference
  feature distribution. The trainset build parallelizes the independent
  hhblits queries (one single-threaded search per thread by default rather
  than one serial search), kills any query that exceeds the per-query
  wall-clock timeout, bounds each MSA (maxseq, maxfilt, Neff and memory) so
  the large homolog sets of curated proteins do not build multi-GB alignments,
  and reuses completed embed/search stages on rerun.
  `server/build_sprot_dates.py` extracts each curated entry's sequence-version
  date for the temporal split. `server/add_eggnog_channel.py` adds the eggnog
  channel to an already-built dataset by invoking emapper from its own conda
  environment and merging the evidence, so a skipped channel does not force
  rerunning the expensive hhblits searches. The rank stage loads a configured model from
  `rank.ltr_model_path` (a relative path resolves against the dbdir) and uses
  its learned ranking in place of the rule fusion, falling back silently to
  group noisy-OR when the model is absent or unloadable;
  `calibration.json` records the active backend.
- A protein-language-model nearest-neighbour evidence channel
  (`search --channels plm_neighbor`). Representative ESM-2 embeddings are
  compared against an annotated reference library (`<dbdir>/plm_ref`); the
  cosine top-k neighbours transfer their GO, EC and Pfam terms, weighted by a
  linear ramp of cosine and aggregated per term with a noisy-OR, so agreement
  across several neighbours accumulates. The channel verifies that the query
  and reference embeddings share the same checkpoint and pooling and disables
  itself on a mismatch. The channel belongs to the independent "sequence"
  evidence group. `server/build_plm_reference.py` builds the library from
  `uniprot_sprot.dat[.gz]` (or a FASTA plus cross-reference index), loads the
  ESM checkpoint once, and writes resumable shards; the embed stage was
  refactored into shared `load_esm_model` / `embed_items` functions so query
  and reference pooling cannot diverge.
- BERT-whitening for the pLM nearest-neighbour channel. ESM-2 residue-mean
  vectors are anisotropic, so unrelated proteins cluster at high cosine and
  the raw channel produces many false positives. The reference build now
  estimates a whitening transform from the library covariance and applies it
  to the reference vectors and each query, followed by L2 normalisation. The
  transform restores family-specific directions: unrelated proteins move to
  cosine near zero while proteins from the same family stay separated at high
  cosine. `server/whiten_existing_library.py` applies the transform to an
  already-built library without re-embedding, and `server/verify_whitening.py`,
  `server/tune_whitening_threshold.py` and `server/eval_whitening_floor.py`
  diagnose the geometry and the cosine floor. On the hold-out set the whitened
  four-channel model reaches Fmax 0.8435 and micro-AUPR 0.8810, the highest
  AUPR of the tracked models; transferred pLM terms fall from about 120000 to
  49856, and pLM-only terms shrink to a small set kept below the high
  confidence threshold.
- Post-hoc confidence calibration for the LTR ranker. The sigmoid applied to
  the LightGBM output maps ranking scores into [0, 1] but does not make them
  probabilities, so a fixed confidence threshold such as 0.80 does not match
  the observed correctness. `server/fit_calibrator.py` rebuilds the same
  temporal split, scores the validation and test proteins with a saved ranker,
  and fits temperature and isotonic calibrators on the validation split only,
  selecting the method by validation negative log-likelihood. The untouched
  test split reports expected calibration error, Brier score and per-bin
  reliability. Fitted calibrators are stored with the calibration helpers and
  loaded at rank time through `rank.calibration_model_path`; an empty path or a
  missing file keeps the identity calibration without error.
- Gene Ontology true-path propagation in the rank stage. A dependency-free OBO
  reader (`pufprobe.ontology_graph`) loads `go-basic.obo` (plain or gzip),
  folds alt ids onto primary terms, excludes obsolete terms, and extends each
  GO evidence entry to its full `is_a` and `part_of` ancestor closure. The
  `regulates` relation is deliberately not traversed. Propagated evidence
  keeps the same channel and evidence group as the direct hit (it is never
  counted as an independent witness), and only the strongest support is kept
  for a given ancestor and channel. The ontology is located from
  `rank.go_obo_path` or by looking for `go-basic.obo` in the rank `--dbdir`;
  when absent the rank stage logs that GO terms are kept at their annotated
  level. `calibration.json` records the propagation state
  (`go_true_path_propagation`, `n_go_propagated`).
- An eggNOG-mapper evidence channel (`search --channels eggnog`). It runs the
  externally installed emapper in its own conda environment against an
  eggNOG v7 database and parses the emapper 3.0 annotation format (OG/COG,
  GO, EC and Pfam fields), feeding the database evidence group.
- A server script `server/12_download_evidence_dbs.sh` that fetches the
  databases backing the fusion evidence channels. It provides STRING v12.5
  COG-level links and per-species interaction networks, the GO `go-basic.obo`
  ontology for true-path propagation, the M-CSA catalytic-site entries and
  residues, and the UniRef30 HHblits build. The script is idempotent and
  resumable, downloads to the data volume (`EVIDENCE_DB_ROOT`), and reuses the
  `hhblits_db_*` symlink convention from `01_download_dbs.sh`.
- Configuration variables for the evidence databases (`EVIDENCE_DB_ROOT`,
  `STRING_VER`, `STRING_TAXIDS`, `STRING_TAXIDS_FILE`, `STRING_WITH_ALIASES`,
  `WITH_GOA`, `EVIDENCE_UNIREF30_RELEASE`). The default STRING panel lists
  verified strain-level reference proteomes for PUFprobe's target pathogens.

### Removed

- The `goatools` pip install from `server/12_download_evidence_dbs.sh`. The
  true-path closure is now computed by the built-in dependency-free
  `pufprobe.ontology_graph` reader, so goatools is no longer required to run
  the pipeline (it remains optional for independent GO enrichment work).

### Changed

- The large UniRef30 HHblits archive is now fetched from the Steinegger lab
  S3 mirror (`opendata.mmseqs.org/uniclust`) first, with the throttled GWDG
  host as a fallback, and is downloaded with `aria2c` multi-connection
  segments. The harmful `--lowest-speed-limit` option was removed: on a host
  that throttles every connection it aborts the whole download. The fetcher
  now uses `--timeout` to drop only fully stalled connections and retries each
  mirror in order.
- The pLM nearest-neighbour channel is now enabled by default in
  `search.channels`, and the rank stage defaults to the whitened four-channel
  model `ltr_fusion_white.pkl`. A missing library or model still falls back to
  the other channels or rule fusion without error, and the earlier
  un-whitened model remains available for comparison.

## [0.6.0] - 2026-10-21

Pfam terms now display their official family descriptions. A compact name
index built from the Pfam-A.hmm.dat metadata replaces PDB template ids and
distant eukaryotic homolog names with the family's own description.

### Added

- A Pfam family name index (`pfam_names.json`), built offline by a companion
  server script (`server/build_pfam_names.py`) that parses the
  `Pfam-A.hmm.dat` metadata and stores each accession with its short
  identifier, official description, and entry type. The trailing version
  number is stripped from accessions.
- A `--dbdir` option on the rank command, which loads the Pfam name index
  and resolves Pfam terms to their official descriptions.

### Changed

- Pfam term names are taken from the official family description rather
  than the strongest evidence label, so a bacterial protein is no longer
  named after a PDB chain or a distant eukaryotic homolog. The index
  carries no scoring weight, and the predicted terms and confidence values
  are unchanged.

## [0.5.0] - 2026-10-21

Structure-based inference is now an active, gated evidence channel. The
pipeline folds family representatives, searches their structures, and
resolves structural hits against Pfam, so a sequence hit and a structural
hit can converge on the same family and provide two independent witnesses.

### Added

- GPU structure modelling with ESMFold for family representatives. A family
  is marked folded when its median pLDDT is at least 70, and low-confidence
  models are retained but not searched. Representatives longer than 1000
  residues are skipped, and precomputed models are accepted.
- Foldseek structure search on folded representatives, with PDB-chain
  targets mapped to Pfam accessions through the SIFTS chain-level summary.
  Structural and sequence evidence therefore converge on the same Pfam
  terms.

### Changed

- Foldseek hits are scored by alignment TM-score, query coverage and E-value
  significance, mirroring the sequence-homology score, instead of Foldseek's
  reported probability, which saturates near 1.0 even for short partial
  alignments. Hits below a TM-score of 0.5 or above an E-value of 1e-2 are
  rejected; coverage is not hard-filtered, because domain-level matches can
  legitimately cover part of a multi-domain template.
- Each direct evidence group must clear a minimum strength gate
  (`rank.group_min_strength`, default 0.40) to count as an independent
  witness. A strong sequence profile paired with a weak structure, or a
  strong structure paired with a twilight-zone sequence, is held at the
  single-group level, while weaker groups still contribute to the fused
  score.
- Genomic context is classified as an indirect group. It can raise
  confidence through the noisy-OR combination, but it does not count as an
  independent witness for the two-group requirement.
- Structure-backed Pfam terms display a protein name rather than a bare PDB
  identifier. The display name is taken from the best non-PDB protein hit,
  falling back to the Pfam accession, while structure-only terms keep the
  PDB name.

## [0.4.0] - 2026-10-20

Term normalization and stronger homology scoring. Independent evidence groups
can only reinforce the same function when database-specific identifiers are
mapped onto a shared ontology.

### Added

- A Swiss-Prot cross-reference index (`search/idmapping.py`) maps accessions to
  Pfam, EC, and GO terms parsed from `uniprot_sprot.dat.gz`. The database
  download script builds `swissprot_xref.json` (`server/build_xref.py`), and
  the search stage supplements homology hits with ontology evidence while
  keeping the original evidence channel. The hit still originates from
  sequence homology, so cross-references from several accessions remain one
  evidence group; a structurally derived hit mapped to the same Pfam later
  supplies a genuinely independent group.
- Gene Ontology terms are emitted as identifiers without true-path
  propagation, because no ontology graph is loaded.

### Changed

- MMseqs2 homology scores now combine sequence identity, coverage, and E-value
  rather than relying on E-value alone. Identity follows a conservative
  transfer curve (around 40% identity is a safe-transfer boundary, below
  roughly 20% transfer is unreliable), and coverage guards against short local
  matches. This prevents weak, low-identity matches from reaching the same
  score as a strong transfer, which also resolves cases where one protein was
  assigned several unrelated functions.
- Hits to targets that carry no function themselves (a Swiss-Prot entry still
  labelled hypothetical or uncharacterized protein) are filtered out, because
  they support homology but not a functional term.

## [0.3.0] - 2026-10-20

Methodological change to evidence fusion. The previous linear weighted sum
could not reach high confidence under independent evidence: the channel
weights summed to 1.0, so even the two strongest independent groups totalled
only 0.42 and a single channel was capped at 0.20, which made the 0.80
annotation threshold mathematically unreachable. Fusion is now organised
around independent evidence groups.

### Changed

- Rule fusion (`rank.fusion_method`, default `group_noisy_or`) now aggregates
  evidence at the level of independent evidence groups (sequence,
  sequence_profile, structure, context, database). Within a group the
  strongest reliability-discounted hit is taken, so channels that share an
  underlying signal (foldseek and pocket both rely on structure) are not
  counted twice. Group strengths are combined with a noisy-OR, which is the
  probability that at least one independent witness supports the term.
- A term backed by fewer than `rank.require_independent_groups` groups is
  capped at `rank.single_group_cap` (default 0.60). A single strong homology
  hit therefore becomes a medium-confidence candidate for manual review,
  while two strong orthogonal groups (for example sequence profile plus
  structure, mapped to the same Pfam/EC/GO term) can reach high confidence.
- Each channel carries a `rank.channel_reliability` prior, the upper
  probability that a perfect hit transfers function. These are conservative
  uncalibrated priors for the rule path; a learned fusion model or held-out
  calibration set overrides them.
- The legacy linear fusion is retained as `rank.fusion_method: weighted_sum`
  for comparison.
- `calibration.json` records the fusion method, single-group cap and required
  group count in addition to the previous fields.

### Note

- Reaching high confidence still requires independent evidence channels to be
  present and to agree on the same function term. The structural channel
  (ESMFold plus Foldseek) and curated databases (CARD/VFDB) must be enabled,
  and hits from different databases must be mapped to a common ontology
  (Pfam/EC/GO) before their evidence can be combined.

## [0.2.7] - 2026-10-20

Performance fix for the report stage on a full cohort.

### Fixed

- The report index was built by scanning every prediction for each CDS record,
  making the stage O(n_cds * n_predictions). On a 70,505-CDS run this was on
  the order of 5 billion comparisons and appeared to hang after the stage
  started. Predictions are now indexed once into a dictionary and each CDS is
  resolved in O(1), reducing the stage to O(n).

## [0.2.6] - 2026-10-20

Fixes for the fold, search, and context stages after a full fallback run
(`FORCE_FALLBACK=1 bash server/03_run_pufprobe.sh`).

### Fixed

- MMseqs2 v18 does not accept a `desc` or `tdesc` format-output column;
  `convertalis` failed with "Format code desc does not exist" and the channel
  produced no evidence. The full target header column is now `theader`, the
  field supported in MMseqs2 v18, and the free-text description is derived
  from it.
- The fold stage ignored the CPU/fallback flag and attempted to load ESMFold
  once per family, failing on every family because of network access and a
  missing dependency. ESMFold is now disabled and all families are marked
  skipped under `--cpu`/fallback unless `fold.allow_cpu_fold` is set.
- The ESMFold model is loaded once before the family loop instead of being
  reloaded for each family; a missing model or dependency is reported a single
  time.
- The context stage aborted with "unbound variable" on older Bash because an
  empty `CONTEXT_OPTS` array was expanded under `set -u`; the array now uses
  the empty-safe expansion.
- MMseqs2 command failures now include the program's stderr in the raised
  error so the cause appears in the log instead of only the exit code.

### Added

- `00_setup_env.sh` installs `einops` together with `fair-esm`; `einops` is
  required by ESMFold but is not pulled in by `fair-esm`.
- Config key `fold.allow_cpu_fold` (default `false`).

## [0.2.5] - 2026-10-20

All code, comments, and documentation shipped on GitHub are now in English; the
project is being prepared for an English-language publication.

### Changed

- Translated the Chinese text in the server configuration, the Tier A hotfix
  script, the summary script, and the summary stage header comments
- Rewrote `server/README_runbook.md` and `docs/server_setup.md` in English
- Translated the text labels in `docs/architecture.svg`; the figure remains a
  valid SVG and now renders only English text
- A repository-wide scan for CJK characters returns no matches

## [0.2.4] - 2026-10-20

Follow-up to v0.2.3: the Bakta protein subset is named
`<prefix>.hypotheticals.faa` (plural), while v0.2.3 only matched the singular
`.hypothetical.faa`, so the subset still overrode the complete set.

### Fixed

- The input collection and sample grouping now skip both
  `.hypotheticals.faa` (the released Bakta name) and `.hypothetical.faa`; the
  regression test builds both spellings and asserts the complete `.faa` is the
  only protein index, with all CDS aligned

## [0.2.3] - 2026-10-20

Hotfix for prepare aligning CDS records after a Bakta run.

### Fixed

- Bakta writes both a complete `<prefix>.faa` and an extra
  `<prefix>.hypothetical.faa` subset holding only the hypothetical proteins.
  When a directory was passed to prepare, both matched the `*.faa` recursion and
  the subset overrode the complete set, so every annotated protein was logged
  as "could not be aligned to any .faa sequence". The input collection and the
  sample grouping now skip `*.hypothetical.faa` and always index the complete
  protein set; a regression test reproduces the Bakta layout

## [0.2.2] - 2026-10-19

Hotfix for database download failures on the first real TIER A run.

### Added

- server/hotfix_tierA.sh: a self-contained, idempotent Tier A database repair
  script that builds Swiss-Prot from the UniProt FASTA, fetches the CARD
  canonical bundle and VFDB, and optionally pulls Foldseek PDB and ESM-2
  weights, writing the same completion markers as 01_download_dbs.sh so the
  workflow can continue without re-transferring the whole release

### Fixed

- The MMseqs2 Swiss-Prot database is requested as `UniProtKB/Swiss-Prot`, the
  valid name in MMseqs2 v18; the old `Swiss-Prot` name aborted with
  "Selected database Swiss-Prot was not found"
- The retired mmseqs.com tarball fallback (which 302-redirects to an HTML
  page) is replaced by a UniProt FASTA fallback: the script downloads
  uniprot_sprot / uniref90 / uniref50 FASTA from ftp.uniprot.org and builds
  the search database with mmseqs createdb + createindex; partial HTML files
  left in the tmp directory are cleaned automatically
- CARD no longer uses the legacy /download/6/*.fasta link (now an HTML page);
  it downloads the canonical https://card.mcmaster.ca/latest/data tar.bz2,
  validates its integrity and extracts protein_fasta_protein_homolog_model.fasta
- ESM weight preloading is non-fatal: a network failure only warns, so the
  smoke run (which uses sequence fallback features) is not blocked
- server_config.yaml gains an overridable UNIPROT_BASE for mirror substitution

## [0.2.1] - 2026-10-19

Hotfix for the first real server run of 00_setup_env.sh.

### Fixed

- The setup script now pins MAMBA_ROOT_PREFIX to the conda base and invokes
  the runner by absolute path, so the environment is created at
  anaconda3/envs/pufprobe instead of being nested inside the currently
  activated environment (e.g. anaconda3/envs/checkm/envs/pufprobe)
- Any incorrectly nested environment left by v0.2.0 is detected and removed
  automatically when the script is re-run
- P2Rank was never available on bioconda and aborted the conda solver; it is
  no longer installed through conda. The script now fetches the portable
  P2Rank 2.5.1 GitHub release into thirdparty/p2rank (with openjdk 17 from
  conda-forge and proxy-mirror fallbacks), and failure is non-blocking
- torch wheels are now tried in order from Aliyun, SJTU and the official
  PyTorch index instead of relying on the official index only
- externals.p2rank() resolves the portable installation through P2RANK_HOME
  and the repo-local thirdparty/p2rank/distro/prank before falling back to
  PATH

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
