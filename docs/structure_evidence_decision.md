# Structure evidence decision

This note records how structural models are obtained for the family
representatives that feed the structure-search channel and the HPDB build. It
answers the question of whether a local AlphaFold2 installation is needed.

## What the structure is used for

Only one representative per protein family is folded, and the other members
inherit the result through family assignment. The model is searched against the
foldseek structural database, and fold-level hits are mapped to Pfam families
through SIFTS. A hit counts when the alignment TM-score reaches the same-fold
boundary of 0.5 and the report E-value passes the configured gate. The model
therefore needs to be accurate at the fold level rather than at atomic detail.

## The three available sources

A representative with a UniProt accession usually has a model already computed
in AlphaFold DB. The prediction API returns the current model URL and mean
pLDDT, so the download stays valid when the model version advances. This path
has no GPU cost and no local sequence search.

A representative absent from AlphaFold DB can be folded with ESMFold on the
local GPU. ESMFold takes a single sequence and needs no multiple-sequence
alignment, so one model is produced in seconds on a current GPU. The B-factor
column of the emitted PDB holds the per-residue pLDDT, which is read with the
same parser used for AlphaFold models.

A local AlphaFold2 installation would build its own multiple-sequence
alignment through large reference databases and run a slower inference
pipeline. The af2 engine is not bundled with the package and is accepted only
through externally supplied precomputed models.

## Decision

A local AlphaFold2 installation is not used. Representatives with a UniProt
accession are taken from AlphaFold DB, and representatives without a model
there are folded with ESMFold. This keeps the GPU workload limited to
sequences that have no precomputed model.

The choice follows from the input set and the intended use. The hypothetical
proteins in the database often have no detectable homolog, so their alignments
are shallow and the alignment-dependent advantage of AlphaFold2 is reduced,
while a single-sequence model such as ESMFold remains applicable. Maintaining
the AlphaFold2 reference databases and its inference stack would add storage
and setup cost without a measurable gain for fold-level matching. ESMFold
models are accurate enough for the TM-score threshold applied downstream.

## Confidence handling

Every model carries a median per-residue pLDDT regardless of its source.
Families below the configured threshold are marked low confidence and their
structure is not used by the structure-search channel, so a weak model cannot
drive a prediction. The source engine and the model version are recorded per
family, which lets the structural evidence be traced when the database is
rebuilt.

## Current execution constraint

The GPU folding step runs on the biflab-1 host through the pufprobe conda
environment. When that host is not reachable from the current network, the
AlphaFold DB download path and the package code can still be prepared, but the
ESMFold fallback cannot be executed until the host is reachable again.
