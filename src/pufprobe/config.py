"""Configuration constants, default thresholds and config-file loading.

Thresholds can be overridden through a JSON config file passed with
``--config``. Unknown keys are rejected so that typos surface early.
"""

import json
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict


DEFAULTS: Dict[str, Any] = {
    # Product strings treated as markers of a PUF record.
    "hypothetical_patterns": [
        "hypothetical protein",
        "uncharacterized protein",
        "unknown protein",
        "protein of unknown function",
        "putative uncharacterized",
        "domain of unknown function",
        "DUF",
    ],
    # prepare / true-gene quality control
    "qc": {
        "min_length_aa": 30,          # Prodigal does not call proteins shorter than ~30 aa
        "smorf_length_aa": 50,        # records below this use the permissive smORF path
        "min_length_smorf_aa": 16,
        "max_low_complexity_fraction": 0.60,
        "max_homopolymer_run": 10,
        "max_gc_skew_abs": 0.30,
        "gc_zscore_cutoff": 3.0,
        "max_overlap_fraction": 0.60,
        "false_orf_score": 0.50,      # score at/above which a record is moved out
        "borderline_score": 0.30,
        "rbs_motifs": ["AGGAGG", "AGGA", "GGAG", "AGGAG", "GAGG"],
    },
    # embedding
    "embed": {
        "default_model": "esm2_t33_650M_UR50D",
        "cpu_model": "esm2_t12_35M_UR50D",
        "batch_size": 8,
        "fallback_feature_dim": 48,
    },
    # clustering
    "cluster": {
        "identity": 0.30,
        "coverage": 0.80,
        "cosine_threshold": 0.80,
        "min_cluster_size": 2,
    },
    # folding
    "fold": {
        "engine": "esmfold",
        "min_plddt": 70.0,
        "max_length_esmfold": 1000,
        "tolerate_disordered": True,
        # When false (default), ESMFold is not loaded under --cpu / fallback.
        "allow_cpu_fold": False,
    },
    # search
    "search": {
        "channels": ["mmseqs", "foldseek", "hhblits"],
        "foldseek_sensitivity": 9.5,
        # Structural hit quality gates: TM-score >= 0.5 is the same-fold
        # boundary, and the foldseek report E-value is kept slightly looser
        # than the sequence channel because structure detects remote homologs.
        "foldseek_min_tmscore": 0.5,
        "foldseek_report_evalue": 1e-2,
        "hhblits_iterations": 2,
        "report_evalue": 1e-3,
    },
    # context
    "context": {
        "neighbour_window": 5,
        "operon_intergenic_bp": 60,
        "min_genomes_for_profile": 4,
        "same_strand_operon": True,
    },
    # ranking / calibration
    "rank": {
        "high_confidence": 0.80,
        "min_confidence_annotate": 0.80,
        "require_independent_groups": 2,
        # Rule-based fusion. "group_noisy_or" (default) aggregates evidence at
        # the level of independent evidence groups: the strongest reliability-
        # discounted hit is taken per group, and group strengths are combined
        # with a noisy-OR. Terms backed by fewer than require_independent_groups
        # are capped at single_group_cap (medium confidence, manual review).
        # "weighted_sum" is the previous linear sum kept for comparison.
        "fusion_method": "group_noisy_or",
        "single_group_cap": 0.60,
        # Reliability prior of a strong hit from each channel, i.e. the upper
        # probability that a perfect hit in that channel transfers function.
        # These are conservative priors for the uncalibrated rule path; a
        # learned fusion model or held-out calibration set overrides them.
        "channel_reliability": {
            "plm_neighbor": 0.75,
            "foldseek": 0.85,
            "remote_homology": 0.95,
            "pocket": 0.60,
            "context": 0.60,
            "database": 0.90,
        },
        # Linear weights used only when fusion_method == "weighted_sum".
        "rule_weights": {
            "plm_neighbor": 0.18,
            "foldseek": 0.22,
            "remote_homology": 0.20,
            "pocket": 0.12,
            "context": 0.16,
            "database": 0.12,
        },
        "priority_weights": {
            "structural_novelty": 0.30,
            "pocket": 0.20,
            "clonal_specificity": 0.20,
            "amr_virulence_neighbourhood": 0.20,
            "small_protein": 0.10,
        },
    },
}


def load_config(path: Path | None = None) -> Dict[str, Any]:
    """Return defaults updated with the optional JSON config file."""
    config = deepcopy(DEFAULTS)
    if path is None:
        return config
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        user = json.load(handle)
    _update_nested(config, user)
    return config


def _update_nested(base: Dict[str, Any], override: Dict[str, Any]) -> None:
    for key, value in override.items():
        if key not in base:
            raise KeyError(f"Unknown config key: {key}")
        if isinstance(value, dict) and isinstance(base[key], dict):
            _update_nested(base[key], value)
        else:
            base[key] = value
