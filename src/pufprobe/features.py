"""Lightweight, dependency-free sequence features.

These descriptors keep the pipeline runnable without protein language models
or GPUs. They are deliberately simple and interpretable; the ESM-2 embeddings
replace them whenever torch and model weights are available.
"""

import math
from collections import Counter
from typing import Dict, List, Sequence

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"

# Kyte-Doolitzer hydropathy, indexed by AA_ORDER.
_HYDROPATHY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
_CHARGE_POS = set("KRH")
_CHARGE_NEG = set("DE")
_TINY = set("ASDG")
_SMALL = set("ACSTP")
_AROMATIC = set("FWY")
_NONPOLAR = set("AVLIPMFW")


def amino_composition(sequence: str) -> List[float]:
    counts = Counter(sequence)
    length = max(len(sequence), 1)
    return [counts.get(aa, 0) / length for aa in AA_ORDER]


def physicochemical(sequence: str) -> List[float]:
    length = max(len(sequence), 1)
    counts = Counter(sequence)
    hydro_values = [_HYDROPATHY.get(aa, 0.0) for aa in sequence]
    mean_hydro = sum(hydro_values) / length
    var_hydro = sum((value - mean_hydro) ** 2 for value in hydro_values) / length
    return [
        math.log(length),
        sum(1 for aa in sequence if aa in _CHARGE_POS) / length,
        sum(1 for aa in sequence if aa in _CHARGE_NEG) / length,
        sum(1 for aa in sequence if aa in _TINY) / length,
        sum(1 for aa in sequence if aa in _SMALL) / length,
        sum(1 for aa in sequence if aa in _AROMATIC) / length,
        sum(1 for aa in sequence if aa in _NONPOLAR) / length,
        mean_hydro,
        math.sqrt(var_hydro),
        # fraction of residues in the most hydrophobic 19-residue window (TM proxy)
        _max_window_hydro(sequence, 19),
    ]


def _max_window_hydro(sequence: str, window: int) -> float:
    if len(sequence) < window:
        return sum(_HYDROPATHY.get(aa, 0.0) for aa in sequence) / max(len(sequence), 1)
    values = [_HYDROPATHY.get(aa, 0.0) for aa in sequence]
    running = sum(values[:window])
    best = running
    for i in range(window, len(values)):
        running += values[i] - values[i - window]
        best = max(best, running)
    return best / window


def low_complexity_fraction(sequence: str, window: int = 12) -> float:
    """Fraction of residues in low-diversity windows."""
    if len(sequence) < window:
        return 1.0 - len(set(sequence)) / max(len(sequence), 1)
    hits = 0
    n_windows = len(sequence) - window + 1
    for i in range(n_windows):
        if len(set(sequence[i : i + window])) <= 4:
            hits += 1
    return hits / n_windows


def longest_homopolymer(sequence: str) -> int:
    best = current = 1
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            current += 1
            best = max(best, current)
        else:
            current = 1
    return best if sequence else 0


def fallback_features(sequence: str) -> List[float]:
    """Fixed-order feature vector: composition (20) + physicochemical (10)."""
    seq = sequence.upper().replace("*", "")
    return amino_composition(seq) + physicochemical(seq)


# ---------------------------------------------------------------------------
# Nucleotide-level coding signals used by true-gene quality control
# ---------------------------------------------------------------------------

def gc_content(sequence: str) -> float:
    seq = sequence.upper()
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def gc_skew(sequence: str) -> float:
    seq = sequence.upper()
    g, c = seq.count("G"), seq.count("C")
    if g + c == 0:
        return 0.0
    return (g - c) / (g + c)


def has_internal_stop(nt_sequence: str, table: int = 11) -> bool:
    stops = {"TAA", "TAG", "TGA"}
    for i in range(0, len(nt_sequence) - 2, 3):
        if nt_sequence[i : i + 3].upper() in stops:
            return True
    return False


def rbs_likelihood(upstream: str, motifs: Sequence[str]) -> bool:
    """True when a Shine-Dalgarno-like motif occurs in the upstream window."""
    if not upstream:
        return False
    window = upstream[-15:].upper().replace("T", "U")
    dna_window = upstream[-15:].upper()
    return any(motif in dna_window for motif in motifs) or any(
        motif.replace("T", "U") in window for motif in motifs
    )
