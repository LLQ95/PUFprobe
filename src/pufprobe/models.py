"""Frozen data contracts shared by every PUFprobe module.

All records passed between pipeline stages are the dataclasses defined here.
They serialize to plain JSON-compatible dictionaries through ``to_dict`` and
can be reconstructed through ``from_dict``. A protein is addressed across
samples by a global key of the form ``sample_id@locus_tag``.
"""

from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Any, Dict, List, Optional


def gkey(sample_id: str, locus_tag: str) -> str:
    """Return the global key for a CDS record."""
    return f"{sample_id}@{locus_tag}"


def split_gkey(key: str) -> tuple[str, str]:
    sample_id, _, locus_tag = key.partition("@")
    return sample_id, locus_tag


class GeneStatus(str, Enum):
    HYPOTHETICAL = "hypothetical"
    CONSERVED_HYPOTHETICAL = "conserved_hypothetical"
    ANNOTATED = "annotated"


class GeneVerdict(str, Enum):
    LIKELY_FALSE_ORF = "likely_false_orf"
    TRUE_PROTEIN_UNKNOWN = "true_protein_unknown"
    BORDERLINE = "borderline"


class EvidenceChannel(str, Enum):
    PLM_NEIGHBOR = "plm_neighbor"
    FOLDSEEK = "foldseek"
    REMOTE_HOMOLOGY = "remote_homology"
    POCKET = "pocket"
    CONTEXT = "context"
    DATABASE = "database"
    EGGNOG = "eggnog"


# Independent evidence groups used by the channel-diversity rule. Channels in
# the same group are not counted twice when deciding high confidence.
#
# eggNOG-mapper is treated as a curated database witness (group "database"),
# independent of the representative profile searches the pipeline runs itself
# (group "sequence_profile"). An eggNOG annotation and a profile hit to the
# same function therefore count as two orthogonal groups.
CHANNEL_GROUPS: Dict[str, str] = {
    EvidenceChannel.PLM_NEIGHBOR.value: "sequence",
    EvidenceChannel.FOLDSEEK.value: "structure",
    EvidenceChannel.REMOTE_HOMOLOGY.value: "sequence_profile",
    EvidenceChannel.POCKET.value: "structure",
    EvidenceChannel.CONTEXT.value: "context",
    EvidenceChannel.DATABASE.value: "database",
    EvidenceChannel.EGGNOG.value: "database",
}


@dataclass
class Sample:
    sample_id: str
    gff_path: Optional[str] = None
    gbk_path: Optional[str] = None
    faa_path: Optional[str] = None
    fna_path: Optional[str] = None
    source: str = "generic"  # bakta, prokka or generic


@dataclass
class CDSRecord:
    sample_id: str
    locus_tag: str
    protein_id: Optional[str]
    sequence_aa: str
    start: int
    end: int
    strand: str
    product: str
    status: GeneStatus
    contig: str = ""
    dbxrefs: List[str] = field(default_factory=list)
    attributes: Dict[str, Any] = field(default_factory=dict)
    sequence_nt: Optional[str] = None
    # true-gene quality control
    verdict: Optional[GeneVerdict] = None
    qc_score: Optional[float] = None
    qc_reasons: List[str] = field(default_factory=list)
    # family assignment
    family_id: Optional[str] = None
    is_representative: bool = False

    @property
    def key(self) -> str:
        return gkey(self.sample_id, self.locus_tag)


@dataclass
class Evidence:
    channel: EvidenceChannel
    target_id: str
    target_label: str
    score: float
    evalue: Optional[float] = None
    coverage: Optional[float] = None
    details: Dict[str, Any] = field(default_factory=dict)

    @property
    def group(self) -> str:
        return CHANNEL_GROUPS[self.channel.value]


@dataclass
class TermPrediction:
    term_type: str  # GO, EC, Pfam, COG, KEGG
    term_id: str
    term_name: str
    raw_score: float
    confidence: float
    evidence: List[Evidence] = field(default_factory=list)
    rationale: str = ""

    @property
    def evidence_groups(self) -> List[str]:
        return sorted({item.group for item in self.evidence})


@dataclass
class ProteinPrediction:
    sample_id: str
    locus_tag: str
    terms: List[TermPrediction] = field(default_factory=list)
    priority_score: Optional[float] = None
    priority_reasons: List[str] = field(default_factory=list)

    @property
    def key(self) -> str:
        return gkey(self.sample_id, self.locus_tag)


@dataclass
class Family:
    family_id: str
    representative_key: str
    members: List[str]
    representative_seq: str
    size: int
    lineage_width: int
    metrics: Dict[str, Any] = field(default_factory=dict)
    structure_path: Optional[str] = None
    plddt: Optional[float] = None
    fold_status: str = "not_folded"  # folded, low_confidence, skipped, failed
    novel_score: Optional[float] = None
    priority_score: Optional[float] = None
    priority_reasons: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Serialization helpers
# ---------------------------------------------------------------------------

def to_dict(obj: Any) -> Any:
    """Recursively convert dataclasses/enums into JSON-compatible values."""
    if hasattr(obj, "__dataclass_fields__"):
        return {key: to_dict(value) for key, value in asdict(obj).items()}
    if isinstance(obj, Enum):
        return obj.value
    if isinstance(obj, dict):
        return {key: to_dict(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_dict(value) for value in obj]
    return obj


def evidence_from_dict(data: dict) -> Evidence:
    data = dict(data)
    data["channel"] = EvidenceChannel(data["channel"])
    return Evidence(**data)


def term_from_dict(data: dict) -> TermPrediction:
    data = dict(data)
    data["evidence"] = [evidence_from_dict(item) for item in data.get("evidence", [])]
    return TermPrediction(**data)


def cds_from_dict(data: dict) -> CDSRecord:
    data = dict(data)
    if data.get("status") is not None:
        data["status"] = GeneStatus(data["status"])
    if data.get("verdict") is not None:
        data["verdict"] = GeneVerdict(data["verdict"])
    return CDSRecord(**data)


def sample_from_dict(data: dict) -> Sample:
    return Sample(**dict(data))


def family_from_dict(data: dict) -> Family:
    return Family(**dict(data))


def prediction_from_dict(data: dict) -> ProteinPrediction:
    data = dict(data)
    data["terms"] = [term_from_dict(item) for item in data.get("terms", [])]
    return ProteinPrediction(**data)


_BUILDERS = {
    Sample: sample_from_dict,
    CDSRecord: cds_from_dict,
    Evidence: evidence_from_dict,
    TermPrediction: term_from_dict,
    ProteinPrediction: prediction_from_dict,
    Family: family_from_dict,
}


def from_dict(data: Any, cls: type) -> Any:
    """Reconstruct an instance of one of the registered model classes."""
    builder = _BUILDERS.get(cls)
    if builder is None:
        raise KeyError(f"No builder registered for {getattr(cls, '__name__', cls)}")
    return builder(data)


# Backwards-compatible alias.
def _evidence_from_dict(data: dict) -> Evidence:
    return evidence_from_dict(data)


def records_from_dicts(rows: List[dict]) -> List[CDSRecord]:
    return [cds_from_dict(row) for row in rows]


def families_from_dicts(rows: List[dict]) -> List[Family]:
    return [family_from_dict(row) for row in rows]


def predictions_from_dicts(rows: List[dict]) -> List[ProteinPrediction]:
    return [prediction_from_dict(row) for row in rows]
