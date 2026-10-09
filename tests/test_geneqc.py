"""geneqc classifier: false ORF vs conserved-hypothetical decisions and reasons."""

from pufprobe.config import load_config
from pufprobe.geneqc.classifier import classify_all
from pufprobe.models import CDSRecord, GeneStatus, GeneVerdict


def _record(locus: str, seq: str, start: int) -> CDSRecord:
    return CDSRecord(
        sample_id="s1",
        locus_tag=locus,
        protein_id=locus,
        sequence_aa=seq,
        start=start,
        end=start + len(seq) * 3,
        strand="+",
        product="hypothetical protein",
        status=GeneStatus.HYPOTHETICAL,
        contig="contig1",
    )


def test_false_orf_is_flagged() -> None:
    config = load_config(None)
    records = [
        _record("FAKE", "AAAAKNAAAAKNAAAAKNAAA", 100),   # 20 aa, low complexity
        _record("KEEP", "MARGKKIGYSGKTTAE" * 10, 1000),    # ~160 aa, normal
    ]
    summary = classify_all(records, config)

    fake, keep = records
    assert fake.verdict == GeneVerdict.LIKELY_FALSE_ORF
    assert fake.qc_score >= config["qc"]["false_orf_score"]
    # Reasons carry the threshold values.
    assert any("(min 30)" in reason for reason in fake.qc_reasons)
    assert any("low_complexity:" in reason for reason in fake.qc_reasons)
    assert any("RBS: not evaluated" in reason for reason in fake.qc_reasons)

    # A normal-length conserved hypothetical must not be called a false ORF.
    assert keep.verdict != GeneVerdict.LIKELY_FALSE_ORF
    assert keep.qc_score < config["qc"]["borderline_score"]
    assert summary["by_verdict"][GeneVerdict.LIKELY_FALSE_ORF.value] == 1


def test_smorf_path_does_not_penalise_length_alone() -> None:
    config = load_config(None)
    seq = "MPTVLISWREQANKDHY" * 2  # 34 aa, normal complexity
    record = _record("SMORF", seq, 1000)
    classify_all([record], config, include_smorf=True)
    assert record.verdict != GeneVerdict.LIKELY_FALSE_ORF
    assert not any("below smORF minimum" in reason for reason in record.qc_reasons)
