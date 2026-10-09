"""One-off generator for the frozen demo genome under examples/demo_genome.

Run once with:  python examples/_make_examples.py
It emits sampleA/sampleB GFF3 + FAA with matched protein lengths.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def rep(motif: str, n: int) -> str:
    return (motif * ((n // len(motif)) + 1))[:n]


HK = "MARGKKIGYSGKTTAE"          # diverse housekeeping motif
CHP = "MPTVLISWREQANKDHY"        # conserved hypothetical motif
DUFM = "MSQDFGPLKVNAWERTY"       # DUF motif
ORPHAN = "MSEVPLINKDWRQTAF"     # sampleA-only orphan
CLONAL = "MQVWPLINKDREQTAF"      # sampleB-only clonal


# (locus, product, length_aa, strand, start, end)
def build(prefix: str, extra_hypothetical: list[str]) -> list[tuple]:
    rows = [
        (f"{prefix}_0001", "DNA gyrase subunit A", 450, "+", 400, 1750),
        (f"{prefix}_0002", "DNA-directed RNA polymerase subunit beta", 600, "+", 1800, 3600),
        (f"{prefix}_0003", "elongation factor Tu", 400, "+", 3650, 4850),
        (f"{prefix}_0004", "30S ribosomal protein S5", 167, "+", 4900, 5400),
        (f"{prefix}_0005", "ABC transporter ATP-binding protein", 250, "+", 5450, 6200),
        (f"{prefix}_0006", "superoxide dismutase", 200, "+", 6250, 6850),
        (f"{prefix}_0007", "hypothetical protein", 150, "+", 6900, 7350),
        (f"{prefix}_0008", "DUF1620 domain-containing protein", 120, "+", 7400, 7760),
        (f"{prefix}_0009", "transposase", 300, "+", 8400, 9300),
        (f"{prefix}_0010", "toxin-antitoxin system toxin RelE", 90, "+", 9350, 9620),
        (f"{prefix}_0011", "DNA polymerase III subunit alpha", 350, "+", 9680, 10730),
    ]
    return rows


def sequence_for(product: str, length: int) -> str:
    low = product.lower()
    if "hypothetical protein" in low and "duf" not in low:
        return rep(CHP, length)
    if "duf" in low:
        return rep(DUFM, length)
    return rep(HK, length)


def write_sample(name: str, rows: list[tuple], special: dict) -> None:
    sdir = ROOT / "demo_genome" / name
    sdir.mkdir(parents=True, exist_ok=True)
    gff = [
        "##gff-version 3",
        "##sequence-region contig1 1 10000",
    ]
    faa = []
    for locus, product, length, strand, start, end in rows:
        gff.append(
            f"contig1\tdemo\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID=gene_{locus};locus_tag={locus}"
        )
        gff.append(
            f"contig1\tdemo\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
            f"ID=cds_{locus};Parent=gene_{locus};locus_tag={locus};product={product}"
        )
        seq = sequence_for(product, length)
        faa.append(f">{locus} {product}")
        for i in range(0, len(seq), 60):
            faa.append(seq[i : i + 60])
    # sample-specific extra CDS
    for locus, product, length, strand, start, end, seq in special:
        gff.append(
            f"contig1\tdemo\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID=gene_{locus};locus_tag={locus}"
        )
        gff.append(
            f"contig1\tdemo\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
            f"ID=cds_{locus};Parent=gene_{locus};locus_tag={locus};product={product}"
        )
        faa.append(f">{locus} {product}")
        for i in range(0, len(seq), 60):
            faa.append(seq[i : i + 60])

    (sdir / f"{name}.gff3").write_text("\n".join(gff) + "\n", encoding="utf-8", newline="\n")
    (sdir / f"{name}.faa").write_text("\n".join(faa) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {name}: {len(rows) + len(special)} CDS")


def main() -> None:
    sampleA = build("SA", [])
    sampleB = build("SB", [])
    write_sample("sampleA", sampleA, special=[
        ("SA_0012", "hypothetical protein", 20, "+", 7800, 7860, "AAAAKNAAAAKNAAAAKNAAA"),
        ("SA_0013", "hypothetical protein", 100, "+", 8000, 8300, rep(ORPHAN, 100)),
    ])
    write_sample("sampleB", sampleB, special=[
        ("SB_0012", "hypothetical protein", 110, "+", 7800, 8130, rep(CLONAL, 110)),
    ])


if __name__ == "__main__":
    main()
