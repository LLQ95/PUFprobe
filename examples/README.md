# Demo genome

A frozen, tiny two-sample toy dataset used by the PUFprobe tests. It contains two
annotated genomes (`sampleA`, `sampleB`) in the layout expected by `pufprobe
prepare`:

```
demo_genome/
  sampleA/sampleA.gff3   # 13 CDS
  sampleA/sampleA.faa
  sampleB/sampleB.gff3   # 12 CDS
  sampleB/sampleB.faa
```

## What each genome contains

Both samples share the same eight annotated housekeeping / context genes with
identical protein sequences:

- DNA gyrase subunit A
- DNA-directed RNA polymerase subunit beta
- elongation factor Tu
- 30S ribosomal protein S5
- ABC transporter ATP-binding protein (efflux keyword)
- superoxide dismutase
- transposase (mobilisation keyword)
- toxin-antitoxin system toxin RelE (toxin keyword)
- DNA polymerase III subunit alpha

They also share two hypothetical proteins with identical sequences but different
locus tags:

- a conserved `hypothetical protein` (SA_0007 / SB_0007)
- a `DUF1620 domain-containing protein` (SA_0008 / SB_0008; treated as
  conserved hypothetical)

Sample-specific records:

- `SA_0012`: a ~20 aa, A-rich low-complexity false ORF (`hypothetical protein`).
- `SA_0013`: an orphan hypothetical protein unique to sampleA.
- `SB_0012`: a clonal-specific hypothetical protein unique to sampleB.

The GFF3 files are valid (`##gff-version 3`, `##sequence-region`, 9 tabbed
columns, `ID`/`Parent`/`locus_tag`/`product` attributes) and the `.faa` headers
are `>locus_tag product`, with protein lengths matching the CDS.

## Demo run

```powershell
# from the repository root, with PYTHONPATH pointing at src
$env:PYTHONPATH = "src"
pufprobe --workdir .\_demo prepare --input examples\demo_genome
pufprobe --workdir .\_demo context
```

The tests under `tests/` exercise exactly these two samples without invoking the
CLI.
