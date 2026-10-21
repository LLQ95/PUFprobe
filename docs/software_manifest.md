# Software and Database Manifest

This document is the bill of materials for PUFprobe. It pins the versions of
the software, Python packages and reference databases used in the tested
environment, and lists the optional tools reserved for later stages. The
machine-readable Python pins are in `requirements-lock.txt`.

Tested environment: host `biflab-1`, NVIDIA RTX A5000, conda env `pufprobe`,
lock date 2026-10-21.

## Core pipeline software

| Component | Version | Role | Source |
|---|---|---|---|
| PUFprobe | 0.6.0 | Second-stage functional inference for uncharacterized CDS | https://github.com/LLQ95/PUFprobe |
| Python | 3.11.17 | Runtime | https://www.python.org |
| MMseqs2 | 18.8cc5c (Release 18) | Sequence search and linclust family clustering | https://github.com/soedinglab/MMseqs2 |
| Foldseek | 10.941cd33 (Release 10) | Structure search against PDB | https://github.com/steineggerlab/foldseek |
| HH-suite | conda bioconda build | Profile HMM (hhblits) remote homology | https://github.com/soedinglab/hh-suite |
| P2Rank | 2.5.1 | Ligand-binding pocket prediction (optional) | https://github.com/rdk/p2rank |

## Python packages

| Package | Version | Role | Source |
|---|---|---|---|
| PyTorch | 2.5.1+cu121 | ESM-2 and ESMFold inference | https://github.com/pytorch/pytorch |
| fair-esm | 2.0.0 | ESM-2 and ESMFold model definitions | https://github.com/facebookresearch/esm |
| NumPy | 2.4.6 | Numerical operations | https://github.com/numpy/numpy |
| Biopython | 1.88 | Sequence and GenBank parsing | https://github.com/biopython/biopython |
| scikit-learn | 1.9.1 | kNN graph, calibration, benchmark metrics | https://github.com/scikit-learn/scikit-learn |
| pandas | 3.0.6 | Tabular aggregation | https://github.com/pandas-dev/pandas |
| python-igraph | 1.0.0 | Graph handling for community detection | https://github.com/igraph/python-igraph |
| leidenalg | 0.12.0 | Leiden community detection | https://github.com/vtraag/leidenalg |
| Click | 8.5.0 | Command-line interface | https://github.com/pallets/click |
| einops | 0.8.2 | Tensor operations for ESMFold | https://github.com/arogozhnikov/einops |
| tqdm | 4.70.1 | Progress reporting | https://github.com/tqdm/tqdm |
| pytest | 9.1.1 | Test suite | https://github.com/pytest-dev/pytest |

## Reference databases

| Database | Version / record count | Refresh date | Role | Source |
|---|---|---|---|---|
| Pfam | Release 38.2, 30,134 family names | 2026-10-21 | Family identifiers and official descriptions | https://www.ebi.ac.uk/interpro/entry/pfam |
| UniProtKB/Swiss-Prot | 791,650 cross-referenced accessions | 2026-10-20 | Sequence search and id-mapping | https://www.uniprot.org |
| SIFTS | 954,958 PDB-chain to Pfam rows | 2026-10-21 | PDB-chain hits to Pfam mapping | https://www.ebi.ac.uk/pdbe/docs/sifts |
| PDB | Foldseek PDB database (EBI/PDBe mirror) | 2026-10-21 | Structure target library | https://www.rcsb.org |
| AlphaFold DB | model v6, 260M+ predicted structures | accessed 2026-10-21 | Precomputed structures for representatives | https://alphafold.ebi.ac.uk |
| CARD | 6,059 homolog proteins | 2026-10-21 | Antimicrobial resistance context | https://card.mcmaster.ca |
| VFDB | set A, 4,755 proteins | 2026-10-21 | Virulence factor context | http://www.mgc.ac.cn/VFs |

## Optional and planned tools

These are not installed in the pinned environment. They are reserved for
benchmark comparisons, the hypothetical protein database (HPDB), and later
evidence channels.

| Tool | Purpose | Source |
|---|---|---|
| InterProScan | Integrated reference annotation for external comparison | https://github.com/ebi-pf-team/interproscan |
| InterPro | Integrated family and domain resource | https://www.ebi.ac.uk/interpro |
| eggNOG-mapper | Orthology and functional annotation comparison | https://github.com/eggnogdb/eggnog-mapper |
| CATH | Protein structure classification | https://www.cathdb.info |
| TED | Protein domain structural resource | https://edwardslab.bmcb.georgetown.edu |
| TM-Vec | Structure-aware vector search without folding | https://github.com/tymor22/tm-vec |
| NetworkX | Additional graph analysis backend | https://github.com/networkx/networkx |
| DeepFRI | Structure-based function prediction comparison | https://github.com/viibrem/DeepFRI |
| ProstT5 | Sequence-structure language model embeddings | https://github.com/mheinzinger/ProstT5 |
| SaProt | Structure-aware protein language model | https://github.com/westlake-repl/SaProt |
| Datasette | Searchable web interface for the HPDB SQLite file | https://github.com/simonw/datasette |
| GTDB | Representative genome taxonomy for sampling | https://gtdb.ecogenomic.org |
| MGnify | Metagenomic protein resource | https://www.ebi.ac.uk/metagenomics |
| ESM Metagenomic Atlas | Metagenomic structure predictions | https://esmatlas.com |

## Reproducibility notes

The conda solver is mamba, and the CUDA driver is 550.120. When a database is
refreshed, its release date and record count are recorded by
`server/05_check_versions.sh`, and every prediction in the output keeps the
channel, method and database version used to derive it. Upgrading fair-esm or
PyTorch is not done without a full regression run, because the ESM packages are
no longer actively maintained.
