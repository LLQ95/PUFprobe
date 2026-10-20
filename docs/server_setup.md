# PUFprobe Server Deployment Checklist

This document describes how to set up the full PUFprobe environment on a fresh
Linux server (Ubuntu 22.04 / CentOS 7+, x86_64, ideally with an NVIDIA GPU).
Work through the sections in order. The CPU fallback path runs without a GPU,
but ESM-2 embedding and ESMFold are noticeably slower.

---

## 1. Getting the code and authenticating with GitHub

### 1.1 With the gh CLI (recommended)

```bash
# Option A: browser device-code flow
gh auth login

# Option B: log in with a token directly (common for CI / non-interactive servers)
echo "<YOUR_PAT>" | gh auth login --with-token
```

Minimum PAT scopes:

- `repo` (read/write private repositories)
- `read:org` (read organization information)

### 1.2 Without the gh CLI

```bash
# After cloning, replace the remote with a token-authenticated URL
git clone https://github.com/LLQ95/PUFprobe.git
cd PUFprobe
git remote set-url origin https://<YOUR_PAT>@github.com/LLQ95/PUFprobe.git
```

> Security notes:
> - Do not put the PAT in any file that is committed (`.git/config` is local and
>   fine; nothing inside the repository is).
> - Prefer `git config credential.helper store` with `~/.git-credentials`
>   (permission 600), or simply use `gh auth login`, to avoid the token appearing
>   in the shell history. You can `export GH_TOKEN=<PAT>` and then clear the
>   plaintext from history.

---

## 2. Conda / Mamba environment

### 2.1 Install micromamba (lightest option)

```bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C ~/ bin/micromamba
~/bin/micromamba shell init -s bash -r ~/micromamba
source ~/.bashrc
```

### 2.2 Or install Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
~/miniconda3/bin/conda init bash
source ~/.bashrc
```

### 2.3 Create the environment

Use the `environment.yml` bundled with the repository:

```bash
cd PUFprobe
micromamba env create -f environment.yml   # or conda env create -f environment.yml
micromamba activate pufprobe
```

If you prefer not to use the yml, the equivalent commands are:

```bash
micromamba create -n pufprobe -c conda-forge -c bioconda \
    python=3.11 click numpy pytest \
    mmseqs2>=13 foldseek>=8 hhsuite scikit-learn
micromamba activate pufprobe
pip install -e .
```

Channel order matters: use `conda-forge` first, then `bioconda`; all bioinformatics
tools come from these two channels. P2Rank is not available on bioconda and is
installed from its GitHub release (see section 3).

---

## 3. External tool versions (conda channels)

| Tool | Install command | Recommended version | Use |
| --- | --- | --- | --- |
| MMseqs2 | `micromamba install -c bioconda mmseqs2` | >=13.45 | Remote homology search, clustering |
| Foldseek | `micromamba install -c bioconda foldseek` | >=8.6 | Structure search |
| HH-suite | `micromamba install -c bioconda hhsuite` | >=3.3.0 | hhblits/hhsearch profile search |
| P2Rank | GitHub release `rdk/p2rank` (not on bioconda) | >=2.4 | Binding-pocket prediction |

Run the self-check after installation:

```bash
pufprobe doctor
```

It reports whether each binary is on the PATH and the torch/CUDA status.

---

## 4. PyTorch and ESM-2

### 4.1 GPU build (choose the index URL for your CUDA version)

```bash
# CUDA 12.1
pip install torch --index-url https://download.pytorch.org/whl/cu121

# CUDA 11.8
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

Check the local CUDA version with `nvidia-smi` (CUDA Version in the top right).

### 4.2 CPU build (when there is no GPU)

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu
```

### 4.3 fair-esm

```bash
pip install fair-esm
```

---

## 5. ESM-2 weights

- The weights are downloaded automatically on first use: the first time the
  `embed` stage runs ESM-2 (default `esm2_t33_650M_UR50D`), fair-esm stores the
  weights in `~/.cache/torch/hub/checkpoints/` (that is,
  `$TORCH_HOME/hub/checkpoints/`).
- For an offline server, download the file on a machine with internet (for
  example `esm2_t33_650M_UR50D.pt`, about 2.5 GB), copy it over with scp, and
  place it as follows:

  ```bash
  mkdir -p ~/.cache/torch/hub/checkpoints/
  cp esm2_t33_650M_UR50D.pt ~/.cache/torch/hub/checkpoints/
  ```

- Alternatively, point an environment variable to a custom directory:

  ```bash
  export TORCH_HOME=/data/models/torch
  mkdir -p $TORCH_HOME/hub/checkpoints
  # Put the weights into $TORCH_HOME/hub/checkpoints/
  ```

CPU mode (`--cpu`) uses the smaller `esm2_t12_35M_UR50D`, whose weights are about
130 MB.

---

## 6. Database download and build

Keep the databases in a single location, e.g. `/data/db/`, and pass it to
`search` with `--dbdir /data/db`.

### 6.1 UniRef90 / UniRef50 (MMseqs sequence databases)

```bash
mkdir -p /data/db/uniref90 && cd /data/db/uniref90
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz
# Build a searchable database from the UniProtKB/Swiss-Prot name (MMseqs v18):
mmseqs databases "UniProtKB/Swiss-Prot" /data/db/mmseqs/swissprot /data/tmp
# Build the uniref90 search database yourself:
mmseqs createdb uniref90.fasta.gz uniref90
mmseqs createindex uniref90 /data/tmp -t 16
# UniRef50 is the same: replace uniref90 with uniref50
```

### 6.2 PDB / AlphaFold DB (Foldseek structure databases)

```bash
# The official one-command download (Foldseek >=8)
mkdir -p /data/db/foldseek && cd /data/db/foldseek
foldseek databases Alphafold /data/db/foldseek/afdb /data/tmp
foldseek databases PDB /data/db/foldseek/pdb /data/tmp
# SwissProt structure database (small, for a quick smoke test)
foldseek databases SwissProt /data/db/foldseek/swissprot /data/tmp
```

### 6.3 BFD / UniClust (hhblits profile databases)

```bash
mkdir -p /data/db/hhblits && cd /data/db/hhblits
# BFD (about 1.8 TB unpacked; use a smaller database when disk is tight)
wget https://bfd.mmseqs.com/bfd_metaclust_clu_complete_id30_c90_final.sseadv.ffdata.tar.gz
# UniClust30 (recommended for routine use, about 60 GB)
wget https://wwwuser.gwdg.de/~compbiol/uniclust/2024_06/UniRef30_2024_06_hhsuite.tar.gz
tar xzf UniRef30_2024_06_hhsuite.tar.gz
```

### 6.4 CARD (antimicrobial resistance database)

```bash
mkdir -p /data/db/card && cd /data/db/card
# Download latest.tar.bz2 from the CARD download page (https://card.mcmaster.ca/download)
# It contains aro_index.json / card_database.fasta. For the RGI database:
pip install rgi
rgi card_annotation -i <json from the latest.tar.bz2>
rgi load --card_json card_json --local
```

### 6.5 VFDB (virulence factor database)

```bash
mkdir -p /data/db/vfdb && cd /data/db/vfdb
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz   # reference protein sequences
gunzip VFDB_setA_pro.fas.gz
```

### 6.6 SPI (Salmonella pathogenicity island reference, optional)

Prepare this from the SPI database or a FASTA directory maintained by your
group, and pass the corresponding genome FASTA directory to the context stage
with `--reference-fna`.

---

## 7. Disk and GPU memory requirements (rough estimates)

| Resource | Size | Notes |
| --- | --- | --- |
| MMseqs UniRef90 database | ~150 GB | including the index |
| Foldseek AlphaFold DB | ~1.2 TB | full database; smoke uses only SwissProt (~2 GB) |
| hhblits UniRef30 database | ~60 GB | BFD needs ~1.8 TB |
| CARD + VFDB | ~2 GB | |
| ESM-2 650M weights | ~2.5 GB | ESMFold weights need another ~3 GB |
| Working directory (tens of genomes) | ~10–30 GB | embeddings plus intermediate files |
| GPU memory: ESM-2 650M inference | ~8 GB | batch 8; reduce the batch on out-of-memory |
| GPU memory: ESMFold | ~16 GB minimum | proteins under 1000 aa ideally use 24 GB; folding is otherwise skipped |

---

## 8. Complete command sequence for dozens of Salmonella genomes

```bash
WORKDIR=/data/runs/salmonella_$(date +%Y%m%d)
DB=/data/db
THREADS=32

# 0) prepare: place all GFF3/FAA files under samples/sampleA/
pufprobe --workdir $WORKDIR --threads $THREADS prepare \
    --input /data/genomes/salmonella/

# 1) embed + cluster (GPU)
pufprobe --workdir $WORKDIR --threads $THREADS embed

# 2) fold (GPU; short sequences are folded automatically with ESMFold)
pufprobe --workdir $WORKDIR --threads $THREADS fold

# 3) search
pufprobe --workdir $WORKDIR --threads $THREADS search \
    --channels mmseqs,foldseek,hhblits --dbdir $DB

# 4) context (neighborhood / phylogenetic profile / microsynteny)
pufprobe --workdir $WORKDIR --threads $THREADS context \
    --card $DB/card --vfdb $DB/vfdb/VFDB_setA_pro.fas \
    --reference-fna /data/genomes/salmonella/

# 5) evidence fusion and ranking (optional --weights fusion.joblib for learned weights)
pufprobe --workdir $WORKDIR --threads $THREADS rank

# 6) report
pufprobe --workdir $WORKDIR --threads $THREADS report --formats json,html,md

# 7) annotation backfill
pufprobe --workdir $WORKDIR --threads $THREADS annotate
```

Estimated runtime (32 cores plus a single GPU on the order of an A100, roughly
30–50 genomes):

| Stage | Time |
| --- | --- |
| prepare | <5 min |
| embed (ESM-2 650M) | 2–4 h (GPU); CPU fallback ~1–2 days |
| fold (ESMFold) | ~1–3 min per representative protein, depending on their number |
| search (full databases) | 3–8 h |
| context / rank / report / annotate | <30 min total |

Resuming: each stage writes a `_DONE` marker, so re-running the same command
after an interruption skips completed stages; force a rerun with `--force`.

---

## 9. Common errors

- CUDA out of memory: reduce the embed batch size (`embed.batch_size` to 2 in the
  config), or let ESMFold skip proteins longer than 1000 aa. As a last resort,
  run `pufprobe --cpu ... embed --fallback-only` to use the feature fallback.
- mmseqs database path errors: make sure the directory given to `--dbdir`
  contains the database prefixes produced by `mmseqs createdb`, not the raw
  FASTA. In `pufprobe doctor`, `mmseqs available=True` only locates the binary;
  the database path must still be set correctly.
- foldseek missing structures: a "no usable structures" message at the
  05_search stage means the family pLDDT was below the threshold or it was
  skipped at 04_fold. The structure search channel is turned off automatically;
  the sequence channel is unaffected.
- GFF loci do not match FAA: if a CDS is not changed by annotate, the
  `locus_tag` in the GFF3 and the first identifier in the FAA header usually
  differ. Compare the ids in the two files with `grep`; prepare already links
  them by locus_tag / protein_id where possible.
- CRLF from Windows: the GFF3/TSV/JSON produced by this workflow are written
  with LF (`newline="\n"`). Whether a Windows editor displays them does not
  affect parsing on Linux; no conversion is needed across platforms.
