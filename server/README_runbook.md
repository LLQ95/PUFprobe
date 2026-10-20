# PUFprobe Server Runbook (biflab-1, user liushiwei)

This runbook is written for a Linux server that cannot be reached directly
from the local machine. After syncing the whole repository to
`/data/liushiwei/PUFprobe`, work through the steps in order starting from step 0.
All scripts live in the `server/` directory and are configured centrally in
`server/configs/server_config.yaml`. They provide logging, are idempotent, and
can resume after an interruption; re-running the same command skips the steps
that are already complete.

The goal is an entirely dry (computational) run covering the full path from
contigs to backfilled functional predictions. Steps that would require wet-lab
work only produce a priority list; no experiments are performed.

## Environment facts and script defaults

| Item | Value |
| --- | --- |
| Host / user | biflab-1 / liushiwei |
| Repository sync root | `/data/liushiwei/PUFprobe` |
| Existing environment bakta_env | bakta 1.12.0, database `/db/student/metagenome/bakta_db/db`, no prokka |
| Existing environment checkm | prokka 1.15.6 |
| New environment | `pufprobe` (python 3.11, standalone; the two environments above are not changed) |
| Invocation | always `conda run -n bakta_env` / `conda run -n checkm`; they are neither activated nor modified |
| GPU | unknown; the script detects it with `nvidia-smi` and provides both a GPU and a CPU path |
| Database root | `/data/liushiwei/PUFprobe/dbs` (overridable with `DB_ROOT`) |
| Threads | defaults to `nproc`; annotation uses 8 threads per strain, 2 strains in parallel |

Every script loads the configuration the same way: it first reads
`server/configs/server_config.yaml` and then allows environment variables to
override it. To process only 3 strains temporarily, for example:
`MAX_SAMPLES=3 bash server/02_annotate.sh`.

## Overall route and expected runtime

| Step | Script | Purpose | Typical time |
| --- | --- | --- | --- |
| Step 0 | none | Sync code | 1–5 min |
| Step 1 | 00_setup_env.sh | Build the standalone environment | 15–40 min (network-dependent) |
| Step 2 | 01_download_dbs.sh (TIER=A) | Smoke databases | 20–60 min |
| Step 3 | 02_annotate.sh | Annotate 2–3 strains | about 20–40 min |
| Step 4 | 03_run_pufprobe.sh | Fallback-mode smoke | under 5 min |
| Step 5 | 04_summary.sh | Check the smoke results | 1 min |
| Step 6 | 01_download_dbs.sh (TIER=B) | Full databases | hours to 1–2 days |
| Step 7 | 02/03/04 | Run the full cohort | hours to days |

---

## Step 0: Sync the code to the server

On a machine that can reach both the local machine and the server, choose one
method.

Method 1, copy the archive with scp (simplest; it does not depend on GitHub
authentication on the server):

```bash
scp PUFprobe-v0.2.5.zip liushiwei@biflab-1:/data/liushiwei/
ssh liushiwei@biflab-1
cd /data/liushiwei
unzip -o PUFprobe-v0.2.5.zip -d PUFprobe
cd PUFprobe
ls server/
```

Under `server/` you should see `00_setup_env.sh` through `04_summary.sh`,
`summarize.py`, `configs/server_config.yaml`, and this manual.

Method 2, clone the public repository directly on the server (requires the
server to reach GitHub; this may fail on some networks, in which case use method 1):

```bash
cd /data/liushiwei
git clone https://github.com/LLQ95/PUFprobe.git
cd PUFprobe
```

Make the scripts executable after syncing (the permission can be lost when
unzipping):

```bash
chmod +x server/*.sh
```

Run a one-time syntax check with `bash -n`, which executes nothing:

```bash
for f in server/*.sh; do bash -n "$f" && echo "OK $f"; done
```

The expected output is five lines of `OK ...`.

Success criterion: both `/data/liushiwei/PUFprobe/server/` and
`src/pufprobe/` exist, as does `pyproject.toml`.

---

## Step 1: Create the standalone pufprobe environment

```bash
cd /data/liushiwei/PUFprobe
bash server/00_setup_env.sh
```

What this step does: it detects conda or mamba and pins MAMBA_ROOT_PREFIX to
the conda root so a nested environment is not created inside the currently
active one. It creates the `pufprobe` environment and installs mmseqs2,
foldseek, hhsuite, scikit-learn, and openjdk 17 from the TUNA conda-forge and
bioconda channels. It runs the `nvidia-smi` check and, according to the CUDA
version supported by the driver, installs the cu121 or cu118 GPU build of
torch; when no GPU is detected it installs the CPU build. Torch wheels are
tried in order from Aliyun, SJTU, and the official source. It installs
fair-esm from the TUNA PyPI mirror. P2Rank is not on bioconda (that channel
does not provide a p2rank package), so the script instead downloads the 2.5.1
portable build from the GitHub release into `thirdparty/p2rank` (it needs Java
17, which is installed with the environment); a download failure only warns
and does not block, and the pocket channel is skipped automatically at runtime.
Finally it installs the repository with `pip install -e .` and runs
`pufprobe doctor`.

Expected output: the log shows `Creating conda environment`,
`Installing bioinformatics tools`, `Installing torch`, `Installing fair-esm`,
`Installing PUFprobe`, and `Tool versions` in order, followed by a JSON doctor
report.

Success criterion:

1. The log ends with `Environment ready`;
2. Manually confirm the binaries are discoverable:

```bash
conda run -n pufprobe mmseqs version
conda run -n pufprobe foldseek version
conda run -n pufprobe hhblits 2>&1 | head -2
conda run -n pufprobe pufprobe doctor
```

In the doctor output mmseqs2, foldseek, and hhblits should be available. With a
GPU the torch section shows `cuda_available=true`; without a GPU it shows false,
which is expected and not an error.

Running in the background: this step depends on the network and may be slower
than expected, so it is best run inside tmux:

```bash
tmux new -s setup
bash server/00_setup_env.sh
# Press Ctrl+b then d to detach; tmux attach -t setup to return
```

The log is at `server/logs/00_setup_env_*.log`.

Common errors:

- conda not found: install Miniconda first. `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`,
  `bash Miniconda3-latest-Linux-x86_64.sh -b`, then log in again and re-run.
- TUNA conda channel timeout: change `CONDA_MIRROR_HOST` in
  `server_config.yaml` to `https://mirrors.bfsu.edu.cn/anaconda` and re-run;
  the existing environment is reused and nothing is installed twice.
- torch wheel download fails or times out: the script tries Aliyun, SJTU, and
  the official source in order. If all three fail, you can skip it for now (the
  installation does not block the other tools) and add it later with
  `conda run -n pufprobe pip install torch --index-url https://mirrors.aliyun.com/pytorch-wheels/cu121/`.
- Driver CUDA version below 11.8: the script installs CPU torch instead and
  notes this in the log; using the GPU requires an administrator to update the
  NVIDIA driver.
- The environment was created at a nested path (for example
  `anaconda3/envs/checkm/envs/pufprobe`): this was caused by an older script
  that did not pin MAMBA_ROOT_PREFIX. v0.2.1 fixes it; a re-run removes the
  nested environment automatically and uses `anaconda3/envs/pufprobe`. The
  manual cleanup command is
  `rm -rf /data/liushiwei/anaconda3/envs/checkm/envs/pufprobe`.
- Installing p2rank reports `p2rank does not exist (perhaps a typo or a missing
  channel)`: bioconda does not provide p2rank, and an older script incorrectly
  listed it in the conda install. v0.2.1 uses the GitHub release portable build;
  if that download fails, the pocket channel is skipped and the rest of the
  pipeline is unaffected.
- P2Rank download fails: the script tries direct GitHub and two proxy mirrors in
  order. If it still fails, download `p2rank_2.5.1.tar.gz` manually
  (https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz),
  unpack it, and rename the directory to
  `/data/liushiwei/PUFprobe/thirdparty/p2rank`, ensuring
  `thirdparty/p2rank/distro/prank` exists.
- A message says the environment already exists: the script reuses it. If the
  environment is half-built, run `conda env remove -n pufprobe` and re-run; this
  only affects the newly created pufprobe environment.

---

## Step 2: Download the smoke databases (tier A)

```bash
cd /data/liushiwei/PUFprobe
df -h /data/liushiwei
TIER=A bash server/01_download_dbs.sh
```

Tier A contents: the MMseqs2 Swiss-Prot sequence database (in the new MMseqs2
its name is `UniProtKB/Swiss-Prot`; when the built-in downloader is unavailable,
it is downloaded automatically from the official UniProt FASTA and built with
createdb/createindex, about 90 MB compressed and several GB after building),
the Foldseek PDB structure database (about 10–15 GB; set `SMOKE_FOLDSEEK=0` in
the config to skip it when disk is tight), the CARD protein homolog model FASTA
(extracted from the official `https://card.mcmaster.ca/latest/data` bundle; the
old `/download/6/...fasta` direct link no longer works), VFDB protein set A,
and the ESM-2 35M small-model weights (a preload failure only warns; the smoke
fallback features do not need it). Everything is at smoke scale and uses about
15–20 GB in total.

Expected output: in order you see `Building MMseqs2 database
'UniProtKB/Swiss-Prot'` (if the built-in downloader is unavailable it then
prints `building from UniProt FASTA` together with createdb/createindex),
`Building Foldseek database 'PDB'`, `Downloading CARD canonical bundle`,
`Downloading VFDB`, and `Pre-loading esm2_t12_35M_UR50D`, followed by the list
of files under `$DBDIR` and `Done`.

Success criterion:

```bash
ls /data/liushiwei/PUFprobe/dbs/dbdir/swissprot.index
ls /data/liushiwei/PUFprobe/dbs/dbdir/foldseek_db.index
grep -c '^>' /data/liushiwei/PUFprobe/dbs/card/protein_fasta_protein_homolog_model.fasta
grep -c '^>' /data/liushiwei/PUFprobe/dbs/vfdb/VFDB_setA_pro.fas
```

The first two `ls` commands do not error; CARD should contain thousands of
records and VFDB hundreds.

This is the step most likely to fail because of the network. MMseqs2 already
has a fallback that builds from the official UniProt FASTA (`ftp.uniprot.org`,
verified reachable), so it no longer depends on the defunct mmseqs.com tarball.
If Foldseek keeps failing, run the same
`foldseek databases PDB <prefix> <tmp>` on a machine with external network
access and copy the whole prefix directory back to `dbs/dbdir/`. The download
script supports resuming (`wget -c` and completion markers), so re-run it after
an interruption. Note that an earlier failed run can leave an HTML placeholder
such as `swissprot.tar.gz` in `dbs/tmp/`; the new script removes it
automatically. If it still looks wrong, run `rm -rf dbs/tmp/*` manually and
re-run.

---

## Step 3: Annotate 2–3 genomes (smoke)

First tell the script where the contigs are. Edit
`server/configs/server_config.yaml` and change `CONTIGS_DIR=` to the real
directory, for example:

```bash
CONTIGS_DIR=/data/liushiwei/salmonella_contigs
```

That directory holds one multi-contig FASTA per strain (extension `.fna`,
`.fa`, `.fasta`, or `.fas`; `.gz` is supported). Process only 3 strains first:

```bash
cd /data/liushiwei/PUFprobe
MAX_SAMPLES=3 bash server/02_annotate.sh
```

For each strain the script calls `conda run -n bakta_env bakta` (outputs
`.gff3`, `.gbff`, `.faa`) and `conda run -n checkm prokka` (outputs `.gff`,
`.gbk`, `.faa`, `.ffn`), collecting results under `annotations/bakta/<strain>/`
and `annotations/prokka/<strain>/`. It writes `.bakta_done` or `.prokka_done`
after each strain, and a re-run skips them.

Expected output: the log starts by listing the 3 genome files found, then shows
`[bakta] <strain> ->` and `[prokka] <strain> ->`, and ends with the two output
directories.

Success criterion:

```bash
ls annotations/bakta/
for d in annotations/bakta/*/; do ls "$d" | grep -E '\.gff3$|\.faa$|\.gbff$'; done
for d in annotations/prokka/*/; do ls "$d" | grep -E '\.gff$|\.faa$|\.gbk$'; done
```

The three products are present in each strain directory. Typical time per
strain: bakta about 10–20 min on 8 threads, prokka about 2–5 min; with two
strains in parallel, 3 strains take about 20–40 min.

Common errors:

- bakta reports a database path error: confirm files such as `light.db` exist
  under `/db/student/metagenome/bakta_db/db`. Do not change bakta_env; if the
  path is wrong, only change the `BAKTA_DB` variable.
- contigs not found: check `CONTIGS_DIR` and the extensions; the script scans
  only the first level of that directory and does not recurse.
- A strain fails: it has no `.bakta_done` marker while the other strains are
  unaffected. Look at the corresponding section in
  `server/logs/02_annotate_*.log`, fix it, and re-run; completed strains are
  skipped.

---

## Step 4: PUFprobe fallback-mode smoke

The smoke stage forces no deep models and verifies that the whole pipeline runs
on pure CPU with no torch calls:

```bash
cd /data/liushiwei/PUFprobe
FORCE_FALLBACK=1 bash server/03_run_pufprobe.sh
```

The bakta annotation results are used by default (`ANNOTATION_SOURCE=bakta`).
The script runs the eight stages prepare, embed, fold, search, context, rank,
report, and annotate in order, and each stage writes a `_DONE` marker. Embed
uses the 48-dimensional sequence feature fallback; fold, with no model, marks
every family as skipped; search automatically enables only the databases
downloaded in step 2 (both mmseqs and foldseek are enabled when their markers
are present, and hhblits is skipped automatically).

Expected output: one line of `=== stage <name> START` for each stage, followed by
`=== stage <name> OK (Ns)` on success. The fold stage logs that families are
skipped, which is expected.

Success criterion:

```bash
RUN=runs/salmonella_bakta
ls $RUN/00_prepare/samples.json $RUN/00_prepare/cds_records.json
ls $RUN/02_embed/embeddings.npz $RUN/02_embed/manifest.json
cat $RUN/02_embed/manifest.json    # source should read fallback
ls $RUN/07_rank/predictions.json $RUN/08_report/report.html $RUN/08_report/report.md
ls $RUN/09_annotate/
find $RUN -maxdepth 2 -name _DONE | sort
```

In `manifest.json` the `source` is fallback, and under `09_annotate/` there are
`.pufprobe.gff3`, `.pufprobe.gbk`, `.pufprobe.faa`, and `changes.tsv` for each
strain. More than ten `_DONE` markers are present. The fallback smoke usually
finishes in under 5 minutes.

Handling failures: the logs are at `server/logs/03_<stage>_*.log`; the script
prints the failed stage and the log path. After fixing, re-run the same command
and completed stages are skipped automatically.

---

## Step 5: Check the smoke results

```bash
bash server/04_summary.sh
```

The results are in `results/salmonella_bakta/`: `summary.md` is the
human-readable overview (per-strain prediction counts, functional term counts,
high-confidence term counts, evidence counts per channel, and wet-lab priority
candidates), and the rest are tab-separated CSV tables
(`per_sample_summary.csv`, `channel_hits.csv`, `family_overview.csv`,
`high_confidence_predictions.csv`, `wetlab_priorities_top30.csv`).

The smoke data cover 3 strains, while the phylogenetic profile needs at least 4
genomes to be informative, so the context channel has few hits during the smoke
run and that is expected. The criterion is that the pipeline runs without
interruption, the products are complete, and the numbers can be checked, rather
than the number of high-confidence predictions.

This verifies the smallest closed loop. The next two steps are the full run.

---

## Step 6: Download the full databases (tier B; check disk first)

First confirm the GPU and disk:

```bash
nvidia-smi
df -h /data/liushiwei
```

Full database disk usage (including unpacking and indexing, approximate):

| Database | Size | Notes |
| --- | --- | --- |
| MMseqs2 UniRef90 | about 150 GB | main database for remote sequence homology |
| Foldseek AlphaFold/Proteome | on the order of 1 TB | structure database; the live size on the download page is authoritative |
| hhblits UniRef30 (2023_02) | about 60–80 GB | profile search database |
| BFD | about 272 GB compressed, 1.8 TB unpacked | not downloaded by default; set `DOWNLOAD_BFD=1` when needed |
| CARD + VFDB | about 2 GB | downloaded in tier A |
| ESM-2 650M weights | about 2.5 GB | ESMFold weights need another ~3 GB |
| Working directory for tens of genomes | about 10–30 GB | embeddings and intermediate products |

When you do not open BFD and do not run the full AlphaFold/UniProt database,
keep at least 300 GB free; to download the full Foldseek Proteome, keep more
than 1.5 TB free. When the free space shown by `df -h` is insufficient, point
`DB_ROOT` at a larger disk first, or leave `FULL_FOLDSEEK_DB` to be downloaded
separately later.

Run it in the background after confirming:

```bash
tmux new -s dbs
cd /data/liushiwei/PUFprobe
TIER=B bash server/01_download_dbs.sh
# Ctrl+b d to detach; tmux attach -t dbs to inspect
```

Tier B contents: MMseqs2 UniRef90, Foldseek AlphaFold/Proteome, UniRef30
hhsuite (after unpacking it creates the `hhblits_db_*` symlinks automatically
under `dbs/dbdir/`), CARD/VFDB, and the ESM-2 650M weights; on a GPU it
preloads the ESMFold weights. UniRef50 can be added by setting
`DOWNLOAD_UNIREF50=1` as needed.

Success criterion:

```bash
ls /data/liushiwei/PUFprobe/dbs/dbdir/uniref90.index
ls /data/liushiwei/PUFprobe/dbs/dbdir/foldseek_db.index
ls /data/liushiwei/PUFprobe/dbs/dbdir/hhblits_db_hhm.ffdata
conda run -n pufprobe mmseqs linfo /data/liushiwei/PUFprobe/dbs/dbdir/uniref90 | head
```

This is the step with the highest network risk in the whole run: the full
Foldseek database and UniRef30 are pulled from sites such as gwdg, are large
and slow, and are most likely to be interrupted during overnight network
fluctuations. All downloads use `wget -c` resuming and completion markers, so
re-run to continue. If a site is unreachable for a long time, download on a
machine with network and copy it back to the corresponding path (after
unpacking UniRef30, make sure the six
`dbs/dbdir/hhblits_db_{hhm,a3m,cs}.ffdata/.ffindex` symlinks exist).

---

## Step 7: Annotate and run the full cohort

First annotate all strains (remove the `MAX_SAMPLES` limit):

```bash
tmux new -s annotate
cd /data/liushiwei/PUFprobe
MAX_SAMPLES=0 bash server/02_annotate.sh
```

Annotating tens of genomes takes a few hours. If you only want one set of bakta
results (you do not need both annotations), set `RUN_PROKKA=0` in the config.

Then run the full PUFprobe pipeline; there are three compute scenarios.

With a GPU and a working driver (`nvidia-smi` is available), use the full
models:

```bash
tmux new -s run
cd /data/liushiwei/PUFprobe
FORCE_FALLBACK=0 bash server/03_run_pufprobe.sh
```

Without a GPU, when you want real embeddings from the small model (slow; the
number of proteins per strain determines it, and it may take hours to a day):

```bash
FORCE_FALLBACK=0 USE_CPU_ESM=1 bash server/03_run_pufprobe.sh
```

Without a GPU and when you only want fast sequence feature results:

```bash
FORCE_FALLBACK=1 bash server/03_run_pufprobe.sh
```

On the GPU path embed uses the 650M model by default (about 8 GB GPU memory,
batch 8); when memory is tight, set `ESM_BATCH_SIZE` to 2 or 4. ESMFold is best
run with 16–24 GB GPU memory; proteins longer than 1000 aa are skipped
automatically, and structures with pLDDT below `MIN_PLDDT` (default 70) do not
enter structure search. To also run ESMFold on CPU (very slow, suitable only for
a few representative proteins), set `FOLD_ON_CPU=1`.

Summarize when the run is complete:

```bash
bash server/04_summary.sh
```

The results are in `results/salmonella_bakta/`. If you also need results from
the prokka annotation, set `ANNOTATION_SOURCE=prokka` and run step 7 once more
(the working and results directories are named separately by source and do not
overwrite each other).

Reference runtime per stage (30–50 strains, 32 cores, a single GPU on the order
of an A100): prepare under 5 min; embed 650M about 2–4 h, CPU fallback about
1–2 days; fold 1–3 min per representative protein depending on the number of
families; search against the full databases about 3–8 h; context, rank, report,
and annotate together under 30 min.

---

## Quick reference for common configuration variables

| Variable | Default | Meaning |
| --- | --- | --- |
| CONTIGS_DIR | empty, required | contigs directory with one FASTA per strain |
| THREADS | nproc | thread count for each tool |
| TIER | A | database tier for the 01 script, A smoke, B full |
| MAX_SAMPLES | 0 | annotate the first N strains, 0 for all; set 3 for smoke |
| RUN_BAKTA / RUN_PROKKA | 1 / 1 | whether to run both annotations |
| ANNOT_JOBS / ANNOT_THREADS | 2 / 8 | parallel strains and threads per strain |
| FORCE_FALLBACK | 1 | when 1, pure sequence features, fastest for smoke |
| USE_CPU_ESM | 0 | set 1 to run the 35M small model without a GPU |
| ESM_BATCH_SIZE | 8 | reduce to 2–4 when GPU memory is tight |
| MIN_PLDDT | 70 | pLDDT gate for structures entering the database |
| HIGH_CONFIDENCE | 0.80 | high-confidence and backfill threshold |
| ANNOTATION_SOURCE | bakta | which annotation set to use, can be prokka |
| WETLAB_TOP_N | 30 | number of wet-lab priority candidates |
| DB_ROOT | /data/liushiwei/PUFprobe/dbs | root for all databases |

## Shortest command sequence starting from contigs

```bash
cd /data/liushiwei/PUFprobe
bash server/00_setup_env.sh
TIER=A bash server/01_download_dbs.sh
# Edit server/configs/server_config.yaml to fill in CONTIGS_DIR
MAX_SAMPLES=3 bash server/02_annotate.sh
FORCE_FALLBACK=1 bash server/03_run_pufprobe.sh
bash server/04_summary.sh
# After the smoke passes
TIER=B bash server/01_download_dbs.sh
MAX_SAMPLES=0 bash server/02_annotate.sh
FORCE_FALLBACK=0 bash server/03_run_pufprobe.sh
bash server/04_summary.sh
```

Annotate 2–3 strains first; from building the environment to checking the smoke
run takes about 1–2 h, with the network as the main variable.

## Steps that may need human intervention

| Step | Reason | Handling |
| --- | --- | --- |
| torch install | access to the pytorch wheel source is unstable on some networks | use the Aliyun or SJTU mirror, see step 1 |
| Full Foldseek download | large data, overseas sites | resume or download elsewhere and copy back |
| UniRef30 download | the gwdg site is unstable on some networks | same; check the six symlinks |
| Full MMseqs2 database | about 150 GB, slow | run in tmux; a re-run resumes automatically |
| GPU availability | driver version and GPU memory unknown | the script detects it; reduce the batch on low memory |
| ESM weights | first download from overseas | pre-download on a connected machine and put them in `$TORCH_HOME/hub/checkpoints/` |
| CARD/VFDB | small but the sites occasionally fluctuate | download in a browser and place them at the configured path |

## Common errors quick reference

- A stage is skipped when re-run: the stage has a `_DONE` marker and this is
  normal resuming; add `--force` to recompute, for example
  `conda run -n pufprobe pufprobe -w runs/xxx -t 8 rank --force`.
- search says a channel has no database: the corresponding prefix file is
  missing under `dbs/dbdir`, and that channel is skipped automatically without
  affecting the others; re-run search after tier B completes.
- fold is all skipped: expected during a smoke (forced fallback) or without a
  GPU; folding representative proteins requires torch and a GPU, or an explicit
  `FOLD_ON_CPU=1`.
- context has very few high-confidence items: the smoke covers only 3 strains and
  the phylogenetic profile needs at least 4 genomes; this improves after the
  full cohort runs.
- GFF and FAA do not match: the script links by locus_tag and the protein header
  identifier; check that the two identifiers in the annotation products agree,
  and do not edit the products by hand.
- Paths with spaces: the scripts quote paths consistently; still, it is best to
  keep databases and working directories under `/data/liushiwei` and avoid
  spaced directory names on mount points.
