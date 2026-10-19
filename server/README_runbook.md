# PUFprobe 服务器运行手册（biflab-1，用户 liushiwei）

本手册面向一台不能从本机直连的 Linux 服务器。把整个仓库同步到 `/data/liushiwei/PUFprobe` 后，从第 0 步开始逐条执行即可。所有脚本都在仓库的 `server/` 目录下，统一在 `server/configs/server_config.yaml` 中集中配置，脚本自带日志、幂等、断点续跑，中断后重跑同一条命令会跳过已完成的步骤。

本手册以纯干性计算为目标，跑通从 contigs 到功能预测回填的全部流程。需要湿实验的环节只产出优先级清单，不执行实验。

## 环境事实与脚本默认值

| 项目 | 取值 |
| --- | --- |
| 主机 / 用户 | biflab-1 / liushiwei |
| 仓库同步根目录 | `/data/liushiwei/PUFprobe` |
| 已有环境 bakta_env | bakta 1.12.0，数据库 `/db/student/metagenome/bakta_db/db`，不含 prokka |
| 已有环境 checkm | prokka 1.15.6 |
| 新建环境 | `pufprobe`（python 3.11，独立，不改动上述两个环境） |
| 调用方式 | 一律 `conda run -n bakta_env` / `conda run -n checkm`，不激活、不修改它们 |
| GPU | 未知，脚本用 `nvidia-smi` 自动检测，给出 GPU 与 CPU 两条路径 |
| 数据库根目录 | `/data/liushiwei/PUFprobe/dbs`（可用 `DB_ROOT` 覆盖） |
| 线程数 | 默认 `nproc`，注释阶段每菌株 8 线程、2 株并行 |

每个脚本都用同样方式加载配置：先读 `server/configs/server_config.yaml`，再允许用环境变量覆盖。例如临时只跑 3 株：`MAX_SAMPLES=3 bash server/02_annotate.sh`。

## 总体路线与预期耗时

| 步骤 | 脚本 | 目的 | 典型耗时 |
| --- | --- | --- | --- |
| 第 0 步 | 无 | 同步代码 | 1–5 分钟 |
| 第 1 步 | 00_setup_env.sh | 建独立环境 | 15–40 分钟（取决于网络） |
| 第 2 步 | 01_download_dbs.sh（TIER=A） | 冒烟小库 | 20–60 分钟 |
| 第 3 步 | 02_annotate.sh | 注释 2–3 株 | 约 20–40 分钟 |
| 第 4 步 | 03_run_pufprobe.sh | 降级模式冒烟 | 5 分钟以内 |
| 第 5 步 | 04_summary.sh | 核对冒烟结果 | 1 分钟 |
| 第 6 步 | 01_download_dbs.sh（TIER=B） | 全量库 | 数小时到 1–2 天 |
| 第 7 步 | 02/03/04 | 全队列运行 | 数小时到数天 |

---

## 第 0 步：同步代码到服务器

在能同时访问本机与服务器的机器上，任选一种方式。

方式一，用 scp 传压缩包（最简单，不依赖服务器上的 GitHub 认证）：

```bash
scp PUFprobe-v0.2.1.zip liushiwei@biflab-1:/data/liushiwei/
ssh liushiwei@biflab-1
cd /data/liushiwei
unzip -o PUFprobe-v0.2.1.zip -d PUFprobe
cd PUFprobe
ls server/
```

`server/` 下应能看到 `00_setup_env.sh` 到 `04_summary.sh`、`summarize.py`、`configs/server_config.yaml` 和本手册。

方式二，直接在服务器上克隆公开仓库（需要服务器能访问 GitHub，国内网络可能失败，失败就用方式一）：

```bash
cd /data/liushiwei
git clone https://github.com/LLQ95/PUFprobe.git
cd PUFprobe
```

同步后给脚本加可执行权限（zip 解压后可能丢失）：

```bash
chmod +x server/*.sh
```

用 `bash -n` 做一次语法自检，不执行任何内容：

```bash
for f in server/*.sh; do bash -n "$f" && echo "OK $f"; done
```

预期屏幕输出五行 `OK ...`。

成功判据：`/data/liushiwei/PUFprobe/server/` 与 `src/pufprobe/` 都存在，`pyproject.toml` 存在。

---

## 第 1 步：创建 pufprobe 独立环境

```bash
cd /data/liushiwei/PUFprobe
bash server/00_setup_env.sh
```

这一步做了什么：检测 conda 或 mamba，把 MAMBA_ROOT_PREFIX 固定到 conda 根目录，避免在当前激活环境下建出嵌套环境；新建环境 `pufprobe`，从清华 TUNA 的 conda-forge 与 bioconda 频道安装 mmseqs2、foldseek、hhsuite、scikit-learn 与 openjdk 17；运行 `nvidia-smi` 检测，按驱动支持的 CUDA 版本装 cu121 或 cu118 的 GPU 版 torch，检测不到 GPU 就装 CPU 版；torch 依次尝试阿里云、上海交大、官方三个轮子源；从清华 PyPI 镜像装 fair-esm；P2Rank 不在 bioconda（该频道没有 p2rank 包），脚本改为从 GitHub Release 下载 2.5.1 便携版到 `thirdparty/p2rank`（需要 Java 17，已随环境安装），下载失败只警告不阻塞，口袋通道会在运行时自动跳过；最后 `pip install -e .` 安装本仓库，并执行 `pufprobe doctor`。

预期屏幕输出：日志依次出现 `Creating conda environment`、`Installing bioinformatics tools`、`Installing torch`、`Installing fair-esm`、`Installing PUFprobe`、`Tool versions`，最后是一段 JSON 形式的 doctor 报告。

成功判据：

1. 日志末尾出现 `Environment ready`；
2. 手动复核二进制可发现：

```bash
conda run -n pufprobe mmseqs version
conda run -n pufprobe foldseek version
conda run -n pufprobe hhblits 2>&1 | head -2
conda run -n pufprobe pufprobe doctor
```

doctor 输出里 mmseqs2、foldseek、hhblits 应为 available；有 GPU 时 torch 段 `cuda_available=true`，无 GPU 时为 false，这是正常的，不算错误。

后台运行：这一步依赖网络，可能比预期慢，建议在 tmux 中执行：

```bash
tmux new -s setup
bash server/00_setup_env.sh
# 按 Ctrl+b 再按 d 脱离；tmux attach -t setup 回到会话
```

日志在 `server/logs/00_setup_env_*.log`。

常见报错：

- conda 找不到：先装 Miniconda。`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`，`bash Miniconda3-latest-Linux-x86_64.sh -b`，重新登录后再跑。
- 清华 conda 频道超时：把 `server_config.yaml` 里的 `CONDA_MIRROR_HOST` 改成 `https://mirrors.bfsu.edu.cn/anaconda` 重跑，已有环境会被复用，不会重复安装。
- torch 轮子下载失败或超时：脚本会自动依次尝试阿里云、上海交大、官方三个源；三个源都失败时可先跳过（安装不阻塞后续工具的部分），事后用 `conda run -n pufprobe pip install torch --index-url https://mirrors.aliyun.com/pytorch-wheels/cu121/` 补装。
- 驱动 CUDA 版本低于 11.8：脚本会改装 CPU torch 并在日志中提示；要用 GPU 需管理员更新 NVIDIA 驱动。
- 环境被建到嵌套路径（如 `anaconda3/envs/checkm/envs/pufprobe`）：这是旧版脚本未固定 MAMBA_ROOT_PREFIX 导致的。v0.2.1 已修复，重跑时会自动删除嵌套环境并改到 `anaconda3/envs/pufprobe`。手工清理命令为 `rm -rf /data/liushiwei/anaconda3/envs/checkm/envs/pufprobe`。
- 安装 p2rank 时报 `p2rank does not exist (perhaps a typo or a missing channel)`：bioconda 不提供 p2rank，旧脚本误把它列入 conda 安装。v0.2.1 已改为 GitHub Release 便携版安装；若便携版下载失败，口袋通道跳过，不影响其余全部流程。
- P2Rank 下载失败：脚本已依次尝试 GitHub 直连与两个代理镜像；仍失败时可手动下载 `p2rank_2.5.1.tar.gz`（https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz），解压后将目录重命名为 `/data/liushiwei/PUFprobe/thirdparty/p2rank`，保证存在 `thirdparty/p2rank/distro/prank`。
- 提示已有环境存在：脚本直接复用；若环境是半成品，可 `conda env remove -n pufprobe` 后重跑，这只影响新建的 pufprobe 环境。

---

## 第 2 步：下载冒烟数据库（档位 A）

```bash
cd /data/liushiwei/PUFprobe
df -h /data/liushiwei
TIER=A bash server/01_download_dbs.sh
```

档位 A 内容：MMseqs2 的 Swiss-Prot 序列库、Foldseek 的 PDB 结构库（约 10–15 GB，磁盘紧张可在配置里设 `SMOKE_FOLDSEEK=0` 跳过）、CARD 蛋白同源模型 FASTA、VFDB 蛋白集 A、ESM-2 35M 小模型权重。全部是冒烟规模，总占用约 15–20 GB。

预期屏幕输出：依次出现 `Building MMseqs2 database 'Swiss-Prot'`、`Building Foldseek database 'PDB'`、`Downloading CARD`、`Downloading VFDB`、`Pre-loading esm2_t12_35M_UR50D`，末尾打印 `$DBDIR` 下的文件清单与 `Done`。

成功判据：

```bash
ls /data/liushiwei/PUFprobe/dbs/dbdir/swissprot.index
ls /data/liushiwei/PUFprobe/dbs/dbdir/foldseek_db.index
grep -c '^>' /data/liushiwei/PUFprobe/dbs/card/protein_fasta_protein_homolog_model.fasta
grep -c '^>' /data/liushiwei/PUFprobe/dbs/vfdb/VFDB_setA_pro.fas
```

前两个 `ls` 不报错；CARD 记录数应为数千条，VFDB 应为数百条。

这一步最容易因网络失败：Foldseek 与 MMseqs2 的库由官方下载器从境外站点拉取，国内网络可能超时。脚本在失败时会给出 MMseqs2 的 mmseqs.com tarball 备选；Foldseek 若持续失败，可在能访问外网的机器上执行同样的 `foldseek databases PDB <prefix> <tmp>` 再把整个前缀目录拷回 `dbs/dbdir/`。下载脚本支持断点续跑（`wget -c` 与完成标记），中断后重跑即可。

---

## 第 3 步：注释 2–3 株基因组（冒烟）

先告诉脚本 contigs 在哪。编辑 `server/configs/server_config.yaml`，把 `CONTIGS_DIR=` 改成真实目录，例如：

```bash
CONTIGS_DIR=/data/liushiwei/salmonella_contigs
```

该目录下每个菌株放一个多 contig FASTA（扩展名 `.fna`、`.fa`、`.fasta` 或 `.fas`，支持 `.gz`）。先只跑 3 株：

```bash
cd /data/liushiwei/PUFprobe
MAX_SAMPLES=3 bash server/02_annotate.sh
```

脚本对每株分别调用 `conda run -n bakta_env bakta`（输出 `.gff3`、`.gbff`、`.faa`）和 `conda run -n checkm prokka`（输出 `.gff`、`.gbk`、`.faa`、`.ffn`），结果归集到 `annotations/bakta/<菌株>/` 与 `annotations/prokka/<菌株>/`。每株完成后写 `.bakta_done` 或 `.prokka_done`，重跑自动跳过。

预期输出：日志开头列出找到的 3 个基因组文件，随后出现 `[bakta] <菌株> ->` 与 `[prokka] <菌株> ->`，末尾给出两个输出目录。

成功判据：

```bash
ls annotations/bakta/
for d in annotations/bakta/*/; do ls "$d" | grep -E '\.gff3$|\.faa$|\.gbff$'; done
for d in annotations/prokka/*/; do ls "$d" | grep -E '\.gff$|\.faa$|\.gbk$'; done
```

每个菌株目录下三种产物齐全。单株典型耗时：bakta 8 线程约 10–20 分钟，prokka 约 2–5 分钟；两株并行时 3 株约 20–40 分钟。

常见报错：

- bakta 报数据库路径错误：核对 `/db/student/metagenome/bakta_db/db` 下存在 `light.db` 等文件；不要改动 bakta_env，路径错了只改 `BAKTA_DB` 变量。
- 找不到 contigs：检查 `CONTIGS_DIR` 与扩展名；脚本只扫描该目录第一层，不递归。
- 某株失败：该株没有 `.bakta_done` 标记，其他株不受影响；看 `server/logs/02_annotate_*.log` 中对应段落，修复后重跑，已完成的株会跳过。

---

## 第 4 步：PUFprobe 降级模式冒烟

冒烟阶段强制不使用任何深度模型，验证整条管线在纯 CPU、无 torch 调用的情况下能跑通：

```bash
cd /data/liushiwei/PUFprobe
FORCE_FALLBACK=1 bash server/03_run_pufprobe.sh
```

默认用 bakta 注释结果（`ANNOTATION_SOURCE=bakta`）。脚本依次执行 prepare、embed、fold、search、context、rank、report、annotate 八个阶段，每个阶段写 `_DONE` 标记；embed 用 48 维序列特征回退，fold 因无模型把所有家族标记为 skipped，search 自动只启用第 2 步已下载的库（mmseqs 与 foldseek 标记齐全则都启用，hhblits 自动跳过）。

预期输出：每个阶段一行 `=== stage <名称> START`，阶段成功后出现 `=== stage <名称> OK (Ns)`；fold 阶段日志里会有家族 skipped 的提示，这是预期行为。

成功判据：

```bash
RUN=runs/salmonella_bakta
ls $RUN/00_prepare/samples.json $RUN/00_prepare/cds_records.json
ls $RUN/02_embed/embeddings.npz $RUN/02_embed/manifest.json
cat $RUN/02_embed/manifest.json    # 应显示 source 为 fallback
ls $RUN/07_rank/predictions.json $RUN/08_report/report.html $RUN/08_report/report.md
ls $RUN/09_annotate/
find $RUN -maxdepth 2 -name _DONE | sort
```

`manifest.json` 里 `source` 为 fallback，`09_annotate/` 下有每个菌株的 `.pufprobe.gff3`、`.pufprobe.gbk`、`.pufprobe.faa` 与 `changes.tsv`。十余个 `_DONE` 标记齐全。降级冒烟全程通常在 5 分钟以内。

失败处理：日志在 `server/logs/03_<阶段>_*.log`，脚本会直接打印失败阶段与日志路径；修复后重跑同一条命令，已完成阶段自动跳过。

---

## 第 5 步：核对冒烟结果

```bash
bash server/04_summary.sh
```

结果在 `results/salmonella_bakta/`：`summary.md` 是中文汇总（每株预测数、功能条目数、高置信条目数、各通道证据数、湿实验优先候选），其余是制表符分隔的 CSV 表（`per_sample_summary.csv`、`channel_hits.csv`、`family_overview.csv`、`high_confidence_predictions.csv`、`wetlab_priorities_top30.csv`）。

冒烟数据是 3 株，系统发育谱需要至少 4 个基因组才会充分发挥作用，因此冒烟时 context 通道命中偏少是正常的。判据是管线无中断、产物齐全、数字可复核，而不是高置信预测的数量。

至此最小闭环已验证。下面两步才是全量运行。

---

## 第 6 步：下载全量数据库（档位 B，务必先核对磁盘）

先确认 GPU 与磁盘：

```bash
nvidia-smi
df -h /data/liushiwei
```

全量库磁盘占用（含解压与索引，约数）：

| 数据库 | 体积 | 说明 |
| --- | --- | --- |
| MMseqs2 UniRef90 | 约 150 GB | 序列远程同源主力库 |
| Foldseek AlphaFold/Proteome | 约 1 TB 量级 | 结构库；以下载页实时大小为准 |
| hhblits UniRef30（2023_02） | 约 60–80 GB | 谱检索库 |
| BFD | 压缩约 272 GB、解压约 1.8 TB | 默认不下载，需要时设 `DOWNLOAD_BFD=1` |
| CARD + VFDB | 约 2 GB | 档位 A 已下载 |
| ESM-2 650M 权重 | 约 2.5 GB | ESMFold 权重另约 3 GB |
| 数十株工作目录 | 约 10–30 GB | 嵌入与中间产物 |

不开 BFD、不跑 AlphaFold/UniProt 全库时，建议预留至少 300 GB 空闲；要下 Foldseek Proteome 全库，建议预留 1.5 TB 以上。`df -h` 显示的可用空间不足时，先把 `DB_ROOT` 改到大容量盘，或把 `FULL_FOLDSEEK_DB` 留到后续再单独下载。

确认后在后台执行：

```bash
tmux new -s dbs
cd /data/liushiwei/PUFprobe
TIER=B bash server/01_download_dbs.sh
# Ctrl+b d 脱离；tmux attach -t dbs 查看
```

档位 B 内容：MMseqs2 UniRef90、Foldseek AlphaFold/Proteome、UniRef30 hhsuite（解压后自动在 `dbs/dbdir/` 建 `hhblits_db_*` 软链接）、CARD/VFDB、ESM-2 650M 权重，有 GPU 时预载 ESMFold 权重。UniRef50 可按需设 `DOWNLOAD_UNIREF50=1`。

成功判据：

```bash
ls /data/liushiwei/PUFprobe/dbs/dbdir/uniref90.index
ls /data/liushiwei/PUFprobe/dbs/dbdir/foldseek_db.index
ls /data/liushiwei/PUFprobe/dbs/dbdir/hhblits_db_hhm.ffdata
conda run -n pufprobe mmseqs linfo /data/liushiwei/PUFprobe/dbs/dbdir/uniref90 | head
```

这一步是全程网络风险最高的一步：Foldseek 全库和 UniRef30 从 gwdg 等境外站点拉取，体量大、耗时长，最可能在夜间网络波动时中断。脚本所有下载都带 `wget -c` 续传和完成标记，直接重跑即可续上；若某站点长期不可达，在有网机器下载后拷回相应路径（UniRef30 解包后保证 `dbs/dbdir/hhblits_db_{hhm,a3m,cs}.ffdata/.ffindex` 六个软链接存在）。

---

## 第 7 步：全队列注释与运行

先注释全部菌株（去掉 `MAX_SAMPLES` 限制）：

```bash
tmux new -s annotate
cd /data/liushiwei/PUFprobe
MAX_SAMPLES=0 bash server/02_annotate.sh
```

数十株规模下注释需要几个小时。若只想要 bakta 一套结果（两种注释不必都跑），把配置里 `RUN_PROKKA=0`。

然后跑 PUFprobe 全流程，分三种算力情形：

有 GPU 且驱动正常（`nvidia-smi` 可用），使用完整模型：

```bash
tmux new -s run
cd /data/liushiwei/PUFprobe
FORCE_FALLBACK=0 bash server/03_run_pufprobe.sh
```

无 GPU 但想用小模型获得真实嵌入（慢，每株蛋白数量决定，可能需要数小时到一天）：

```bash
FORCE_FALLBACK=0 USE_CPU_ESM=1 bash server/03_run_pufprobe.sh
```

无 GPU 且只求快速序列特征结果：

```bash
FORCE_FALLBACK=1 bash server/03_run_pufprobe.sh
```

GPU 路径下 embed 默认使用 650M 模型（约 8 GB 显存，batch 8）；显存不足时把 `ESM_BATCH_SIZE` 调到 2 或 4。ESMFold 建议 16–24 GB 显存，超过 1000 aa 的蛋白会被自动跳过，pLDDT 低于 `MIN_PLDDT`（默认 70）的结构不进入结构检索。若要在 CPU 上也跑 ESMFold（非常慢，仅适合少量代表蛋白），设 `FOLD_ON_CPU=1`。

运行完成后汇总：

```bash
bash server/04_summary.sh
```

结果在 `results/salmonella_bakta/`。若同时需要 prokka 注释口径的结果，把 `ANNOTATION_SOURCE=prokka` 再跑一次第 7 步（工作目录与结果目录按来源分开命名，不会覆盖）。

各阶段参考耗时（30–50 株，32 核、单卡 A100 量级）：prepare 少于 5 分钟；embed 650M 约 2–4 小时，CPU 回退约 1–2 天；fold 每个代表蛋白 1–3 分钟，取决于家族数；search 全库约 3–8 小时；context、rank、report、annotate 合计少于 30 分钟。

---

## 常用配置变量速查

| 变量 | 默认值 | 含义 |
| --- | --- | --- |
| CONTIGS_DIR | 空，必填 | 每个菌株一个 FASTA 的 contigs 目录 |
| THREADS | nproc | 各工具线程数 |
| TIER | A | 01 脚本数据库档位，A 冒烟、B 全量 |
| MAX_SAMPLES | 0 | 注释前 N 株，0 为全部；冒烟设 3 |
| RUN_BAKTA / RUN_PROKKA | 1 / 1 | 是否跑两套注释 |
| ANNOT_JOBS / ANNOT_THREADS | 2 / 8 | 并行菌株数与每菌株线程 |
| FORCE_FALLBACK | 1 | 1 时纯序列特征，冒烟最快 |
| USE_CPU_ESM | 0 | 无 GPU 时设 1 跑 35M 小模型 |
| ESM_BATCH_SIZE | 8 | GPU 显存紧张时调到 2–4 |
| MIN_PLDDT | 70 | 结构入库的 pLDDT 门限 |
| HIGH_CONFIDENCE | 0.80 | 高置信与回填阈值 |
| ANNOTATION_SOURCE | bakta | 使用哪套注释，可换 prokka |
| WETLAB_TOP_N | 30 | 湿实验优先候选条数 |
| DB_ROOT | /data/liushiwei/PUFprobe/dbs | 全部数据库根目录 |

## 从 contigs 开始的最短操作序列

```bash
cd /data/liushiwei/PUFprobe
bash server/00_setup_env.sh
TIER=A bash server/01_download_dbs.sh
# 编辑 server/configs/server_config.yaml 填好 CONTIGS_DIR
MAX_SAMPLES=3 bash server/02_annotate.sh
FORCE_FALLBACK=1 bash server/03_run_pufprobe.sh
bash server/04_summary.sh
# 冒烟通过后
TIER=B bash server/01_download_dbs.sh
MAX_SAMPLES=0 bash server/02_annotate.sh
FORCE_FALLBACK=0 bash server/03_run_pufprobe.sh
bash server/04_summary.sh
```

先注释 2–3 株，从建环境到冒烟核对约 1–2 小时，主要变量是网络。

## 可能需要人工介入的环节

| 环节 | 原因 | 处理 |
| --- | --- | --- |
| torch 安装 | 国内访问 pytorch 轮子源不稳定 | 改用阿里云或上海交大镜像，见第 1 步 |
| Foldseek 全库下载 | 数据量大、境外站点 | 断点续跑或异地下载拷回 |
| UniRef30 下载 | gwdg 站点国内访问不稳定 | 同上，核对六个软链接 |
| MMseqs2 全库 | 约 150 GB，耗时长 | 放 tmux 后台，重跑自动续传 |
| GPU 可用性 | 驱动版本、显存未知 | 脚本自动检测；显存不足调小 batch |
| ESM 权重 | 首次从境外下载 | 可在有网机器预下载后放入 `$TORCH_HOME/hub/checkpoints/` |
| CARD/VFDB | 小库但站点偶尔波动 | 浏览器下载后放到配置指定路径 |

## 常见报错速查

- 阶段重跑却跳过：阶段有 `_DONE` 标记，属于正常断点续跑；需要重算加 `--force`，例如 `conda run -n pufprobe pufprobe -w runs/xxx -t 8 rank --force`。
- search 提示某通道没有库：`dbs/dbdir` 下缺少对应前缀文件，该通道自动跳过，不影响其他通道；下完档位 B 后重跑 search。
- fold 全部 skipped：冒烟（强制回退）或无 GPU 环境下这是预期；要折叠代表蛋白需装 torch 且有 GPU，或显式设 `FOLD_ON_CPU=1`。
- context 高置信条目很少：冒烟只有 3 株，系统发育谱需要至少 4 个基因组，全队列运行后会改善。
- GFF 与 FAA 对不上：脚本按 locus_tag 与蛋白头标记关联；检查注释产物中两类标识是否一致，不要手工改产物。
- 中文或带空格路径：脚本对路径统一加引号；仍建议数据库与工作目录放在 `/data/liushiwei` 下，避免在挂载点用空格目录名。
