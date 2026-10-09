# PUFprobe 服务器部署清单

本文档面向一台全新的 Linux 服务器（Ubuntu 22.04 / CentOS 7+，x86_64，建议带 NVIDIA GPU），照本清单逐条执行即可搭起 PUFprobe 全流程环境。无 GPU 也能跑通 CPU 回退路径，但 ESM-2 嵌入与 ESMFold 会明显变慢。

---

## 1. 获取代码与 GitHub 认证

### 1.1 有 gh CLI（推荐）

```bash
# 方式一：浏览器设备码流程
gh auth login

# 方式二：直接用 Token 登录（CI / 无交互服务器常用）
echo "<YOUR_PAT>" | gh auth login --with-token
```

PAT 所需最小 scope：

- `repo`（私有仓库读写）
- `read:org`（读取组织信息）

### 1.2 没有 gh CLI

```bash
# 克隆后把 remote 换成带 PAT 的地址
git clone https://github.com/LLQ95/PUFprobe.git
cd PUFprobe
git remote set-url origin https://<YOUR_PAT>@github.com/LLQ95/PUFprobe.git
```

> 安全提醒：
> - 不要把 PAT 写进任何会提交的文件（`.git/config` 是本地的，可以；仓库内任何文件都不行）。
> - 推荐用 `git config credential.helper store` 配合 `~/.git-credentials`（权限 600），或直接用 `gh auth login`，避免 PAT 出现在 shell 历史里：`export GH_TOKEN=<PAT>` 后 `unset` 掉历史中的明文。

---

## 2. Conda / Mamba 环境

### 2.1 安装 micromamba（最轻量）

```bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C ~/ bin/micromamba
~/bin/micromamba shell init -s bash -r ~/micromamba
source ~/.bashrc
```

### 2.2 或安装 Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
~/miniconda3/bin/conda init bash
source ~/.bashrc
```

### 2.3 创建环境

用仓库自带的 `environment.yml`：

```bash
cd PUFprobe
micromamba env create -f environment.yml   # 或 conda env create -f environment.yml
micromamba activate pufprobe
```

如果不想用 yml，逐条安装等价命令：

```bash
micromamba create -n pufprobe -c conda-forge -c bioconda \
    python=3.11 click numpy pytest \
    mmseqs2>=13 foldseek>=8 hhsuite p2rank scikit-learn
micromamba activate pufprobe
pip install -e .
```

channel 顺序：先 `conda-forge`，再 `bioconda`，所有生物工具都从这两个频道出。

---

## 3. 外部工具版本（conda 渠道）

| 工具 | 安装命令 | 建议版本 | 用途 |
| --- | --- | --- | --- |
| MMseqs2 | `micromamba install -c bioconda mmseqs2` | >=13.45 | 序列远程同源检索、聚类 |
| Foldseek | `micromamba install -c bioconda foldseek` | >=8.6 | 结构检索 |
| HH-suite | `micromamba install -c bioconda hhsuite` | >=3.3.0 | hhblits/hhsearch 谱检索 |
| P2Rank | `micromamba install -c bioconda p2rank` | >=2.4 | 结合口袋预测 |

安装后自检：

```bash
pufprobe doctor
```

会列出各二进制是否在 PATH、torch/CUDA 状态。

---

## 4. PyTorch 与 ESM-2

### 4.1 GPU 版（按 CUDA 版本选 index-url）

```bash
# CUDA 12.1
pip install torch --index-url https://download.pytorch.org/whl/cu121

# CUDA 11.8
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

查看本机 CUDA：`nvidia-smi` 右上角的 CUDA Version。

### 4.2 CPU 版备选（无 GPU 时）

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu
```

### 4.3 fair-esm

```bash
pip install fair-esm
```

---

## 5. ESM-2 权重

- **首次自动下载**：`embed` 阶段第一次跑 ESM-2（默认 `esm2_t33_650M_UR50D`）时，fair-esm 会自动把权重下载到 `~/.cache/torch/hub/checkpoints/`（即 `$TORCH_HOME/hub/checkpoints/`）。
- **离线服务器手动下载**：在有网机器上下载对应文件（如 `esm2_t33_650M_UR50D.pt`，约 2.5 GB），scp 到服务器后放到：

  ```bash
  mkdir -p ~/.cache/torch/hub/checkpoints/
  cp esm2_t33_650M_UR50D.pt ~/.cache/torch/hub/checkpoints/
  ```

- 或设置环境变量指向自定义目录：

  ```bash
  export TORCH_HOME=/data/models/torch
  mkdir -p $TORCH_HOME/hub/checkpoints
  # 把权重放进 $TORCH_HOME/hub/checkpoints/
  ```

CPU 模式（`--cpu`）默认使用更小的 `esm2_t12_35M_UR50D`，权重约 130 MB。

---

## 6. 数据库下载与建库

建议统一放在 `/data/db/`，用 `--dbdir /data/db` 传给 `search`。

### 6.1 UniRef90 / UniRef50（MMseqs 序列库）

```bash
mkdir -p /data/db/uniref90 && cd /data/db/uniref90
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz
mmseqs databases --tmpdir /data/tmp UniProtKB/Swiss-Prot swissprot /data/db/mmseqs
# 自建 uniref90 搜索库：
mmseqs createdb uniref90.fasta.gz uniref90
mmseqs createindex uniref90 /data/tmp -t 16
# UniRef50 同理：把 uniref90 换成 uniref50
```

### 6.2 PDB / AlphaFold DB（Foldseek 结构库）

```bash
# 官方脚本一键下载（Foldseek >=8）
mkdir -p /data/db/foldseek && cd /data/db/foldseek
foldseek databases Alphafold /data/db/foldseek/afdb /data/tmp
foldseek databases PDB /data/db/foldseek/pdb /data/tmp
# SwissProt 结构库（小，快速冒烟用）
foldseek databases SwissProt /data/db/foldseek/swissprot /data/tmp
```

### 6.3 BFD / UniClust（hhblits 谱库）

```bash
mkdir -p /data/db/hhblits && cd /data/db/hhblits
# BFD（约 1.8 TB 解包后，磁盘紧张可用更小编制）
wget https://bfd.mmseqs.com/bfd_metaclust_clu_complete_id30_c90_final.sseadv.ffdata.tar.gz
# UniClust30（推荐日常用，约 60 GB）
wget https://wwwuser.gwdg.de/~compbiol/uniclust/2024_06/UniRef30_2024_06_hhsuite.tar.gz
tar xzf UniRef30_2024_06_hhsuite.tar.gz
```

### 6.4 CARD（耐药库）

```bash
mkdir -p /data/db/card && cd /data/db/card
# 从 CARD Download 页（https://card.mcmaster.ca/download）下载 latest.tar.bz2
# 解出 aro_index.json / card_database.fasta；RGI 库：
pip install rgi
rgi card_annotation -i <latest.tar.bz2 解压后的 json>
rgi load --card_json card_json --local
```

### 6.5 VFDB（毒力因子库）

```bash
mkdir -p /data/db/vfdb && cd /data/db/vfdb
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz   # 参考蛋白序列
gunzip VFDB_setA_pro.fas.gz
```

### 6.6 SPI（沙门氏菌致病岛参考，可选）

从 SPI 数据库 / 课题组维护的 fasta 目录准备，context 阶段用 `--reference-fna` 传入对应基因组 fna 目录。

---

## 7. 磁盘与显存需求（约略估计）

| 资源 | 规模 | 说明 |
| --- | --- | --- |
| MMseqs UniRef90 库 | ~150 GB | 含索引 |
| Foldseek AlphaFold DB | ~1.2 TB | 全库；冒烟只用 SwissProt ~2 GB |
| hhblits UniRef30 库 | ~60 GB | BFD 需 ~1.8 TB |
| CARD + VFDB | ~2 GB | |
| ESM-2 650M 权重 | ~2.5 GB | ESMFold 权重另需 ~3 GB |
| 工作目录（数十株菌） | ~10–30 GB | 嵌入向量 + 中间结果 |
| **显存：ESM-2 650M 推理** | ~8 GB | batch 8；OOM 就把 batch 调小 |
| **显存：ESMFold** | ~16 GB 起步 | 1000 aa 以内蛋白建议 24 GB；否则自动跳过折叠 |

---

## 8. 数十株沙门氏菌完整命令串

```bash
WORKDIR=/data/runs/salmonella_$(date +%Y%m%d)
DB=/data/db
THREADS=32

# 0) 准备：所有样本 GFF3/FAA 放在 samples/sampleA/ 下
pufprobe --workdir $WORKDIR --threads $THREADS prepare \
    --input /data/genomes/salmonella/

# 1) 嵌入 + 聚类（GPU）
pufprobe --workdir $WORKDIR --threads $THREADS embed

# 2) 折叠（GPU；短序列自动跑 ESMFold）
pufprobe --workdir $WORKDIR --threads $THREADS fold

# 3) 检索
pufprobe --workdir $WORKDIR --threads $THREADS search \
    --channels mmseqs,foldseek,hhblits --dbdir $DB

# 4) 上下文（邻域 / 系统发育谱 / 微共线性）
pufprobe --workdir $WORKDIR --threads $THREADS context \
    --card $DB/card --vfdb $DB/vfdb/VFDB_setA_pro.fas \
    --reference-fna /data/genomes/salmonella/

# 5) 证据融合排序（可选 --weights fusion.joblib 用学习权重）
pufprobe --workdir $WORKDIR --threads $THREADS rank

# 6) 报告
pufprobe --workdir $WORKDIR --threads $THREADS report --formats json,html,md

# 7) 注释回填
pufprobe --workdir $WORKDIR --threads $THREADS annotate
```

预计耗时（32 核 + 单卡 A100 量级，约 30–50 株）：

| 阶段 | 耗时 |
| --- | --- |
| prepare | <5 min |
| embed（ESM-2 650M） | 2–4 h（GPU）；CPU 回退则 ~1–2 天 |
| fold（ESMFold） | 每代表蛋白 ~1–3 min，视代表数而定 |
| search（全库） | 3–8 h |
| context / rank / report / annotate | <30 min 合计 |

断点续跑：每个阶段写 `_DONE` 标记，中断后重跑同一命令自动跳过已完成阶段；强制重跑加 `--force`。

---

## 9. 常见报错排查

- **CUDA out of memory**：先把 embed 的 batch 调小（config 里 `embed.batch_size` 改成 2）；或对 >1000 aa 的蛋白 ESMFold 会自动跳过；实在不行用 `pufprobe --cpu ... embed --fallback-only` 走特征回退。
- **mmseqs database path errors**：确认 `--dbdir` 指向的目录下用的是 `mmseqs createdb` 产出的库前缀（不是原始 fasta）；`pufprobe doctor` 里 `mmseqs available=True` 只是找到二进制，库路径要自己给对。
- **foldseek 缺结构**：05_search 阶段日志提示 "no usable structures" —— 说明 04_fold 里该家族 pLDDT 低于阈值或被 skip，结构检索通道自动关闭，不影响序列通道。
- **GFF 中 locus 与 FAA 对不上**：annotate 没有改到某条 CDS，多半是 GFF3 里 `locus_tag` 与 FAA header 首标记不一致；用 `grep` 对一下两个文件里的 id，prepare 阶段已尽量按 locus_tag / protein_id 关联。
- **Windows 生成文件的 CRLF**：本流程产出的 GFF3/TSV/JSON 都以 LF 写出（`newline="\n"`），Windows 编辑器显示与否不影响 Linux 端解析；跨平台运行无需转换。
