# hameth_v2

`hameth_v2` 是 HamDis 项目中基于 read-pattern Hamdist 的深度平衡版本。它把
BAM 转 PAT、无监督 batch 计算、等信息量 bin 筛选和等 CpG 权重的 cell x cell
矩阵整理为可安装的命令行软件。

与旧方法相比，`hameth_v2` 主要解决两个问题：

1. 选 bin 时，不再让低覆盖 bin 因少数比较的随机波动而天然获得较高 SD。
2. 构建细胞距离矩阵时，不再让 CpG-rich bin 因包含更多 CpGs 而天然拥有更高权重。

## 安装

推荐在已有的 `new_libis` 环境中安装：

```bash
conda activate new_libis
cd /path/to/HamDis/hameth_v2
pip install -e . --no-deps
hameth_v2 --version
```

Python 依赖：

- Python >= 3.7
- numpy
- pandas
- pysam
- tqdm

外部程序：

- `wgbstools`：`prepare` 的 BAM 转 PAT 步骤需要。
- `samtools`：缺少 BAM index 时自动建立 `.bai`。

## 三步完整流程

### 1. BAM 转 PAT

```bash
hameth_v2 prepare \
  --bam-dir /path/to/bam \
  --pat-dir /path/to/pat \
  --workers 30 \
  --bam-batch-num 3 \
  --genome hg38 \
  --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
                chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
                chr20 chr21 chr22 chrX chrY
```

`--bam-batch-num 3` 表示把 BAM 列表分成 3 组并发处理，而不是依次处理 3 组。
总 worker 数仍受 `--workers` 控制。默认还会在 PAT 根目录生成每条染色体的
`chr*_allcells_combined.gz`；后续步骤不依赖这些合并文件，可用
`--skip-merge-allcells` 节省空间。

输出结构：

```text
pat/
  chr1/sample_001.pat.gz
  chr1/sample_001.pat.gz.csi
  chr2/sample_001.pat.gz
  ...
```

### 2. 无监督选 bin

```bash
hameth_v2 hamdist \
  --pat-dir /path/to/pat \
  --output-dir /path/to/hameth_v2_run \
  --batch-size 50 \
  --batch-num 40 \
  --workers 30 \
  --top-n 20000 \
  --region-dir /path/to/region_files \
  --cpg-dir /path/to/CpG_reference \
  --genome hg38 \
  --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
                chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
                chr20 chr21 chr22 chrX chrY \
  --bin-size 2000 \
  --batch-seed 20261008 \
  --selection-seed 20261008 \
  --rarefaction-k 100 \
  --repeats 20
```

这条命令自动完成两部分：

1. 按旧版经过验证的 read-pattern Hamdist 规则生成 40 个无监督 batch 的 merged
   PAT 和每个候选 bin 的 metric。
2. 对这些 batch-bin 做 K=100 等信息量稀释、20 次重复和 Q10 稳定性排序。

默认要求至少 90% 的 batches 有效；`bn=40` 时即 36 批。可通过
`--min-valid-batches` 显式修改。由于当前 batch 内的细胞状态使用 `uint64` mask，
`--batch-size` 必须小于等于 64。

主要结果：

```text
hameth_v2_run/
  middle_res/metrics/                         # 原始 batch-bin metrics
  middle_res/mergepats/                       # 每批合并 PAT
  balanced_selection/K100/
    SD_rarefied_pat_hamdist_top20000_ranked.tsv
    SD_rarefied_pat_hamdist_top20000.bed
    rarefied_pat_hamdist_all_eligible.tsv.gz
    score_depth_correlations.tsv
    run_parameters.json
```

若已经有旧版 Hammeth/HamDis 的完整 `middle_res/metrics` 和
`middle_res/mergepats`，无需重跑 batch 计算，可直接运行：

```bash
hameth_v2 select-bins \
  --run-dir /path/to/existing_run \
  --cpg-dir /path/to/CpG_reference \
  --output-dir /path/to/existing_run/balanced_selection/K100 \
  --batch-count 40 \
  --rarefaction-k 100 \
  --min-valid-batches 36 \
  --repeats 20 \
  --top-n 20000 \
  --selection-seed 20261008 \
  --workers 20
```

### 3. 构建 cell x cell Hamdist 矩阵

元数据至少需要 `cell_id` 列；可选的 `MajorType` 仅用于保存注释，不参与选
CpG 或计算距离。如果 PAT 文件名与 `cell_id` 不一致，可增加 `pat_stem` 列。

```bash
hameth_v2 matrix \
  --meta /path/to/meta.tsv \
  --pat-dir /path/to/pat \
  --cpg-dir /path/to/CpG_reference \
  --bins /path/to/hameth_v2_run/balanced_selection/K100/SD_rarefied_pat_hamdist_top20000_ranked.tsv \
  --output-dir /path/to/hameth_v2_run/fixed20_matrix \
  --cpgs-per-bin 20 \
  --min-observed-cells 2 \
  --workers 30
```

每个入选 bin 固定选择 20 个 CpGs。若某个 bin 无法提供足够 CpGs，则明确丢弃并
记录原因。两个细胞仅在双方共同观测的 selected CpGs 上计算：

```text
Hamdist = 状态不一致的 CpG 数 / 双方共同观测的 CpG 数
```

主要输出：

- `all_cells_hamdist_matrix.npz`
- `all_cells_hamdist_matrix.tsv.gz`
- `all_cells_shared_cpg_count.tsv.gz`
- `fixed_cpg_bin_support.tsv`
- `matrix_target_cpg_index.tsv.gz`
- `cell_manifest.tsv`
- `matrix_summary.tsv`
- `run_parameters.json`

只需要 NPZ 时可增加 `--skip-tsv`。若确认覆盖已有结果，增加 `--force`。

## 参考目录格式

`--cpg-dir` 中每条染色体对应一个 `CpG_chr*.bed`，至少包含位置和连续 CpG
编号。`--region-dir` 中包含旧 batch metric 生成步骤所需的
`region_chr*.bed`。

本流程沿用现有 HamDis reference 的坐标约定：区间和 CpG position 按 1-based、
两端包含进行匹配。输出虽然使用 `.bed` 后缀，但在与标准 0-based half-open BED
工具组合前应先转换坐标。

## 方法详情

稀释统计、Q10、fixed-CpG 矩阵以及缓存规则详见 [METHOD_CN.md](METHOD_CN.md)。
完整三步模板见 [examples/run_full_pipeline.sh](examples/run_full_pipeline.sh)。
