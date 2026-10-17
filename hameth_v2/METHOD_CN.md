# hameth_v2 方法说明

## 1. 为什么需要深度平衡

原始流程把基因组切成固定长度 bin，然后在多个随机细胞 batch 中计算 Hamdist。
如果一个 bin 只有很少的有效 CpG/细胞配对，它的 batch 间数值容易产生很大的随机
波动；直接按 SD 排序时，这些低信息量 bin 会被优先选中。相反，CpG 多、覆盖深的
bin 因大量比较被平均，SD 往往更小。这会把测序深度与 CpG 密度混入生物差异。

`hameth_v2` 分别在选 bin 和构建矩阵时控制信息量：

- 选 bin：每个 batch-bin 固定抽取 K 个有效比较。
- 构建矩阵：每个保留 bin 固定贡献相同数量的 CpGs。

## 2. 候选 bin 与无监督 batch

用户可通过 `--bin-size` 指定固定长度候选区间，例如 2 kb 或 5 kb。随后从全部
单细胞中随机抽取 `batch-size` 个细胞，重复 `batch-num` 次。细胞类型标签不参与
抽样、选 bin 或打分。

每个 batch 的单细胞 PAT 会按染色体合并。PAT pattern 只有在同一候选 bin 内至少
包含两个有效 `C/T` CpGs 时，才参与该 bin 的支持统计。

## 3. 每个 batch-bin 的充分统计量

对某个 CpG，设观测到有效甲基化状态的细胞数为 `n`，其中甲基化和非甲基化
细胞数分别为 `m` 和 `u`。同一细胞在同一 CpG 同时出现 C 与 T 时视为冲突并排除。

在一个 batch-bin 内累加：

```text
N = sum[n * (n - 1) / 2]
M = sum[m * u]
```

`N` 是可比较的 cell-pair/CpG 总数，`M` 是其中状态不一致的比较数。仍沿用旧流程
的有效性门槛，默认要求该 batch-bin 的 `HammingCount >= 5`。

## 4. K 次等信息量稀释

只有 `N >= K` 的 batch-bin 才能进入下一步。默认 `K=100`，从 N 个比较中无放回
抽取固定 K 个：

```text
X ~ Hypergeometric(M, N-M, K)
Z = 2 * asin(sqrt((X + 0.375) / (K + 0.75)))
```

固定 K 后，不同 bin 和不同 batch 的打分都基于相同数量的有效比较。角度变换 Z
用于减弱二项比例方差随均值改变的问题。

## 5. 跨 batch SD、重复稀释与 Q10

对同一个 bin，在所有有效 batches 的 Z 值之间计算样本 SD。默认要求至少 90%
batches 有效，例如 40 批至少 36 批。

由于超几何稀释本身有随机性，流程默认独立重复 20 次，每次得到一个跨 batch SD。
最终排序分数不是最大值或均值，而是 20 个 SD 的 10% 分位数 Q10。

Q10 是保守稳定性分数：只有在绝大多数稀释重复中都保持较高 SD 的 bin 才能获得
高分；偶然一次抽到极端组成不会把一个 bin 推到前列。随机种子固定后结果可复现。

## 6. 缓存

扫描 merged PAT 生成 M、N 的过程最耗时，结果按染色体保存在
`sufficient_statistics_cache/`。仅修改 K、重复次数、随机种子或 Top-N 时会自动
复用缓存。只有 merged PAT、metrics、CpG reference 或支持统计规则变化时才使用
`--force-cache`。

## 7. 固定 CpG 数构建细胞矩阵

对 ranked bins 中每个候选 bin：

1. 从 indexed single-cell PAT 统计每个 CpG 的 C_count 和 T_count。
2. `C_count / (C_count + T_count) > 0.5` 记为甲基化；等于 0.5 时记为非甲基化。
3. 去掉被少于 `--min-observed-cells` 个细胞观测到的 CpGs。
4. 按观测细胞数从高到低选择固定 `--cpgs-per-bin` 个 CpGs；并列时按 CpG ID。
5. 无法选满固定数量 CpGs 的 bin 被丢弃并写入 QC 表。

对细胞 i、j，仅在双方共同观测的 selected CpGs 上计算：

```text
Hamdist(i,j) = mismatch_count(i,j) / shared_observed_cpg_count(i,j)
```

因此每个保留 bin 提供同样数量的候选特征，CpG-rich bin 不会仅因绝对 CpG 数更多
而被赋予更高的预设权重。同时输出 pairwise shared-CpG count，便于识别证据不足的
细胞对。

## 8. 两个阶段不要混淆

选 bin 阶段保留 read/PAT pattern 层面的规则，用于寻找 batch 间稳定变化区域；
矩阵阶段从单细胞 PAT 累积每个 CpG 的状态，并控制每个 bin 的 CpG 数，用于比较
细胞。两阶段目标不同，因此统计汇总方式也不同，但共同核心仍是甲基化状态之间的
Hamdist。
