#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=100G
set -euo pipefail

# 1. 提取 VMR 内 CpGs
mkdir -p tmp

echo "[Step 1] Extract CpGs in VMRs ..."
for f in *.bed; do
    sample=$(basename $f .bed)
    if [[ ! -s tmp/${sample}_VMR.bed ]]; then
        bedtools intersect -a "$f" -b VMRs_5%_chr1-22_3col_fixed.txt -wa > tmp/${sample}_VMR.bed
    fi
done

echo "[Step 2] Computing Hamming distance with memory optimization..."
python <<'PYCODE'

import os
import glob
import numpy as np
import pandas as pd
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# ----------------------------
# 参数设置
# ----------------------------
tmp_dir = "tmp"
out_long = "hamming_long_v4.tsv"
out_matrix = "hamming_distance_matrix_v4.tsv"
nproc = 16
batch_size = 5000  # 减少批量大小

# 使用内存映射文件
memmap_file = "X_memmap.dat"

# ----------------------------
# 1. 构建全局 CpG 列索引
# ----------------------------
print("Step 1: Collect all CpG positions...")
all_positions = set()
file_list = sorted(glob.glob(os.path.join(tmp_dir, "*_VMR.bed")))
sample_names = [os.path.basename(f).replace("_VMR.bed","") for f in file_list]

for f in file_list:
    with open(f) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom, pos = parts[0], int(parts[1])
            all_positions.add((chrom,pos))

all_positions = sorted(all_positions)
pos2idx = {pos:i for i,pos in enumerate(all_positions)}
n_cpgs = len(all_positions)
n_cells = len(sample_names)
print(f"Total CpG positions: {n_cpgs}")

# ----------------------------
# 2. 构建 memory-mapped 矩阵
# ----------------------------
print("Step 2: Build memory-mapped matrix...")
# 使用 uint8 类型节省空间 (0: 未甲基化, 1: 甲基化, 255: 缺失值)
X = np.memmap(memmap_file, dtype=np.uint8, mode='w+', shape=(n_cells, n_cpgs))
X[:] = 255  # 初始化为缺失值

for i, f in enumerate(tqdm(file_list)):
    with open(f) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom, pos, ratio = parts[0], int(parts[1]), float(parts[3])
            idx = pos2idx[(chrom,pos)]
            # 将甲基化比率转换为二进制值
            binary_val = 1 if ratio > 0.5 else 0
            X[i, idx] = binary_val

# 确保数据写入磁盘
X.flush()
del X  # 释放内存映射

# ----------------------------
# 3. 定义 Hamming distance 函数
# ----------------------------
def hamming_row(idx_pair):
    # 重新打开内存映射文件为只读模式
    X_readonly = np.memmap(memmap_file, dtype=np.uint8, mode='r', shape=(n_cells, n_cpgs))
    i, j = idx_pair
    v1 = X_readonly[i]
    v2 = X_readonly[j]
    
    # 找出两个样本都有数据的位置
    valid_mask = (v1 != 255) & (v2 != 255)
    if np.sum(valid_mask) == 0:
        return sample_names[i], sample_names[j], np.nan
    
    # 计算不同值的比例
    dist = np.mean(v1[valid_mask] != v2[valid_mask])
    return sample_names[i], sample_names[j], dist

# ----------------------------
# 4. 分批并行计算 Hamming distance
# ----------------------------
print("Step 3: Compute Hamming distance...")
pairs = list(combinations(range(n_cells),2))
total_batches = (len(pairs) + batch_size - 1) // batch_size

with open(out_long, "w") as out:
    out.write("sample1\tsample2\thamming\n")

    # 分批处理
    for batch_num in range(0, len(pairs), batch_size):
        batch_pairs = pairs[batch_num:batch_num+batch_size]
        current_batch = batch_num // batch_size + 1
        print(f"Processing batch {current_batch}/{total_batches} ({len(batch_pairs)} pairs)")

        with ProcessPoolExecutor(max_workers=nproc) as executor:
            futures = {executor.submit(hamming_row, p): p for p in batch_pairs}
            for fut in tqdm(as_completed(futures), total=len(batch_pairs), desc=f"Batch {current_batch}"):
                s1, s2, dist = fut.result()
                out.write(f"{s1}\t{s2}\t{dist}\n")
                out.flush()

# ----------------------------
# 5. 转换为对称矩阵
# ----------------------------
print("Step 4: Convert to symmetric matrix...")
df = pd.read_csv(out_long, sep="\t")
matrix = df.pivot(index="sample1", columns="sample2", values="hamming")
matrix = matrix.combine_first(matrix.T)
for s in matrix.index:
    matrix.loc[s, s] = 0
matrix.to_csv(out_matrix, sep="\t")
print("Done! Hamming distance matrix saved to:", out_matrix)

# 清理临时文件
if os.path.exists(memmap_file):
    os.remove(memmap_file)
PYCODE
