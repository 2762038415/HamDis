#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import subprocess
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from tqdm import tqdm


def safe_makedirs(path):
    os.makedirs(path, exist_ok=True)


def run_cmd(cmd):
    print(f"[CMD] {cmd}")
    subprocess.run(cmd, shell=True, check=True)


# ----------------------------
# Step 1. *_dedup.bed -> *.with_ratio.bed
# 输出列：chr start end ratio meth total
# ----------------------------
def convert_bed_to_ratio(input_dir, pattern="*_dedup.bed", overwrite=False):
    bed_files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not bed_files:
        raise FileNotFoundError(f"在 {input_dir} 下没有找到 {pattern}")

    out_files = []
    print("[Step 1] Converting *_dedup.bed to *.with_ratio.bed ...")
    for f in bed_files:
        out = os.path.splitext(f)[0] + ".with_ratio.bed"

        if (not overwrite) and os.path.exists(out) and os.path.getsize(out) > 0:
            print(f"[Skip] exists: {out}")
            out_files.append(out)
            continue

        with open(f) as fin, open(out, "w") as fout:
            for line in fin:
                line = line.rstrip("\n")
                if not line:
                    continue

                parts = line.split("\t")
                if len(parts) < 5:
                    continue

                chrom = parts[0]
                start = parts[1]
                end = parts[2]
                meth = parts[3]
                total = parts[4]

                try:
                    meth_f = float(meth)
                    total_f = float(total)
                except ValueError:
                    continue

                ratio = "NA" if total_f == 0 else str(meth_f / total_f)
                fout.write(f"{chrom}\t{start}\t{end}\t{ratio}\t{meth}\t{total}\n")

        print(f"[Done] {os.path.basename(f)} -> {os.path.basename(out)}")
        out_files.append(out)

    return out_files


# ----------------------------
# Step 2. 提取目标区域内 CpGs
# bedtools intersect -wa
# ----------------------------
def extract_vmrs(ratio_bed_files, region_bed, tmp_dir, overwrite=False):
    safe_makedirs(tmp_dir)
    out_files = []

    print("[Step 2] Extract CpGs in target regions ...")
    for f in ratio_bed_files:
        sample = os.path.basename(f).replace(".with_ratio.bed", "")
        out = os.path.join(tmp_dir, f"{sample}_VMR.bed")

        if (not overwrite) and os.path.exists(out) and os.path.getsize(out) > 0:
            print(f"[Skip] exists: {out}")
            out_files.append(out)
            continue

        cmd = f'bedtools intersect -a "{f}" -b "{region_bed}" -wa > "{out}"'
        run_cmd(cmd)
        out_files.append(out)

    return out_files


# ----------------------------
# Step 3. 收集全局 CpG 列索引
# 用 (chrom, start)
# ----------------------------
def collect_all_positions(file_list):
    print("[Step 3] Collect all CpG positions ...")
    all_positions = set()

    for f in file_list:
        with open(f) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                chrom, pos = parts[0], int(parts[1])
                all_positions.add((chrom, pos))

    all_positions = sorted(all_positions)
    pos2idx = {pos: i for i, pos in enumerate(all_positions)}
    print(f"[Info] Total CpG positions: {len(all_positions)}")
    return all_positions, pos2idx


# ----------------------------
# Step 4. 构建 memmap 矩阵
# 255 = missing
# ratio > 0.5 -> 1
# else -> 0
# ----------------------------
def build_memmap_matrix(file_list, sample_names, pos2idx, memmap_file, dtype=np.uint8):
    n_cells = len(sample_names)
    n_cpgs = len(pos2idx)

    print("[Step 4] Build memory-mapped matrix ...")
    X = np.memmap(memmap_file, dtype=dtype, mode="w+", shape=(n_cells, n_cpgs))
    X[:] = 255

    for i, f in enumerate(tqdm(file_list, desc="Loading VMR beds")):
        with open(f) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue

                chrom = parts[0]
                pos = int(parts[1])
                ratio_raw = parts[3]

                if ratio_raw == "NA":
                    continue

                try:
                    ratio = float(ratio_raw)
                except ValueError:
                    continue

                idx = pos2idx[(chrom, pos)]
                binary_val = 1 if ratio > 0.5 else 0
                X[i, idx] = binary_val

    X.flush()
    del X
    return n_cells, n_cpgs


# ----------------------------
# Step 5. Hamming distance
# ----------------------------
def _hamming_row(idx_pair, memmap_file, n_cells, n_cpgs, sample_names):
    X_readonly = np.memmap(memmap_file, dtype=np.uint8, mode="r", shape=(n_cells, n_cpgs))
    i, j = idx_pair
    v1 = X_readonly[i]
    v2 = X_readonly[j]

    valid_mask = (v1 != 255) & (v2 != 255)
    if np.sum(valid_mask) == 0:
        return sample_names[i], sample_names[j], np.nan

    dist = np.mean(v1[valid_mask] != v2[valid_mask])
    return sample_names[i], sample_names[j], dist


def compute_hamming(memmap_file, n_cells, n_cpgs, sample_names, out_long, out_matrix,
                    nproc=20, batch_size=5000):
    print("[Step 5] Compute Hamming distance ...")
    pairs = list(combinations(range(n_cells), 2))
    total_batches = (len(pairs) + batch_size - 1) // batch_size

    with open(out_long, "w") as out:
        out.write("sample1\tsample2\thamming\n")

        for batch_num in range(0, len(pairs), batch_size):
            batch_pairs = pairs[batch_num:batch_num + batch_size]
            current_batch = batch_num // batch_size + 1
            print(f"[Info] Processing batch {current_batch}/{total_batches} ({len(batch_pairs)} pairs)")

            with ProcessPoolExecutor(max_workers=nproc) as executor:
                futures = {
                    executor.submit(
                        _hamming_row,
                        p,
                        memmap_file,
                        n_cells,
                        n_cpgs,
                        sample_names,
                    ): p for p in batch_pairs
                }

                for fut in tqdm(as_completed(futures), total=len(batch_pairs), desc=f"Batch {current_batch}"):
                    s1, s2, dist = fut.result()
                    out.write(f"{s1}\t{s2}\t{dist}\n")
                    out.flush()

    print("[Step 6] Convert to symmetric matrix ...")
    df = pd.read_csv(out_long, sep="\t")
    matrix = df.pivot(index="sample1", columns="sample2", values="hamming")
    matrix = matrix.combine_first(matrix.T)

    all_samples = sorted(set(sample_names))
    matrix = matrix.reindex(index=all_samples, columns=all_samples)
    for s in all_samples:
        matrix.loc[s, s] = 0

    matrix.to_csv(out_matrix, sep="\t")
    print(f"[Done] Hamming distance matrix saved to: {out_matrix}")


# ----------------------------
# 主入口，供 commands/matrix.py 调用
# ----------------------------
def run_matrix_pipeline(
    input_dir,
    region_bed,
    tmp_dir="tmp",
    out_long="tutorial_hamming_long_v4.tsv",
    out_matrix="tutorial_hamming_distance_matrix_v4.tsv",
    memmap_file="X_memmap.dat",
    nproc=20,
    batch_size=5000,
    overwrite=False,
    keep_tmp=False,
):
    input_dir = os.path.abspath(input_dir)
    region_bed = os.path.abspath(region_bed)
    tmp_dir = os.path.abspath(tmp_dir)
    out_long = os.path.abspath(out_long)
    out_matrix = os.path.abspath(out_matrix)
    memmap_file = os.path.abspath(memmap_file)

    if not os.path.exists(region_bed):
        raise FileNotFoundError(f"region bed 不存在: {region_bed}")

    ratio_bed_files = convert_bed_to_ratio(
        input_dir=input_dir,
        pattern="*_dedup.bed",
        overwrite=overwrite
    )

    extract_vmrs(
        ratio_bed_files=ratio_bed_files,
        region_bed=region_bed,
        tmp_dir=tmp_dir,
        overwrite=overwrite
    )

    file_list = sorted(glob.glob(os.path.join(tmp_dir, "*_VMR.bed")))
    if not file_list:
        raise RuntimeError(f"{tmp_dir} 下没有生成 *_VMR.bed")

    sample_names = [os.path.basename(f).replace("_VMR.bed", "") for f in file_list]

    _, pos2idx = collect_all_positions(file_list)
    if len(pos2idx) == 0:
        raise RuntimeError("目标区域内没有任何可用 CpG 位点，无法构建矩阵")

    n_cells, n_cpgs = build_memmap_matrix(
        file_list=file_list,
        sample_names=sample_names,
        pos2idx=pos2idx,
        memmap_file=memmap_file
    )

    compute_hamming(
        memmap_file=memmap_file,
        n_cells=n_cells,
        n_cpgs=n_cpgs,
        sample_names=sample_names,
        out_long=out_long,
        out_matrix=out_matrix,
        nproc=nproc,
        batch_size=batch_size
    )

    if not keep_tmp:
        if os.path.exists(memmap_file):
            os.remove(memmap_file)
            print(f"[Clean] Removed {memmap_file}")

    print("\nAll done.")
    print(f"with_ratio beds: {input_dir}/*.with_ratio.bed")
    print(f"VMR beds       : {tmp_dir}/*_VMR.bed")
    print(f"long table     : {out_long}")
    print(f"matrix         : {out_matrix}")
