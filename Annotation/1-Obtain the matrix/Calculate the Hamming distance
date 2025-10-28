#stpe1. Extract CpGs from the VMRs
mkdir -p tmp

echo "[Step 1] Extract CpGs in VMRs ..."
for f in *.bed; do
    sample=$(basename $f .bed)
    if [[ ! -s tmp/${sample}_VMR.bed ]]; then
        bedtools intersect -a "$f" -b VMRs.txt -wa > tmp/${sample}_VMR.bed
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
# parameter setting
# ----------------------------
tmp_dir = "tmp"
out_long = "hamming_long_v4.tsv"
out_matrix = "hamming_distance_matrix_v4.tsv"
nproc = 16
batch_size = 5000  # Control the batch size

# Use memory-mapped files
memmap_file = "X_memmap.dat"

# ----------------------------
# 1.Build a global CpGs column index
# ----------------------------
print("Step 1: Collect all CpGs positions...")
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
# 2. Build a memory-mapped matrix
# ----------------------------
print("Step 2: Build memory-mapped matrix...")
# Using the uint8 type saves space (0: unmethylated, 1: methylated, 255: missing value)
X = np.memmap(memmap_file, dtype=np.uint8, mode='w+', shape=(n_cells, n_cpgs))
X[:] = 255  # Initialize as missing values

for i, f in enumerate(tqdm(file_list)):
    with open(f) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom, pos, ratio = parts[0], int(parts[1]), float(parts[3])
            idx = pos2idx[(chrom,pos)]
            # Convert the methylation ratio to a binary value and initialize it as a missing value
            binary_val = 1 if ratio > 0.5 else 0
            X[i, idx] = binary_val

# Ensure that the data is written to the disk.
X.flush()
del X  # Release memory mapping

# ----------------------------
# 3. Define the Hamming distance function
# ----------------------------
def hamming_row(idx_pair):
    # Reopen the memory-mapped file in read-only mode
    X_readonly = np.memmap(memmap_file, dtype=np.uint8, mode='r', shape=(n_cells, n_cpgs))
    i, j = idx_pair
    v1 = X_readonly[i]
    v2 = X_readonly[j]
    
    # Identify the positions where data is available for both samples
    valid_mask = (v1 != 255) & (v2 != 255)
    if np.sum(valid_mask) == 0:
        return sample_names[i], sample_names[j], np.nan
    
    # Calculate the proportion of different values
    dist = np.mean(v1[valid_mask] != v2[valid_mask])
    return sample_names[i], sample_names[j], dist

# ----------------------------
# 4. Batch parallel computation of Hamming distance
# ----------------------------
print("Step 3: Compute Hamming distance...")
pairs = list(combinations(range(n_cells),2))
total_batches = (len(pairs) + batch_size - 1) // batch_size

with open(out_long, "w") as out:
    out.write("sample1\tsample2\thamming\n")

    # batch handling
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
# 5. Convert to a symmetric matrix
# ----------------------------
print("Step 4: Convert to symmetric matrix...")
df = pd.read_csv(out_long, sep="\t")
matrix = df.pivot(index="sample1", columns="sample2", values="hamming")
matrix = matrix.combine_first(matrix.T)
for s in matrix.index:
    matrix.loc[s, s] = 0
matrix.to_csv(out_matrix, sep="\t")
print("Done! Hamming distance matrix saved to:", out_matrix)

# Delete temporary files
if os.path.exists(memmap_file):
    os.remove(memmap_file)
PYCODE
