#!/usr/bin/env bash

set -euo pipefail

BAM_DIR=/path/to/bam
PAT_DIR=/path/to/output/pat
RUN_DIR=/path/to/output/hameth_v2_run
REGION_DIR=/path/to/region_files
CPG_DIR=/path/to/CpG_reference
META=/path/to/meta.tsv
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

hameth_v2 prepare \
    --bam-dir "${BAM_DIR}" \
    --pat-dir "${PAT_DIR}" \
    --workers 30 \
    --bam-batch-num 3 \
    --genome hg38 \
    --chromosomes "${CHROMS[@]}"

hameth_v2 hamdist \
    --pat-dir "${PAT_DIR}" \
    --output-dir "${RUN_DIR}" \
    --batch-size 50 \
    --batch-num 40 \
    --workers 30 \
    --top-n 20000 \
    --region-dir "${REGION_DIR}" \
    --cpg-dir "${CPG_DIR}" \
    --genome hg38 \
    --chromosomes "${CHROMS[@]}" \
    --bin-size 2000 \
    --rarefaction-k 100 \
    --repeats 20

hameth_v2 matrix \
    --meta "${META}" \
    --pat-dir "${PAT_DIR}" \
    --cpg-dir "${CPG_DIR}" \
    --bins "${RUN_DIR}/balanced_selection/K100/SD_rarefied_pat_hamdist_top20000_ranked.tsv" \
    --output-dir "${RUN_DIR}/fixed20_matrix" \
    --cpgs-per-bin 20 \
    --min-observed-cells 2 \
    --workers 30
