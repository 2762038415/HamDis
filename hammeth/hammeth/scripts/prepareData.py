#!/usr/bin/env python3

import os
import glob
import argparse
import subprocess
import multiprocessing
import gzip
import pandas as pd


############################################
# 工具函数
############################################

def check_bam_index(bam_file):
    """Ensure BAM index exists."""
    bai = bam_file + ".bai"

    if not os.path.exists(bai):
        print(f"[INFO] Creating BAM index for {bam_file}")

        cmd = ["samtools", "index", bam_file]
        res = subprocess.run(cmd, capture_output=True, text=True)

        if res.returncode != 0:
            raise RuntimeError(f"samtools index failed:\n{res.stderr}")


############################################
# BAM → PAT
############################################

def bam_to_pat_worker(args):

    bam_file, pat_dir, chr_list, region_dir, genome = args

    bam_name = os.path.basename(bam_file)
    bam_prefix = os.path.splitext(bam_name)[0]

    check_bam_index(bam_file)

    for chrom in chr_list:

        out_dir = os.path.join(pat_dir, chrom)
        os.makedirs(out_dir, exist_ok=True)

        pat_file = os.path.join(out_dir, f"{bam_prefix}.pat.gz")

        if os.path.exists(pat_file):
            print(f"[SKIP] {pat_file} exists")
            continue

        region_file = os.path.join(region_dir, f"region_{chrom}.bed")

        if not os.path.exists(region_file):
            raise FileNotFoundError(region_file)

        cmd = [
            "wgbstools",
            "bam2pat",
            "--genome", genome,
            "-f",
            "-r", chrom,
            bam_file,
            "-o", out_dir,
            "--no_beta"
        ]

        print("[RUN]", " ".join(cmd))

        res = subprocess.run(cmd)

        if res.returncode != 0:
            raise RuntimeError("bam2pat failed")


############################################
# merge PAT
############################################

def merge_gz_by_chr(pat_dir, chr_list):

    print("\n[INFO] Start merging PAT files")

    for chrom in chr_list:

        chr_dir = os.path.join(pat_dir, chrom)

        if not os.path.isdir(chr_dir):
            print(f"[WARN] {chr_dir} not found")
            continue

        gz_files = glob.glob(os.path.join(chr_dir, "*.gz"))

        if not gz_files:
            print(f"[WARN] no gz files in {chr_dir}")
            continue

        out_file = os.path.join(pat_dir, f"{chrom}_allcells_combined.gz")

        if os.path.exists(out_file):
            print(f"[SKIP] {out_file} exists")
            continue

        print(f"[MERGE] {chrom} ({len(gz_files)} files)")

        with gzip.open(out_file, "wt") as outfile:

            for f in gz_files:

                cell = os.path.basename(f).split(".")[0]

                with gzip.open(f, "rt") as infile:

                    for line in infile:

                        line = line.rstrip()

                        if not line:
                            continue

                        cols = line.split()

                        if len(cols) == 4:
                            cols.append(cell)
                        else:
                            cols[4] = cell

                        outfile.write("\t".join(cols) + "\n")


############################################
# 主程序
############################################

def main():

    parser = argparse.ArgumentParser(
        description="Convert BAM to PAT and merge by chromosome"
    )

    parser.add_argument(
        "-i",
        dest="bam_dir",
        required=True,
        help="input BAM directory"
    )

    parser.add_argument(
        "-o",
        dest="pat_dir",
        default=None,
        help="output PAT directory"
    )

    parser.add_argument(
        "-m",
        dest="threads",
        type=int,
        default=8,
        help="number of processes"
    )

    parser.add_argument(
        "-cl",
        dest="chr_list",
        nargs="+",
        required=True,
        help="chromosome list"
    )

    parser.add_argument(
        "-r",
        dest="region_dir",
        required=True,
        help="region directory"
    )

    parser.add_argument(
        "-cp",
        dest="cpg_dir",
        required=True,
        help="CpG directory (reserved)"
    )

    parser.add_argument(
        "-g",
        dest="genome",
        required=True,
        help="reference genome (hg38/mm10)"
    )

    args = parser.parse_args()

    bam_dir = args.bam_dir
    chr_list = args.chr_list
    region_dir = args.region_dir
    genome = args.genome
    threads = args.threads

    if args.pat_dir is None:
        pat_dir = os.path.join(os.path.dirname(bam_dir), "pat")
    else:
        pat_dir = args.pat_dir

    os.makedirs(pat_dir, exist_ok=True)

    bam_files = glob.glob(os.path.join(bam_dir, "*.bam"))

    if not bam_files:
        raise RuntimeError("No BAM files found")

    print(f"[INFO] BAM files: {len(bam_files)}")
    print(f"[INFO] Threads: {threads}")

    worker_args = [
        (bam, pat_dir, chr_list, region_dir, genome)
        for bam in bam_files
    ]

    with multiprocessing.Pool(threads) as pool:
        pool.map(bam_to_pat_worker, worker_args)

    merge_gz_by_chr(pat_dir, chr_list)

    print("\n[FINISHED]")


############################################

if __name__ == "__main__":
    main()
