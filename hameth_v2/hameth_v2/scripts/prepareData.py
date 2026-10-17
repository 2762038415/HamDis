#!/usr/bin/env python3

import argparse
import glob
import gzip
import multiprocessing
import os
import subprocess


def check_bam_index(bam_file):
    bai_candidates = [bam_file + ".bai", os.path.splitext(bam_file)[0] + ".bai"]
    if any(os.path.exists(path) for path in bai_candidates):
        return

    print("[INFO] Creating BAM index for {}".format(bam_file))
    cmd = ["samtools", "index", bam_file]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError("samtools index failed:\n{}".format(res.stderr))


def split_into_batches(items, batch_num):
    items = list(items)
    if not items:
        return []

    batch_num = max(1, int(batch_num))
    batch_num = min(batch_num, len(items))
    base_size, remainder = divmod(len(items), batch_num)

    batches = []
    start = 0
    for idx in range(batch_num):
        size = base_size + (1 if idx < remainder else 0)
        end = start + size
        if size > 0:
            batches.append(items[start:end])
        start = end
    return batches


def split_worker_counts(total_workers, group_num):
    total_workers = max(1, int(total_workers))
    group_num = max(1, int(group_num))
    base_size, remainder = divmod(total_workers, group_num)
    counts = []
    for idx in range(group_num):
        count = base_size + (1 if idx < remainder else 0)
        counts.append(max(1, count))
    return counts


def normalize_group_count(requested_group_num, item_count, total_workers):
    if item_count <= 0:
        return 0

    requested_group_num = max(1, int(requested_group_num))
    total_workers = max(1, int(total_workers))
    return min(requested_group_num, item_count, total_workers)


def bam_to_pat_worker(args):
    bam_file, pat_dir, chr_list, genome = args

    bam_name = os.path.basename(bam_file)
    bam_prefix = os.path.splitext(bam_name)[0]

    check_bam_index(bam_file)

    for chrom in chr_list:
        out_dir = os.path.join(pat_dir, chrom)
        os.makedirs(out_dir, exist_ok=True)

        pat_file = os.path.join(out_dir, "{}.pat.gz".format(bam_prefix))
        if os.path.exists(pat_file) and os.path.getsize(pat_file) > 0:
            print("[SKIP] {}".format(pat_file))
            continue

        cmd = [
            "wgbstools",
            "bam2pat",
            "--genome",
            genome,
            "-f",
            "-r",
            chrom,
            "-@",
            "1",
            bam_file,
            "-o",
            out_dir,
            "--no_beta",
        ]
        print("[RUN] {}".format(" ".join(cmd)))
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(
                "bam2pat failed for {} {}\nstdout:\n{}\nstderr:\n{}".format(
                    bam_file,
                    chrom,
                    (res.stdout or "").strip(),
                    (res.stderr or "").strip(),
                )
            )
        if (not os.path.exists(pat_file)) or os.path.getsize(pat_file) == 0:
            raise RuntimeError(
                "bam2pat produced no PAT output for {} {}\nexpected: {}\nstdout:\n{}\nstderr:\n{}\n"
                "note: wgbstools may print argument/reference errors while still exiting with status 0".format(
                    bam_file,
                    chrom,
                    pat_file,
                    (res.stdout or "").strip(),
                    (res.stderr or "").strip(),
                )
            )


def run_bam_batch_group(bam_batch, pat_dir, chr_list, genome, batch_threads, batch_idx, total_batches):
    print(
        "[INFO] Start BAM group {}/{}: {} BAMs, workers={}".format(
            batch_idx,
            total_batches,
            len(bam_batch),
            batch_threads,
        )
    )

    worker_args = [
        (bam, pat_dir, chr_list, genome)
        for bam in bam_batch
    ]

    with multiprocessing.Pool(batch_threads) as pool:
        pool.map(bam_to_pat_worker, worker_args)

    print(
        "[INFO] Finished BAM group {}/{}: {} BAMs".format(
            batch_idx,
            total_batches,
            len(bam_batch),
        )
    )


def process_bam_batches(bam_files, pat_dir, chr_list, genome, threads, bam_batch_num):
    actual_group_num = normalize_group_count(bam_batch_num, len(bam_files), threads)
    bam_batches = split_into_batches(sorted(bam_files), actual_group_num)
    total_batches = len(bam_batches)
    worker_counts = split_worker_counts(threads, total_batches)

    print("[INFO] Requested BAM groups: {}".format(bam_batch_num))
    print("[INFO] Concurrent BAM groups: {}, total workers: {}".format(total_batches, threads))
    print("[INFO] Worker split across BAM groups: {}".format(worker_counts))

    if total_batches == 1:
        run_bam_batch_group(
            bam_batch=bam_batches[0],
            pat_dir=pat_dir,
            chr_list=chr_list,
            genome=genome,
            batch_threads=worker_counts[0],
            batch_idx=1,
            total_batches=1,
        )
        return

    processes = []
    for batch_idx, (bam_batch, batch_threads) in enumerate(zip(bam_batches, worker_counts), 1):
        proc = multiprocessing.Process(
            target=run_bam_batch_group,
            args=(bam_batch, pat_dir, chr_list, genome, batch_threads, batch_idx, total_batches),
        )
        proc.start()
        processes.append((batch_idx, proc))

    failed_groups = []
    for batch_idx, proc in processes:
        proc.join()
        if proc.exitcode != 0:
            failed_groups.append((batch_idx, proc.exitcode))

    if failed_groups:
        raise RuntimeError("BAM groups failed: {}".format(failed_groups))


def merge_gz_by_chr(pat_dir, chr_list):
    print("\n[INFO] Start merging PAT files")

    for chrom in chr_list:
        chr_dir = os.path.join(pat_dir, chrom)
        if not os.path.isdir(chr_dir):
            print("[WARN] {} not found".format(chr_dir))
            continue

        gz_files = sorted(glob.glob(os.path.join(chr_dir, "*.pat.gz")))
        if not gz_files:
            print("[WARN] no gz files in {}".format(chr_dir))
            continue

        out_file = os.path.join(pat_dir, "{}_allcells_combined.gz".format(chrom))
        if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
            print("[SKIP] {}".format(out_file))
            continue

        print("[MERGE] {} ({} files)".format(chrom, len(gz_files)))

        with gzip.open(out_file, "wt") as outfile:
            for gz_file in gz_files:
                cell = os.path.basename(gz_file).split(".")[0]
                with gzip.open(gz_file, "rt") as infile:
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


def main():
    parser = argparse.ArgumentParser(description="Convert BAM to PAT and merge by chromosome")
    parser.add_argument("-i", dest="bam_dir", required=True, help="input BAM directory")
    parser.add_argument("-o", dest="pat_dir", default=None, help="output PAT directory")
    parser.add_argument("-m", dest="threads", type=int, default=8, help="number of processes")
    parser.add_argument("-cl", dest="chr_list", nargs="+", required=True, help="chromosome list")
    parser.add_argument("-r", dest="region_dir", default=None, help="deprecated; retained for compatibility")
    parser.add_argument("-cp", dest="cpg_dir", default=None, help="deprecated; retained for compatibility")
    parser.add_argument("-g", dest="genome", required=True, help="reference genome (hg38/mm10)")
    parser.add_argument(
        "--bam_batch_num",
        dest="bam_batch_num",
        type=int,
        default=1,
        help="split BAM files into N concurrent groups; 1 means a single group",
    )
    parser.add_argument(
        "--skip-merge-allcells",
        action="store_true",
        help="skip optional chromosome-wide combined PAT files",
    )

    args = parser.parse_args()

    bam_dir = args.bam_dir
    chr_list = args.chr_list
    genome = args.genome
    threads = max(1, int(args.threads))
    bam_batch_num = max(1, int(args.bam_batch_num))

    if args.pat_dir is None:
        pat_dir = os.path.join(os.path.dirname(bam_dir), "pat")
    else:
        pat_dir = args.pat_dir

    os.makedirs(pat_dir, exist_ok=True)

    bam_files = sorted(glob.glob(os.path.join(bam_dir, "*.bam")))
    if not bam_files:
        raise RuntimeError("No BAM files found")

    print("[INFO] BAM files: {}".format(len(bam_files)))
    print("[INFO] Threads: {}".format(threads))
    print("[INFO] BAM group parameter: {}".format(bam_batch_num))

    process_bam_batches(
        bam_files=bam_files,
        pat_dir=pat_dir,
        chr_list=chr_list,
        genome=genome,
        threads=threads,
        bam_batch_num=bam_batch_num,
    )

    if args.skip_merge_allcells:
        print("[INFO] Skipping optional all-cell PAT merge")
    else:
        merge_gz_by_chr(pat_dir, chr_list)
    print("\n[FINISHED]")


if __name__ == "__main__":
    main()
