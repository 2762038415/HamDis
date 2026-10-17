#!/usr/bin/env python3

"""Command-line interface for hameth_v2."""

from __future__ import print_function

import argparse
import math
import os
import subprocess
import sys
from pathlib import Path

from . import __version__


def implementation_script(name):
    return Path(__file__).resolve().parent / "scripts" / name


def run_script(name, arguments):
    command = [sys.executable, str(implementation_script(name))] + list(arguments)
    print("[COMMAND] {}".format(" ".join(command)), flush=True)
    subprocess.run(command, check=True)


def min_valid_batches(batch_count, requested=None):
    if requested is None:
        return int(math.ceil(0.9 * batch_count))
    return int(requested)


def selection_command(args, run_dir, output_dir, batch_count):
    required_batches = min_valid_batches(batch_count, args.min_valid_batches)
    if required_batches < 2 or required_batches > batch_count:
        raise ValueError("--min-valid-batches must be between 2 and batch count")
    command = [
        "--run-dir",
        os.path.abspath(run_dir),
        "--cpg-dir",
        os.path.abspath(args.cpg_dir),
        "--output-dir",
        os.path.abspath(output_dir),
        "--batch-count",
        str(batch_count),
        "--min-hamming-count",
        str(args.min_hamming_count),
        "--rarefaction-k",
        str(args.rarefaction_k),
        "--min-valid-batches",
        str(required_batches),
        "--repeats",
        str(args.repeats),
        "--top-n",
        str(args.top_n),
        "--seed",
        str(args.selection_seed),
        "--workers",
        str(args.workers),
    ]
    if args.chromosomes:
        command.extend(["--chromosomes"] + list(args.chromosomes))
    if args.force_cache:
        command.append("--force-cache")
    return command


def run_prepare(args):
    command = [
        "-i",
        args.bam_dir,
        "-o",
        args.pat_dir,
        "-m",
        str(args.workers),
        "-g",
        args.genome,
        "-cl",
    ] + list(args.chromosomes)
    command.extend(["--bam_batch_num", str(args.bam_batch_num)])
    if args.region_dir:
        command.extend(["-r", args.region_dir])
    if args.cpg_dir:
        command.extend(["-cp", args.cpg_dir])
    if args.skip_merge_allcells:
        command.append("--skip-merge-allcells")
    run_script("prepareData.py", command)


def run_select_bins(args):
    command = selection_command(
        args,
        run_dir=args.run_dir,
        output_dir=args.output_dir,
        batch_count=args.batch_count,
    )
    run_script("select_balanced_hamdist_bins.py", command)


def run_hamdist(args):
    if args.batch_size > 64:
        raise ValueError(
            "--batch-size must be <=64 because PAT support uses a uint64 cell mask"
        )
    if args.batch_size < 2 or args.batch_num < 2:
        raise ValueError("batch size and batch count must both be at least 2")
    if args.rarefaction_k < 2:
        raise ValueError("--rarefaction-k must be at least 2")

    from .scripts.legacy_hamdist_core import run_hamdist as run_legacy_hamdist

    run_dir = os.path.abspath(args.output_dir)
    selection_dir = args.selection_dir or os.path.join(
        run_dir, "balanced_selection", "K{}".format(args.rarefaction_k)
    )
    print("[1/2] Building unsupervised batch PAT metrics", flush=True)
    run_legacy_hamdist(
        pat_dir=args.pat_dir,
        output_dir=run_dir,
        batch_size=args.batch_size,
        batch_num=args.batch_num,
        m_process=args.workers,
        bin_num=args.top_n,
        region_dir=args.region_dir,
        cpg_file_dir=args.cpg_dir,
        genome=args.genome,
        chr_list=args.chromosomes,
        binsize=args.bin_size,
        cpgs_per_bin=args.candidate_cpgs_per_bin,
        seed=args.batch_seed,
        min_hm_count=args.min_hamming_count,
        min_cpg=args.min_pattern_cpg,
        max_dis=args.max_read_distance,
        max_possible_read=args.max_reads_per_bin,
        basedon0=args.basedon0,
        vmr="hmtop",
        vmr_bed=None,
        pat2bed_batch_num=1,
        pat2bed_workers=1,
        pat2beta_threads=1,
        overlap_pat2bed=False,
        skip_pat2bed=True,
    )

    print("[2/2] Selecting depth-balanced K/Q10 bins", flush=True)
    command = selection_command(
        args,
        run_dir=run_dir,
        output_dir=selection_dir,
        batch_count=args.batch_num,
    )
    run_script("select_balanced_hamdist_bins.py", command)
    print(
        "[RESULT] ranked bins: {}".format(
            os.path.join(
                os.path.abspath(selection_dir),
                "SD_rarefied_pat_hamdist_top{}_ranked.tsv".format(args.top_n),
            )
        ),
        flush=True,
    )


def run_matrix(args):
    command = [
        "--meta",
        args.meta,
        "--pat-dir",
        args.pat_dir,
        "--cpg-dir",
        args.cpg_dir,
        "--bins",
        args.bins,
        "--output-dir",
        args.output_dir,
        "--workers",
        str(args.workers),
        "--cpgs-per-bin",
        str(args.cpgs_per_bin),
        "--min-observed-cells",
        str(args.min_observed_cells),
        "--fetch-padding-cpg",
        str(args.fetch_padding_cpg),
        "--cell-id-column",
        args.cell_id_column,
        "--label-column",
        args.label_column,
        "--pat-stem-column",
        args.pat_stem_column,
        "--max-cells",
        str(args.max_cells),
    ]
    if args.skip_tsv:
        command.append("--skip-tsv")
    if args.force:
        command.append("--force")
    run_script("build_fixed_cpg_cell_matrix.py", command)


def add_selection_options(parser, include_run_paths):
    if include_run_paths:
        parser.add_argument("--run-dir", required=True)
        parser.add_argument("--output-dir", required=True)
        parser.add_argument("--batch-count", type=int, default=40)
    parser.add_argument("--cpg-dir", required=True)
    parser.add_argument("--chromosomes", nargs="+", default=None)
    parser.add_argument("--min-hamming-count", type=int, default=5)
    parser.add_argument("--rarefaction-k", type=int, default=100)
    parser.add_argument("--min-valid-batches", type=int, default=None)
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--top-n", type=int, default=20000)
    parser.add_argument("--selection-seed", type=int, default=20261008)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--force-cache", action="store_true")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="hameth_v2",
        description=(
            "Depth-balanced PAT Hamdist bin selection and equal-CpG cell matrix"
        ),
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command")

    prepare = subparsers.add_parser("prepare", help="convert BAM files to indexed PAT")
    prepare.add_argument("-i", "--bam-dir", required=True)
    prepare.add_argument("-o", "--pat-dir", required=True)
    prepare.add_argument("-m", "--workers", type=int, default=8)
    prepare.add_argument("-g", "--genome", required=True)
    prepare.add_argument("-cl", "--chromosomes", nargs="+", required=True)
    prepare.add_argument("--bam-batch-num", type=int, default=1)
    prepare.add_argument(
        "--skip-merge-allcells",
        action="store_true",
        help="do not create optional chromosome-wide combined PAT files",
    )
    prepare.add_argument(
        "-r", "--region-dir", default=None, help="deprecated compatibility option"
    )
    prepare.add_argument(
        "-cp", "--cpg-dir", default=None, help="deprecated compatibility option"
    )
    prepare.set_defaults(func=run_prepare)

    hamdist = subparsers.add_parser(
        "hamdist", help="build batch metrics and select depth-balanced bins"
    )
    hamdist.add_argument("-i", "--pat-dir", required=True)
    hamdist.add_argument("-o", "--output-dir", required=True)
    hamdist.add_argument("-bs", "--batch-size", type=int, default=50)
    hamdist.add_argument("-bn", "--batch-num", type=int, default=40)
    hamdist.add_argument("-m", "--workers", type=int, default=20)
    hamdist.add_argument("-n", "--top-n", type=int, default=20000)
    hamdist.add_argument("-r", "--region-dir", required=True)
    hamdist.add_argument("-cp", "--cpg-dir", required=True)
    hamdist.add_argument("-g", "--genome", required=True)
    hamdist.add_argument("-cl", "--chromosomes", nargs="+", required=True)
    hamdist.add_argument("--bin-size", type=int, default=2000)
    hamdist.add_argument("--candidate-cpgs-per-bin", type=int, default=0)
    hamdist.add_argument("--batch-seed", type=int, default=20261008)
    hamdist.add_argument("--selection-seed", type=int, default=20261008)
    hamdist.add_argument("--min-hamming-count", type=int, default=5)
    hamdist.add_argument("--min-pattern-cpg", type=int, default=2)
    hamdist.add_argument("--max-read-distance", type=int, default=600)
    hamdist.add_argument("--max-reads-per-bin", type=int, default=3000)
    hamdist.add_argument("--basedon0", type=int, choices=[0, 1], default=0)
    hamdist.add_argument("--rarefaction-k", type=int, default=100)
    hamdist.add_argument("--min-valid-batches", type=int, default=None)
    hamdist.add_argument("--repeats", type=int, default=20)
    hamdist.add_argument("--selection-dir", default=None)
    hamdist.add_argument("--force-cache", action="store_true")
    hamdist.set_defaults(func=run_hamdist)

    select_bins = subparsers.add_parser(
        "select-bins", help="apply K/Q10 selection to an existing Hammeth run"
    )
    add_selection_options(select_bins, include_run_paths=True)
    select_bins.set_defaults(func=run_select_bins)

    matrix = subparsers.add_parser(
        "matrix", help="build an equal-CpG cell-by-cell Hamdist matrix"
    )
    matrix.add_argument("--meta", required=True)
    matrix.add_argument("--pat-dir", required=True)
    matrix.add_argument("--cpg-dir", required=True)
    matrix.add_argument("--bins", required=True)
    matrix.add_argument("--output-dir", required=True)
    matrix.add_argument("--workers", type=int, default=20)
    matrix.add_argument("--cpgs-per-bin", type=int, default=20)
    matrix.add_argument("--min-observed-cells", type=int, default=2)
    matrix.add_argument("--fetch-padding-cpg", type=int, default=512)
    matrix.add_argument("--cell-id-column", default="cell_id")
    matrix.add_argument("--label-column", default="MajorType")
    matrix.add_argument("--pat-stem-column", default="pat_stem")
    matrix.add_argument("--max-cells", type=int, default=0)
    matrix.add_argument("--skip-tsv", action="store_true")
    matrix.add_argument("--force", action="store_true")
    matrix.set_defaults(func=run_matrix)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    args.func(args)
    return 0


if __name__ == "__main__":
    main()
