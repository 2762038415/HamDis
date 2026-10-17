#!/usr/bin/env python

"""Build a cell-by-cell Hamdist matrix with equal CpG weight per selected bin.

The ranked bins come from step1_select_balanced_hamdist_bins.py.  For every bin,
this program selects exactly a configurable number of reference CpGs (20 by
default), prioritizing CpGs observed in more cells.  Bins that cannot provide
the full number are dropped.  Consequently, every retained bin contributes the
same number of CpG features to the cell-by-cell matrix.

For a cell pair, Hamdist is:

    discordant methylation states / jointly observed selected CpGs

Single-cell states are reconstructed from indexed PAT files.  A CpG is called
methylated when C_count / (C_count + T_count) > 0.5; ties are nonmethylated.
Cell-type labels are retained for downstream plotting but never used to select
CpGs or calculate distances.
"""

from __future__ import print_function

import argparse
import json
import multiprocessing
import os
import re
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import pysam


_PROFILE_TARGET_BY_CHROM = None
_PROFILE_FETCH_RANGES = None
_PROFILE_PAT_DIR = None
_PROFILE_CHROMOSOMES = None
_PROFILE_OBSERVED = None
_PROFILE_METHYLATED = None

_PAIR_OBSERVED = None
_PAIR_METHYLATED = None
_PAIR_POPCOUNT = None
_PAIR_CELL_COUNT = None


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meta", required=True, help="CSV/TSV metadata table.")
    parser.add_argument(
        "--pat-dir",
        required=True,
        help="Single-cell PAT root: PAT_DIR/chr*/<cell>.pat.gz plus .csi.",
    )
    parser.add_argument(
        "--cpg-dir",
        required=True,
        help="Reference directory containing CpG_chr*.bed files.",
    )
    parser.add_argument(
        "--bins",
        required=True,
        help="Ranked TSV produced by step 1 (must contain rank/chrom/start/end/SD).",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", "20")),
    )
    parser.add_argument("--cpgs-per-bin", type=int, default=20)
    parser.add_argument(
        "--min-observed-cells",
        type=int,
        default=2,
        help="A selectable CpG must be observed in at least this many cells.",
    )
    parser.add_argument(
        "--fetch-padding-cpg",
        type=int,
        default=512,
        help="Fetch this many CpG indices before a target range to recover spanning PATs.",
    )
    parser.add_argument("--cell-id-column", default="cell_id")
    parser.add_argument("--label-column", default="MajorType")
    parser.add_argument(
        "--pat-stem-column",
        default="pat_stem",
        help="Optional metadata column containing PAT filename without .pat.gz.",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=0,
        help="Use the first N cells for a smoke test; 0 means all cells.",
    )
    parser.add_argument(
        "--skip-tsv",
        action="store_true",
        help="Write only compressed NPZ matrices, not the human-readable matrix TSV.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow replacing outputs in an existing output directory.",
    )
    return parser.parse_args()


def chromosome_key(value):
    match = re.fullmatch(r"chr(\d+)", str(value))
    if match:
        return (0, int(match.group(1)))
    if str(value) in ("chrX", "chrY"):
        return (1, {"chrX": 23, "chrY": 24}[str(value)])
    return (2, str(value))


def ensure_directory(path):
    os.makedirs(path, exist_ok=True)


def atomic_to_csv(frame, path, compression=None, **kwargs):
    directory = os.path.dirname(os.path.abspath(path))
    ensure_directory(directory)
    suffix = ".tsv.gz" if compression == "gzip" else ".tsv"
    handle = tempfile.NamedTemporaryFile(
        prefix=".tmp.", suffix=suffix, dir=directory, delete=False
    )
    temporary = handle.name
    handle.close()
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            compression=compression,
            **kwargs
        )
        os.replace(temporary, path)
    except Exception:
        if os.path.exists(temporary):
            os.unlink(temporary)
        raise


def atomic_savez(path, **arrays):
    directory = os.path.dirname(os.path.abspath(path))
    ensure_directory(directory)
    handle = tempfile.NamedTemporaryFile(
        prefix=".tmp.", suffix=".npz", dir=directory, delete=False
    )
    temporary = handle.name
    handle.close()
    try:
        np.savez_compressed(temporary, **arrays)
        os.replace(temporary, path)
    except Exception:
        if os.path.exists(temporary):
            os.unlink(temporary)
        raise


def read_metadata(path):
    extension = os.path.splitext(path)[1].lower()
    separator = "\t" if extension in (".tsv", ".txt") else ","
    return pd.read_csv(path, sep=separator, dtype=str)


def discover_pat_stems(pat_dir, chromosome):
    chromosome_dir = os.path.join(pat_dir, chromosome)
    if not os.path.isdir(chromosome_dir):
        raise FileNotFoundError(chromosome_dir)
    stems = {}
    duplicate_ids = set()
    for name in os.listdir(chromosome_dir):
        if not name.endswith(".pat.gz"):
            continue
        stem = name[:-7]
        cell_id = name.split(".", 1)[0]
        if cell_id in stems and stems[cell_id] != stem:
            duplicate_ids.add(cell_id)
        stems[cell_id] = stem
        stems.setdefault(stem, stem)
    if duplicate_ids:
        raise ValueError(
            "Ambiguous PAT stems for cell IDs: {}".format(
                ", ".join(sorted(duplicate_ids)[:5])
            )
        )
    return stems


def load_manifest(args, first_chromosome):
    metadata = read_metadata(args.meta)
    if args.cell_id_column not in metadata.columns:
        raise ValueError(
            "Metadata lacks cell ID column {!r}".format(args.cell_id_column)
        )
    keep_columns = [args.cell_id_column]
    if args.label_column in metadata.columns:
        keep_columns.append(args.label_column)
    if args.pat_stem_column in metadata.columns:
        keep_columns.append(args.pat_stem_column)
    manifest = metadata[keep_columns].copy()
    manifest = manifest.rename(columns={args.cell_id_column: "cell_id"})
    manifest["cell_id"] = manifest["cell_id"].str.replace(r"\.bam$", "", regex=True)
    if manifest["cell_id"].isna().any() or manifest["cell_id"].duplicated().any():
        raise ValueError("Metadata contains missing or duplicate cell IDs")

    stem_lookup = discover_pat_stems(args.pat_dir, first_chromosome)
    if args.pat_stem_column in manifest.columns:
        explicit_stems = manifest[args.pat_stem_column].str.replace(
            r"\.pat\.gz$", "", regex=True
        )
        if args.pat_stem_column != "pat_stem":
            manifest = manifest.drop(columns=[args.pat_stem_column])
        manifest["pat_stem"] = explicit_stems
    else:
        manifest["pat_stem"] = manifest["cell_id"].map(stem_lookup)
    missing = manifest.loc[manifest["pat_stem"].isna(), "cell_id"].tolist()
    if missing:
        raise FileNotFoundError(
            "{} metadata cells lack {} PAT files; examples: {}".format(
                len(missing), first_chromosome, ", ".join(missing[:5])
            )
        )
    if args.label_column in manifest.columns:
        manifest = manifest.rename(columns={args.label_column: "cell_type"})
        manifest["cell_type"] = manifest["cell_type"].fillna("NA").str.strip()
        sort_columns = ["cell_type", "cell_id"]
    else:
        manifest["cell_type"] = "NA"
        sort_columns = ["cell_id"]
    manifest = manifest.sort_values(sort_columns).reset_index(drop=True)
    if args.max_cells:
        manifest = manifest.head(args.max_cells).copy()
    manifest.insert(0, "row_index", np.arange(len(manifest), dtype=int))
    return manifest


def load_bins(path):
    frame = pd.read_csv(path, sep="\t")
    rank_column = "SD_rank" if "SD_rank" in frame.columns else "rank"
    score_column = "SD" if "SD" in frame.columns else "rarefied_hamdist_sd_q10"
    required = {rank_column, score_column, "chrom", "start", "end"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError("Ranked-bin table lacks {}".format(sorted(missing)))
    frame = frame[[rank_column, score_column, "chrom", "start", "end"]].copy()
    frame.columns = ["SD_rank", "SD", "chrom", "start", "end"]
    frame["SD_rank"] = pd.to_numeric(frame["SD_rank"], errors="raise").astype(int)
    frame["SD"] = pd.to_numeric(frame["SD"], errors="raise")
    frame["start"] = pd.to_numeric(frame["start"], errors="raise").astype(int)
    frame["end"] = pd.to_numeric(frame["end"], errors="raise").astype(int)
    frame = frame.sort_values("SD_rank").reset_index(drop=True)
    if frame.empty or frame["SD_rank"].tolist() != list(range(1, len(frame) + 1)):
        raise ValueError("Ranked bins must have contiguous ranks 1..N")
    if frame.duplicated(["chrom", "start", "end"]).any():
        raise ValueError("Ranked-bin table contains duplicate intervals")
    if (frame["end"] < frame["start"]).any():
        raise ValueError("Ranked-bin table contains end < start")
    return frame


def load_reference_chromosome(path):
    positions = []
    identifiers = []
    with open(path, "r") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                positions.append(int(fields[1]))
                identifiers.append(int(fields[2]))
            except ValueError:
                if line_number == 1:
                    continue
                raise
    if not positions:
        raise ValueError("No CpGs in {}".format(path))
    positions = np.asarray(positions, dtype=np.int64)
    identifiers = np.asarray(identifiers, dtype=np.int64)
    if np.any(positions[1:] < positions[:-1]):
        raise ValueError("CpG reference is not position sorted: {}".format(path))
    return positions, identifiers


def merge_ranges(ranges):
    merged = []
    for start, end in sorted(ranges):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [tuple(values) for values in merged]


def build_target_index(bins, cpg_dir, padding, chromosomes):
    target_rows = []
    fetch_ranges = {}
    bin_counts = {}
    for chromosome in chromosomes:
        chromosome_bins = bins.loc[bins["chrom"].eq(chromosome)].sort_values("start")
        if chromosome_bins.empty:
            continue
        reference_path = os.path.join(cpg_dir, "CpG_{}.bed".format(chromosome))
        if not os.path.isfile(reference_path):
            raise FileNotFoundError(reference_path)
        positions, identifiers = load_reference_chromosome(reference_path)
        raw_ranges = []
        chromosome_cpg_count = 0
        for row in chromosome_bins.itertuples(index=False):
            left = int(np.searchsorted(positions, row.start, side="left"))
            right = int(np.searchsorted(positions, row.end, side="right"))
            count = right - left
            if count < 1:
                raise ValueError(
                    "Selected bin {}:{}-{} has no reference CpGs".format(
                        chromosome, row.start, row.end
                    )
                )
            bin_counts[int(row.SD_rank)] = count
            chromosome_cpg_count += count
            for position, identifier in zip(positions[left:right], identifiers[left:right]):
                target_rows.append(
                    {
                        "chrom": chromosome,
                        "position": int(position),
                        "reference_cpg_id": int(identifier),
                        "bin_SD_rank": int(row.SD_rank),
                        "bin_start": int(row.start),
                        "bin_end": int(row.end),
                        "bin_SD": float(row.SD),
                    }
                )
            raw_ranges.append(
                (
                    max(1, int(identifiers[left]) - padding),
                    int(identifiers[right - 1]),
                )
            )
        fetch_ranges[chromosome] = merge_ranges(raw_ranges)
        print(
            "[TARGET] {} bins={:,} CpGs={:,} indexed_ranges={:,}".format(
                chromosome,
                len(chromosome_bins),
                chromosome_cpg_count,
                len(fetch_ranges[chromosome]),
            ),
            flush=True,
        )

    target = pd.DataFrame(target_rows)
    chromosome_order = {chromosome: index for index, chromosome in enumerate(chromosomes)}
    target["_chrom_order"] = target["chrom"].map(chromosome_order)
    target = target.sort_values(
        ["_chrom_order", "position", "reference_cpg_id"]
    ).drop(columns="_chrom_order").reset_index(drop=True)
    if target.duplicated(["chrom", "position"]).any():
        raise ValueError("Selected bins overlap or produce duplicate CpG positions")
    if target.duplicated(["chrom", "reference_cpg_id"]).any():
        raise ValueError("Selected bins produce duplicate PAT CpG identifiers")
    target.insert(0, "target_cpg_index", np.arange(len(target), dtype=int))
    target_by_chrom = {
        chromosome: dict(
            zip(
                group["reference_cpg_id"].astype(int),
                group["target_cpg_index"].astype(int),
            )
        )
        for chromosome, group in target.groupby("chrom", sort=False)
    }
    return target, target_by_chrom, fetch_ranges, bin_counts


def init_profile_worker(observed_path, methylated_path, shape):
    global _PROFILE_OBSERVED, _PROFILE_METHYLATED
    _PROFILE_OBSERVED = np.memmap(
        observed_path, dtype=np.uint8, mode="r+", shape=shape
    )
    _PROFILE_METHYLATED = np.memmap(
        methylated_path, dtype=np.uint8, mode="r+", shape=shape
    )


def build_cell_profile(job):
    row_index, cell_id, pat_stem = job
    target_count = sum(len(values) for values in _PROFILE_TARGET_BY_CHROM.values())
    methylated_counts = np.zeros(target_count, dtype=np.uint32)
    total_counts = np.zeros(target_count, dtype=np.uint32)
    parsed_records = 0
    maximum_pattern_length = 0
    for chromosome in _PROFILE_CHROMOSOMES:
        target_lookup = _PROFILE_TARGET_BY_CHROM.get(chromosome)
        if not target_lookup:
            continue
        pat_path = os.path.join(_PROFILE_PAT_DIR, chromosome, pat_stem + ".pat.gz")
        index_path = pat_path + ".csi"
        if not os.path.isfile(pat_path) or not os.path.isfile(index_path):
            raise FileNotFoundError(
                "Missing PAT or CSI for {} {}: {}".format(cell_id, chromosome, pat_path)
            )
        with pysam.TabixFile(pat_path, index=index_path) as tabix:
            for start_id, end_id in _PROFILE_FETCH_RANGES[chromosome]:
                for line in tabix.fetch(chromosome, max(0, start_id - 1), end_id):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 4:
                        continue
                    first_cpg_id = int(fields[1])
                    pattern = fields[2]
                    count = int(fields[3])
                    if count < 1:
                        continue
                    parsed_records += 1
                    maximum_pattern_length = max(maximum_pattern_length, len(pattern))
                    for offset, symbol in enumerate(pattern):
                        if symbol == ".":
                            continue
                        target_index = target_lookup.get(first_cpg_id + offset)
                        if target_index is None:
                            continue
                        if symbol not in ("C", "T"):
                            raise ValueError(
                                "Unexpected PAT symbol {!r} in {}".format(symbol, pat_path)
                            )
                        total_counts[target_index] += count
                        if symbol == "C":
                            methylated_counts[target_index] += count
    observed = total_counts > 0
    methylated = observed & (methylated_counts * 2 > total_counts)
    _PROFILE_OBSERVED[row_index, :] = np.packbits(observed)
    _PROFILE_METHYLATED[row_index, :] = np.packbits(methylated)
    return (
        row_index,
        cell_id,
        int(np.count_nonzero(observed)),
        parsed_records,
        maximum_pattern_length,
    )


def calculate_cpg_cell_coverage(observed_path, packed_shape, target_count):
    packed = np.memmap(observed_path, dtype=np.uint8, mode="r", shape=packed_shape)
    coverage = np.zeros(target_count, dtype=np.int32)
    for row_index in range(packed_shape[0]):
        coverage += np.unpackbits(packed[row_index])[:target_count].astype(np.int32)
        if row_index == 0 or (row_index + 1) % 250 == 0 or row_index + 1 == packed_shape[0]:
            print(
                "[COVERAGE] {}/{} cell profiles".format(row_index + 1, packed_shape[0]),
                flush=True,
            )
    return coverage


def select_fixed_cpgs(target, coverage, bin_counts, cpgs_per_bin, minimum_cells):
    annotated = target.copy()
    annotated["observed_cell_count"] = coverage
    annotated["selected_for_matrix"] = False
    selected_indices = []
    support_rows = []
    for rank, frame in annotated.groupby("bin_SD_rank", sort=True):
        supported = frame.loc[frame["observed_cell_count"].ge(minimum_cells)].copy()
        supported = supported.sort_values(
            ["observed_cell_count", "reference_cpg_id"], ascending=[False, True]
        )
        chosen = supported.head(cpgs_per_bin)
        support_rows.append(
            {
                "SD_rank": int(rank),
                "reference_cpg_count": int(bin_counts[int(rank)]),
                "supported_cpg_count": int(len(supported)),
                "selected_cpg_count": int(len(chosen)),
                "selected_min_observed_cell_count": (
                    int(chosen["observed_cell_count"].min()) if len(chosen) else 0
                ),
                "selected_median_observed_cell_count": (
                    float(chosen["observed_cell_count"].median()) if len(chosen) else 0.0
                ),
            }
        )
        if len(chosen) == cpgs_per_bin:
            selected_indices.extend(chosen["target_cpg_index"].astype(int).tolist())

    support = pd.DataFrame(support_rows).sort_values("SD_rank")
    support["included_in_matrix"] = support["selected_cpg_count"].eq(cpgs_per_bin)
    failed = support.loc[~support["included_in_matrix"]]
    if not failed.empty:
        print(
            "[WARNING] dropping {} of {} bins lacking {} CpGs observed in >= {} cells".format(
                len(failed), len(support), cpgs_per_bin, minimum_cells
            ),
            flush=True,
        )
    selected_indices = np.asarray(sorted(selected_indices), dtype=np.int64)
    expected = int(support["included_in_matrix"].sum()) * cpgs_per_bin
    if expected < 1:
        raise ValueError("No bins have enough supported CpGs for the matrix")
    if len(selected_indices) != expected or len(np.unique(selected_indices)) != expected:
        raise ValueError("Fixed-CpG selection is not one-to-one")
    annotated.loc[selected_indices, "selected_for_matrix"] = True
    selected = annotated.loc[selected_indices].copy().reset_index(drop=True)
    selected.insert(0, "matrix_cpg_index", np.arange(len(selected), dtype=int))
    return annotated, selected, support, selected_indices


def write_reduced_profiles(
    full_observed_path,
    full_methylated_path,
    full_shape,
    target_count,
    selected_indices,
    output_dir,
):
    selected_count = len(selected_indices)
    reduced_shape = (full_shape[0], (selected_count + 7) // 8)
    observed_path = os.path.join(output_dir, "observed_fixed_cpg_bits.uint8.memmap")
    methylated_path = os.path.join(output_dir, "methylated_fixed_cpg_bits.uint8.memmap")
    full_observed = np.memmap(
        full_observed_path, dtype=np.uint8, mode="r", shape=full_shape
    )
    full_methylated = np.memmap(
        full_methylated_path, dtype=np.uint8, mode="r", shape=full_shape
    )
    reduced_observed = np.memmap(
        observed_path, dtype=np.uint8, mode="w+", shape=reduced_shape
    )
    reduced_methylated = np.memmap(
        methylated_path, dtype=np.uint8, mode="w+", shape=reduced_shape
    )
    reduced_observed[:] = 0
    reduced_methylated[:] = 0
    observed_counts = np.zeros(full_shape[0], dtype=np.int32)
    for row_index in range(full_shape[0]):
        observed = np.unpackbits(full_observed[row_index])[:target_count][selected_indices]
        methylated = np.unpackbits(full_methylated[row_index])[:target_count][selected_indices]
        reduced_observed[row_index] = np.packbits(observed)
        reduced_methylated[row_index] = np.packbits(methylated)
        observed_counts[row_index] = int(np.count_nonzero(observed))
        if row_index == 0 or (row_index + 1) % 250 == 0 or row_index + 1 == full_shape[0]:
            print(
                "[REDUCE] {}/{} cell profiles".format(row_index + 1, full_shape[0]),
                flush=True,
            )
    reduced_observed.flush()
    reduced_methylated.flush()
    del reduced_observed
    del reduced_methylated
    return observed_path, methylated_path, reduced_shape, observed_counts


def init_pair_worker(observed_path, methylated_path, shape, cell_count):
    global _PAIR_OBSERVED, _PAIR_METHYLATED, _PAIR_POPCOUNT, _PAIR_CELL_COUNT
    _PAIR_OBSERVED = np.memmap(observed_path, dtype=np.uint8, mode="r", shape=shape)
    _PAIR_METHYLATED = np.memmap(
        methylated_path, dtype=np.uint8, mode="r", shape=shape
    )
    values = np.arange(256, dtype=np.uint8).reshape(256, 1)
    _PAIR_POPCOUNT = np.unpackbits(values, axis=1).sum(axis=1)
    _PAIR_CELL_COUNT = cell_count


def calculate_distance_row(first):
    pair_count = _PAIR_CELL_COUNT - first - 1
    distances = np.full(pair_count, np.nan, dtype=np.float32)
    shared_counts = np.zeros(pair_count, dtype=np.int32)
    first_observed = _PAIR_OBSERVED[first]
    first_methylated = _PAIR_METHYLATED[first]
    for offset, second in enumerate(range(first + 1, _PAIR_CELL_COUNT)):
        valid_bits = np.bitwise_and(first_observed, _PAIR_OBSERVED[second])
        shared = int(_PAIR_POPCOUNT[valid_bits].sum(dtype=np.uint64))
        shared_counts[offset] = shared
        if shared == 0:
            continue
        mismatch_bits = np.bitwise_and(
            np.bitwise_xor(first_methylated, _PAIR_METHYLATED[second]), valid_bits
        )
        mismatches = int(_PAIR_POPCOUNT[mismatch_bits].sum(dtype=np.uint64))
        distances[offset] = mismatches / float(shared)
    return first, distances, shared_counts


def calculate_pair_matrix(
    observed_path, methylated_path, shape, observed_counts, workers
):
    cell_count = shape[0]
    distance_matrix = np.full((cell_count, cell_count), np.nan, dtype=np.float32)
    shared_matrix = np.zeros((cell_count, cell_count), dtype=np.int32)
    np.fill_diagonal(distance_matrix, 0.0)
    np.fill_diagonal(shared_matrix, observed_counts)
    context = multiprocessing.get_context("fork")
    with ProcessPoolExecutor(
        max_workers=min(workers, max(1, cell_count - 1)),
        initializer=init_pair_worker,
        initargs=(observed_path, methylated_path, shape, cell_count),
        mp_context=context,
    ) as executor:
        for completed, (first, distances, shared_counts) in enumerate(
            executor.map(calculate_distance_row, range(cell_count - 1), chunksize=1), 1
        ):
            second_indices = np.arange(first + 1, cell_count)
            distance_matrix[first, second_indices] = distances
            distance_matrix[second_indices, first] = distances
            shared_matrix[first, second_indices] = shared_counts
            shared_matrix[second_indices, first] = shared_counts
            if completed == 1 or completed % 25 == 0 or completed == cell_count - 1:
                print(
                    "[PAIRS] completed first rows {}/{}".format(
                        completed, cell_count - 1
                    ),
                    flush=True,
                )
    return distance_matrix, shared_matrix


def validate_args(args):
    if args.workers < 1 or args.cpgs_per_bin < 1 or args.min_observed_cells < 1:
        raise ValueError("workers, cpgs-per-bin, and min-observed-cells must be positive")
    if args.fetch_padding_cpg < 1 or args.max_cells < 0:
        raise ValueError("fetch-padding-cpg must be positive; max-cells cannot be negative")
    for path in (args.meta, args.bins):
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    for path in (args.pat_dir, args.cpg_dir):
        if not os.path.isdir(path):
            raise FileNotFoundError(path)
    final_path = os.path.join(args.output_dir, "all_cells_hamdist_matrix.npz")
    if os.path.isfile(final_path) and not args.force:
        raise FileExistsError(
            "{} already exists; use --force to replace outputs".format(final_path)
        )


def main():
    args = parse_args()
    validate_args(args)
    ensure_directory(args.output_dir)
    started = time.time()
    bins = load_bins(args.bins)
    chromosomes = sorted(bins["chrom"].astype(str).unique(), key=chromosome_key)
    manifest = load_manifest(args, chromosomes[0])
    if len(manifest) < 2:
        raise ValueError("At least two cells are required")
    target, target_by_chrom, fetch_ranges, bin_counts = build_target_index(
        bins, args.cpg_dir, args.fetch_padding_cpg, chromosomes
    )
    atomic_to_csv(
        target,
        os.path.join(args.output_dir, "all_target_cpg_index.tsv.gz"),
        compression="gzip",
        index=False,
    )
    atomic_to_csv(
        manifest, os.path.join(args.output_dir, "cell_manifest.tsv"), index=False
    )

    global _PROFILE_TARGET_BY_CHROM, _PROFILE_FETCH_RANGES
    global _PROFILE_PAT_DIR, _PROFILE_CHROMOSOMES
    _PROFILE_TARGET_BY_CHROM = target_by_chrom
    _PROFILE_FETCH_RANGES = fetch_ranges
    _PROFILE_PAT_DIR = args.pat_dir
    _PROFILE_CHROMOSOMES = chromosomes

    target_count = len(target)
    packed_shape = (len(manifest), (target_count + 7) // 8)
    observed_path = os.path.join(args.output_dir, "observed_all_cpg_bits.uint8.memmap")
    methylated_path = os.path.join(
        args.output_dir, "methylated_all_cpg_bits.uint8.memmap"
    )
    for path in (observed_path, methylated_path):
        matrix = np.memmap(path, dtype=np.uint8, mode="w+", shape=packed_shape)
        matrix[:] = 0
        matrix.flush()
        del matrix

    print(
        "[INFO] cells={:,} cell_types={} bins={:,} reference_CpGs={:,} workers={}".format(
            len(manifest),
            manifest["cell_type"].nunique(),
            len(bins),
            target_count,
            args.workers,
        ),
        flush=True,
    )
    profile_jobs = list(
        manifest[["row_index", "cell_id", "pat_stem"]].itertuples(
            index=False, name=None
        )
    )
    profile_rows = []
    context = multiprocessing.get_context("fork")
    with ProcessPoolExecutor(
        max_workers=min(args.workers, len(profile_jobs)),
        initializer=init_profile_worker,
        initargs=(observed_path, methylated_path, packed_shape),
        mp_context=context,
    ) as executor:
        for completed, result in enumerate(
            executor.map(build_cell_profile, profile_jobs, chunksize=1), 1
        ):
            profile_rows.append(result)
            if completed == 1 or completed % 25 == 0 or completed == len(profile_jobs):
                print(
                    "[PROFILE] {}/{} latest={} observed_CpGs={:,} max_pattern={}".format(
                        completed, len(profile_jobs), result[1], result[2], result[4]
                    ),
                    flush=True,
                )
    profile_summary = pd.DataFrame(
        profile_rows,
        columns=[
            "row_index",
            "cell_id",
            "observed_all_target_cpg_count",
            "parsed_pat_record_count",
            "maximum_pattern_length",
        ],
    ).sort_values("row_index")
    if profile_summary["observed_all_target_cpg_count"].le(0).any():
        raise ValueError("At least one cell has no observed target CpGs")
    if profile_summary["maximum_pattern_length"].max() > args.fetch_padding_cpg:
        raise ValueError(
            "A PAT pattern is longer than --fetch-padding-cpg; rerun with a larger padding"
        )

    coverage = calculate_cpg_cell_coverage(observed_path, packed_shape, target_count)
    annotated, selected, support, selected_indices = select_fixed_cpgs(
        target,
        coverage,
        bin_counts,
        args.cpgs_per_bin,
        args.min_observed_cells,
    )
    atomic_to_csv(
        annotated,
        os.path.join(args.output_dir, "all_target_cpg_coverage_and_selection.tsv.gz"),
        compression="gzip",
        index=False,
    )
    atomic_to_csv(
        selected,
        os.path.join(args.output_dir, "matrix_target_cpg_index.tsv.gz"),
        compression="gzip",
        index=False,
    )
    atomic_to_csv(
        support,
        os.path.join(args.output_dir, "fixed_cpg_bin_support.tsv"),
        index=False,
    )

    reduced_observed, reduced_methylated, reduced_shape, observed_counts = (
        write_reduced_profiles(
            observed_path,
            methylated_path,
            packed_shape,
            target_count,
            selected_indices,
            args.output_dir,
        )
    )
    profile_summary["observed_fixed_target_cpg_count"] = observed_counts
    atomic_to_csv(
        manifest.merge(
            profile_summary, on=["row_index", "cell_id"], validate="one_to_one"
        ),
        os.path.join(args.output_dir, "cell_profile_summary.tsv"),
        index=False,
    )

    distance_matrix, shared_matrix = calculate_pair_matrix(
        reduced_observed,
        reduced_methylated,
        reduced_shape,
        observed_counts,
        args.workers,
    )
    cell_count = len(manifest)
    upper = np.triu_indices(cell_count, k=1)
    distances = distance_matrix[upper]
    shared = shared_matrix[upper]
    nonfinite = int((~np.isfinite(distances)).sum())
    if nonfinite:
        raise ValueError(
            "Complete matrix contains {} cell pairs with no shared selected CpGs".format(
                nonfinite
            )
        )
    if not np.allclose(distance_matrix, distance_matrix.T, atol=1e-7):
        raise ValueError("Hamdist matrix is not symmetric")
    if float(distances.min()) < 0 or float(distances.max()) > 1:
        raise ValueError("Hamdist values fall outside [0,1]")

    cell_ids = np.asarray(manifest["cell_id"].astype(str).tolist(), dtype=str)
    atomic_savez(
        os.path.join(args.output_dir, "all_cells_hamdist_matrix.npz"),
        matrix=distance_matrix,
        shared_cpg_count=shared_matrix,
        cell_ids=cell_ids,
    )
    if not args.skip_tsv:
        matrix_frame = pd.DataFrame(distance_matrix, index=cell_ids, columns=cell_ids)
        matrix_frame.index.name = "cell_id"
        atomic_to_csv(
            matrix_frame,
            os.path.join(args.output_dir, "all_cells_hamdist_matrix.tsv.gz"),
            compression="gzip",
            index=True,
            na_rep="NA",
            float_format="%.8g",
        )
        shared_frame = pd.DataFrame(shared_matrix, index=cell_ids, columns=cell_ids)
        shared_frame.index.name = "cell_id"
        atomic_to_csv(
            shared_frame,
            os.path.join(args.output_dir, "all_cells_shared_cpg_count.tsv.gz"),
            compression="gzip",
            index=True,
        )

    elapsed = time.time() - started
    summary = pd.DataFrame(
        [
            ("cell_count", cell_count),
            ("cell_type_count", manifest["cell_type"].nunique()),
            ("input_selected_bin_count", len(bins)),
            ("retained_equal_weight_bin_count", int(support["included_in_matrix"].sum())),
            ("dropped_bin_count", int((~support["included_in_matrix"]).sum())),
            ("all_reference_cpg_count", target_count),
            ("fixed_cpgs_per_bin", args.cpgs_per_bin),
            ("matrix_target_cpg_count", len(selected_indices)),
            ("cell_pair_count", len(distances)),
            ("hamdist_min", float(distances.min())),
            ("hamdist_median", float(np.median(distances))),
            ("hamdist_mean", float(np.mean(distances))),
            ("hamdist_max", float(distances.max())),
            ("shared_cpg_min", int(shared.min())),
            ("shared_cpg_median", float(np.median(shared))),
            ("shared_cpg_mean", float(np.mean(shared))),
            ("shared_cpg_max", int(shared.max())),
            ("maximum_observed_pattern_length", int(profile_summary["maximum_pattern_length"].max())),
            ("elapsed_seconds", elapsed),
        ],
        columns=["metric", "value"],
    )
    atomic_to_csv(summary, os.path.join(args.output_dir, "matrix_summary.tsv"), index=False)
    parameters = vars(args).copy()
    parameters.update(
        {
            "chromosomes": chromosomes,
            "elapsed_seconds": elapsed,
            "state_rule": "methylated iff C_count/(C_count+T_count) > 0.5",
            "distance_rule": "discordant jointly observed CpGs / jointly observed CpGs",
        }
    )
    with open(os.path.join(args.output_dir, "run_parameters.json"), "w") as handle:
        json.dump(parameters, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        "[FINISHED] cells={:,} retained_bins={:,} target_CpGs={:,} elapsed={:.1f}s".format(
            cell_count,
            int(support["included_in_matrix"].sum()),
            len(selected_indices),
            elapsed,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
