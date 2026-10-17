#!/usr/bin/env python

"""Select depth-balanced Hamdist bins directly from a completed Hammeth run.

For each genomic bin and each Hammeth sampling batch, this program reconstructs
PAT-level sufficient statistics:

    N = number of valid cell-pair/CpG comparisons
    M = number of methylation-state mismatches among those comparisons

It then samples exactly K comparisons without replacement, applies a
variance-stabilizing arcsine transform, and calculates the SD across batches.
The rarefaction is repeated, and the 10th percentile (Q10) of the repeated SDs
is used as the conservative ranking score.  Cell-type labels are never used.

The PAT scan is expensive, so sufficient statistics are cached per chromosome.
Changing K, repeat count, seed, or Top-N reuses that cache.  Use --force-cache
only when the source metrics/PAT files or the support-counting rules changed.
"""

from __future__ import print_function

import argparse
import glob
import json
import os
import re
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd


METRIC_COLUMNS = ["chrom", "start", "end", "HammingCount"]
BATCH_RE = re.compile(r"_merge_(\d+)_MeConcord\.txt$")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Completed Hammeth hamdist run containing middle_res/metrics and mergepats.",
    )
    parser.add_argument(
        "--cpg-dir",
        required=True,
        help="Reference directory containing CpG_chr*.bed files.",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--chromosomes",
        nargs="+",
        default=None,
        help="Default: discover chromosomes under middle_res/metrics.",
    )
    parser.add_argument("--batch-count", type=int, default=40)
    parser.add_argument("--min-hamming-count", type=int, default=5)
    parser.add_argument("--rarefaction-k", type=int, default=100)
    parser.add_argument("--min-valid-batches", type=int, default=36)
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--top-n", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20261008)
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", "20")),
    )
    parser.add_argument(
        "--force-cache",
        action="store_true",
        help="Rebuild per-chromosome PAT sufficient-statistics caches.",
    )
    return parser.parse_args()


def chromosome_key(value):
    match = re.fullmatch(r"chr(\d+)", str(value))
    if match:
        return (0, int(match.group(1)))
    if str(value) in ("chrX", "chrY"):
        return (1, {"chrX": 23, "chrY": 24}[str(value)])
    return (2, str(value))


def chromosome_seed_index(value):
    match = re.fullmatch(r"chr(\d+)", str(value))
    if match:
        return int(match.group(1))
    return {"chrX": 23, "chrY": 24}.get(str(value), 2000)


def ensure_directory(path):
    os.makedirs(path, exist_ok=True)


def atomic_to_csv(frame, path, compression=None, header=True):
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
            index=False,
            header=header,
            compression=compression,
            na_rep="NA",
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


def discover_chromosomes(run_dir):
    metric_root = os.path.join(run_dir, "middle_res", "metrics")
    chromosomes = [
        os.path.basename(path)
        for path in glob.glob(os.path.join(metric_root, "chr*"))
        if os.path.isdir(path)
    ]
    if not chromosomes:
        raise FileNotFoundError("No chromosome directories under {}".format(metric_root))
    return sorted(chromosomes, key=chromosome_key)


def batch_id_from_metric(path):
    match = BATCH_RE.search(os.path.basename(path))
    if not match:
        raise ValueError("Cannot parse batch ID from {}".format(path))
    return int(match.group(1))


def discover_batches(run_dir, chromosome, expected_count):
    metric_dir = os.path.join(run_dir, "middle_res", "metrics", chromosome)
    metric_paths = glob.glob(os.path.join(metric_dir, "*_MeConcord.txt"))
    metrics = {batch_id_from_metric(path): path for path in metric_paths}
    expected = set(range(expected_count))
    if set(metrics) != expected:
        raise ValueError(
            "{} metric batches differ from 0..{}; missing={}, extra={}".format(
                chromosome,
                expected_count - 1,
                sorted(expected - set(metrics)),
                sorted(set(metrics) - expected),
            )
        )
    merge_dir = os.path.join(run_dir, "middle_res", "mergepats")
    pats = {
        batch_id: os.path.join(
            merge_dir, "{}_merge_{}.pat".format(chromosome, batch_id)
        )
        for batch_id in range(expected_count)
    }
    missing = [path for path in pats.values() if not os.path.isfile(path)]
    if missing:
        raise FileNotFoundError("Missing merged PAT: {}".format(missing[0]))
    return metrics, pats


def load_metric(path, expected_chromosome):
    frame = pd.read_csv(path, sep="\t", usecols=METRIC_COLUMNS)
    if frame.empty:
        raise ValueError("Empty metric file: {}".format(path))
    frame["chrom"] = frame["chrom"].astype(str)
    frame["start"] = pd.to_numeric(frame["start"], errors="raise").astype(np.int64)
    frame["end"] = pd.to_numeric(frame["end"], errors="raise").astype(np.int64)
    if set(frame["chrom"]) != {expected_chromosome}:
        raise ValueError("Unexpected chromosome in {}".format(path))
    if frame.duplicated(["chrom", "start", "end"]).any():
        raise ValueError("Duplicate bins in {}".format(path))
    return frame.sort_values(["start", "end"]).reset_index(drop=True)


def build_cpg_to_bin(cpg_path, coordinates):
    starts = coordinates["start"].to_numpy(dtype=np.int64)
    ends = coordinates["end"].to_numpy(dtype=np.int64)
    if len(starts) > 1 and np.any(starts[1:] <= ends[:-1]):
        raise ValueError("Overlapping or unsorted bins in metric coordinates")

    orders = []
    positions = []
    with open(cpg_path, "r") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                positions.append(int(fields[1]))
                orders.append(int(fields[2]))
            except ValueError:
                if line_number == 1:
                    continue
                raise
    if not orders:
        raise ValueError("No CpGs in {}".format(cpg_path))

    orders = np.asarray(orders, dtype=np.int64)
    positions = np.asarray(positions, dtype=np.int64)
    if len(positions) > 1 and np.any(positions[1:] < positions[:-1]):
        raise ValueError("CpG reference is not position sorted: {}".format(cpg_path))
    cpg_to_bin = np.full(int(orders.max()) + 1, -1, dtype=np.int32)
    candidate = np.searchsorted(ends, positions, side="left")
    inside = candidate < len(starts)
    inside_indices = np.flatnonzero(inside)
    candidate_bins = candidate[inside]
    inside[inside_indices] = positions[inside_indices] >= starts[candidate_bins]
    cpg_to_bin[orders[inside]] = candidate[inside].astype(np.int32)
    return cpg_to_bin


def make_popcount16():
    values = np.arange(65536, dtype=np.uint32)
    counts = np.zeros(65536, dtype=np.uint8)
    for _ in range(16):
        counts += (values & 1).astype(np.uint8)
        values >>= 1
    return counts


POPCOUNT16 = make_popcount16()


def popcount_uint64(values):
    result = np.zeros(values.shape, dtype=np.uint8)
    for shift in (0, 16, 32, 48):
        chunks = ((values >> np.uint64(shift)) & np.uint64(65535)).astype(np.uint16)
        result += POPCOUNT16[chunks]
    return result


def scan_pat_sufficient_statistics(pat_path, cpg_to_bin, bin_count):
    """Return effective CpGs, valid comparisons, and mismatches for one batch.

    A PAT record contributes only when at least two C/T calls from that record
    fall in the same bin, matching Hammeth's read-level comparison semantics.
    For one cell/CpG, conflicting C and T evidence is treated as ambiguous and
    excluded rather than counted toward either state.
    """
    methylated_cells = np.zeros(len(cpg_to_bin), dtype=np.uint64)
    unmethylated_cells = np.zeros(len(cpg_to_bin), dtype=np.uint64)
    cell_codes = {}

    with open(pat_path, "r") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            try:
                first_cpg = int(fields[1])
            except ValueError:
                if line_number == 1:
                    continue
                raise
            pattern = fields[2]
            cell_id = fields[4]
            if cell_id not in cell_codes:
                cell_codes[cell_id] = len(cell_codes)
                if len(cell_codes) > 64:
                    raise ValueError(
                        "{} contains >64 cells; this implementation uses uint64 masks. "
                        "Use Hammeth batch size <=64.".format(pat_path)
                    )
            cell_bit = np.uint64(1) << np.uint64(cell_codes[cell_id])

            hits_by_bin = {}
            for offset, state in enumerate(pattern):
                if state not in ("C", "T"):
                    continue
                cpg_order = first_cpg + offset
                if cpg_order < 0 or cpg_order >= len(cpg_to_bin):
                    continue
                bin_index = int(cpg_to_bin[cpg_order])
                if bin_index >= 0:
                    hits_by_bin.setdefault(bin_index, []).append((cpg_order, state))

            for hits in hits_by_bin.values():
                if len(hits) < 2:
                    continue
                for cpg_order, state in hits:
                    if state == "C":
                        methylated_cells[cpg_order] |= cell_bit
                    else:
                        unmethylated_cells[cpg_order] |= cell_bit

    conflicts = methylated_cells & unmethylated_cells
    methylated_count = popcount_uint64(methylated_cells & ~conflicts).astype(np.int64)
    unmethylated_count = popcount_uint64(unmethylated_cells & ~conflicts).astype(np.int64)
    observed_count = methylated_count + unmethylated_count
    mapped = cpg_to_bin >= 0
    mapped_bins = cpg_to_bin[mapped]

    effective = observed_count[mapped] >= 2
    effective_count = np.bincount(
        mapped_bins[effective], minlength=bin_count
    ).astype(np.int32)
    common_per_cpg = observed_count * (observed_count - 1) // 2
    mismatch_per_cpg = methylated_count * unmethylated_count
    common = np.bincount(
        mapped_bins,
        weights=common_per_cpg[mapped],
        minlength=bin_count,
    ).astype(np.int64)
    mismatch = np.bincount(
        mapped_bins,
        weights=mismatch_per_cpg[mapped],
        minlength=bin_count,
    ).astype(np.int64)
    return effective_count, common, mismatch, len(cell_codes)


def cache_path(config, chromosome):
    return os.path.join(
        config["output_dir"],
        "sufficient_statistics_cache",
        "{}_batch{}_minHC{}.npz".format(
            chromosome, config["batch_count"], config["min_hamming_count"]
        ),
    )


def build_or_load_sufficient_statistics(config, chromosome):
    path = cache_path(config, chromosome)
    if os.path.isfile(path) and not config["force_cache"]:
        with np.load(path) as cached:
            return {
                "starts": cached["starts"],
                "ends": cached["ends"],
                "pipeline_valid": cached["pipeline_valid"].astype(bool),
                "effective": cached["effective"],
                "common": cached["common"],
                "mismatch": cached["mismatch"],
                "cell_counts": cached["cell_counts"],
                "cache_path": path,
                "cached": True,
            }

    metrics, pats = discover_batches(
        config["run_dir"], chromosome, config["batch_count"]
    )
    first_metric = load_metric(metrics[0], chromosome)
    coordinates = first_metric[["chrom", "start", "end"]].copy()
    coordinate_values = coordinates[["start", "end"]].to_numpy()
    bin_count = len(coordinates)
    cpg_path = os.path.join(config["cpg_dir"], "CpG_{}.bed".format(chromosome))
    if not os.path.isfile(cpg_path):
        raise FileNotFoundError(cpg_path)
    cpg_to_bin = build_cpg_to_bin(cpg_path, coordinates)

    shape = (bin_count, config["batch_count"])
    pipeline_valid = np.zeros(shape, dtype=bool)
    effective = np.zeros(shape, dtype=np.int32)
    common = np.zeros(shape, dtype=np.int64)
    mismatch = np.zeros(shape, dtype=np.int64)
    cell_counts = np.zeros(config["batch_count"], dtype=np.int16)
    for batch_id in range(config["batch_count"]):
        metric = first_metric if batch_id == 0 else load_metric(metrics[batch_id], chromosome)
        if not np.array_equal(metric[["start", "end"]].to_numpy(), coordinate_values):
            raise ValueError(
                "Metric coordinates differ for {} batch {}".format(chromosome, batch_id)
            )
        hamming_count = pd.to_numeric(
            metric["HammingCount"], errors="coerce"
        ).fillna(0).to_numpy(dtype=np.int64)
        pipeline_valid[:, batch_id] = hamming_count >= config["min_hamming_count"]
        batch_effective, batch_common, batch_mismatch, cell_count = (
            scan_pat_sufficient_statistics(pats[batch_id], cpg_to_bin, bin_count)
        )
        effective[:, batch_id] = batch_effective
        common[:, batch_id] = batch_common
        mismatch[:, batch_id] = batch_mismatch
        cell_counts[batch_id] = cell_count
        print(
            "[{}] batch {}/{} cells={} metric_valid_bins={:,}".format(
                chromosome,
                batch_id + 1,
                config["batch_count"],
                cell_count,
                int(np.count_nonzero(pipeline_valid[:, batch_id])),
            ),
            flush=True,
        )

    if int(cell_counts.min()) < 1 or int(cell_counts.max()) > 64:
        raise ValueError("Unexpected batch cell counts: {}".format(cell_counts.tolist()))
    atomic_savez(
        path,
        starts=coordinates["start"].to_numpy(dtype=np.int64),
        ends=coordinates["end"].to_numpy(dtype=np.int64),
        pipeline_valid=pipeline_valid,
        effective=effective,
        common=common,
        mismatch=mismatch,
        cell_counts=cell_counts,
    )
    return {
        "starts": coordinates["start"].to_numpy(dtype=np.int64),
        "ends": coordinates["end"].to_numpy(dtype=np.int64),
        "pipeline_valid": pipeline_valid,
        "effective": effective,
        "common": common,
        "mismatch": mismatch,
        "cell_counts": cell_counts,
        "cache_path": path,
        "cached": False,
    }


def stable_standard_deviation(values, valid):
    count = valid.sum(axis=1).astype(np.int16)
    safe = np.where(valid, values, 0.0)
    means = safe.sum(axis=1) / count.astype(np.float64)
    centered = np.where(valid, values - means[:, None], 0.0)
    variance = np.square(centered).sum(axis=1) / (count - 1).astype(np.float64)
    return np.sqrt(np.maximum(variance, 0.0))


def process_chromosome(config, chromosome):
    started = time.time()
    statistics = build_or_load_sufficient_statistics(config, chromosome)
    common = statistics["common"]
    mismatch = statistics["mismatch"]
    effective = statistics["effective"]
    valid = statistics["pipeline_valid"] & (common >= config["rarefaction_k"])
    valid_batch_n = valid.sum(axis=1)
    eligible = valid_batch_n >= config["min_valid_batches"]
    if not np.any(eligible):
        return chromosome, None, None, 0, time.time() - started, statistics["cached"]

    starts = statistics["starts"][eligible]
    ends = statistics["ends"][eligible]
    common = common[eligible]
    mismatch = mismatch[eligible]
    effective = effective[eligible]
    valid = valid[eligible]
    valid_batch_n = valid_batch_n[eligible]
    if np.any(mismatch[valid] < 0) or np.any(mismatch[valid] > common[valid]):
        raise ValueError("Invalid mismatch/common counts for {}".format(chromosome))

    repeat_scores = np.empty((len(starts), config["repeats"]), dtype=np.float32)
    chromosome_seed = (
        config["seed"]
        + 100003 * chromosome_seed_index(chromosome)
        + 1009 * config["rarefaction_k"]
    )
    for repeat in range(config["repeats"]):
        rng = np.random.RandomState(chromosome_seed + repeat)
        sampled = np.zeros(common.shape, dtype=np.float64)
        sampled[valid] = rng.hypergeometric(
            mismatch[valid],
            common[valid] - mismatch[valid],
            config["rarefaction_k"],
        )
        transformed = np.full(common.shape, np.nan, dtype=np.float64)
        transformed[valid] = 2.0 * np.arcsin(
            np.sqrt(
                (sampled[valid] + 0.375)
                / (config["rarefaction_k"] + 0.75)
            )
        )
        repeat_scores[:, repeat] = stable_standard_deviation(
            transformed, valid
        ).astype(np.float32)

    valid_common = np.where(valid, common.astype(np.float64), np.nan)
    valid_effective = np.where(valid, effective.astype(np.float64), np.nan)
    summary = pd.DataFrame(
        {
            "chrom": chromosome,
            "start": starts,
            "end": ends,
            "K": config["rarefaction_k"],
            "valid_batch_n": valid_batch_n.astype(int),
            "valid_batch_fraction": valid_batch_n / float(config["batch_count"]),
            "median_total_common_comparisons": np.nanmedian(valid_common, axis=1),
            "mean_total_common_comparisons": np.nanmean(valid_common, axis=1),
            "median_unique_effective_cpg": np.nanmedian(valid_effective, axis=1),
            "mean_unique_effective_cpg": np.nanmean(valid_effective, axis=1),
            "rarefied_hamdist_sd_q10": np.quantile(repeat_scores, 0.10, axis=1),
            "rarefied_hamdist_sd_median": np.median(repeat_scores, axis=1),
            "rarefied_hamdist_sd_q90": np.quantile(repeat_scores, 0.90, axis=1),
            "rarefied_hamdist_sd_repeat_sd": np.std(repeat_scores, axis=1, ddof=1),
        }
    )
    summary["rarefied_hamdist_excess_sd_q10"] = np.sqrt(
        np.maximum(
            np.square(summary["rarefied_hamdist_sd_q10"])
            - 1.0 / config["rarefaction_k"],
            0.0,
        )
    )
    per_chrom_dir = os.path.join(config["output_dir"], "per_chromosome")
    ensure_directory(per_chrom_dir)
    summary_path = os.path.join(
        per_chrom_dir, "{}_rarefied_summary.tsv.gz".format(chromosome)
    )
    score_path = os.path.join(
        per_chrom_dir, "{}_repeat_scores.npy".format(chromosome)
    )
    atomic_to_csv(summary, summary_path, compression="gzip")
    temporary_score = score_path + ".tmp.npy"
    np.save(temporary_score, repeat_scores)
    os.replace(temporary_score, score_path)
    return (
        chromosome,
        summary_path,
        score_path,
        len(summary),
        time.time() - started,
        statistics["cached"],
    )


def add_global_stability(summary, repeat_scores, top_n):
    eligible_n, repeats = repeat_scores.shape
    if eligible_n < top_n:
        raise ValueError(
            "Only {} eligible bins are available for Top{}".format(eligible_n, top_n)
        )
    rank_matrix = np.empty((eligible_n, repeats), dtype=np.int32)
    top_hits = np.zeros(eligible_n, dtype=np.int16)
    for repeat in range(repeats):
        order = np.argsort(-repeat_scores[:, repeat], kind="mergesort")
        ranks = np.empty(eligible_n, dtype=np.int32)
        ranks[order] = np.arange(1, eligible_n + 1, dtype=np.int32)
        rank_matrix[:, repeat] = ranks
        top_hits[order[:top_n]] += 1
    result = summary.copy()
    result["top{}_repeat_frequency".format(top_n)] = top_hits / float(repeats)
    result["median_repeat_rank"] = np.median(rank_matrix, axis=1)
    result["repeat_rank_q90"] = np.quantile(rank_matrix, 0.90, axis=1)
    return result


def correlation_rows(frame):
    score_columns = [
        "rarefied_hamdist_sd_q10",
        "rarefied_hamdist_sd_median",
        "rarefied_hamdist_sd_repeat_sd",
    ]
    depth_columns = [
        "median_total_common_comparisons",
        "mean_total_common_comparisons",
        "median_unique_effective_cpg",
        "mean_unique_effective_cpg",
    ]
    rows = []
    for score in score_columns:
        for depth in depth_columns:
            selected = frame[[score, depth]].replace([np.inf, -np.inf], np.nan).dropna()
            if (
                len(selected) < 2
                or selected[score].nunique() < 2
                or selected[depth].nunique() < 2
            ):
                pearson = np.nan
                spearman = np.nan
            else:
                pearson = selected[score].corr(selected[depth])
                spearman = selected[score].rank().corr(selected[depth].rank())
            rows.append(
                {
                    "score": score,
                    "depth": depth,
                    "n": len(selected),
                    "pearson_r": pearson,
                    "spearman_r": spearman,
                }
            )
    return pd.DataFrame(rows)


def validate_args(args):
    if args.workers < 1 or args.batch_count < 2:
        raise ValueError("--workers must be >=1 and --batch-count must be >=2")
    if args.min_hamming_count < 0 or args.rarefaction_k < 2:
        raise ValueError("--min-hamming-count must be >=0 and --rarefaction-k >=2")
    if not 2 <= args.min_valid_batches <= args.batch_count:
        raise ValueError("--min-valid-batches must be in [2, batch-count]")
    if args.repeats < 2 or args.top_n < 1:
        raise ValueError("--repeats must be >=2 and --top-n must be positive")
    for path in (args.run_dir, args.cpg_dir):
        if not os.path.isdir(path):
            raise FileNotFoundError(path)


def main():
    args = parse_args()
    validate_args(args)
    ensure_directory(args.output_dir)
    chromosomes = (
        sorted(args.chromosomes, key=chromosome_key)
        if args.chromosomes
        else discover_chromosomes(args.run_dir)
    )
    config = vars(args).copy()
    config["chromosomes"] = chromosomes
    started = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(chromosomes))) as executor:
        futures = {
            executor.submit(process_chromosome, config, chromosome): chromosome
            for chromosome in chromosomes
        }
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(
                "[DONE] {} eligible={:,} elapsed={:.1f}s cached_statistics={}".format(
                    result[0], result[3], result[4], result[5]
                ),
                flush=True,
            )

    results.sort(key=lambda item: chromosome_key(item[0]))
    nonempty = [result for result in results if result[1] is not None]
    if not nonempty:
        raise ValueError("No bins passed the valid-batch requirement")
    summaries = [pd.read_csv(result[1], sep="\t") for result in nonempty]
    scores = [np.load(result[2]) for result in nonempty]
    summary = pd.concat(summaries, ignore_index=True)
    repeat_scores = np.concatenate(scores, axis=0)
    summary = add_global_stability(summary, repeat_scores, args.top_n)
    summary["_chrom_order"] = summary["chrom"].map(
        {chromosome: index for index, chromosome in enumerate(chromosomes)}
    )
    summary = summary.sort_values(["_chrom_order", "start", "end"]).drop(
        columns="_chrom_order"
    ).reset_index(drop=True)
    all_path = os.path.join(args.output_dir, "rarefied_pat_hamdist_all_eligible.tsv.gz")
    atomic_to_csv(summary, all_path, compression="gzip")

    frequency_column = "top{}_repeat_frequency".format(args.top_n)
    ranked = summary.sort_values(
        [
            "rarefied_hamdist_sd_q10",
            frequency_column,
            "rarefied_hamdist_sd_median",
            "rarefied_hamdist_sd_repeat_sd",
            "chrom",
            "start",
        ],
        ascending=[False, False, False, True, True, True],
    ).head(args.top_n).copy()
    ranked.insert(0, "SD_rank", np.arange(1, len(ranked) + 1, dtype=int))
    ranked.insert(1, "SD", ranked["rarefied_hamdist_sd_q10"].to_numpy())
    ranked_path = os.path.join(
        args.output_dir,
        "SD_rarefied_pat_hamdist_top{}_ranked.tsv".format(args.top_n),
    )
    atomic_to_csv(ranked, ranked_path)

    bed = ranked[["chrom", "start", "end"]].copy()
    bed["_chrom_order"] = bed["chrom"].map(
        {chromosome: index for index, chromosome in enumerate(chromosomes)}
    )
    bed = bed.sort_values(["_chrom_order", "start", "end"]).drop(
        columns="_chrom_order"
    )
    bed_path = os.path.join(
        args.output_dir, "SD_rarefied_pat_hamdist_top{}.bed".format(args.top_n)
    )
    atomic_to_csv(bed, bed_path, header=False)
    atomic_to_csv(
        correlation_rows(summary),
        os.path.join(args.output_dir, "score_depth_correlations.tsv"),
    )

    elapsed = time.time() - started
    parameters = config.copy()
    parameters.update(
        {
            "elapsed_seconds": elapsed,
            "eligible_bin_count": len(summary),
            "selected_bin_count": len(ranked),
            "ranked_output": ranked_path,
            "bed_output": bed_path,
            "score_definition": "Q10 of repeated across-batch SD after K-comparison rarefaction",
        }
    )
    with open(os.path.join(args.output_dir, "run_parameters.json"), "w") as handle:
        json.dump(parameters, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        "[FINISHED] eligible={:,} selected={:,} elapsed={:.1f}s output={}".format(
            len(summary), len(ranked), elapsed, ranked_path
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
