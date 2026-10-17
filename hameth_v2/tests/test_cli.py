import os
import subprocess
import sys
import tempfile
import unittest

from hameth_v2 import __version__
from hameth_v2.cli import build_parser, min_valid_batches


class CliTests(unittest.TestCase):
    def test_version(self):
        self.assertEqual(__version__, "2.0.0")

    def test_default_valid_batch_rule(self):
        self.assertEqual(min_valid_batches(40), 36)
        self.assertEqual(min_valid_batches(20), 18)

    def test_subcommands_exist(self):
        parser = build_parser()
        choices = next(
            action.choices
            for action in parser._actions
            if isinstance(action, __import__("argparse")._SubParsersAction)
        )
        self.assertEqual(
            set(choices), {"prepare", "hamdist", "select-bins", "matrix"}
        )

    def test_balanced_selection_smoke(self):
        with tempfile.TemporaryDirectory() as root:
            run_dir = os.path.join(root, "run")
            cpg_dir = os.path.join(root, "cpg")
            output_dir = os.path.join(root, "selection")
            metric_dir = os.path.join(run_dir, "middle_res", "metrics", "chr1")
            merge_dir = os.path.join(run_dir, "middle_res", "mergepats")
            os.makedirs(metric_dir)
            os.makedirs(merge_dir)
            os.makedirs(cpg_dir)

            with open(os.path.join(cpg_dir, "CpG_chr1.bed"), "w") as handle:
                handle.write("chr1\t110\t1\nchr1\t120\t2\n")
                handle.write("chr1\t210\t3\nchr1\t220\t4\n")

            patterns = [
                ["CC", "CC", "TT", "TT"],
                ["CC", "CT", "TC", "TT"],
                ["CC", "CC", "CC", "TT"],
                ["CC", "CT", "CT", "TT"],
            ]
            cells = ["cellA", "cellB", "cellC", "cellD"]
            for batch_id, batch_patterns in enumerate(patterns):
                metric_path = os.path.join(
                    metric_dir,
                    "chr1_merge_{}_MeConcord.txt".format(batch_id),
                )
                with open(metric_path, "w") as handle:
                    handle.write("chrom\tstart\tend\tHammingCount\n")
                    handle.write("chr1\t100\t199\t8\n")
                    handle.write("chr1\t200\t299\t8\n")
                pat_path = os.path.join(
                    merge_dir, "chr1_merge_{}.pat".format(batch_id)
                )
                with open(pat_path, "w") as handle:
                    for cell, pattern in zip(cells, batch_patterns):
                        handle.write("chr1\t1\t{}\t1\t{}\n".format(pattern, cell))
                        handle.write("chr1\t3\tCT\t1\t{}\n".format(cell))

            command = [
                sys.executable,
                "-m",
                "hameth_v2",
                "select-bins",
                "--run-dir",
                run_dir,
                "--cpg-dir",
                cpg_dir,
                "--output-dir",
                output_dir,
                "--batch-count",
                "4",
                "--min-hamming-count",
                "5",
                "--rarefaction-k",
                "4",
                "--min-valid-batches",
                "3",
                "--repeats",
                "5",
                "--top-n",
                "2",
                "--selection-seed",
                "17",
                "--workers",
                "1",
                "--chromosomes",
                "chr1",
            ]
            completed = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            self.assertEqual(completed.returncode, 0, msg=completed.stderr)
            self.assertTrue(
                os.path.isfile(
                    os.path.join(
                        output_dir,
                        "SD_rarefied_pat_hamdist_top2_ranked.tsv",
                    )
                )
            )


if __name__ == "__main__":
    unittest.main()
