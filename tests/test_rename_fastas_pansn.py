import csv
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "rename_fastas_pansn.sh"


class RenameFastasTests(unittest.TestCase):
    def test_fastas_run_concurrently_and_worker_failure_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "good.fa").write_text(">contig\nACGT\n")
            (root / "bad.fa").write_text("not a fasta\n")
            design = root / "samples.csv"
            design.write_text("sample,hap1,hap2\nsample,good.fa,bad.fa\n")
            tools = root / "tools"
            tools.mkdir()
            wrapper = tools / "gzip"
            wrapper.write_text(
                f"#!{sys.executable}\n"
                "import os, pathlib, sys, time\n"
                f"root = pathlib.Path({str(root)!r})\n"
                "if '-cdf' in sys.argv:\n"
                "    (root / (pathlib.Path(sys.argv[-1]).name + '.started')).touch()\n"
                "    deadline = time.monotonic() + 10\n"
                "    while len(list(root.glob('*.started'))) < 2:\n"
                "        if time.monotonic() > deadline:\n"
                "            sys.exit('FASTA jobs did not start concurrently')\n"
                "        time.sleep(0.02)\n"
                f"os.execv({shutil.which('gzip')!r}, ['gzip', *sys.argv[1:]])\n"
            )
            wrapper.chmod(0o755)
            output = root / "output"
            result = subprocess.run(
                ["bash", str(SCRIPT), str(design), str(output)],
                env={**os.environ, "PATH": str(tools) + os.pathsep + os.environ["PATH"]},
                capture_output=True, text=True, timeout=20,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Input contains no FASTA headers", result.stderr)
            self.assertNotIn("did not start concurrently", result.stderr)
            self.assertEqual([path.name for path in output.iterdir()], ["sample_hap1.fa.gz"])
            with gzip.open(output / "sample_hap1.fa.gz", "rt") as handle:
                self.assertEqual(handle.read(), ">sample#1#contig\nACGT\n")

    def test_plain_and_compressed_csv_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sequence = ">ctg1 description\nACGT\nTGCA\n>ctg2\nNNAA\n"
            (root / "plain, reads.fa").write_text(sequence)
            with gzip.open(root / "compressed.fa.gz", "wt") as handle:
                handle.write(sequence)
            design = root / "samples.csv"
            with design.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["sample", "hap1", "hap2"])
                writer.writerow(["child", "plain, reads.fa", "compressed.fa.gz"])
            output = root / "output assemblies"
            result = subprocess.run(
                ["bash", str(SCRIPT), str(design), str(output)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            for hap in (1, 2):
                with gzip.open(output / f"child_hap{hap}.fa.gz", "rt") as handle:
                    self.assertEqual(
                        handle.read(),
                        f">child#{hap}#ctg1 description\nACGT\nTGCA\n>child#{hap}#ctg2\nNNAA\n",
                    )
            self.assertEqual((root / "plain, reads.fa").read_text(), sequence)
            repeat = subprocess.run(
                ["bash", str(SCRIPT), str(design), str(output)],
                capture_output=True, text=True,
            )
            self.assertNotEqual(repeat.returncode, 0)
            self.assertIn("Refusing to overwrite", repeat.stderr)
            self.assertEqual(len(list(output.iterdir())), 2)

    def test_rejects_invalid_input(self):
        for contents in (">sample#1#contig\nACGT\n", "not a fasta\n"):
            with self.subTest(contents=contents), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                (root / "input.fa").write_text(contents)
                (root / "samples.csv").write_text("sample,hap1,hap2\nsample,input.fa,input.fa\n")
                output = root / "output"
                result = subprocess.run(
                    ["bash", str(SCRIPT), str(root / "samples.csv"), str(output)],
                    capture_output=True, text=True,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertTrue(result.stderr)
                self.assertEqual(list(output.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
