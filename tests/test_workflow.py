import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = os.environ.get("NEXTFLOW", "nextflow")

MOCK_TOOL = """#!/usr/bin/env python3
import gzip
import json
from pathlib import Path
import sys

args = sys.argv[1:]
output = Path(args[args.index("-o") + 1])
with gzip.open(args[-1], "rt") as handle:
    assert handle.read().startswith("@read")
if Path(sys.argv[0]).name == "yak":
    output.write_text(json.dumps({"reads": args[-1]}))
else:
    if "-1" in args:
        for option in ["-1", "-2"]:
            data = json.loads(Path(args[args.index(option) + 1]).read_text())
            assert data["reads"].endswith(".fastq.gz")
        mode = "dip"
    else:
        mode = "bp"
    for hap in [1, 2]:
        Path(str(output) + f".{mode}.hap{hap}.p_ctg.gfa").write_text("S\\tcontig\\tACGT\\n")
print(json.dumps({"tool": Path(sys.argv[0]).name, "args": args}))
"""


@unittest.skipUnless(shutil.which(NEXTFLOW), "Nextflow is required")
@unittest.skipUnless(shutil.which("samtools") and shutil.which("bgzip"), "samtools and bgzip are required")
class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="hifi-trios-test-")
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.tools = self.directory / "tools"
        self.tools.mkdir()
        for tool in ["hifiasm", "yak"]:
            executable = self.tools / tool
            executable.write_text(MOCK_TOOL)
            executable.chmod(0o755)
        self.config = self.directory / "test.config"
        self.config.write_text("apptainer.enabled = false\ndocker.enabled = false\n")
        self.bam = self.directory / "input reads.bam"
        subprocess.run(
            ["samtools", "view", "-b", "-o", str(self.bam), "-"],
            input="@HD\tVN:1.6\tSO:unknown\nread\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
            text=True,
            check=True,
            capture_output=True,
        )

    def run_workflow(self, rows=None, arguments=(), success=True, header=None, before_tasks=True):
        run_dir = self.directory / f"run-{len(list(self.directory.glob('run-*')))}"
        run_dir.mkdir()
        command = [
            NEXTFLOW, "run", str(ROOT / "main.nf"), "-ansi-log", "false",
            "-c", str(self.config), "-work-dir", str(run_dir / "work"),
            "-with-trace", str(run_dir / "trace.tsv"), "--cpus", "1", "--memory", "1G",
        ]
        if rows is not None:
            design = run_dir / "design.csv"
            with design.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(header or ["BAM", "sample_id", "sample_family", "role"])
                writer.writerows(rows)
            command += ["--design", str(design)]
        command += list(arguments)
        env = {**os.environ, "PATH": str(self.tools) + os.pathsep + os.environ["PATH"]}
        env["NXF_OFFLINE"] = "true"
        env.setdefault("NXF_SYNTAX_PARSER", "v2")
        result = subprocess.run(command, cwd=run_dir, env=env, capture_output=True, text=True, timeout=180)
        output = result.stdout + result.stderr
        if success:
            self.assertEqual(result.returncode, 0, output)
        else:
            self.assertNotEqual(result.returncode, 0, output)
            if before_tasks:
                self.assertFalse(list((run_dir / "work").glob("*/*/.command.sh")), output)
        return run_dir, output

    def rows(self, samples):
        return [[str(self.bam), sample, family, role] for sample, family, role in samples]

    def test_families_siblings_parent_order_and_sample_outputs(self):
        samples = [
            ("child_b", "family2", "proband"),
            ("mother_a", "family1", "mat"),
            ("sibling_a", "family1", "proband"),
            ("father_b", "family2", "pat"),
            ("child_a", "family1", "proband"),
            ("solo_mother", "family3", "mat"),
            ("mother_b", "family2", "mat"),
            ("father_a", "family1", "pat"),
        ]
        run_dir, _ = self.run_workflow(self.rows(samples))
        outputs = list((run_dir / "assemblies").glob("*/*.fa.gz"))
        self.assertEqual(len(outputs), 2 * len(samples))
        for sample, _, _ in samples:
            for hap in [1, 2]:
                path = run_dir / "assemblies" / f"{sample}_diploid" / f"{sample}_hap{hap}.fa.gz"
                with gzip.open(path, "rt") as handle:
                    self.assertEqual(handle.read(), f">{sample}#{hap}#contig\nACGT\n")
        calls = []
        for log in (run_dir / "work").glob("*/*/.command.log"):
            calls.extend(json.loads(line) for line in log.read_text().splitlines() if line.startswith('{"tool":'))
        self.assertEqual(sum(call["tool"] == "yak" for call in calls), 5)
        assemblies = [call["args"] for call in calls if call["tool"] == "hifiasm"]
        self.assertEqual(len(assemblies), len(samples))
        expected = {
            "child_a": ("father_a.yak", "mother_a.yak"),
            "sibling_a": ("father_a.yak", "mother_a.yak"),
            "child_b": ("father_b.yak", "mother_b.yak"),
        }
        for args in assemblies:
            sample = args[-1].removesuffix(".fastq.gz")
            if sample in expected:
                self.assertEqual((args[args.index("-1") + 1], args[args.index("-2") + 1]), expected[sample])
            else:
                self.assertNotIn("-1", args)
                self.assertNotIn("-2", args)
        fastqs = list((run_dir / "work").glob("*/*/*.fastq.gz"))
        actual_fastqs = [path for path in fastqs if not path.is_symlink()]
        self.assertEqual(len(actual_fastqs), len(samples))
        for path in actual_fastqs:
            with gzip.open(path, "rt") as handle:
                self.assertEqual(handle.read(), "@read\nACGT\n+\nIIII\n")

    def test_parent_only(self):
        run_dir, _ = self.run_workflow(self.rows([("mother", "family1", "mat")]))
        self.assertEqual(len(list((run_dir / "assemblies").glob("*/*.fa.gz"))), 2)

    def test_pansn_existing_assemblies(self):
        hap1 = self.directory / "hap1.fa"
        hap2 = self.directory / "hap2.fa.gz"
        sequence = ">contig1 description\nACGT\nTGCA\n>contig2\nNNAA\n"
        hap1.write_text(sequence)
        with gzip.open(hap2, "wt") as handle:
            handle.write(sequence)
        assemblies = self.directory / "assemblies.csv"
        with assemblies.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sample", "hap1", "hap2"])
            writer.writerow(["existing", hap1, hap2])
        run_dir, _ = self.run_workflow(arguments=["--assemblies", str(assemblies)])
        for hap in [1, 2]:
            output = run_dir / "assemblies" / "existing_diploid" / f"existing_hap{hap}.fa.gz"
            with gzip.open(output, "rt") as handle:
                self.assertEqual(
                    handle.read(),
                    f">existing#{hap}#contig1 description\nACGT\nTGCA\n>existing#{hap}#contig2\nNNAA\n",
                )
        self.assertEqual(hap1.read_text(), sequence)
        with gzip.open(hap2, "rt") as handle:
            self.assertEqual(handle.read(), sequence)
        hap1.write_text(">already#1#contig\nACGT\n")
        run_dir, output = self.run_workflow(
            arguments=["--assemblies", str(assemblies)], success=False, before_tasks=False,
        )
        self.assertIn("Contig identifiers must not contain #", output)
        self.assertFalse(list((run_dir / "assemblies").glob("*/*.fa.gz")))

    def test_invalid_bam_conversion_fails(self):
        bad_bam = self.directory / "invalid.bam"
        bad_bam.write_text("not a BAM")
        run_dir, output = self.run_workflow(
            [[str(bad_bam), "mother", "family1", "mat"]], success=False, before_tasks=False,
        )
        self.assertIn("bam_to_fastq", output)
        self.assertIn("terminated with an error exit status", output)
        self.assertFalse(list((run_dir / "assemblies").glob("*/*.fa.gz")))

    def test_invalid_designs_fail_before_tasks(self):
        valid = [("mother", "family1", "mat"), ("father", "family1", "pat"), ("child", "family1", "proband")]
        cases = [
            ([], "at least one sample"),
            ([("child", "family1", "proband")], "exactly one 'mat'"),
            ([valid[0], valid[2]], "exactly one 'pat'"),
            (valid + [("mother2", "family1", "mat")], "exactly one 'mat'"),
            (valid + [("child", "family1", "proband")], "Duplicate sample_id"),
            ([("sample", "family1", "other")], "Invalid role"),
            ([("../sample", "family1", "mat")], "Invalid sample or family identifier"),
            ([("sample", "", "mat")], "non-empty 'sample_family'"),
        ]
        for samples, message in cases:
            with self.subTest(samples=samples):
                _, output = self.run_workflow(self.rows(samples), success=False)
                self.assertIn(message, output)
        _, output = self.run_workflow(
            [[str(self.directory / "missing.bam"), "sample", "family1", "mat"]], success=False,
        )
        self.assertIn("missing.bam", output)
        _, output = self.run_workflow(self.rows(valid), header=["BAM", "sample_id", "sample_family", "wrong"], success=False)
        self.assertIn("non-empty 'role'", output)

    def test_variant_routing_and_existing_assemblies(self):
        reference = self.directory / "ref.fa"
        reference.write_text(">chr1\nACGT\n")
        samples = [("mother", "family1", "mat"), ("father", "family1", "pat"), ("child", "family1", "proband")]
        for mode, suffixes in [("dipcall", [".dip.vcf.gz"]), ("svim_asm", [".dip.vcf.gz", ".snvs.dip.vcf.gz"])]:
            with self.subTest(mode=mode):
                run_dir, _ = self.run_workflow(
                    self.rows(samples), arguments=["-stub-run", f"--{mode}", "true", "--ref", str(reference)],
                )
                self.assertEqual(
                    {path.name for path in (run_dir / "variants").iterdir()},
                    {sample + suffix for sample, _, _ in samples for suffix in suffixes},
                )
        assemblies = self.directory / "assemblies.csv"
        hap1 = self.directory / "hap1.fa"
        hap2 = self.directory / "hap2.fa"
        hap1.write_text(">contig\nACGT\n")
        hap2.write_text(">contig\nACGT\n")
        with assemblies.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sample", "hap1", "hap2"])
            writer.writerow(["existing", hap1, hap2])
        run_dir, _ = self.run_workflow(arguments=[
            "-stub-run", "--assemblies", str(assemblies), "--dipcall", "true", "--ref", str(reference),
        ])
        self.assertTrue((run_dir / "variants" / "existing.dip.vcf.gz").exists())


if __name__ == "__main__":
    unittest.main()
