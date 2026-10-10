"""Terminal execution must leave notebook algorithms and production guards intact."""

import ast
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import run_regression_bootstrap as runner


class RunnerConfigurationTests(unittest.TestCase):
    def test_formal_defaults_and_smoke_are_explicit(self):
        formal = runner.parse_args([])
        self.assertEqual((formal.n, formal.p, formal.N, formal.B), (1000, 1000, 1000, 500))
        self.assertEqual(formal.backend, "gpu")
        self.assertEqual(formal.gpu_devices, (0,))
        self.assertEqual(formal.methods[0], "Ridge")
        smoke = runner.parse_args(["--smoke", "--backend", "cpu"])
        self.assertEqual((smoke.n, smoke.p, smoke.N, smoke.B), (40, 40, 2, 4))
        invalid = (["--n", "40"], ["--methods", "SC", "SC"], ["--resume-from", "."],
                   ["--gpu-devices", "0,0"], ["--N", "0"], ["--smoke", "--p", "19"],
                   ["--smoke", "--n", "1000", "--p", "1000", "--N", "1000", "--B", "500"])
        for argv in invalid:
            with self.subTest(argv=argv), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                runner.parse_args(argv)

    def test_only_whitelisted_module_assignments_change(self):
        notebook = {"cells": [{"cell_type": "code", "source": "N = 1000\ndef f():\n    N = 99\n    return N\n"}]}
        cells = runner.prepare_cells(notebook, {"N": 2})
        self.assertEqual(ast.literal_eval(cells[0][1].body[0].value), 2)
        original = ast.parse(notebook["cells"][0]["source"]).body[1]
        self.assertEqual(ast.dump(cells[0][1].body[1]), ast.dump(original))
        with self.assertRaises(ValueError):
            runner.prepare_cells(notebook, {"kappa": 1.0})
        with self.assertRaises(ValueError):
            runner.prepare_cells(notebook, {"B": 2})
        notebook["cells"][0]["source"] += "N = 8\n"
        with self.assertRaises(ValueError):
            runner.prepare_cells(notebook, {"N": 2})

    def test_all_notebook_nonconfiguration_statements_are_identical(self):
        overrides = {"N": 2, "B": 4, "n": 40, "p": 40, "backend": "cpu",
                     "gpu_devices": (0,), "cpu_workers": 1, "blas_threads": 1,
                     "batch_size": 2, "output_root": "example", "resume_from": None}
        for method, name in runner.NOTEBOOKS.items():
            with self.subTest(method=method):
                notebook = json.loads((runner.ROOT / name).read_text(encoding="utf-8"))
                for index, transformed in runner.prepare_cells(notebook, overrides):
                    source = notebook["cells"][index]["source"]
                    original = ast.parse("".join(source) if isinstance(source, list) else source)
                    self.assertEqual(len(original.body), len(transformed.body))
                    for before, after in zip(original.body, transformed.body):
                        changed = (isinstance(before, ast.Assign) and len(before.targets) == 1
                                   and isinstance(before.targets[0], ast.Name)
                                   and before.targets[0].id in overrides)
                        if not changed:
                            self.assertEqual(ast.dump(before), ast.dump(after))
                    compile(transformed, name, "exec")


class RunnerIntegrationTests(unittest.TestCase):
    def test_interrupted_worker_resumes_only_uncommitted_replications(self):
        from regression_results import BootstrapRun

        with tempfile.TemporaryDirectory() as temporary:
            batch = Path(temporary)
            arguments = ["--smoke", "--backend", "cpu", "--methods", "Ridge",
                         "--cpu-workers", "1", "--blas-threads", "1", "--batch-size", "2",
                         "--_worker-method", "Ridge", "--_batch-directory", temporary]
            original_record = BootstrapRun.record

            def interrupted_record(instance, *args, **kwargs):
                original_record(instance, *args, **kwargs)
                raise RuntimeError("test interruption after a committed replication")

            with contextlib.redirect_stdout(io.StringIO()), mock.patch.object(BootstrapRun, "record", interrupted_record):
                with self.assertRaisesRegex(RuntimeError, "test interruption"):
                    runner.run_worker(runner.parse_args(arguments))
            status = json.loads((batch / "Ridge.json").read_text(encoding="utf-8"))
            directory = Path(status["run_directory"])
            saved = directory / "replications" / "rep_000000.npz"
            first_bytes = saved.read_bytes()
            metadata = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["status"], "failed")
            self.assertEqual(metadata["completed_replications"], 1)
            with contextlib.redirect_stdout(io.StringIO()):
                result = runner.run_worker(runner.parse_args(arguments + ["--resume-from", str(directory)]))
            self.assertEqual(result, 0)
            self.assertEqual(first_bytes, saved.read_bytes())
            self.assertEqual(len(list((directory / "replications").glob("*.npz"))), 2)
            metadata = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["status"], "complete")

    def test_isolated_cpu_smoke_archive_and_completed_resume(self):
        with tempfile.TemporaryDirectory() as temporary:
            command = [sys.executable, str(runner.ROOT / "run_regression_bootstrap.py"),
                       "--smoke", "--backend", "cpu", "--methods", "Ridge", "SC",
                       "--cpu-workers", "1", "--blas-threads", "1", "--batch-size", "2",
                       "--output-root", temporary]
            first = subprocess.run(command, cwd=temporary, capture_output=True, text=True, timeout=90)
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            manifests = list(Path(temporary).glob("*/manifest.json"))
            self.assertEqual(len(manifests), 1)
            manifest = json.loads(manifests[0].read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(set(manifest["runs"]), {"Ridge", "SC"})
            self.assertNotEqual(manifest["runs"]["Ridge"]["pid"], manifest["runs"]["SC"]["pid"])
            ridge_directory = Path(manifest["runs"]["Ridge"]["run_directory"])
            for report in manifest["runs"].values():
                self.assertEqual(report["completed_replications"], 2)
                directory = Path(report["run_directory"])
                metadata = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
                self.assertFalse(metadata["is_manuscript_configuration"])
                self.assertEqual(metadata["config"]["execution"]["backend"], "cpu")
                self.assertEqual(metadata["config"]["terminal_runner"]["version"], runner.VERSION)
                archive = directory / "terminal_runner"
                self.assertTrue((archive / "effective_code.py").exists())
                self.assertTrue((archive / "runner.py").exists())
                self.assertIn("checkpoint 2/2 committed", Path(report["log"]).read_text(encoding="utf-8"))
            replications = {path.name: path.read_bytes() for path in (ridge_directory / "replications").glob("*.npz")}
            resume = command.copy()
            resume.remove("SC")
            resume += ["--resume-from", str(ridge_directory)]
            second = subprocess.run(resume, cwd=temporary, capture_output=True, text=True, timeout=90)
            self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
            self.assertEqual(replications, {path.name: path.read_bytes() for path in (ridge_directory / "replications").glob("*.npz")})
            changed = resume.copy()
            changed[changed.index("--batch-size") + 1] = "3"
            rejected = subprocess.run(changed, cwd=temporary, capture_output=True, text=True, timeout=90)
            self.assertNotEqual(rejected.returncode, 0)
            metadata = json.loads((ridge_directory / "metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["status"], "complete")


if __name__ == "__main__":
    unittest.main()
