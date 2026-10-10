"""Independent scan provenance, per-outer records and validation-plan guards."""

import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

import run_regression_threshold_scan as runner
from regression_results import _array_hash, _digest


class ScanFixtures(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        rng = np.random.default_rng(42)
        self.X = rng.normal(size=(24, 24))
        self.beta = np.r_[np.repeat([2.0, -2.0, 1.0, -1.0], 5), np.zeros(4)]
        self.source = {"files": {"fixture.py": "# test fixture\n"}}
        self.config = {"stage": "calibration", "method": "SC", "percentiles": [0.0, 50.0, 100.0],
                       "n": 24, "p": 24, "N": 2, "B": 6, "sigma": 2.0, "confidence": 0.95,
                       "seed_root": 20261007, "design_seed": 24, "noise_model": "normal",
                       "threshold_grid": [0.001, 0.001, 1.0], "preprocess_blas_threads": 1,
                       "method_settings": {"selection": "oracle_rank", "rank_max": 24}}

    def values(self, index=0):
        rng = np.random.default_rng(index + 11)
        estimates = np.broadcast_to(self.beta[None, :, None], (3, 24, 2)).copy()
        scan = np.column_stack((0.001 * np.arange(1, 1001), np.ones((1000, 2))))
        return dict(estimates=estimates, sigma=np.ones((3, 2)), thresholds=np.full((3, 2), 0.2),
                    bootstrap_errors=rng.random((3, 6, 2)), selected_parameter=20,
                    diagnostics={"selected_count": 20, "at_iteration_cap": False, "tdp_satisfied": None},
                    threshold_scan=scan, percentiles=self.config["percentiles"])

    def old_run(self):
        path = self.root / "old-run"
        path.mkdir()
        np.savez_compressed(path / "data.npz", X=self.X, beta=self.beta)
        metadata = {"status": "complete", "identity": "fixture-input", "config": {"design_seed": 24},
                    "fingerprints": {"X": _array_hash(self.X), "beta": _array_hash(self.beta)}}
        (path / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
        return path


class ScanPersistenceTests(ScanFixtures):
    def test_resume_preserves_records_and_rejects_configuration_or_source_drift(self):
        run = runner.ScanRun(self.config, self.X, self.beta, self.source, self.root)
        with run:
            run.record(0, **self.values())
        self.assertEqual(run.metadata["status"], "incomplete")
        first = run.directory / "replications" / "rep_000000.npz"
        original = first.read_bytes()
        for config, source in ((dict(self.config, B=7), self.source),
                               (self.config, {"files": {"fixture.py": "changed\n"}})):
            with self.assertRaisesRegex(ValueError, "Resume requires identical"):
                runner.ScanRun(config, self.X, self.beta, source, self.root, run.directory)
        resumed = runner.ScanRun(self.config, self.X, self.beta, self.source, self.root, run.directory)
        with resumed:
            self.assertEqual(resumed.start_index, 1)
            resumed.record(1, **self.values(1))
        self.assertEqual(first.read_bytes(), original)
        metadata = runner.validate_run(run.directory)
        self.assertFalse(metadata["is_manuscript_configuration"])
        records = list(runner.read_scan_records(run.directory, metadata))
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["bootstrap_errors"].shape, (3, 6, 2))
        summary = runner.load_json(run.directory / "summary.json")
        self.assertEqual(summary["status"], "complete")
        self.assertEqual(len(summary["results"]), 3)
        self.assertEqual(summary["results"][0]["coverage_debiased"], 1.0)
        self.assertEqual(summary["results"][0]["mcse_coverage_debiased"], 0.0)

    def test_interrupted_record_commit_is_recovered_and_active_lock_is_respected(self):
        run = runner.ScanRun(self.config, self.X, self.beta, self.source, self.root)
        with self.assertRaisesRegex(RuntimeError, "test interruption"):
            with run:
                run.record(0, **self.values())
                with self.assertRaises(RuntimeError):
                    with runner.ScanRun(self.config, self.X, self.beta, self.source, self.root, run.directory):
                        pass
                raise RuntimeError("test interruption")
        self.assertEqual(run.metadata["status"], "failed")
        self.assertFalse((run.directory / ".active.lock").exists())
        resumed = runner.ScanRun(self.config, self.X, self.beta, self.source, self.root, run.directory)
        with resumed:
            self.assertEqual(resumed.start_index, 1)
            resumed.record(1, **self.values(1))
        self.assertEqual(runner.validate_run(run.directory)["status"], "complete")

    def test_corrupt_metrics_missing_index_and_nonfinite_data_are_rejected(self):
        run = runner.ScanRun(self.config, self.X, self.beta, self.source, self.root)
        with run:
            bad = self.values()
            bad["bootstrap_errors"][0, 0, 0] = np.nan
            with self.assertRaises(ValueError):
                run.record(0, **bad)
            run.record(0, **self.values())
            run.record(1, **self.values(1))
        path = run.directory / "replications" / "rep_000001.npz"
        with np.load(path, allow_pickle=False) as saved:
            arrays = {name: saved[name].copy() for name in saved.files}
        arrays["metrics"][0, 0] = 0.0
        np.savez_compressed(path, **arrays)
        with self.assertRaisesRegex(ValueError, "Corrupt derived metrics"):
            list(runner.read_scan_records(run.directory))
        path.rename(path.with_name("rep_000002.npz"))
        with self.assertRaisesRegex(ValueError, "contiguous prefix"):
            list(runner.read_scan_records(run.directory))

    def test_validation_plan_is_frozen_independent_and_fingerprinted(self):
        plan = {"version": "regression-threshold-selection-v1", "status": "ready_for_validation",
                "calibration_seed_root": 20261007, "validation_seed_root": 20261008,
                "validation_percentiles": {"SC": [0.0, 50.0, 100.0]}}
        runner.validate_plan(plan, "SC", self.config["percentiles"], 20261008)
        for seed, percentiles in ((20261007, self.config["percentiles"]), (20261008, [50.0])):
            with self.assertRaises(ValueError):
                runner.validate_plan(plan, "SC", percentiles, seed)
        config = dict(self.config, stage="validation", seed_root=20261008, selection_plan_sha256=_digest(plan))
        run = runner.ScanRun(config, self.X, self.beta, self.source, self.root, selection_plan=plan)
        with run:
            run.record(0, **self.values())
            run.record(1, **self.values(1))
        runner.validate_run(run.directory)
        plan["validation_seed_root"] = 20261009
        (run.directory / "selection_plan.json").write_text(json.dumps(plan), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Selection plan fingerprint"):
            runner.validate_run(run.directory)


class ScanRunnerTests(ScanFixtures):
    def test_ridge_spacing_guard_and_child_arguments(self):
        common = ["--stage", "calibration", "--data-from", str(self.root)]
        default = runner.parse_args([*common, "--methods", "Ridge"])
        self.assertIsNone(default.ridge_alpha_step)
        self.assertNotIn("--ridge-alpha-step", runner.child_arguments(default, "Ridge", self.root))
        requested = runner.parse_args([*common, "--methods", "Ridge", "--ridge-alpha-step", "0.5"])
        child = runner.parse_args(runner.child_arguments(requested, "Ridge", self.root)[3:])
        self.assertEqual(child.ridge_alpha_step, 0.5)
        self.assertEqual(child._worker_method, "Ridge")
        invalid = [["--methods", "Ridge", "--ridge-alpha-step", value]
                   for value in ("0", "-0.5", "nan", "inf")]
        invalid += [["--methods", "SC", "--ridge-alpha-step", "0.5"],
                    ["--methods", "Ridge", "SC", "--ridge-alpha-step", "0.5"],
                    ["--methods", "Ridge", "--ridge-alpha-step", "0.5",
                     "--_worker-method", "SC", "--_batch-directory", str(self.root)]]
        for extra in invalid:
            with self.subTest(extra=extra), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                runner.parse_args([*common, *extra])
        settings = dict(C=200000, m=100, alpha_divisor=1000)
        self.assertEqual(runner.ridge_grid_stride(0.5, settings), 500)
        self.assertEqual(settings["m"], 100)
        for value in (0.0005, 0.3, 201.0, float("inf")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                runner.ridge_grid_stride(value, settings)

    def test_terminal_ridge_spacing_reaches_worker_and_archive(self):
        # Preserve n=1000's production alpha scale while keeping the smoke design narrow.
        self.X = np.random.default_rng(43).normal(size=(1000, 24))
        old = self.old_run()
        output = self.root / "ridge-output"
        command = [sys.executable, str(runner.ROOT / "run_regression_threshold_scan.py"),
                   "--stage", "calibration", "--methods", "Ridge", "--percentiles", "0",
                   "--ridge-alpha-step", "0.5", "--data-from", str(old), "--output-root", str(output),
                   "--N", "1", "--B", "4", "--smoke-test", "--backend", "cpu",
                   "--cpu-workers", "1", "--blas-threads", "1", "--batch-size", "2"]
        completed = subprocess.run(command, cwd=self.root, capture_output=True, text=True, timeout=120)
        logs = "\n".join(path.read_text(encoding="utf-8") for path in output.glob("*/*.log"))
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr + logs)
        manifest_path, = output.glob("*/manifest.json")
        manifest = runner.load_json(manifest_path)
        self.assertEqual(manifest["status"], "complete")
        self.assertEqual(manifest["configuration"]["ridge_alpha_step"], 0.5)
        directory = Path(manifest["runs"]["Ridge"]["run_directory"])
        metadata = runner.validate_run(directory)
        self.assertEqual(metadata["config"]["method_settings"],
                         dict(selection="oracle_grid", C=200000, m=500, alpha_divisor=1000))
        record, = runner.read_scan_records(directory, metadata)
        self.assertEqual(float(record["selected_parameter"]) % 0.5, 0.0)

    def test_cli_defaults_and_rejections(self):
        common = ["--data-from", str(self.root)]
        baseline = runner.parse_args(["--stage", "baseline", *common])
        self.assertEqual(baseline.methods, ["HBF", "FAR"])
        self.assertEqual(baseline.percentiles, [50.0])
        self.assertEqual((baseline.N, baseline.B, baseline.seed_root), (1000, 500, 20261006))
        coarse = runner.parse_args(["--stage", "calibration", *common])
        self.assertEqual(coarse.percentiles, list(range(0, 101, 10)))
        fine = runner.parse_args(["--stage", "refinement", *common])
        self.assertEqual(fine.percentiles, list(range(5, 100, 10)))
        self.assertEqual(coarse.seed_root, fine.seed_root)
        for extra in (["--percentiles", "50", "50"], ["--methods", "SC", "SC"],
                      ["--percentiles", "101"], ["--seed-root", "-1"]):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                runner.parse_args(["--stage", "calibration", *common, *extra])
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            runner.parse_args(["--stage", "validation", *common])

    def test_worker_selects_once_per_outer_and_saves_all_percentiles(self):
        from regression_spectral import SpectralDesign

        old = self.old_run()
        batch = self.root / "worker"
        batch.mkdir()
        args = runner.parse_args(["--stage", "calibration", "--methods", "Ridge", "--percentiles", "0", "50", "100",
                                  "--N", "2", "--B", "6", "--data-from", str(old), "--smoke-test", "--backend", "cpu",
                                  "--cpu-workers", "1", "--blas-threads", "1", "--batch-size", "3",
                                  "--_worker-method", "Ridge", "--_batch-directory", str(batch)])
        original = SpectralDesign.select
        calls = []

        def select(design, *args, **kwargs):
            calls.append(1)
            return original(design, *args, **kwargs)

        with contextlib.redirect_stdout(io.StringIO()), mock.patch.object(SpectralDesign, "select", select):
            self.assertEqual(runner.run_worker(args), 0)
        self.assertEqual(len(calls), 2)
        report = runner.load_json(batch / "Ridge.json")
        self.assertEqual(report["status"], "complete")
        directory = Path(report["run_directory"])
        metadata = runner.validate_run(directory)
        self.assertEqual(metadata["config"]["percentiles"], [0.0, 50.0, 100.0])
        self.assertEqual(metadata["config"]["method_settings"],
                         dict(selection="oracle_grid", C=200000, m=100, alpha_divisor=24))
        self.assertIn("regression_threshold_scan.py", runner.load_json(directory / "source.json")["files"])
        record = next(runner.read_scan_records(directory))
        self.assertEqual(record["bootstrap_errors"].shape, (3, 6, 2))

    def test_terminal_baseline_runs_new_methods_in_separate_processes(self):
        old = self.old_run()
        output = self.root / "output"
        command = [sys.executable, str(runner.ROOT / "run_regression_threshold_scan.py"),
                   "--stage", "baseline", "--data-from", str(old), "--output-root", str(output),
                   "--N", "1", "--B", "4", "--smoke-test", "--backend", "cpu",
                   "--cpu-workers", "1", "--blas-threads", "1", "--batch-size", "2"]
        completed = subprocess.run(command, cwd=self.root, capture_output=True, text=True, timeout=120)
        logs = "\n".join(path.read_text(encoding="utf-8") for path in output.glob("*/*.log"))
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr + logs)
        manifest_path, = output.glob("*/manifest.json")
        manifest = runner.load_json(manifest_path)
        self.assertEqual(manifest["status"], "complete")
        self.assertEqual(set(manifest["runs"]), {"HBF", "FAR"})
        self.assertNotEqual(manifest["runs"]["HBF"]["pid"], manifest["runs"]["FAR"]["pid"])
        for report in manifest["runs"].values():
            metadata = runner.validate_run(report["run_directory"])
            self.assertEqual(metadata["config"]["percentiles"], [50.0])
            self.assertEqual(len(list(runner.read_scan_records(report["run_directory"]))), 1)


if __name__ == "__main__":
    unittest.main()
