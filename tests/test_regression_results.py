"""Small persistence tests; no notebook Monte Carlo cells are executed."""

import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from regression_results import BootstrapRun, METHODS, METRICS, VERSION, export_regression_table, is_production_config


class RegressionResultsTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.notebook = self.root / "test.ipynb"
        self.code = "resume_from = None\noutput_root = 'outputs'\nN = 2\n"
        self.write_notebook()
        self.X = np.eye(3)
        self.beta = np.array([1., -1., 0.])
        self.config = {"method": "SC", "n": 3, "p": 3, "N": 2, "B": 3,
                       "sigma": 2., "confidence": .95, "method_settings": {"selection": "oracle_rank", "rank_max": 3},
                       "threshold_percentile": 50, "seed_root": 20261006, "design_seed": 24,
                       "noise_model": "normal", "threshold_grid": [.001, .001, 1.]}

    def write_notebook(self, output=None):
        notebook = {"cells": [{"cell_type": "code", "source": self.code.splitlines(keepends=True),
                               "outputs": [] if output is None else [output]}], "metadata": {}}
        self.notebook.write_text(json.dumps(notebook), encoding="utf-8")

    def create(self, config=None, resume_from=None, X=None, beta=None):
        return BootstrapRun(config or self.config, self.X if X is None else X,
                            self.beta if beta is None else beta, self.notebook,
                            output_root=self.root / "runs", resume_from=resume_from)

    def values(self):
        return {"selected_parameter": 2., "thresholds": [.01, .02],
                "estimates": np.column_stack((self.beta + np.array([.3, .4, 0.]), self.beta)),
                "sigma": [1., 2.], "bootstrap_errors": np.array([[.1, 0.], [.2, 0.], [.5, 0.]]),
                "threshold_scan": [[.001, 1., 2.], [.002, .5, 1.]],
                "threshold_columns": ["bn", "ordinary", "debiased"],
                "diagnostics": {"cap_hit": False, "max_steps": 10}}

    def test_complete_run_stores_metrics_and_all_raw_inputs(self):
        with self.create() as run:
            run.record(0, **self.values())
            run.record(1, **self.values())
        self.assertEqual(run.start_index, 2)
        self.assertEqual(run.summary["coverage_debiased"], 1.)
        self.assertEqual(run.summary["coverage_original"], 1.)
        self.assertAlmostEqual(run.summary["l2_original"], .5)
        self.assertAlmostEqual(run.summary["linf_original"], .4)
        self.assertAlmostEqual(run.summary["variance_absolute_error_original"], 3.)
        metadata = json.loads((run.directory / "metadata.json").read_text())
        self.assertEqual(metadata["status"], "complete")
        self.assertFalse(metadata["is_manuscript_configuration"])
        with np.load(run.directory / "replications/rep_000000.npz", allow_pickle=False) as saved:
            np.testing.assert_array_equal(saved["variance"], [1., 4.])
            np.testing.assert_allclose(saved["quantiles"], [.47, 0.])
            self.assertEqual(saved["bootstrap_errors"].shape, (3, 2))
            self.assertEqual(saved["threshold_scan"].shape, (2, 3))
            self.assertEqual(saved["metric_names"].tolist(), list(METRICS))
        self.assertTrue((run.directory / "summary.csv").exists())
        with (run.directory / "summary.csv").open(newline="") as stream:
            summary_rows = list(csv.DictReader(stream))
        self.assertEqual(summary_rows[0]["status"], "complete")
        self.assertEqual(summary_rows[0]["target_replications"], "2")
        self.assertTrue(set(("scipy", "scikit-learn", "numba", "joblib")) <= set(metadata["runtime"]))
        with (run.directory / "replications.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 2)

    def test_keyboard_interrupt_preserves_prefix_and_resume_finishes(self):
        run = self.create()
        with self.assertRaises(KeyboardInterrupt):
            with run:
                run.record(0, **self.values())
                raise KeyboardInterrupt("test interruption")
        metadata = json.loads((run.directory / "metadata.json").read_text())
        self.assertEqual(metadata["status"], "failed")
        self.assertEqual(metadata["completed_replications"], 1)
        self.assertFalse((run.directory / ".active.lock").exists())
        with self.create(resume_from=run.directory) as resumed:
            self.assertEqual(resumed.start_index, 1)
            resumed.record(1, **self.values())
        with self.create() as uninterrupted:
            for x in range(2):
                uninterrupted.record(x, **self.values())
        self.assertEqual(resumed.summary, uninterrupted.summary)

    def test_normal_early_exit_is_incomplete_and_new_run_never_overwrites(self):
        with self.create() as first:
            first.record(0, **self.values())
        with self.create() as second:
            pass
        self.assertNotEqual(first.directory, second.directory)
        self.assertEqual(json.loads((first.directory / "metadata.json").read_text())["status"], "incomplete")

    def test_resume_reads_committed_files_not_stale_checkpoint(self):
        with self.create() as run:
            run.record(0, **self.values())
        checkpoint = run.directory / "checkpoint.json"
        checkpoint.write_text('{"completed_replications":0}', encoding="utf-8")
        metadata_path = run.directory / "metadata.json"
        metadata = json.loads(metadata_path.read_text())
        metadata["completed_replications"] = 0
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
        (run.directory / "replications/.unfinished.tmp").write_bytes(b"incomplete")
        with self.create(resume_from=run.directory) as resumed:
            self.assertEqual(resumed.start_index, 1)
            resumed.record(1, **self.values())

    def test_resume_accepts_output_and_io_path_changes_but_rejects_numeric_code(self):
        with self.create() as run:
            run.record(0, **self.values())
        self.code = "resume_from = 'different/path'\noutput_root = 'new/output'\nN = 2\n"
        self.write_notebook({"output_type": "stream", "text": "new output"})
        with self.create(resume_from=run.directory) as resumed:
            resumed.record(1, **self.values())
        self.code += "N = 3\n"
        self.write_notebook()
        with self.assertRaisesRegex(ValueError, "identical"):
            self.create(resume_from=run.directory)

    def test_resume_rejects_changed_config_or_data(self):
        with self.create() as run:
            run.record(0, **self.values())
        for configuration in ({**self.config, "B": 4}, {**self.config, "threshold_percentile": 60}):
            with self.assertRaisesRegex(ValueError, "identical"):
                self.create(config=configuration, resume_from=run.directory)
        with self.assertRaisesRegex(ValueError, "identical"):
            self.create(resume_from=run.directory, X=2*self.X)
        with self.assertRaisesRegex(ValueError, "identical"):
            self.create(resume_from=run.directory, beta=self.beta+.1)

    def test_acceleration_sources_are_archived_and_invalidate_resume_when_changed(self):
        from regression_results import _source

        repository = Path(__file__).resolve().parents[1]
        snapshot = _source(repository / 'BootstrapARk.ipynb')
        self.assertEqual(VERSION, 'regression-bootstrap-v2')
        self.assertTrue({'regression_common.py', 'regression_spectral.py',
                         'regression_bootstrap.py', 'regression_execution.py',
                         'regression_results.py'} <= snapshot['helpers'].keys())
        helper = self.root / 'regression_spectral.py'
        helper.write_text('version = 1\n', encoding='utf-8')
        with self.create() as run:
            run.record(0, **self.values())
        helper.write_text('version = 2\n', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'identical'):
            self.create(resume_from=run.directory)

    def test_resume_rejects_changed_execution_settings(self):
        config = {**self.config, 'execution': {'backend': 'cpu', 'batch_size': 64}}
        with self.create(config=config) as run:
            run.record(0, **self.values())
        config['execution']['batch_size'] = 32
        with self.assertRaisesRegex(ValueError, 'identical'):
            self.create(config=config, resume_from=run.directory)

    def test_record_validates_index_shapes_nonfinite_and_no_overwrite(self):
        with self.create() as run:
            with self.assertRaises(ValueError):
                run.record(1, **self.values())
            for name, value in (("bootstrap_errors", np.zeros((2, 2))), ("sigma", [np.nan, 1.]),
                                ("estimates", np.zeros((2, 3))), ("sigma", [-1., 1.])):
                arguments = self.values()
                arguments[name] = value
                with self.assertRaises(ValueError):
                    run.record(0, **arguments)
            run.record(0, **self.values())
            with self.assertRaises(ValueError):
                run.record(0, **self.values())

    def test_failed_atomic_replication_write_does_not_advance_checkpoint(self):
        run = self.create()
        original = np.savez_compressed
        def failing(stream, **arrays):
            if "index" in arrays:
                stream.write(b"partial")
                raise OSError("simulated full disk")
            return original(stream, **arrays)
        with self.assertRaises(OSError):
            with run:
                with mock.patch("regression_results.np.savez_compressed", side_effect=failing):
                    run.record(0, **self.values())
        self.assertEqual(list((run.directory / "replications").glob("*.npz")), [])
        self.assertEqual(run.start_index, 0)

    def test_derived_overflow_is_rejected_before_committing_replication(self):
        with self.create() as run:
            for name, value in (("estimates", np.full((3, 2), 1e200)), ("sigma", [1e200, 1e200])):
                arguments = self.values()
                arguments[name] = value
                with self.assertRaisesRegex(ValueError, "statistics must be finite"):
                    run.record(0, **arguments)
                self.assertEqual(run.start_index, 0)
                self.assertEqual(list((run.directory / "replications").glob("*.npz")), [])
                checkpoint = json.loads((run.directory / "checkpoint.json").read_text())
                self.assertEqual(checkpoint["completed_replications"], 0)
            run.record(0, **self.values())

    def test_finite_large_variance_metrics_do_not_overflow_during_averaging(self):
        arguments = self.values()
        arguments["sigma"] = [1e154, 1e154]
        with self.create() as run:
            for x in range(2):
                run.record(x, **arguments)
        self.assertTrue(np.isfinite(run.summary["variance_absolute_error_original"]))
        self.assertAlmostEqual(run.summary["variance_absolute_error_original"] / 1e308, 1.)

    def test_active_run_cannot_be_opened_twice(self):
        with self.create() as first:
            second = self.create(resume_from=first.directory)
            with self.assertRaisesRegex(RuntimeError, "locked"):
                with second:
                    pass

    def complete_methods(self):
        directories = []
        for method in METHODS:
            with self.create(config={**self.config, "method": method}) as run:
                for x in range(2):
                    run.record(x, **self.values())
            directories.append(run.directory)
        return directories

    def test_export_requires_complete_unique_methods_and_explicit_nonproduction(self):
        directories = self.complete_methods()
        with self.assertRaisesRegex(ValueError, "Nonproduction"):
            export_regression_table(directories, self.root / "export")
        for invalid in (directories[:-1], directories + [directories[0]]):
            with self.assertRaises(ValueError):
                export_regression_table(invalid, self.root / "export", allow_nonproduction=True)
        paths = export_regression_table(directories, self.root / "export", allow_nonproduction=True)
        self.assertEqual(paths["csv"].name, "regression_table.csv")
        self.assertEqual(paths["tex"].name, "regression_table.tex")
        self.assertEqual(paths["manifest"].name, "regression_table_manifest.json")
        self.assertNotIn("NON_MANUSCRIPT", paths["tex"].read_text())
        with paths["csv"].open(newline="") as stream:
            rows = list(csv.reader(stream))
        self.assertEqual(rows[0], ["metric", *METHODS])
        self.assertEqual(rows[1][0], "threshold_percentile")
        self.assertIn(r"\begin{tabular}", paths["tex"].read_text())
        self.assertIn(r"Coverage I & 100.0\%", paths["tex"].read_text())
        self.assertIn(r"$b_n^\circ$", paths["tex"].read_text())
        self.assertEqual(rows[2][1], "1.0")
        manifest = json.loads(paths["manifest"].read_text())
        self.assertTrue(manifest["non_manuscript_export"])
        self.assertEqual(manifest["common_settings"]["N"], self.config["N"])
        self.assertEqual(manifest["common_settings"]["B"], self.config["B"])
        self.assertEqual(set(manifest["runs"]), set(METHODS))
        self.assertEqual(manifest["runs"]["SC"]["directory"], str(directories[0].resolve()))
        self.assertEqual(manifest["runs"]["SC"]["config"], self.config)
        self.assertIn("identity", manifest["runs"]["SC"])
        self.assertIn("runtime", manifest["runs"]["SC"])
        with self.assertRaises(FileExistsError):
            export_regression_table(directories, self.root / "export", allow_nonproduction=True)

    def test_export_recomputes_summary_and_rejects_old_version_or_incomplete(self):
        directories = self.complete_methods()
        metadata_path = directories[0] / "metadata.json"
        metadata = json.loads(metadata_path.read_text())
        metadata["summary"] = {key: 12345 for key in METRICS}
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
        paths = export_regression_table(directories, self.root / "export", allow_nonproduction=True)
        self.assertNotIn("12345", paths["csv"].read_text())
        for key, value in (("version", "legacy"), ("status", "incomplete"), ("completed_replications", 1)):
            invalid = dict(metadata)
            invalid[key] = value
            metadata_path.write_text(json.dumps(invalid), encoding="utf-8")
            with self.assertRaises(ValueError):
                export_regression_table(directories, self.root / ("invalid-"+key), allow_nonproduction=True)

    def test_export_rejects_different_common_design_and_changed_artifacts(self):
        directories = self.complete_methods()
        with self.create(config={**self.config, "method": "SC"}, X=2*self.X) as altered:
            for x in range(2):
                altered.record(x, **self.values())
        with self.assertRaisesRegex(ValueError, "share"):
            export_regression_table([altered.directory, *directories[1:]], self.root / "export", allow_nonproduction=True)
        np.savez_compressed(directories[0] / "data.npz", X=self.X+1, beta=self.beta)
        with self.assertRaisesRegex(ValueError, "fingerprints"):
            export_regression_table(directories, self.root / "export", allow_nonproduction=True)

    def test_production_flag_cannot_override_actual_configuration(self):
        config = {**self.config, "is_manuscript_configuration": True}
        self.assertFalse(is_production_config(config, self.X))
        with self.create(config=config) as run:
            self.assertFalse(run.metadata["is_manuscript_configuration"])
        self.assertEqual(run.metadata["version"], VERSION)

    def test_export_rejects_different_runtime(self):
        directories = self.complete_methods()
        path = directories[0] / "metadata.json"
        metadata = json.loads(path.read_text())
        metadata["runtime"]["numpy"] = "different"
        path.write_text(json.dumps(metadata), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "runtime"):
            export_regression_table(directories, self.root / "export", allow_nonproduction=True)

    def test_export_rejects_changed_shared_helper_code(self):
        directories = []
        helper = self.root / "regression_common.py"
        for index, method in enumerate(METHODS):
            helper.write_text("VALUE = " + str(1 if index == 0 else 2) + "\n", encoding="utf-8")
            with self.create(config={**self.config, "method": method}) as run:
                for x in range(2):
                    run.record(x, **self.values())
            directories.append(run.directory)
        with self.assertRaisesRegex(ValueError, "helper code"):
            export_regression_table(directories, self.root / "export", allow_nonproduction=True)

    def test_missing_prefix_record_is_not_silently_resumed(self):
        with self.create() as run:
            run.record(0, **self.values())
            run.record(1, **self.values())
        (run.directory / "replications/rep_000000.npz").unlink()
        with self.assertRaisesRegex(ValueError, "contiguous"):
            self.create(resume_from=run.directory)


if __name__ == "__main__":
    unittest.main()
