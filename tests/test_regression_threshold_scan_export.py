import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

import export_regression_threshold_scan as export
import regression_results as rr


def row(method, percentile, counts=(950, 950), radii=(1, 1), stage="calibration"):
    result = {"stage": stage, "seed_root": export.SEEDS[stage], "method": method,
              "percentile": percentile, "N": 1000, "B": 500,
              "coverage_count_original": counts[0], "coverage_count_debiased": counts[1]}
    result.update(dict.fromkeys(export.STATISTICS, 1.0))
    result.update(dict.fromkeys(export.MCSE_NAMES, 0.01))
    result.update(coverage_original=counts[0] / 1000, coverage_debiased=counts[1] / 1000,
                  mean_radius_original=radii[0], mean_radius_debiased=radii[1])
    return result


def coarse_rows():
    return [row(method, q) for method in export.METHODS for q in export.COARSE]


def common():
    return {"n": 1000, "p": 1000, "N": 1000, "B": 500, "sigma": 2,
            "confidence": 0.95, "seed_root": 20261007, "design_seed": 24,
            "noise_model": "normal", "threshold_grid": [0.001, 0.001, 1.0],
            "preprocess_blas_threads": 1, "X": "design", "beta": "truth"}


class ThresholdSelectionTests(unittest.TestCase):
    def test_sample_aggregation_uses_mean_errors_radii_and_mcse(self):
        config = {"N": 4, "B": 500, "stage": "calibration", "seed_root": 20261007, "method": "Ridge"}
        samples = np.ones((4, 1, 10))
        samples[:, 0, 0] = [1, 1, 0, 0]
        samples[:, 0, 1] = [1, 1, 1, 0]
        samples[:, 0, 2] = [1, 2, 3, 4]
        samples[:, 0, 8] = [2, 4, 6, 8]
        result = export.summarize_samples(config, [30], samples)[0]
        self.assertEqual(result["coverage_original"], 0.5)
        self.assertEqual(result["coverage_count_debiased"], 3)
        self.assertEqual(result["l2_original"], 2.5)
        self.assertEqual(result["mean_radius_original"], 5)
        self.assertAlmostEqual(result["mcse_coverage_original"], 0.25)
        self.assertAlmostEqual(result["mcse_l2_original"], np.std([1, 2, 3, 4], ddof=1) / 2)
        samples[0, 0, 0] = 0.5
        with self.assertRaisesRegex(ValueError, "binary"):
            export.summarize_samples(config, [30], samples)

    def test_selection_uses_counts_to_avoid_floating_pseudo_ties(self):
        rows = coarse_rows()
        for item in rows:
            if item["method"] == "HBF":
                item.update(row("HBF", item["percentile"], (700, 700)))
        next(item for item in rows if item["method"] == "HBF" and item["percentile"] == 10).update(row("HBF", 10, (940, 940)))
        next(item for item in rows if item["method"] == "HBF" and item["percentile"] == 20).update(row("HBF", 20, (960, 960)))
        choices, refinement, _ = export.choose_percentiles(rows)
        self.assertEqual(choices["HBF"]["bestII"]["percentile"], 10)
        self.assertEqual(choices["HBF"]["bestBoth"]["percentile"], 10)
        self.assertEqual(refinement, {})

    def test_radii_tiebreakers_are_branch_specific_then_choose_smaller_q(self):
        rows = coarse_rows()
        for item in rows:
            if item["method"] == "Ridge":
                item.update(row("Ridge", item["percentile"], (500, 500)))
        for q, radii in ((10, (10, 1)), (20, (0.5, 2)), (30, (0.5, 2))):
            next(item for item in rows if item["method"] == "Ridge" and item["percentile"] == q).update(row("Ridge", q, radii=radii))
        choices, _, _ = export.choose_percentiles(rows)
        self.assertEqual(choices["Ridge"]["bestII"]["percentile"], 10)
        self.assertEqual(choices["Ridge"]["bestBoth"]["percentile"], 20)

    def test_refinement_requires_both_objectives_and_reports_failure_honestly(self):
        rows = coarse_rows()
        for item in rows:
            if item["method"] == "FAR":
                item.update(row("FAR", item["percentile"], (500, 950)))
        choices, refinement, _ = export.choose_percentiles(rows)
        self.assertTrue(choices["FAR"]["bestII"]["meets_target"])
        self.assertFalse(choices["FAR"]["bestBoth"]["meets_target"])
        self.assertEqual(refinement, {"FAR": list(export.MIDPOINTS)})
        rows.extend(row("FAR", q, (600, 950), stage="refinement") for q in export.MIDPOINTS)
        plan = export._selection_plan(rows, [], common())
        self.assertEqual(plan["status"], "ready_for_validation")
        self.assertFalse(plan["selections"]["FAR"]["bestBoth"]["meets_target"])
        self.assertEqual(plan["validation_percentiles"]["FAR"], [5])

    def test_nonmonotonicity_is_reported_without_interpolation(self):
        rows = coarse_rows()
        for q, counts in ((0, (900, 900)), (10, (950, 970)), (20, (920, 940))):
            next(item for item in rows if item["method"] == "SC" and item["percentile"] == q).update(row("SC", q, counts))
        _, _, monotonicity = export.choose_percentiles(rows)
        self.assertEqual(monotonicity["SC"], {"original": "nonmonotone", "debiased": "nonmonotone"})
        rows[0]["seed_root"] = 20261008
        with self.assertRaisesRegex(ValueError, "Only the calibration seed"):
            export.choose_percentiles(rows)

    def test_threshold_includes_exactly_one_point_five_percentage_points(self):
        rows = coarse_rows()
        for item in rows:
            if item["method"] == "FAR":
                item.update(row("FAR", item["percentile"], (935, 965)))
        choices, refinement, _ = export.choose_percentiles(rows)
        self.assertTrue(choices["FAR"]["bestBoth"]["meets_target"])
        self.assertEqual(refinement, {})


class ThresholdExportTests(unittest.TestCase):
    def test_actual_scan_storage_roundtrip_for_all_nine_methods(self):
        from regression_threshold_scan import scan_method_settings
        from run_regression_threshold_scan import ScanRun

        with tempfile.TemporaryDirectory() as temporary:
            paths = []
            X, beta = np.eye(2), np.ones(2)
            for method in export.METHODS:
                config = {**common(), "stage": "calibration", "method": method, "n": 2, "p": 2,
                          "N": 2, "B": 3, "percentiles": list(export.COARSE),
                          "largest_singular_value": 1.0,
                          "method_settings": scan_method_settings(method, 2, 1.0)}
                config.pop("X")
                config.pop("beta")
                with ScanRun(config, X, beta, {"files": {"fixture.py": "# same synthetic source"}}, temporary) as run:
                    for index in range(2):
                        scan = np.zeros((1000, 3))
                        scan[:, 0] = 0.001 * np.arange(1, 1001)
                        run.record(index, selected_parameter=1, diagnostics={},
                                   estimates=np.zeros((11, 2, 2)), sigma=np.ones((11, 2)),
                                   thresholds=np.ones((11, 2)), bootstrap_errors=np.ones((11, 3, 2)),
                                   threshold_scan=scan, percentiles=config["percentiles"])
                paths.append(run.directory)
            with self.assertRaisesRegex(ValueError, "Nonformal"):
                export._load_runs(paths, "scan")
            # Only the test bypasses production dimensions, never the CLI or export API.
            with patch.object(export, "_require_formal", return_value=None):
                rows, loaded, shared = export._load_runs(paths, "scan")
            self.assertEqual(len(rows), 99)
            self.assertEqual(len(loaded), 9)
            self.assertEqual(shared["N"], 2)
            self.assertTrue(all(item["coverage_original"] == 1 for item in rows))
            self.assertTrue(all(np.isclose(item["l2_original"], np.sqrt(2)) for item in rows))

    def test_full_matrix_and_ready_plan_export_without_plots_and_without_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary:
            with patch.object(export, "_load_runs", return_value=(coarse_rows(), [], common())):
                first = export.export_threshold_scan([], temporary, plots=False)
                second = export.export_threshold_scan([], temporary, plots=False)
            self.assertNotEqual(first["directory"], second["directory"])
            with first["aggregate_csv"].open(newline="", encoding="utf-8") as stream:
                records = list(csv.DictReader(stream))
            self.assertEqual(len(records), 99)
            self.assertEqual(set(record["method"] for record in records), set(export.METHODS))
            plan = json.loads(first["selection_plan"].read_text())
            self.assertEqual(plan["validation_seed_root"], 20261008)
            self.assertEqual(plan["status"], "ready_for_validation")
            self.assertFalse(plan["interpolation_used"])
            self.assertIn(r"\begin{longtable}", first["aggregate_tex"].read_text())
            self.assertFalse(list(first["directory"].glob("*.pdf")))

    def test_pending_refinement_has_report_but_no_validation_plan(self):
        rows = coarse_rows()
        for item in rows:
            if item["method"] == "FAR":
                item.update(row("FAR", item["percentile"], (0, 0)))
        with tempfile.TemporaryDirectory() as temporary:
            with patch.object(export, "_load_runs", return_value=(rows, [], common())):
                result = export.export_threshold_scan([], temporary, plots=False)
            self.assertIsNone(result["selection_plan"])
            self.assertEqual(result["refinement_percentiles"], {"FAR": list(export.MIDPOINTS)})
            self.assertTrue((result["directory"] / "selection_report.json").exists())

    def test_formal_guards_reject_seed_and_dimension_changes(self):
        config = {**common(), "method": "Ridge"}
        export._require_formal(config, "calibration")
        for key, value in (("seed_root", 20261006), ("N", 50), ("B", 20), ("p", 1500)):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "Nonformal"):
                export._require_formal({**config, key: value}, "calibration")

    def test_validation_requires_independent_seed_selected_points_and_same_plan(self):
        inputs = [{"fingerprints": {"source": "frozen-source"}}]
        plan = export._selection_plan(coarse_rows(), inputs, common())
        validation_rows = [row(method, 0, (930, 940), stage="validation") for method in export.METHODS]
        independent = {**common(), "seed_root": 20261008}
        with tempfile.TemporaryDirectory() as temporary:
            loaded = []
            for method in export.METHODS:
                directory = Path(temporary) / method
                directory.mkdir()
                (directory / "selection_plan.json").write_text(json.dumps(plan), encoding="utf-8")
                loaded.append({"directory": str(directory), "method": method,
                               "config": {"selection_plan_sha256": rr._digest(plan)},
                               "fingerprints": {"source": "frozen-source"}})
            result = export._check_validation_plan(plan, validation_rows, loaded, independent)
            self.assertEqual(len(result), 18)
            self.assertFalse(next(item for item in result if item["selection_rule"] == "bestBoth")["validation_meets_target"])
            self.assertTrue(next(item for item in result if item["selection_rule"] == "bestII")["validation_meets_target"])
            with self.assertRaisesRegex(ValueError, "independently seeded"):
                export._check_validation_plan({**plan, "validation_seed_root": 20261007}, validation_rows, loaded, independent)
            with self.assertRaisesRegex(ValueError, "independently seeded"):
                export._check_validation_plan(plan, validation_rows, loaded, common())
            with self.assertRaisesRegex(ValueError, "prespecified"):
                export._check_validation_plan(plan, validation_rows + [row("FAR", 5, stage="validation")], loaded, independent)
            loaded[0]["config"]["selection_plan_sha256"] = "changed"
            with self.assertRaisesRegex(ValueError, "provenance"):
                export._check_validation_plan(plan, validation_rows, loaded, independent)


if __name__ == "__main__":
    unittest.main()
