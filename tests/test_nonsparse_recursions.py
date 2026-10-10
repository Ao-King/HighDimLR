import ast
import copy
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np
from threadpoolctl import threadpool_limits

import verify_nonsparse_recursions as verifier


ROOT = Path(__file__).resolve().parents[1]


class NonsparseRecursionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.threads = threadpool_limits(limits=1)
        cls.threads.__enter__()
        rng = np.random.default_rng(617)
        cls.X = rng.normal(size=(4, 4)) + 3 * np.eye(4)
        cls.truth = rng.normal(size=4)
        cls.e = rng.normal(size=4)
        cls.Y = cls.X @ cls.truth + cls.e
        cls.U, cls.s, cls.Vt = np.linalg.svd(cls.X)
        cls.delta = .6 * np.linalg.norm(cls.e)

    @classmethod
    def tearDownClass(cls):
        cls.threads.__exit__(None, None, None)

    def direct_pair(self, notebook, method, cap=7, delta=None):
        delta = self.delta if delta is None else delta
        selected, fixed, options = verifier.load_kernels(notebook, method, self.X,
                                                        self.truth, self.s, delta, cap=cap)
        ordinary, count = verifier.invoke(selected, self.X, self.Y, options, cap, self.truth, delta)
        if method == "Landweber":
            debiased = verifier.invoke(fixed, self.X, self.Y, options, 2 * count, self.truth, delta)
        else:
            refit = verifier.invoke(fixed, self.X, self.X @ ordinary, options, count, self.truth, delta)
            debiased = 2 * ordinary - refit
        return ordinary, debiased, count

    def test_original_dense_and_coordinate_selection_and_debias_all_methods(self):
        for rule, filename in verifier.NOTEBOOKS.items():
            for method in verifier.METHODS:
                with self.subTest(rule=rule, method=method):
                    ordinary, debiased, count = self.direct_pair(ROOT / filename, method)
                    report = verifier.verify_method(ROOT / filename, rule, method,
                        self.X, self.Y, self.truth, self.U, self.s, self.Vt,
                        ordinary, debiased, count, self.delta, cap=7)
                    self.assertTrue(report["passed"], report)

    def test_optional_full_dense_three_methods_and_default_skips_selectors(self):
        for rule, filename in verifier.NOTEBOOKS.items():
            for method in ("HBF", "SOAR", "Nesterov"):
                ordinary, debiased, count = self.direct_pair(ROOT / filename, method)
                for enabled in (False, True):
                    with self.subTest(rule=rule, method=method, enabled=enabled):
                        with mock.patch.object(verifier, "invoke", wraps=verifier.invoke) as calls:
                            report = verifier.verify_method(ROOT / filename, rule, method,
                                self.X, self.Y, self.truth, self.U, self.s, self.Vt,
                                ordinary, debiased, count, self.delta, cap=7,
                                full_dense_methods=("HBF", "SOAR", "Nesterov") if enabled else ())
                        self.assertTrue(report["passed"], report)
                        full = report["full_dense_reference"]
                        self.assertEqual(full["enabled"], enabled)
                        dense_selections = [call for call in calls.call_args_list
                            if isinstance(call.args[1], np.ndarray)
                            and call.args[0].__name__ == verifier.METHODS[method]]
                        self.assertEqual(len(dense_selections), int(enabled))
                        if enabled:
                            self.assertEqual(dense_selections[0].args[4], 7)
                            self.assertEqual(full["selected_count"], count)
                            self.assertEqual(full["iteration_cap"], 7)
                            self.assertTrue(full["passed"])
                            self.assertIn("full_dense_fixed_parameter_debias", report["timings_seconds"])

    def test_tdp_first_update_and_count_mismatch_are_not_hidden(self):
        for method in verifier.METHODS:
            with self.subTest(method=method):
                ordinary, debiased, count = self.direct_pair(ROOT / "NonsparseDP.ipynb", method, delta=1e9)
                self.assertEqual(count, 1)
                report = verifier.verify_method(ROOT / "NonsparseDP.ipynb", "TDP", method,
                    self.X, self.Y, self.truth, self.U, self.s, self.Vt,
                    ordinary, debiased, count + 1, 1e9, cap=7, full_dense_methods=(method,))
                self.assertFalse(report["passed"])
                self.assertFalse(report["checks"]["selected_count"]["close"])
                self.assertTrue(report["checks"]["ordinary"]["close"])
                self.assertFalse(report["full_dense_reference"]["checks"]["selected_count_archived"]["close"])

    def test_ast_transform_changes_only_hbf_square_zero_matrix(self):
        original = ast.parse("A = np.zeros((m,m))\nb = np.zeros((m,k_max+1))\nc = np.zeros(m)\n")
        transformed = verifier.MatrixZeros().visit(copy.deepcopy(original))
        self.assertEqual(ast.unparse(transformed),
                         "A = _zero_matrix(m)\nb = np.zeros((m, k_max + 1))\nc = np.zeros(m)")

    def test_original_aosr_floor_and_explicit_smoke_floor(self):
        for method in verifier.METHODS:
            with self.subTest(method=method):
                for minimum, expected in ((None, 12), (5, 5)):
                    selected, _, options = verifier.load_kernels(ROOT / "Nonsparse.ipynb", method,
                        self.X, self.truth, self.s, self.delta, cap=12, minimum_steps=minimum)
                    _, count = verifier.invoke(selected, self.X, np.zeros(4), options, 12, self.truth, self.delta)
                    self.assertEqual(count, expected)

    def test_ark_inverse_and_far_configuration_are_original(self):
        for filename in verifier.NOTEBOOKS.values():
            _, fixed, options = verifier.load_kernels(ROOT / filename, "ARk", self.X,
                                                      self.truth, self.s, self.delta, cap=5)
            self.assertEqual(options["kappa"], .5)
            first = verifier.invoke(fixed, self.X, self.Y, options, 1, self.truth, self.delta)
            np.testing.assert_array_equal(first, np.zeros(4))
            _, _, far = verifier.load_kernels(ROOT / filename, "FAR", self.X,
                                               self.truth, self.s, self.delta, cap=5)
            self.assertEqual(far["theta"], 1.8)
            self.assertEqual(far["t"], 5e-4)

    def make_fixture(self, directory):
        sources = directory / "sources"
        sources.mkdir()
        digests, rows, estimates = {}, [], {}
        for rule, filename in verifier.NOTEBOOKS.items():
            shutil.copyfile(ROOT / filename, sources / filename)
            digests["sources/" + filename] = hashlib.sha256((sources / filename).read_bytes()).hexdigest()
            for method in verifier.METHODS:
                ordinary, debiased, count = self.direct_pair(sources / filename, method)
                rows.append(dict(rule=rule, method=method, k0=count, parameter=count))
                estimates[f"{rule}__{method}__ordinary"] = ordinary
                estimates[f"{rule}__{method}__debiased"] = debiased
        np.savez(directory / "data.npz", X=self.X, Y=self.Y, beta=self.truth, e=self.e,
                 U=self.U, S=self.s, Vt=self.Vt)
        np.savez(directory / "estimates.npz", **estimates)
        summary = dict(status="complete", configuration=dict(n=4, p=4, smoke_test=True),
                       sources=digests, rows=rows)
        (directory / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
        return summary

    def test_archived_roundtrip_and_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.make_fixture(directory)
            args = ["--run-directory", str(directory), "--allow-nonformal", "--max-steps", "7",
                    "--full-dense-methods", "HBF", "SOAR", "Nesterov"]
            self.assertEqual(verifier.main(args), 0)
            report = json.loads((directory / "recursion-verification.json").read_text())
            self.assertEqual(len(report["methods"]), 14)
            self.assertTrue(report["passed"])
            self.assertEqual(report["full_dense_methods"], ["HBF", "SOAR", "Nesterov"])
            self.assertIn("not replay", report["limitations"])
            with self.assertRaises(SystemExit):
                verifier.main(args)

    def test_hash_and_formal_dimension_guards(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.make_fixture(directory)
            with self.assertRaisesRegex(ValueError, "Nonformal dimensions"):
                verifier.verify_run(directory)
            summary = json.loads((directory / "summary.json").read_text())
            summary["configuration"]["smoke_test"] = False
            (directory / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "AOSR floor"):
                verifier.verify_run(directory, allow_nonformal=True, minimum_steps=5)
            notebook = directory / "sources" / "Nonsparse.ipynb"
            notebook.write_bytes(notebook.read_bytes() + b"\n")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                verifier.verify_run(directory, allow_nonformal=True, cap=7)


if __name__ == "__main__":
    unittest.main()
