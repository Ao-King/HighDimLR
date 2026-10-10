"""HBF/FAR spectral paths against the sparse notebooks' discrete recurrences."""

import ast
import copy
import json
from pathlib import Path
import sys
import unittest
from unittest import mock

import numpy as np
from scipy.special import gamma

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from regression_spectral import SpectralDesign


def notebook_namespace(method, count, options, delta=0.0, filename="Sparse.ipynb"):
    notebook = json.loads((ROOT / filename).read_text(encoding="utf-8"))
    function_name = "hbf" if method == "HBF" else "Frac"
    namespace = {"np": np, "gamma": gamma, "delta": delta}
    for cell in notebook["cells"]:
        if cell["cell_type"] != "code":
            continue
        tree = ast.parse("".join(cell["source"]))
        if not any(isinstance(node, ast.FunctionDef) and node.name == function_name for node in tree.body):
            continue
        if method == "FAR":
            # Execute the notebook's actual A/B coefficient construction, with
            # only the test dimensions/order/step substituted before the loop.
            prefix = []
            for node in tree.body:
                if isinstance(node, ast.FunctionDef):
                    break
                node = copy.deepcopy(node)
                if isinstance(node, ast.Assign) and len(node.targets) == 1:
                    target = node.targets[0]
                    overrides = {"max_iter": count, "t": options["step"], "theta": options.get("theta", 1.5)}
                    if isinstance(target, ast.Name) and target.id in overrides:
                        node.value = ast.copy_location(ast.Constant(overrides[target.id]), node.value)
                prefix.append(node)
            exec(compile(ast.fix_missing_locations(ast.Module(body=prefix, type_ignores=[])),
                         "<notebook FAR weights>", "exec"), namespace)
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name in (function_name, function_name + "_debias"):
                node = copy.deepcopy(node)
                if method == "HBF":
                    for statement in node.body:
                        if (isinstance(statement, ast.Assign) and len(statement.targets) == 1
                                and isinstance(statement.targets[0], ast.Name) and statement.targets[0].id == "eta"):
                            statement.value = ast.copy_location(ast.Constant(options.get("eta", 5)), statement.value)
                exec(compile(ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[])),
                             "<notebook discrete recursion>", "exec"), namespace)
        return namespace, function_name
    raise AssertionError("Missing notebook method")


def direct(method, X, Y, count, options, *, fixed=True, delta=0.0, filename="Sparse.ipynb"):
    namespace, name = notebook_namespace(method, count, options, delta, filename)
    args = (X, Y, options["step"], count) if method == "HBF" else (X, Y, count)
    return namespace[name + ("_debias" if fixed else "")](*args)


class ExtendedSpectralFilterTests(unittest.TestCase):
    def test_fixed_counts_match_notebook_ordinary_and_two_fit_debias(self):
        rng = np.random.default_rng(859)
        for X in (rng.normal(size=(6, 4))*.2, rng.normal(size=(4, 6))*.2,
                  np.array([[.5, 0., 0.], [0., .3, 0.], [0., 0., 0.], [0., 0., 0.]])):
            design = SpectralDesign(X)
            Y = rng.normal(size=X.shape[0])
            for method in ("HBF", "FAR"):
                for count in (0, 1, 2, 7):
                    for filename in ("Sparse.ipynb", "sparse15.ipynb"):
                        with self.subTest(shape=X.shape, method=method, count=count, notebook=filename):
                            options = {"step": .05, "eta": 5, "theta": 1.5}
                            expected = direct(method, X, Y, count, options, filename=filename)
                            refit = direct(method, X, X@expected, count, options, filename=filename)
                            ordinary, debiased = design.pair(Y, method, count, options)
                            np.testing.assert_allclose(ordinary, expected, rtol=3e-12, atol=3e-14)
                            np.testing.assert_allclose(debiased, 2*expected-refit, rtol=3e-12, atol=3e-14)

    def test_batched_responses_and_nondefault_parameters(self):
        rng = np.random.default_rng(362)
        X, Y = rng.normal(size=(4, 6))*.2, rng.normal(size=(4, 3))
        design = SpectralDesign(X)
        for method, options in (("HBF", {"step": .03, "eta": .7}),
                                ("FAR", {"step": .03, "theta": 1.2})):
            ordinary, debiased = design.pair(Y, method, 8, options)
            for j in range(Y.shape[1]):
                expected = direct(method, X, Y[:, j], 8, options)
                refit = direct(method, X, X@expected, 8, options)
                np.testing.assert_allclose(ordinary[:, j], expected, rtol=3e-12, atol=3e-14)
                np.testing.assert_allclose(debiased[:, j], 2*expected-refit, rtol=3e-12, atol=3e-14)

    def test_formal_steps_and_caps_match_large_spectral_mode_recursions(self):
        X = np.diag(np.sqrt([537050.2834, 1000., 1.]))
        Y = np.array([100., -2., .1])
        design = SpectralDesign(X)
        for method, options, counts in (
                ("HBF", {"step": 5e-4, "eta": 5}, (1, 20, 10000)),
                ("FAR", {"step": 5e-5, "theta": 1.5}, (1, 20, 500))):
            for count in counts:
                with self.subTest(method=method, count=count):
                    expected = direct(method, X, Y, count, options)
                    refit = direct(method, X, X@expected, count, options)
                    ordinary, debiased = design.pair(Y, method, count, options)
                    np.testing.assert_allclose(ordinary, expected, rtol=5e-10, atol=2e-12)
                    np.testing.assert_allclose(debiased, 2*expected-refit, rtol=5e-10, atol=2e-12)

    def test_trajectory_cache_reuses_and_extends_once_per_parameter_tuple(self):
        design = SpectralDesign(np.diag([.5, .3]))
        Y = np.ones(2)
        for method, name, extra, alternate in (("HBF", "_advance_hbf", "eta", 7),
                                               ("FAR", "_advance_far", "theta", 1.2)):
            options = {"step": .05, "max_steps": 8, "delta": 0.0}
            with mock.patch.object(design, name, wraps=getattr(design, name)) as advance:
                design.select(Y, method, options, np.zeros(2))
                self.assertEqual(advance.call_count, 8)
                design.select(2*Y, method, options, np.zeros(2))
                design.pair(Y, method, 4, options)
                design.pair(Y, method, 8, options)
                self.assertEqual(advance.call_count, 8)
                gain = design.filter(method, 10, options)
                self.assertEqual(advance.call_count, 10)
                gain[:] = 0
                self.assertGreater(np.linalg.norm(design.filter(method, 10, options)), 0)
                design.filter(method, 4, dict(options, **{extra: alternate}))
                self.assertEqual(advance.call_count, 14)
                design.filter(method, 4, dict(options, step=.01))
                self.assertEqual(advance.call_count, 18)

    def test_zero_design_aliases_and_invalid_parameters(self):
        design = SpectralDesign(np.zeros((2, 4)))
        for method, alias in (("HBF", "hbf"), ("FAR", "frac")):
            options = {"step": .05}
            for count in (0, 1, 7):
                pair = design.pair(np.ones(2), alias, count, options)
                for value in pair:
                    np.testing.assert_array_equal(value, np.zeros(4))
                np.testing.assert_array_equal(design.filter(method, count, options), np.zeros(2))
        for method, options in (("HBF", {"step": .1, "eta": -.1}),
                                ("HBF", {"step": .1, "eta": np.inf}),
                                ("FAR", {"step": .1, "theta": 0}),
                                ("FAR", {"step": .1, "theta": np.nan}),
                                ("FAR", {"step": -.1}),
                                ("HBF", {"step": .1, "eta": True})):
            with self.subTest(method=method, options=options), self.assertRaises(ValueError):
                design.filter(method, 1, options)
        for method in ("HBF", "FAR"):
            with self.assertRaises(ValueError):
                design.filter(method, .5, {"step": .1})


class ExtendedSpectralSelectionTests(unittest.TestCase):
    def test_first_hit_zero_cap_nonhit_cap_and_fixed_selected_debias(self):
        X = np.array([[.4, .1], [.1, .3], [.1, -.1]])
        Y, truth = np.array([1., -.5, .2]), np.array([.2, .4])
        design = SpectralDesign(X)
        for method in ("HBF", "FAR"):
            options = {"step": .05, "eta": 5, "theta": 1.5}
            at_two = direct(method, X, Y, 2, options)
            at_three = direct(method, X, Y, 3, options)
            threshold = (np.linalg.norm(Y-X@at_two)+np.linalg.norm(Y-X@at_three))/2
            for cap, delta in ((0, threshold), (1, 0.), (12, threshold), (7, 0.), (12, 100.)):
                with self.subTest(method=method, cap=cap, delta=delta):
                    expected, count = direct(method, X, Y, cap, options, fixed=False, delta=delta)
                    ordinary, debiased, selected, diagnostic = design.select(
                        Y, method, dict(options, max_steps=cap, delta=delta), truth)
                    self.assertEqual(selected, count)
                    if cap == 12:
                        self.assertEqual(count, 1 if delta == 100. else 3)
                    self.assertEqual(diagnostic["at_iteration_cap"], count == cap)
                    self.assertEqual(diagnostic["tdp_satisfied"], np.linalg.norm(Y-X@ordinary) <= delta)
                    np.testing.assert_allclose(ordinary, expected, rtol=3e-12, atol=3e-14)
                    refit = direct(method, X, X@expected, count, options)
                    np.testing.assert_allclose(debiased, 2*expected-refit, rtol=3e-12, atol=3e-14)

    def test_large_mode_formal_tdp_stops_at_notebook_count(self):
        X = np.diag(np.sqrt([537050.2834, 1000., 1.]))
        truth, noise = np.array([.4, -.2, 1.]), np.array([2., -2., 2.])
        Y = X@truth + noise
        design = SpectralDesign(X)
        for method, step, cap in (("HBF", 5e-4, 10000), ("FAR", 5e-5, 500)):
            options = {"step": step, "max_steps": cap, "delta": np.linalg.norm(noise)}
            expected, count = direct(method, X, Y, cap, options, fixed=False, delta=options["delta"])
            ordinary, _, parameter, diagnostic = design.select(Y, method, options, truth)
            self.assertEqual(parameter, count)
            np.testing.assert_allclose(ordinary, expected, rtol=5e-10, atol=2e-12)
            self.assertTrue(diagnostic["tdp_satisfied"] or diagnostic["at_iteration_cap"])

    def test_first_hit_only_extends_cache_to_selected_count(self):
        design = SpectralDesign(np.eye(2))
        for method, name in (("HBF", "_advance_hbf"), ("FAR", "_advance_far")):
            with mock.patch.object(design, name, wraps=getattr(design, name)) as advance:
                _, _, count, _ = design.select(np.ones(2), method,
                    {"step": .1, "max_steps": 1000, "delta": 100.}, np.zeros(2))
                self.assertEqual(count, 1)
                self.assertEqual(advance.call_count, 1)


if __name__ == "__main__":
    unittest.main()
