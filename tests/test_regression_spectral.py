"""Compare spectral filters with the unchanged notebook discrete kernels."""

import ast
import copy
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from regression_spectral import SpectralDesign


METHODS = ("Landweber", "Showalter", "ARk", "SOAR", "Nesterov")
NOTEBOOKS = {"Landweber": ("Bootstraplandweber", "landweber"),
             "Showalter": ("Bootstrapshowalter", "showalter"),
             "ARk": ("BootstrapARk", "ARk"),
             "SOAR": ("Bootstrapsecond", "second"),
             "Nesterov": ("BootstrapNesterov", "Nesterov")}


def settings(method, **changes):
    values = {"step": .05, "max_steps": 12, "minimum_steps": 500,
              "kappa": .5, "s_star": .5, "omega": 3}
    values.update(changes)
    return values


def legacy(method, X, Y, count, options, fixed=True, truth=None, delta=0):
    name, function_name = NOTEBOOKS[method]
    notebook = json.loads((ROOT / (name + ".ipynb")).read_text(encoding="utf-8"))
    definitions = {}
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            for node in ast.parse("".join(cell["source"])).body:
                if isinstance(node, ast.FunctionDef):
                    definitions[node.name] = node
    env = {"np": np, "p": X.shape[1], "X": X, "I": np.eye(X.shape[1]),
           "delta": delta, "t4": options["step"], "kappa": options["kappa"]}
    names = (function_name, function_name + "_debias")
    if method == "ARk":
        names += ("Dk",)
    for name in names:
        node = copy.deepcopy(definitions[name])
        node.decorator_list = []
        exec(compile(ast.Module(body=[node], type_ignores=[]), "<legacy notebook>", "exec"), env)
    if method == "ARk":
        env["D"] = np.asarray([env["Dk"](k) for k in range(count + 1)])
    args = [X, Y, options["step"], count]
    if fixed:
        if method == "Nesterov":
            args.append(options["omega"])
        return env[function_name + "_debias"](*args)
    if method == "Landweber":
        args.append(truth)
    elif method in ("Showalter", "SOAR", "Nesterov"):
        args.append(delta)
        if method == "Nesterov":
            args.append(options["omega"])
    return env[function_name](*args)


def direct_ridge(X, Y, alpha):
    if alpha == 0:
        ordinary = np.linalg.lstsq(X, Y, rcond=None)[0]
        return ordinary, ordinary.copy()
    system = X.T @ X + alpha*np.eye(X.shape[1])
    ordinary = np.linalg.solve(system, X.T @ Y)
    refit = np.linalg.solve(system, X.T @ (X @ ordinary))
    return ordinary, 2*ordinary-refit


class SpectralFilterTests(unittest.TestCase):
    def test_discrete_filters_and_two_fits_match_legacy_at_all_counts(self):
        rng = np.random.default_rng(817)
        designs = [rng.normal(size=(6, 4))*.2, rng.normal(size=(4, 6))*.2,
                   np.array([[.5, 0., 0.], [0., .3, 0.], [0., 0., 0.], [0., 0., 0.]])]
        for X in designs:
            design = SpectralDesign(X)
            Y = rng.normal(size=X.shape[0])
            for method in METHODS:
                for count in (0, 1, 2, 7):
                    with self.subTest(shape=X.shape, method=method, count=count):
                        options = settings(method)
                        expected = legacy(method, X, Y, count, options)
                        second = legacy(method, X, X@expected, count, options)
                        original, debiased = design.pair(Y, method, count, options)
                        np.testing.assert_allclose(original, expected, rtol=2e-12, atol=2e-14)
                        np.testing.assert_allclose(debiased, 2*expected-second, rtol=2e-12, atol=2e-14)

    def test_sc_and_ridge_include_zero_rank_and_singular_zero_alpha(self):
        rng = np.random.default_rng(22)
        for X in (rng.normal(size=(7, 4)), rng.normal(size=(4, 7)),
                  np.array([[1., 0., 0.], [0., .5, 0.], [0., 0., 0.], [0., 0., 0.]])):
            Y = rng.normal(size=X.shape[0])
            U, s, Vt = np.linalg.svd(X, full_matrices=True)
            design = SpectralDesign(X, (U, s, Vt))
            for rank in (0, 1, len(s)):
                truncated = np.zeros_like(s)
                truncated[:rank] = s[:rank]
                sigma = np.zeros((X.shape[0], X.shape[1]))
                np.fill_diagonal(sigma, truncated)
                square = sigma.T@sigma
                inverse = np.zeros_like(square)
                np.divide(1., square, out=inverse, where=square != 0)
                expected = Vt.T@inverse@Vt@X.T@Y
                original, debiased = design.pair(Y, "SC", rank, {})
                np.testing.assert_allclose(original, expected, rtol=1e-12, atol=1e-13)
                np.testing.assert_allclose(debiased, expected, rtol=1e-12, atol=1e-13)
            for alpha in (0., .01, 10.):
                expected = direct_ridge(X, Y, alpha)
                actual = design.pair(Y, "Ridge", alpha, {})
                for a, b in zip(actual, expected):
                    np.testing.assert_allclose(a, b, rtol=3e-11, atol=1e-12)

        X = np.diag([1., np.finfo(float).eps/8, 0.])
        Y = np.ones(3)
        actual = SpectralDesign(X).pair(Y, "Ridge", 0, {})
        expected = np.linalg.lstsq(X, Y, rcond=None)[0]
        np.testing.assert_array_equal(actual[0], expected)
        np.testing.assert_array_equal(actual[1], expected)

    def test_batched_pair_matches_individual_legacy_fits(self):
        rng = np.random.default_rng(17)
        X = rng.normal(size=(4, 6))*.2
        Y = rng.normal(size=(4, 3))
        design = SpectralDesign(X)
        for method in METHODS:
            options = settings(method)
            original, debiased = design.pair(Y, method, 6, options)
            for j in range(3):
                expected = legacy(method, X, Y[:, j], 6, options)
                second = legacy(method, X, X@expected, 6, options)
                np.testing.assert_allclose(original[:, j], expected, rtol=2e-12, atol=2e-14)
                np.testing.assert_allclose(debiased[:, j], 2*expected-second, rtol=2e-12, atol=2e-14)

    def test_sc_pair_and_selection_are_exact_independent_copies(self):
        rng = np.random.default_rng(761)
        for X in (rng.normal(size=(7, 4)), rng.normal(size=(4, 7)),
                  np.diag([1., .1, 0.])):
            design = SpectralDesign(X)
            responses = rng.normal(size=(X.shape[0], 3))
            for Y in (responses[:, 0], responses):
                for rank in (0, 1, design.rank_size):
                    ordinary, debiased = design.pair(Y, "SC", rank, {})
                    np.testing.assert_array_equal(ordinary, debiased)
                    self.assertFalse(np.shares_memory(ordinary, debiased))
                    before = ordinary.copy()
                    debiased[...] = 0
                    np.testing.assert_array_equal(ordinary, before)
            truth = rng.normal(size=X.shape[1])
            ordinary, debiased, _, _ = design.select(responses[:, 0], "SC", {}, truth)
            np.testing.assert_array_equal(ordinary, debiased)
            self.assertFalse(np.shares_memory(ordinary, debiased))
            before = ordinary.copy()
            debiased[...] = 0
            np.testing.assert_array_equal(ordinary, before)

    def test_full_svd_is_trimmed_and_input_cache_is_not_mutable(self):
        for shape in ((5, 3), (3, 5)):
            X = np.random.default_rng(9).normal(size=shape)
            design = SpectralDesign(X, np.linalg.svd(X, full_matrices=True))
            self.assertEqual(design.U.shape, (shape[0], min(shape)))
            self.assertEqual(design.Vt.shape, (min(shape), shape[1]))
            np.testing.assert_allclose(design.U@np.diag(design.s)@design.Vt, X, atol=1e-13)
            X[:] = 0
            self.assertGreater(np.linalg.norm(design.X), 0)
            gain = design.filter("Ridge", .1, {})
            gain[:] = 0
            self.assertGreater(np.linalg.norm(design.filter("Ridge", .1, {})), 0)
            self.assertEqual(hash(design), hash(design))

    def test_default_steps_at_large_spectral_modes_and_nondefault_omega_kappa(self):
        X = np.diag(np.sqrt([537050.2834, 1000., 1.]))
        Y = np.array([100., -2., .1])
        design = SpectralDesign(X)
        for method, step in (("Landweber", 2.5e-6), ("Showalter", 2.5e-6),
                             ("ARk", 5e-5), ("SOAR", 2.5e-4),
                             ("Nesterov", .99/537050.2834)):
            for count in (1, 20, 750):
                options = settings(method, step=step, omega=5, kappa=1.5)
                expected = legacy(method, X, Y, count, options)
                second = legacy(method, X, X@expected, count, options)
                original, debiased = design.pair(Y, method, count, options)
                np.testing.assert_allclose(original, expected, rtol=3e-11, atol=2e-13)
                np.testing.assert_allclose(debiased, 2*expected-second, rtol=3e-11, atol=2e-13)

    def test_invalid_and_nonfinite_inputs_are_rejected(self):
        with self.assertRaises(FloatingPointError):
            SpectralDesign([[np.inf]])
        design = SpectralDesign(np.eye(2))
        for method, parameter, options in (("unknown", 1, {}), ("SC", 3, {}),
                                           ("Ridge", -.1, {}), ("Nesterov", 1.5, {"step": .1}),
                                           ("Landweber", 2, {"step": np.inf})):
            with self.subTest(method=method), self.assertRaises(ValueError):
                design.filter(method, parameter, options)
        with self.assertRaises(FloatingPointError):
            design.pair(np.array([np.nan, 1.]), "SC", 1, {})
        with self.assertRaises(FloatingPointError):
            design.filter("Landweber", 20, {"step": 1e30})


class SpectralSelectionTests(unittest.TestCase):
    def test_tdp_first_hit_and_cap_match_legacy_returned_iterates(self):
        X = np.array([[.4, .1], [.1, .3], [.1, -.1]])
        Y = np.array([1., -.5, .2])
        truth = np.array([.2, .4])
        design = SpectralDesign(X)
        for method in METHODS[1:]:
            options = settings(method)
            at_two = legacy(method, X, Y, 2, options)
            at_three = legacy(method, X, Y, 3, options)
            delta = (np.linalg.norm(Y-X@at_two)+np.linalg.norm(Y-X@at_three))/2
            for cap, threshold in ((0, delta), (1, 0.), (12, delta), (7, 0.)):
                with self.subTest(method=method, cap=cap, delta=threshold):
                    chosen = dict(options, max_steps=cap, delta=threshold)
                    expected, count = legacy(method, X, Y, cap, options, fixed=False, delta=threshold)
                    original, debiased, selected, diagnostics = design.select(Y, method, chosen, truth)
                    self.assertEqual(selected, count)
                    self.assertEqual(diagnostics["selected_count"], count)
                    self.assertEqual(diagnostics["at_iteration_cap"], count == cap)
                    self.assertEqual(diagnostics["tdp_satisfied"], np.linalg.norm(Y-X@original) <= threshold)
                    np.testing.assert_allclose(original, expected, rtol=2e-12, atol=2e-14)
                    second = legacy(method, X, X@expected, count, options)
                    np.testing.assert_allclose(debiased, 2*expected-second, rtol=2e-12, atol=2e-14)

    def test_landweber_aosr_first_turn_minimum_and_cap_match_legacy(self):
        X = np.ones((1, 1))
        Y = np.array([2.])
        options = settings("Landweber", step=.001, max_steps=700)
        design = SpectralDesign(X)
        for cap, reference_count in ((0, 10), (1, 10), (7, 600), (700, 10), (700, 600)):
            truth = legacy("Landweber", X, Y, reference_count, options)
            expected, count = legacy("Landweber", X, Y, cap, options, fixed=False, truth=truth)
            original, debiased, selected, diagnostics = design.select(
                Y, "Landweber", dict(options, max_steps=cap), truth)
            self.assertEqual(selected, count)
            np.testing.assert_allclose(original, expected, rtol=1e-12, atol=1e-14)
            expected_debiased = legacy("Landweber", X, Y, 2*count, options)
            np.testing.assert_allclose(debiased, expected_debiased, rtol=1e-12, atol=1e-14)
            self.assertIsNone(diagnostics["tdp_satisfied"])
        _, _, count, _ = design.select(np.zeros(1), "Landweber", options, np.zeros(1))
        self.assertEqual(count, 500)

    def test_oracle_grids_and_first_tie_conventions(self):
        rng = np.random.default_rng(451)
        X = rng.normal(size=(4, 6))
        Y = rng.normal(size=4)
        truth = rng.normal(size=6)
        design = SpectralDesign(X)
        sc_errors = [(np.linalg.norm(design.pair(Y, "SC", k, {})[0]-truth), k)
                     for k in range(4, 0, -1)]
        _, expected_rank = min(sc_errors, key=lambda item: item[0])
        _, _, rank, _ = design.select(Y, "SC", {"rank_max": 4}, truth)
        self.assertEqual(rank, expected_rank)
        ridge_options = {"C": 8, "m": 2, "alpha_divisor": X.shape[0]}
        errors = [(np.linalg.norm(direct_ridge(X, Y, c/X.shape[0])[0]-truth), c/X.shape[0])
                  for c in range(0, 9, 2)]
        _, expected_alpha = min(errors, key=lambda item: item[0])
        original, debiased, alpha, _ = design.select(Y, "Ridge", ridge_options, truth)
        self.assertEqual(alpha, expected_alpha)
        for actual, expected in zip((original, debiased), direct_ridge(X, Y, alpha)):
            np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-12)
        zeros = np.zeros(4)
        zero_truth = np.zeros(6)
        self.assertEqual(design.select(zeros, "SC", {}, zero_truth)[2], 4)
        self.assertEqual(design.select(zeros, "Ridge", ridge_options, zero_truth)[2], 0)

    def test_tdp_does_not_stop_before_first_update(self):
        design = SpectralDesign(np.eye(2))
        Y = np.ones(2)
        for method in METHODS[1:]:
            _, _, count, diagnostic = design.select(
                Y, method, settings(method, delta=100), np.zeros(2))
            self.assertEqual(count, 1)
            self.assertTrue(diagnostic["tdp_satisfied"])

    def test_default_scale_tdp_decisions_match_direct_discrete_recursions(self):
        X = np.diag(np.sqrt([537050.2834, 1000., 1.]))
        truth = np.array([.4, -.2, 1.])
        noise = np.array([2., -2., 2.])
        Y = X@truth+noise
        design = SpectralDesign(X)
        for method, step, tau in (("Showalter", 2.5e-6, .9), ("ARk", 5e-5, .9),
                                  ("SOAR", 2.5e-4, 1.01),
                                  ("Nesterov", .99/537050.2834, 1.01)):
            options = settings(method, step=step, max_steps=1000, delta=tau*np.linalg.norm(noise))
            expected, count = legacy(method, X, Y, 1000, options, fixed=False, delta=options["delta"])
            original, _, parameter, diagnostics = design.select(Y, method, options, truth)
            self.assertEqual(parameter, count)
            np.testing.assert_allclose(original, expected, rtol=3e-11, atol=2e-13)
            self.assertTrue(diagnostics["tdp_satisfied"] or diagnostics["at_iteration_cap"])


if __name__ == "__main__":
    unittest.main()
