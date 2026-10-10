"""Landweber TDP retains the discrete estimator and fixed-count bootstrap."""

from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from regression_bootstrap import bootstrap_errors
from regression_spectral import SpectralDesign


def direct_landweber(X, Y, step, count):
    estimate = np.zeros((X.shape[1],) + Y.shape[1:])
    for _ in range(count):
        estimate = estimate + step * (X.T @ (Y - X @ estimate))
    return estimate


def direct_pair(X, Y, step, count):
    ordinary = direct_landweber(X, Y, step, count)
    refit = direct_landweber(X, X @ ordinary, step, count)
    return ordinary, 2 * ordinary - refit


def direct_tdp(X, Y, step, cap, delta):
    estimate = np.zeros(X.shape[1])
    count = 0
    for count in range(1, cap + 1):
        estimate = estimate + step * (X.T @ (Y - X @ estimate))
        if np.linalg.norm(Y - X @ estimate) <= delta:
            break
    return estimate, count


def small_designs():
    rng = np.random.default_rng(608)
    U, _ = np.linalg.qr(rng.normal(size=(6, 4)))
    V, _ = np.linalg.qr(rng.normal(size=(4, 4)))
    return (np.diag([2., 1., .25, 0.]),
            U @ np.diag([20., 2., .01, 1e-7]) @ V.T,
            rng.normal(size=(3, 5)))


class LandweberTDPTests(unittest.TestCase):
    def test_fixed_ordinary_and_debiased_fits_match_discrete_recurrence(self):
        rng = np.random.default_rng(771)
        for X in small_designs():
            design = SpectralDesign(X)
            step = .8 / design.s[0]**2
            responses = rng.normal(size=(X.shape[0], 3))
            for Y in (responses[:, 0], responses):
                for count in (0, 1, 2, 9, 120):
                    with self.subTest(shape=X.shape, response_shape=Y.shape, count=count):
                        actual = design.pair(Y, "Landweber", count, {"step": step})
                        expected = direct_pair(X, Y, step, count)
                        for fit, reference in zip(actual, expected):
                            np.testing.assert_allclose(fit, reference, rtol=3e-11, atol=3e-13)

    def test_tdp_first_hit_matches_direct_discrete_recurrence(self):
        rng = np.random.default_rng(377)
        for X in small_designs():
            design = SpectralDesign(X)
            step = .8 / design.s[0]**2
            Y = rng.normal(size=X.shape[0])
            residuals = [np.linalg.norm(Y - X @ direct_landweber(X, Y, step, count))
                         for count in (3, 4)]
            delta = sum(residuals) / 2
            options = {"step": step, "max_steps": 12, "selection": "TDP",
                       "delta": delta, "minimum_steps": 500}
            ordinary, debiased, count, diagnostics = design.select(
                Y, "Landweber", options, np.zeros(X.shape[1]))
            expected, expected_count = direct_tdp(X, Y, step, 12, delta)
            self.assertEqual(count, 4)
            self.assertEqual(count, expected_count)
            self.assertEqual(diagnostics["selected_count"], count)
            self.assertFalse(diagnostics["at_iteration_cap"])
            self.assertTrue(diagnostics["tdp_satisfied"])
            self.assertEqual(diagnostics["stopping_threshold"], delta)
            self.assertEqual(diagnostics["parameter_kind"], "steps")
            np.testing.assert_allclose(ordinary, expected, rtol=3e-11, atol=3e-13)
            fixed = direct_pair(X, Y, step, count)
            np.testing.assert_allclose(debiased, fixed[1], rtol=3e-11, atol=3e-13)
            changed_truth = design.select(Y, "Landweber", options, np.full(X.shape[1], 100.))
            self.assertEqual(changed_truth[2], count)
            np.testing.assert_array_equal(changed_truth[0], ordinary)
            np.testing.assert_array_equal(changed_truth[1], debiased)

    def test_tdp_starts_at_update_one_and_supports_zero_caps_and_zero_response(self):
        X = np.eye(2)
        design = SpectralDesign(X)
        for Y, cap, delta, expected_count, satisfied in (
                (np.ones(2), 12, 100., 1, True),
                (np.zeros(2), 12, 0., 1, True),
                (np.ones(2), 0, 100., 0, True),
                (np.ones(2), 0, 0., 0, False),
                (np.zeros(2), 0, 0., 0, True),
                (np.ones(2), 1, 0., 1, False),
                (np.ones(2), 7, 0., 7, False)):
            with self.subTest(Y=Y, cap=cap, delta=delta):
                options = {"selection": "TDP", "step": .25,
                           "max_steps": cap, "delta": delta}
                ordinary, debiased, count, diagnostics = design.select(
                    Y, "Landweber", options, np.zeros(2))
                expected, expected_direct_count = direct_tdp(X, Y, .25, cap, delta)
                self.assertEqual(count, expected_count)
                self.assertEqual(count, expected_direct_count)
                self.assertEqual(diagnostics["at_iteration_cap"], count == cap)
                self.assertEqual(diagnostics["tdp_satisfied"], satisfied)
                np.testing.assert_allclose(ordinary, expected, rtol=1e-13, atol=1e-14)
                np.testing.assert_allclose(debiased, direct_pair(X, Y, .25, count)[1],
                                           rtol=1e-13, atol=1e-14)

    def test_tdp_requires_valid_delta_and_selection(self):
        design = SpectralDesign(np.eye(2))
        options = {"selection": "TDP", "step": .25, "max_steps": 0}
        with self.assertRaisesRegex(ValueError, "requires delta"):
            design.select(np.ones(2), "Landweber", options, np.zeros(2))
        for delta in (None, True, -1., np.nan, np.inf, "1", 1j):
            with self.subTest(delta=delta), self.assertRaises(ValueError):
                design.select(np.ones(2), "Landweber", dict(options, delta=delta), np.zeros(2))
        for selection in (None, "tdp", "oracle_rank", "unknown"):
            with self.subTest(selection=selection), self.assertRaisesRegex(ValueError, "AOSR or TDP"):
                design.select(np.ones(2), "Landweber", dict(options, selection=selection), np.zeros(2))
        for method, selection in (("SC", "TDP"), ("Ridge", "AOSR"),
                                  ("Showalter", "AOSR"), ("Nesterov", "AOSR")):
            with self.subTest(method=method), self.assertRaisesRegex(ValueError, "selection must be"):
                design.select(np.ones(2), method, dict(options, selection=selection), np.zeros(2))

    def test_implicit_and_explicit_aosr_preserve_minimum_and_diagnostics(self):
        design = SpectralDesign(np.ones((1, 1)))
        options = {"step": .001, "max_steps": 700}
        default = design.select(np.zeros(1), "Landweber", options, np.zeros(1))
        explicit = design.select(np.zeros(1), "Landweber", dict(options, selection="AOSR"), np.zeros(1))
        self.assertEqual(default[2], 500)
        self.assertEqual(default[2:], explicit[2:])
        np.testing.assert_array_equal(default[0], explicit[0])
        np.testing.assert_array_equal(default[1], explicit[1])
        self.assertIsNone(default[3]["tdp_satisfied"])
        self.assertNotIn("stopping_threshold", default[3])

    def test_bootstrap_uses_selected_fixed_count_for_both_fits(self):
        X = np.array([[2., .1], [.1, 1.], [.2, -.3]])
        Y = np.array([1., -.5, .1])
        design = SpectralDesign(X)
        step = .8 / design.s[0]**2
        residuals = [np.linalg.norm(Y - X @ direct_landweber(X, Y, step, count))
                     for count in (2, 3)]
        options = {"selection": "TDP", "step": step, "max_steps": 30,
                   "delta": sum(residuals) / 2}
        ordinary, debiased, count, _ = design.select(Y, "Landweber", options, np.zeros(2))
        self.assertEqual(count, 3)
        centers = np.column_stack((ordinary, debiased))
        sigma = np.array([.2, .3])
        thresholds = np.array([.02, .04])
        B, seed_root, outer_index = 7, 20261007, 2
        actual = bootstrap_errors(design, "Landweber", count, options, centers,
                                  sigma, thresholds, B, seed_root, outer_index,
                                  backend="cpu", batch_size=2, cpu_workers=2, blas_threads=1)
        expected = np.empty((B, 2))
        for replicate in range(B):
            rng = np.random.default_rng((seed_root, outer_index, 1, replicate))
            for branch in range(2):
                response = X @ centers[:, branch] + rng.normal(0., sigma[branch], X.shape[0])
                fitted = direct_pair(X, response, step, count)[branch]
                fitted[np.abs(fitted) <= thresholds[branch]] = 0
                expected[replicate, branch] = np.max(np.abs(fitted - centers[:, branch]))
        np.testing.assert_allclose(actual, expected, rtol=3e-12, atol=3e-14)


if __name__ == "__main__":
    unittest.main()
