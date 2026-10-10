"""Optional HBF Stormer--Verlet against the paper's matrix scheme (45)."""

from pathlib import Path
import sys
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from regression_spectral import SpectralDesign


def matrix_path(X, Y, count, h, eta):
    """Equation (45), directly in coefficient space, with zero initial data."""
    shape = (X.shape[1],) + Y.shape[1:]
    position = np.zeros(shape)
    velocity = np.zeros(shape)
    positions = [position.copy()]
    for _ in range(count):
        half = (velocity + h/2 * (X.T @ (Y-X @ position))) / (1+eta*h/2)
        position = position + h*half
        velocity = half + h/2 * (X.T @ (Y-X @ position)-eta*half)
        positions.append(position.copy())
    return positions


def three_term_path(X, Y, count, h, eta):
    """Eliminate velocity; derive the special first step from zero velocity."""
    previous = np.zeros((X.shape[1],) + Y.shape[1:])
    positions = [previous.copy()]
    if count == 0:
        return positions
    current = h*h/(2+eta*h) * (X.T @ Y)
    positions.append(current.copy())
    alpha = 2*h*h/(2+eta*h)
    momentum = (2-eta*h)/(2+eta*h)
    for _ in range(1, count):
        following = current + momentum*(current-previous) + alpha*(X.T @ (Y-X @ current))
        previous, current = current, following
        positions.append(current.copy())
    return positions


class HbfVerletTests(unittest.TestCase):
    def test_fixed_counts_and_debias_match_matrix_scheme(self):
        rng = np.random.default_rng(729)
        designs = (rng.normal(size=(7, 4))*.3, rng.normal(size=(4, 7))*.3,
                   np.array([[.8, 0., 0.], [0., .4, 0.], [0., 0., 0.], [0., 0., 0.]]))
        options = {"step": .05, "eta": 2., "integrator": "stormer-verlet"}
        for X in designs:
            design = SpectralDesign(X)
            batch = rng.normal(size=(X.shape[0], 3))
            for Y in (batch[:, 0], batch):
                expected_path = matrix_path(X, Y, 20, options["step"], options["eta"])
                for count in (0, 1, 2, 7, 20):
                    with self.subTest(shape=X.shape, response_ndim=Y.ndim, count=count):
                        expected = expected_path[count]
                        refit = matrix_path(X, X@expected, count, options["step"], options["eta"])[-1]
                        ordinary, debiased = design.pair(Y, "HBF", count, options)
                        np.testing.assert_allclose(ordinary, expected, rtol=3e-12, atol=3e-14)
                        np.testing.assert_allclose(debiased, 2*expected-refit, rtol=3e-12, atol=3e-14)

    def test_eliminated_three_term_and_special_first_step(self):
        X = np.array([[1., .3], [.2, .8], [.1, -.1]])
        Y = np.array([.5, -.2, 1.])
        h, eta, count = .1, 2., 40
        options = {"step": h, "eta": eta, "integrator": "stormer-verlet"}
        design = SpectralDesign(X)
        expected = three_term_path(X, Y, count, h, eta)
        direct = matrix_path(X, Y, count, h, eta)
        np.testing.assert_allclose(expected, direct, rtol=2e-12, atol=2e-14)
        np.testing.assert_allclose(design.filter("HBF", 1, options),
                                   h*h/(2+eta*h)*design.s, rtol=1e-15)
        for k in (1, 2, 10, count):
            np.testing.assert_allclose(design.pair(Y, "HBF", k, options)[0],
                                       expected[k], rtol=2e-12, atol=2e-14)

    def test_first_tdp_hit_and_fixed_selected_debias(self):
        X = np.array([[.4, .1], [.1, .3], [.1, -.1]])
        Y = np.array([1., -.5, .2])
        h, eta = .05, 5.
        path = matrix_path(X, Y, 12, h, eta)
        residuals = [np.linalg.norm(Y-X@beta) for beta in path]
        threshold = (residuals[2]+residuals[3])/2
        design = SpectralDesign(X)
        for cap, delta in ((0, threshold), (1, 0.), (12, threshold), (7, 0.), (12, 100.)):
            with self.subTest(cap=cap, delta=delta):
                expected_count = next((k for k in range(1, cap+1) if residuals[k] <= delta), cap)
                options = {"step": h, "eta": eta, "integrator": "stormer-verlet",
                           "max_steps": cap, "delta": delta, "selection": "TDP"}
                ordinary, debiased, count, diagnostics = design.select(Y, "HBF", options, np.zeros(2))
                self.assertEqual(count, expected_count)
                if cap == 12:
                    self.assertEqual(count, 1 if delta == 100. else 3)
                self.assertEqual(diagnostics["at_iteration_cap"], count == cap)
                self.assertEqual(diagnostics["tdp_satisfied"], np.linalg.norm(Y-X@ordinary) <= delta)
                expected = path[count]
                refit = matrix_path(X, X@expected, count, h, eta)[-1]
                np.testing.assert_allclose(ordinary, expected, rtol=3e-12, atol=3e-14)
                np.testing.assert_allclose(debiased, 2*expected-refit, rtol=3e-12, atol=3e-14)

    def test_caches_are_isolated_and_extend_only_as_needed(self):
        design = SpectralDesign(np.diag([.8, .3]))
        options = {"step": .1, "eta": 2.}
        verlet = dict(options, integrator="stormer-verlet")
        with mock.patch.object(design, "_advance_hbf", wraps=design._advance_hbf) as rk4, \
                mock.patch.object(design, "_advance_hbf_verlet", wraps=design._advance_hbf_verlet) as sv:
            original = design.filter("HBF", 6, options)
            np.testing.assert_array_equal(original, design.filter("HBF", 6, dict(options, integrator="rk4")))
            self.assertEqual(rk4.call_count, 6)
            self.assertEqual(sv.call_count, 0)
            different = design.filter("HBF", 6, verlet)
            self.assertGreater(np.linalg.norm(original-different), 1e-7)
            self.assertEqual(sv.call_count, 6)
            design.filter("HBF", 3, verlet)
            design.filter("HBF", 8, verlet)
            self.assertEqual(sv.call_count, 8)
            self.assertEqual(rk4.call_count, 6)
            self.assertEqual(len(design._trajectories), 2)
            np.testing.assert_array_equal(original, design.filter("HBF", 6, options))
        fresh = SpectralDesign(np.eye(2))
        with mock.patch.object(fresh, "_advance_hbf_verlet", wraps=fresh._advance_hbf_verlet) as sv:
            _, _, count, _ = fresh.select(np.ones(2), "HBF", dict(verlet, delta=100., max_steps=1000), np.zeros(2))
            self.assertEqual(count, 1)
            self.assertEqual(sv.call_count, 1)

    def test_stability_bound_and_invalid_integrators(self):
        design = SpectralDesign(np.diag([4., .5]))
        for eta, bound in ((1., np.sqrt(2)/4), (10., .2)):
            for step in (bound, np.nextafter(bound, np.inf)):
                result = design.filter("HBF", 3, {"step": step, "eta": eta, "integrator": "stormer-verlet"})
                self.assertTrue(np.isfinite(result).all())
            with self.assertRaisesRegex(ValueError, "step must satisfy"):
                design.filter("HBF", 3, {"step": bound*(1+1e-8), "eta": eta, "integrator": "stormer-verlet"})
        for eta in (0., -1.):
            with self.assertRaises(ValueError):
                design.filter("HBF", 1, {"step": .1, "eta": eta, "integrator": "stormer-verlet"})
        for integrator in (None, "unknown", "verlet", 1):
            with self.assertRaisesRegex(ValueError, "integrator"):
                design.filter("HBF", 1, {"step": .1, "integrator": integrator})
        zero = SpectralDesign(np.zeros((2, 3)))
        np.testing.assert_array_equal(zero.filter("HBF", 5,
            {"step": 2., "eta": 1., "integrator": "stormer-verlet"}), np.zeros(2))
        with self.assertRaisesRegex(ValueError, "step must satisfy"):
            zero.filter("HBF", 1, {"step": 2.1, "eta": 1., "integrator": "stormer-verlet"})

    def test_default_rk4_keeps_legacy_keys_and_eta_zero(self):
        design = SpectralDesign(np.diag([2., .5]))
        options = {"step": .1, "eta": 0.}
        self.assertEqual(design._key("HBF", 3, options), ("HBF", 3, .1, 0.))
        self.assertEqual(design._key("HBF", 3, dict(options, integrator="rk4")), ("HBF", 3, .1, 0.))
        gain = design.filter("HBF", 3, options)
        explicit = SpectralDesign(design.X).filter("HBF", 3, dict(options, integrator="rk4"))
        np.testing.assert_array_equal(gain, explicit)


if __name__ == "__main__":
    unittest.main()
