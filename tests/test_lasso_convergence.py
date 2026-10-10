"""Convergence certificates, bounded retries, and the default-design regression."""

from pathlib import Path
import sys
import time
import unittest
from unittest.mock import patch
import warnings

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import Lasso
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from regression_common import lasso_coefficients


class FakeLasso:
    def __init__(self, outcomes, **parameters):
        self.outcomes = iter(outcomes)
        self.parameters = parameters
        self.budgets = []
        self.starting_coefficients = []

    def set_params(self, **parameters):
        self.parameters.update(parameters)
        return self

    def fit(self, X, Y):
        self.budgets.append(self.parameters["max_iter"])
        self.starting_coefficients.append(getattr(self, "coef_", None))
        gap, coefficients, warning = next(self.outcomes)
        self.coef_ = np.asarray(coefficients, dtype=float)
        self.dual_gap_ = gap
        self.n_iter_ = self.parameters["max_iter"]
        if warning:
            warnings.warn("mock convergence warning", ConvergenceWarning)
        return self


class LassoConvergenceTests(unittest.TestCase):
    def setUp(self):
        self.X = np.array([[1.], [-1.]])
        self.Y = np.array([2., -2.])

    def fake(self, outcomes):
        model = FakeLasso(outcomes)

        def factory(**parameters):
            model.parameters.update(parameters)
            return model

        return model, patch("sklearn.linear_model.Lasso", side_effect=factory)

    def test_warning_retries_from_previous_fit_without_changing_model(self):
        model, factory = self.fake([(1.0, [1.], True), (1e-7, [1.5], False)])
        with factory as constructor, warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = lasso_coefficients(self.X, self.Y, 0.5)
        self.assertEqual(constructor.call_count, 1)
        self.assertEqual(model.budgets, [10000, 30000])
        np.testing.assert_array_equal(model.starting_coefficients[1], [1.])
        np.testing.assert_array_equal(result, [1.5])
        self.assertFalse(caught)
        self.assertTrue(model.parameters["warm_start"])
        self.assertFalse(model.parameters["fit_intercept"])
        self.assertEqual(model.parameters["alpha"], 0.5)
        self.assertEqual(model.parameters["tol"], 1e-6)
        self.assertEqual(model.parameters["selection"], "cyclic")

    def test_large_gap_is_rejected_even_without_warning(self):
        model, factory = self.fake([(1.0, [1.], False)] * 3)
        with factory, self.assertRaisesRegex(RuntimeError, "140000.*No unconverged"):
            lasso_coefficients(self.X, self.Y, 0.5)
        self.assertEqual(model.budgets, [10000, 30000, 100000])

    def test_warning_is_rejected_even_when_reported_gap_is_small(self):
        model, factory = self.fake([(0.0, [1.], True)] * 3)
        with factory, self.assertRaisesRegex(RuntimeError, "did not converge"):
            lasso_coefficients(self.X, self.Y, 0.5)
        self.assertEqual(len(model.budgets), 3)

    def test_nonfinite_or_invalid_certificates_fail_immediately(self):
        for gap, coefficients in ((np.nan, [1.]), (np.inf, [1.]),
                                  (0.0, [np.nan]), (0.0, [np.inf]), (-1.0, [1.])):
            with self.subTest(gap=gap, coefficients=coefficients):
                model, factory = self.fake([(gap, coefficients, False)])
                with factory, self.assertRaises(FloatingPointError):
                    lasso_coefficients(self.X, self.Y, 0.5)
                self.assertEqual(model.budgets, [10000])

    def test_roundoff_negative_gap_and_zero_response_are_accepted(self):
        model, factory = self.fake([(-np.finfo(float).eps, [0.], False)])
        with factory:
            np.testing.assert_array_equal(lasso_coefficients(self.X, np.zeros(2), .5), [0.])
        self.assertEqual(model.budgets, [10000])
        np.testing.assert_array_equal(lasso_coefficients(self.X, np.zeros(2), .5), [0.])

    def test_alpha_zero_preserves_minimum_norm_least_squares(self):
        X = np.array([[1., 1.], [-1., -1.]])
        with patch("sklearn.linear_model.Lasso") as constructor:
            np.testing.assert_allclose(lasso_coefficients(X, self.Y, 0), [1., 1.])
            constructor.assert_not_called()
        with patch("numpy.linalg.lstsq", return_value=(np.array([np.nan]),)):
            with self.assertRaises(FloatingPointError):
                lasso_coefficients(self.X, self.Y, 0)

    def test_default_size_design_converges_after_old_budget_fails(self):
        with threadpool_limits(limits=1):
            rng = np.random.RandomState(24)
            covariance = np.full((1000, 1000), 0.5)
            np.fill_diagonal(covariance, 2.0)
            X = rng.multivariate_normal(np.zeros(1000), covariance, 1000)
            beta = np.random.RandomState(25).uniform(-2, 2, 1000)
            Y = X @ beta + np.random.RandomState(26).normal(0, 2, 1000)
            legacy = Lasso(alpha=.018, fit_intercept=False, max_iter=10000, tol=1e-6)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", ConvergenceWarning)
                legacy.fit(X, Y)
            self.assertTrue(any(issubclass(w.category, ConvergenceWarning) for w in caught))
            threshold = 1e-6 * np.dot(Y, Y) / Y.size
            self.assertGreater(legacy.dual_gap_, threshold)

            models = []

            def factory(**parameters):
                model = Lasso(**parameters)
                models.append(model)
                return model

            started = time.monotonic()
            with patch("sklearn.linear_model.Lasso", side_effect=factory), \
                    warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = lasso_coefficients(X, Y, .018)
            self.assertEqual(len(models), 1)
            self.assertFalse(any(issubclass(w.category, ConvergenceWarning) for w in caught))
            self.assertGreater(models[0].max_iter, 10000)
            self.assertLessEqual(models[0].dual_gap_, threshold)
            self.assertTrue(np.all(np.isfinite(result)))
            np.testing.assert_array_equal(result, models[0].coef_)
            self.assertGreater(np.linalg.norm(result-legacy.coef_), .01)
            print(f"Default Lasso: legacy gap={legacy.dual_gap_:.9g}, "
                  f"converged gap={models[0].dual_gap_:.9g}, threshold={threshold:.9g}, "
                  f"last budget={models[0].max_iter}, "
                  f"last iterations={models[0].n_iter_}, "
                  f"seconds={time.monotonic()-started:.2f}")


if __name__ == "__main__":
    unittest.main()
