import unittest
from unittest.mock import patch

import numpy as np

import regression_threshold_scan as grid_module
from regression_bootstrap import bootstrap_errors
from regression_execution import threshold_pair
from regression_spectral import SpectralDesign
from regression_threshold_scan import bootstrap_errors_grid, threshold_grid_pair


class ThresholdGridTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(75)
        self.design = SpectralDesign(rng.normal(size=(8, 11)))
        self.truth = np.zeros(11)
        self.truth[:2] = (1, -1)
        self.estimates = rng.normal(size=(11, 2)) * .3
        self.percentiles = [0, 25, 50, 75, 100]

    def test_thresholds_curves_and_centers_match_individual_oracle_calls(self):
        original = self.estimates.copy()
        centers, thresholds, scan = threshold_grid_pair(self.estimates, self.truth, self.percentiles)
        for index, percentile in enumerate(self.percentiles):
            expected = threshold_pair(self.estimates, self.truth, percentile)
            np.testing.assert_array_equal(centers[index], expected[0])
            np.testing.assert_array_equal(thresholds[index], expected[1])
            np.testing.assert_array_equal(scan, expected[2])
        np.testing.assert_array_equal(self.estimates, original)

    def test_complete_B_errors_match_each_existing_bootstrap_call(self):
        centers = np.broadcast_to(self.estimates, (3, 11, 2)).copy()
        centers[2, 0, 0] += .1
        sigma = np.array([[.5, .6], [.5, .6], [.7, .8]])
        thresholds = np.array([[.1, .2], [.3, .4], [.2, .3]])
        methods = dict(SC=4, Ridge=.3, Landweber=6, Showalter=6, ARk=6,
                       SOAR=6, Nesterov=6, HBF=6, FAR=6)
        for method, parameter in methods.items():
            settings = dict(step=1e-4, eta=5, theta=1.5, kappa=.5, s_star=.5, omega=3)
            expected = np.stack([bootstrap_errors(self.design, method, parameter, settings,
                center, scale, threshold, 7, 53, 4, batch_size=3, cpu_workers=1, blas_threads=1)
                for center,scale,threshold in zip(centers,sigma,thresholds)])
            for workers in (1, 2):
                actual = bootstrap_errors_grid(self.design, method, parameter, settings,
                    centers,sigma,thresholds,7,53,4,batch_size=3,cpu_workers=workers,blas_threads=1)
                np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)

    def test_equal_centers_and_noise_reuse_unthresholded_fits(self):
        centers = np.broadcast_to(self.estimates, (5, 11, 2)).copy()
        sigma = np.full((5, 2), .5)
        thresholds = np.arange(10).reshape(5, 2) * .05
        with patch.object(grid_module, "_cpu_grid_chunk", wraps=grid_module._cpu_grid_chunk) as calls:
            result = bootstrap_errors_grid(self.design,"Ridge",.3,{},centers,sigma,thresholds,
                                           7,53,4,batch_size=3,cpu_workers=1,blas_threads=1)
        self.assertEqual(calls.call_count, 3)
        self.assertEqual(result.shape, (5, 7, 2))

    def test_zero_noise_and_zero_iteration_are_finite_and_do_not_mutate_inputs(self):
        centers = np.broadcast_to(self.estimates, (2, 11, 2)).copy()
        initial = centers.copy()
        sigma = np.zeros((2, 2))
        thresholds = np.zeros((2, 2))
        actual = bootstrap_errors_grid(self.design,"HBF",0,{"step":1e-4,"eta":5},
                                      centers,sigma,thresholds,4,5,0,blas_threads=1)
        expected = np.max(np.abs(centers), axis=1)
        np.testing.assert_array_equal(actual, np.broadcast_to(expected[:,None,:],(2,4,2)))
        np.testing.assert_array_equal(centers, initial)

    def test_invalid_percentiles_shapes_and_nonfinite_inputs_are_rejected(self):
        for percentiles in ([], [0, 0], [-1], [101], [np.nan], [50+1j]):
            with self.assertRaises(ValueError):
                threshold_grid_pair(self.estimates, self.truth, percentiles)
        centers = np.broadcast_to(self.estimates, (2, 11, 2)).copy()
        for sigma, thresholds in ((np.zeros((3,2)),np.zeros((2,2))),
                                   (np.full((2,2),-1.),np.zeros((2,2))),
                                   (np.zeros((2,2)),np.full((2,2),np.nan))):
            with self.assertRaises(ValueError):
                bootstrap_errors_grid(self.design,"Ridge",.3,{},centers,sigma,thresholds,4,5,0)


if __name__ == "__main__":
    unittest.main()
