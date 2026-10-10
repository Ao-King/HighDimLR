"""Batched bootstrap equivalence, RNG preservation, and optional real CUDA tests."""

from pathlib import Path
import sys
import threading
import types
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from regression_bootstrap import bootstrap_errors, validate_backend
from regression_spectral import SpectralDesign
from test_notebook_regressions import BOOTSTRAPS, bootstrap_environment


METHODS = {"BootstrapARk": "ARk", "BootstrapNesterov": "Nesterov",
           "Bootstrapsecond": "SOAR", "Bootstraplandweber": "Landweber",
           "Bootstrapshowalter": "Showalter", "Bootstrapcut": "SC", "Bootstrapridge": "Ridge"}


class NumpyCupyStub:
    """Exercise device dispatch/caching only; this is not CUDA hardware testing."""

    __version__ = "numpy-control-flow-stub"
    float64 = np.float64
    ones = staticmethod(np.ones)
    abs = staticmethod(np.abs)
    max = staticmethod(np.max)
    where = staticmethod(np.where)
    column_stack = staticmethod(np.column_stack)

    def __init__(self):
        self.local = threading.local()
        self.lock = threading.Lock()
        self.uploads = []
        self.downloads = []
        self.synchronizations = []
        owner = self

        class Device:
            def __init__(self, index):
                self.index = index

            def __enter__(self):
                self.previous = getattr(owner.local, "device", None)
                owner.local.device = self.index
                return self

            def __exit__(self, *args):
                owner.local.device = self.previous

        runtime = types.SimpleNamespace(getDeviceCount=lambda: 2,
                                        getDeviceProperties=lambda index: {"name": f"SIMULATED {index}"},
                                        runtimeGetVersion=lambda: 0)
        stream = types.SimpleNamespace(synchronize=self.synchronize)
        self.cuda = types.SimpleNamespace(Device=Device, runtime=runtime, get_current_stream=lambda: stream)

    def asarray(self, value):
        with self.lock:
            self.uploads.append((self.local.device, np.shape(value)))
        return np.array(value, copy=True)

    def asnumpy(self, value):
        with self.lock:
            self.downloads.append((self.local.device, np.shape(value)))
        return np.array(value, copy=True)

    def synchronize(self):
        with self.lock:
            self.synchronizations.append(self.local.device)


def legacy_case(name, count=4, B=9, outer_index=2):
    env, _ = bootstrap_environment(name, count=count)
    label, sigma_label, replay, step_name, _, _ = BOOTSTRAPS[name]
    thresholds = np.array([.11, .21])
    centers = np.column_stack((env["beta_" + label], env["beta_" + label + "_debias"]))
    centers[np.abs(centers) <= thresholds] = 0
    sigma = np.array([.17, .41])
    for index, suffix in enumerate(("", "_debias")):
        env["beta_" + label + suffix] = centers[:, index].copy()
        env["sigma_" + sigma_label + suffix] = sigma[index]
        env["bn_" + label + suffix] = thresholds[index]
    method = METHODS[name]
    if method == "SC":
        parameter, settings = env["best_k"], {}
    elif method == "Ridge":
        parameter, settings = env["best_c"] / env["n"], {}
    else:
        parameter, settings = count, {"step": env[step_name]}
        if method == "ARk":
            settings["kappa"] = env["kappa"]
        elif method == "Nesterov":
            settings["omega"] = env["omega"]
        elif method == "SOAR":
            settings["s_star"] = .5
    reference = env["compute_norminf_" + label]
    expected = np.asarray([reference((20261006, outer_index, 1, seed)) for seed in range(B)])
    design = SpectralDesign(env["X"], svd=(env["U"], env["S"], env["Vt"]))
    return design, method, parameter, settings, centers, sigma, thresholds, expected


class RegressionBootstrapTests(unittest.TestCase):
    def test_all_seven_methods_match_legacy_per_seed_fits(self):
        for name in METHODS:
            case = legacy_case(name)
            design, method, parameter, settings, centers, sigma, thresholds, expected = case
            for batch, workers in ((1, 1), (4, 3), (64, 2)):
                with self.subTest(method=method, batch=batch, workers=workers):
                    actual = bootstrap_errors(design, method, parameter, settings, centers, sigma,
                                              thresholds, 9, 20261006, 2, batch_size=batch,
                                              cpu_workers=workers, blas_threads=1)
                    self.assertEqual(actual.dtype, np.float64)
                    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-13)

    def test_zero_count_fixed_fit_matches_legacy(self):
        for name in METHODS:
            if METHODS[name] in ("SC", "Ridge"):
                continue
            design, method, parameter, settings, centers, sigma, thresholds, expected = legacy_case(name, count=0, B=3)
            with self.subTest(method=method):
                actual = bootstrap_errors(design, method, parameter, settings, centers, sigma,
                                          thresholds, 3, 20261006, 2, batch_size=2, blas_threads=1)
                np.testing.assert_array_equal(actual, expected)

    def test_seed_partition_preserves_random_stream_and_inputs(self):
        design, method, parameter, settings, centers, sigma, thresholds, expected = legacy_case("Bootstrapridge")
        snapshots = [value.copy() for value in (design.X, design.U, design.s, design.Vt, centers, sigma, thresholds)]
        state = np.random.get_state()
        actual = bootstrap_errors(design, method, parameter, settings, centers, sigma,
                                  thresholds, 9, 20261006, 2, batch_size=3, blas_threads=1)
        current = np.random.get_state()
        self.assertEqual(state[0], current[0])
        np.testing.assert_array_equal(state[1], current[1])
        self.assertEqual(state[2:], current[2:])
        for original, snapshot in zip((design.X, design.U, design.s, design.Vt, centers, sigma, thresholds), snapshots):
            np.testing.assert_array_equal(original, snapshot)
        np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-13)
        next_outer = bootstrap_errors(design, method, parameter, settings, centers, sigma,
                                      thresholds, 9, 20261006, 3, batch_size=3, blas_threads=1)
        self.assertFalse(np.array_equal(actual, next_outer))

    def test_threshold_closed_boundary_is_in_original_coordinates(self):
        design = SpectralDesign(np.eye(2))
        centers = np.array([[.2, .4], [.7, .8]])
        result = bootstrap_errors(design, "SC", 2, {}, centers, [0., 0.], [.2, .4],
                                  5, 123, 0, batch_size=2, blas_threads=1)
        np.testing.assert_array_equal(result, np.tile([.2, .4], (5, 1)))

    def test_sc_branches_are_exactly_equal_for_identical_responses(self):
        rng = np.random.default_rng(631)
        design = SpectralDesign(rng.normal(size=(13, 8)))
        center = rng.normal(size=8)
        centers = np.column_stack((center, center))
        shared_noise = rng.normal(size=(13, 7))
        with mock.patch("regression_bootstrap._noise", return_value=np.stack((shared_noise, shared_noise))):
            for method in ("SC", "cut"):
                with self.subTest(method=method):
                    errors = bootstrap_errors(design, method, 5, {}, centers, [.2, .2], [.1, .1],
                                              7, 123, 0, batch_size=3, cpu_workers=2, blas_threads=1)
                    np.testing.assert_array_equal(errors[:, 0], errors[:, 1])

    def test_cpu_backend_does_not_import_cupy(self):
        with mock.patch("regression_bootstrap.importlib.import_module", side_effect=ModuleNotFoundError("not installed")):
            self.assertEqual(validate_backend()["backend"], "cpu")
            with self.assertRaisesRegex(RuntimeError, "CuPy.*no CPU fallback"):
                validate_backend("gpu")

    def test_invalid_devices_fail_explicitly(self):
        for devices in ((), (0, 0), (-1,), (False,)):
            with self.subTest(devices=devices), self.assertRaises(ValueError):
                validate_backend("gpu", devices)
        fake = mock.Mock()
        fake.cuda.runtime.getDeviceCount.return_value = 1
        with mock.patch("regression_bootstrap._cupy", return_value=fake):
            with self.assertRaisesRegex(RuntimeError, "only 1 device"):
                validate_backend("gpu", (2,))

    def test_invalid_inputs_are_rejected(self):
        design = SpectralDesign(np.eye(2))
        arguments = dict(design=design, method="Ridge", parameter=.2, settings={},
                         centers=np.zeros((2, 2)), sigma=[1., 1.], thresholds=[0., 0.],
                         B=3, seed_root=1, outer_index=0, blas_threads=1)
        for key, value in (("centers", np.zeros((2, 1))), ("sigma", [np.nan, 1.]),
                           ("sigma", [-1., 1.]), ("thresholds", [-1., 0.]),
                           ("B", 0), ("batch_size", 0), ("cpu_workers", False),
                           ("blas_threads", 0), ("outer_index", -1), ("backend", "auto")):
            with self.subTest(key=key), self.assertRaises(ValueError):
                bootstrap_errors(**{**arguments, key: value})

    def test_nonfinite_filter_or_arithmetic_is_not_silently_returned(self):
        design = SpectralDesign(np.eye(2))
        with mock.patch.object(design, "filter", return_value=np.array([np.nan, 1.])):
            with self.assertRaises(ValueError):
                bootstrap_errors(design, "Ridge", .1, {}, np.zeros((2, 2)), [1., 1.], [0., 0.], 2, 1, 0)
        with mock.patch.object(design, "filter", return_value=np.full(2, 1e200)):
            with self.assertRaises(FloatingPointError):
                bootstrap_errors(design, "Ridge", .1, {}, np.zeros((2, 2)), [1., 1.], [0., 0.], 2, 1, 0)

    def test_gpu_worker_error_does_not_fall_back_to_cpu(self):
        design = SpectralDesign(np.eye(2))
        with mock.patch("regression_bootstrap.validate_backend", return_value={"backend": "gpu"}), \
                mock.patch("regression_bootstrap._cupy", return_value=mock.Mock()), \
                mock.patch("regression_bootstrap._gpu_worker", side_effect=RuntimeError("simulated device failure")), \
                mock.patch("regression_bootstrap._cpu_chunk") as cpu:
            with self.assertRaisesRegex(RuntimeError, "device failure.*no CPU fallback"):
                bootstrap_errors(design, "Ridge", .1, {}, np.zeros((2, 2)), [1., 1.], [0., 0.],
                                 4, 1, 0, backend="gpu", gpu_devices=(0, 1), batch_size=1)
            cpu.assert_not_called()

    def test_numpy_stub_exercises_two_device_batches_and_fresh_centers_and_filters(self):
        rng = np.random.default_rng(955)
        design = SpectralDesign(rng.normal(size=(13, 8)))
        centers = rng.normal(size=(8, 2))
        stub = NumpyCupyStub()
        for call, alpha in enumerate((.4, 2.2)):
            current_centers = centers if call == 0 else -2*centers + .3
            arguments = (design, "Ridge", alpha, {}, current_centers, [.2, .4], [.12, .21],
                         11, 20261006, call)
            expected = bootstrap_errors(*arguments, backend="cpu", batch_size=2,
                                        cpu_workers=2, blas_threads=1)
            with mock.patch("regression_bootstrap._cupy", return_value=stub):
                actual = bootstrap_errors(*arguments, backend="gpu", gpu_devices=(0, 1),
                                          batch_size=2, blas_threads=1)
            np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-14)
        for device in (0, 1):
            # Static U/V uploads occur once per device across both outer calls.
            self.assertEqual(stub.uploads.count((device, (13, 8))), 1)
            self.assertEqual(stub.uploads.count((device, (8, 8))), 1)
            self.assertEqual(stub.uploads.count((device, (8,))), 4)
            # Response means and coefficient centers must be refreshed each time.
            self.assertGreaterEqual(stub.uploads.count((device, (13, 2))), 2)
            self.assertEqual(stub.uploads.count((device, (8, 2))), 2)
            self.assertIn(device, stub.synchronizations)
        batches = [(device, shape) for device, shape in stub.downloads if len(shape) == 2]
        self.assertEqual(len(batches), 12)
        self.assertEqual(sum(shape[0] for _, shape in batches), 22)
        self.assertEqual({device for device, _ in batches}, {0, 1})
        self.assertTrue(all(shape[1] == 2 for _, shape in batches))


class RealGPURegressionBootstrapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            validate_backend("gpu", (0,))
        except RuntimeError as error:
            raise unittest.SkipTest(f"Real CUDA bootstrap not tested: {error}")

    def test_real_gpu_matches_cpu_and_legacy_for_seven_methods(self):
        for name in METHODS:
            design, method, parameter, settings, centers, sigma, thresholds, expected = legacy_case(name)
            with self.subTest(method=method):
                actual = bootstrap_errors(design, method, parameter, settings, centers, sigma,
                                          thresholds, 9, 20261006, 2, backend="gpu", batch_size=4,
                                          gpu_devices=(0,), blas_threads=1)
                np.testing.assert_allclose(actual, expected, rtol=5e-10, atol=5e-12)
                again = bootstrap_errors(design, method, parameter, settings, centers, sigma,
                                         thresholds, 9, 20261006, 2, backend="gpu", batch_size=2,
                                         gpu_devices=(0,), blas_threads=1)
                np.testing.assert_allclose(again, actual, rtol=5e-10, atol=5e-12)

    def test_two_real_gpus_preserve_batch_order(self):
        from regression_bootstrap import _cupy
        if _cupy().cuda.runtime.getDeviceCount() < 2:
            self.skipTest("Two-device CUDA test requires at least two GPUs")
        design, method, parameter, settings, centers, sigma, thresholds, expected = legacy_case("Bootstrapridge")
        actual = bootstrap_errors(design, method, parameter, settings, centers, sigma,
                                  thresholds, 9, 20261006, 2, backend="gpu", batch_size=2,
                                  gpu_devices=(0, 1), blas_threads=1)
        np.testing.assert_allclose(actual, expected, rtol=5e-10, atol=5e-12)


if __name__ == "__main__":
    unittest.main()
