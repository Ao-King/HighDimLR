"""Test iterative bootstrap debiasing without running notebook experiments."""

import ast
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
METHODS = (
    ("BootstrapARk", "ARk", "t4"),
    ("Bootstrapsecond", "second", "t5"),
    ("BootstrapNesterov", "Nesterov", "t1"),
)


def environment(notebook, method, step_name):
    rng = np.random.default_rng(23)
    X = rng.normal(size=(7, 4)) / 2
    env = dict(np=np, X=X, n=7, p=4, I=np.eye(4), kappa=1.1,
               omega=3, delta=-1.0, k_max=7, max_iter=7)
    env[step_name] = 0.12
    env["k_" + method] = 6
    env["beta_" + method] = np.array([0.3, -0.2, 0.1, 0.4])
    env["beta_" + method + "_debias"] = np.array([-0.2, 0.5, 0.3, -0.1])
    for suffix in ("", "_debias"):
        env["sigma_" + method + suffix] = 0.2
        env["bn_" + method + suffix] = 0.0
    names = {method, method + "_debias", "compute_norminf_" + method, "Dk"}
    with (ROOT / (notebook + ".ipynb")).open(encoding="utf-8") as handle:
        cells = json.load(handle)["cells"]
    definitions = []
    for cell in cells:
        if cell["cell_type"] != "code":
            continue
        for node in ast.parse("".join(cell["source"])).body:
            if isinstance(node, ast.FunctionDef) and node.name in names:
                node.decorator_list = []
                definitions.append(node)
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(ast.fix_missing_locations(module), notebook, "exec"), env)
    if method == "ARk":
        env["D"] = np.array([env["Dk"](k) for k in range(9)])
    return env


def linear_operator(env, method, step_name, steps):
    """Build L from first-fit basis responses, independently of the replay."""
    previous_delta = env["delta"]
    env["delta"] = -1.0
    try:
        return np.column_stack([
            first_fit(env, method, step_name, response, steps)[0]
            for response in np.eye(env["n"])
        ])
    finally:
        env["delta"] = previous_delta


def first_fit(env, method, step_name, response, steps):
    args = [env["X"], response, env[step_name]]
    if method == "Nesterov":
        return env[method](*args, steps, env["delta"], env["omega"])
    if method == "second":
        return env[method](*args, steps, env["delta"])
    return env[method](*args, steps)


def replay_fit(env, method, step_name, response, steps):
    args = [env["X"], response, env[step_name], steps]
    if method == "Nesterov":
        args.append(env["omega"])
    return env[method + "_debias"](*args)


class BootstrapDebiasTests(unittest.TestCase):
    def test_bootstrap_matches_fixed_operator_and_thresholding(self):
        for notebook, method, step_name in METHODS:
            for thresholded in (False, True):
                with self.subTest(method=method, thresholded=thresholded):
                    env = environment(notebook, method, step_name)
                    X = env["X"]
                    center = env["beta_" + method + "_debias"]
                    noise = np.linspace(-0.2, 0.3, env["n"])
                    noise_debias = np.linspace(0.4, -0.1, env["n"])
                    response = X @ center + noise_debias
                    steps = 6
                    L = linear_operator(env, method, step_name, steps)
                    first_expected = L @ response
                    corrected = (2 * L - L @ X @ L) @ response
                    self.assertGreater(np.linalg.norm(corrected - first_expected), 1e-5)
                    threshold = float(np.median(np.abs(corrected))) if thresholded else 0
                    env["bn_" + method + "_debias"] = threshold
                    expected = corrected.copy()
                    expected[np.abs(expected) <= threshold] = 0
                    if thresholded:
                        self.assertEqual(np.count_nonzero(expected), 2)
                    first_calls, replay_calls = [], []
                    fit, replay = env[method], env[method + "_debias"]

                    def record_fit(*args):
                        first_calls.append(args)
                        raise AssertionError("Bootstrap must not reselect its parameter")

                    def record_replay(X_arg, Y_arg, step, count, *args):
                        result = replay(X_arg, Y_arg, step, count, *args)
                        replay_calls.append((Y_arg.copy(), result.copy(), count))
                        return result

                    env[method] = record_fit
                    env[method + "_debias"] = record_replay
                    snapshots = {name: value.copy() for name, value in env.items()
                                 if isinstance(value, np.ndarray)}
                    with patch.object(np.random, "default_rng") as generator:
                        generator.return_value.normal.side_effect = [noise, noise_debias]
                        actual = env["compute_norminf_" + method](5)
                    self.assertEqual(len(first_calls), 0)
                    self.assertEqual(len(replay_calls), 3)
                    np.testing.assert_allclose(replay_calls[1][0], response)
                    np.testing.assert_allclose(replay_calls[1][1], first_expected, atol=1e-12)
                    np.testing.assert_allclose(replay_calls[2][0], X @ first_expected, atol=1e-12)
                    self.assertTrue(all(call[2] == steps for call in replay_calls))
                    np.testing.assert_allclose(replay_calls[2][1], L @ X @ first_expected,
                                               atol=1e-12)
                    ordinary = L @ (X @ env["beta_" + method] + noise)
                    ordinary[np.abs(ordinary) <= env["bn_" + method]] = 0
                    self.assertAlmostEqual(actual[0], np.linalg.norm(
                        ordinary - env["beta_" + method], np.inf))
                    self.assertAlmostEqual(actual[1], np.linalg.norm(expected - center, np.inf))
                    if not thresholded:
                        old_error = np.linalg.norm(first_expected - center, np.inf)
                        self.assertGreater(abs(actual[1] - old_error), 1e-5)
                    for name, original in snapshots.items():
                        np.testing.assert_array_equal(env[name], original)

    def test_replay_starts_from_zero_and_runs_exact_selected_count(self):
        for notebook, method, step_name in METHODS:
            for steps in (1, 3, 6):
                with self.subTest(method=method, steps=steps):
                    env = environment(notebook, method, step_name)
                    L = linear_operator(env, method, step_name, steps)
                    response = np.linspace(-0.3, 0.5, env["n"])
                    # Even an immediately satisfied discrepancy rule cannot stop replay.
                    env["delta"] = np.inf
                    expected = L @ response
                    np.testing.assert_allclose(
                        replay_fit(env, method, step_name, response, steps), expected,
                        atol=1e-12)
                    np.testing.assert_array_equal(
                        replay_fit(env, method, step_name, np.zeros(env["n"]), steps),
                        np.zeros(env["p"]))
                    np.testing.assert_allclose(
                        replay_fit(env, method, step_name, response, steps), expected,
                        atol=1e-12)

    def test_bootstrap_preserves_originally_selected_single_step(self):
        for notebook, method, step_name in METHODS:
            with self.subTest(method=method):
                env = environment(notebook, method, step_name)
                L = linear_operator(env, method, step_name, 1)
                env["k_" + method] = 1
                env["delta"] = np.inf
                center = env["beta_" + method + "_debias"]
                noise = np.linspace(-0.2, 0.3, env["n"])
                expected = (2 * L - L @ env["X"] @ L) @ (env["X"] @ center + noise)
                with patch.object(np.random, "default_rng") as generator:
                    generator.return_value.normal.side_effect = [noise, noise]
                    actual = env["compute_norminf_" + method](6)
                self.assertAlmostEqual(actual[1], np.linalg.norm(expected - center, np.inf))

    def test_ark_replay_is_not_truncated_by_global_max_iter(self):
        env = environment(*METHODS[0])
        L = linear_operator(env, "ARk", "t4", 6)
        response = np.linspace(-0.3, 0.5, env["n"])
        env["max_iter"] = 2
        result = env["ARk_debias"](env["X"], response, env["t4"], 6)
        np.testing.assert_allclose(result, L @ response, atol=1e-12)
        short = env["ARk_debias"](env["X"], response, env["t4"], 2)
        self.assertGreater(np.linalg.norm(result - short), 1e-5)


if __name__ == "__main__":
    unittest.main()
