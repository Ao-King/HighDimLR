"""Run notebook functions and statements without their large simulations."""

import ast
import contextlib
import importlib.util
import io
import json
import math
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BOOTSTRAPS = {
    "BootstrapARk": ("ARk", "ARk", "ARk_debias", "t4", "k_ARk", "ARk"),
    "BootstrapNesterov": ("Nesterov", "Nesterov", "Nesterov_debias", "t1", "k_Nesterov", "Nesterov"),
    "Bootstrapsecond": ("second", "second", "second_debias", "t5", "k_second", "second"),
    "Bootstraplandweber": ("Landweber", "Landweber", "landweber_debias", "t1", "k_Landweber", "landweber"),
    "Bootstrapshowalter": ("showalter", "showalter", "showalter_debias", "t2", "k_showalter", "showalter"),
    "Bootstrapcut": ("cut", "cut", None, None, None, "compute_norm_cut"),
    "Bootstrapridge": ("Ridge", "Ridge", None, None, None, "compute_ridge"),
}
NESTEROV = ("Nonsparse", "NonsparseDP", "Sparse", "normal", "BootstrapNesterov")


def trees(name):
    with (ROOT / (name + ".ipynb")).open(encoding="utf-8") as stream:
        notebook = json.load(stream)
    return [ast.parse("".join(cell["source"])) for cell in notebook["cells"]
            if cell["cell_type"] == "code"]


def execute(nodes, env):
    module = ast.Module(body=nodes, type_ignores=[])
    exec(compile(ast.fix_missing_locations(module), "<notebook regression>", "exec"), env)


def functions(name, env):
    definitions = []
    for tree in trees(name):
        for node in tree.body:
            if isinstance(node, ast.FunctionDef):
                node.decorator_list = []
                definitions.append(node)
    execute(definitions, env)
    return env


def nesterov_reference(X, Y, step, count, omega):
    previous = np.zeros(X.shape[1])
    if count == 0:
        return previous
    current = step * X.T @ Y
    for k in range(1, count):
        extrapolated = current + (k - 1) / (k + omega) * (current - previous)
        previous, current = current, extrapolated + step * X.T @ (Y - X @ extrapolated)
    return current


def bootstrap_environment(name, count=4):
    label, sigma_label, replay, step_name, count_name, selector = BOOTSTRAPS[name]
    rng = np.random.default_rng(921)
    X = rng.normal(size=(7, 4)) / 3
    U, singular, Vt = np.linalg.svd(X, full_matrices=True)
    env = dict(np=np, n=7, p=4, X=X, U=U, S=singular, Vt=Vt, I=np.eye(4),
               omega=3, kappa=1.1, delta=-1.0, max_iter=20, k_max=20,
               best_c=2.3, best_k=2, t4=0.12)
    env["beta_" + label] = np.array([0.3, -0.2, 0.1, 0.4])
    env["beta_" + label + "_debias"] = np.array([-0.2, 0.5, 0.3, -0.1])
    for suffix in ("", "_debias"):
        env["sigma_" + sigma_label + suffix] = 0.2
        env["bn_" + label + suffix] = 0.0
    inverse = np.zeros(4)
    inverse[:2] = 1 / singular[:2] ** 2
    env["S_inv"] = np.diag(inverse)
    if replay:
        env[step_name] = 0.12
        env[count_name] = count
    functions(name, env)
    if name == "BootstrapARk":
        env["D"] = np.array([env["Dk"](k) for k in range(30)])
    if name == "Bootstrapcut":
        operator = Vt.T @ env["S_inv"] @ Vt @ X.T
    elif name == "Bootstrapridge":
        operator = np.linalg.solve(X.T @ X + env["best_c"] / env["n"] * env["I"], X.T)
    elif name in ("Bootstraplandweber", "Bootstrapshowalter"):
        H = X.T @ X
        step = env[step_name]
        residual = np.eye(4) - step * H
        if name == "Bootstrapshowalter":
            residual += step**2 * H @ H / 2 - step**3 * H @ H @ H / 6
            residual += step**4 * H @ H @ H @ H / 24
        operator = np.linalg.solve(H, np.eye(4) - np.linalg.matrix_power(residual, count)) @ X.T
    else:
        extra = [env["omega"]] if name == "BootstrapNesterov" else []
        operator = np.column_stack([
            env[replay](X, response, env[step_name], count, *extra)
            for response in np.eye(env["n"])
        ])
    return env, operator


class NotebookRegressionTests(unittest.TestCase):
    def test_nesterov_matches_manuscript_for_all_notebooks(self):
        X = np.array([[1.0, 0.2], [0.0, 0.7], [0.3, 0.4]])
        Y = np.array([1.0, -0.5, 0.3])
        for name in NESTEROV:
            env = functions(name, {"np": np})
            for count in (0, 1, 2, 3, 6):
                for omega in (1.0, 5.0):
                    with self.subTest(notebook=name, count=count, omega=omega):
                        expected = nesterov_reference(X, Y, 0.1, count, omega)
                        choice = np.linalg.lstsq(X, Y, rcond=None)[0] if name == "Nonsparse" else -1.0
                        actual, selected = env["Nesterov"](X, Y, 0.1, count, choice, omega)
                        self.assertEqual(selected, count)
                        np.testing.assert_allclose(actual, expected, atol=1e-13)
                        replay = env["Nesterov_debias"](X, Y, 0.1, selected, omega)
                        np.testing.assert_allclose(replay, expected, atol=1e-13)
        scalar = functions("BootstrapNesterov", {"np": np})
        actual, _ = scalar["Nesterov"](np.eye(1), np.ones(1), 0.1, 3, -1.0, 3.0)
        np.testing.assert_allclose(actual, [0.2872])
        self.assertNotAlmostEqual(actual[0], 1 - 0.9**3)

    def test_discrepancy_count_is_replayed_without_extra_step(self):
        for name in NESTEROV:
            if name == "Nonsparse":
                continue
            with self.subTest(notebook=name):
                env = functions(name, {"np": np})
                X, Y = np.eye(1), np.ones(1)
                first, count = env["Nesterov"](X, Y, 0.1, 8, np.inf, 3)
                self.assertEqual(count, 1)
                second = env["Nesterov_debias"](X, X @ first, 0.1, count, 3)
                np.testing.assert_allclose(2 * first - second, [0.19])

    def test_original_iterative_bootstrap_fits_honor_zero_and_positive_caps(self):
        for name, config in BOOTSTRAPS.items():
            _, _, replay, step_name, _, selector = config
            if replay is None:
                continue
            for cap in (0, 1, 4):
                with self.subTest(notebook=name, cap=cap):
                    env, _ = bootstrap_environment(name)
                    X, Y, step = env["X"], np.linspace(-0.3, 0.5, 7), env[step_name]
                    args = [X, Y, step, cap]
                    if name == "BootstrapNesterov":
                        args.extend([-1.0, env["omega"]])
                    elif name in ("Bootstrapsecond", "Bootstrapshowalter"):
                        args.append(-1.0)
                    elif name == "Bootstraplandweber":
                        args.append(np.zeros(4))
                    fitted, count = env[selector](*args)
                    self.assertEqual(count, cap)
                    replay_args = [X, Y, step, count]
                    if name == "BootstrapNesterov":
                        replay_args.append(env["omega"])
                    np.testing.assert_allclose(fitted, env[replay](*replay_args), atol=1e-13)

    def test_other_notebooks_replay_every_selected_update(self):
        rng = np.random.default_rng(18)
        X, Y = rng.normal(size=(5, 3)) / 3, rng.normal(size=5)
        wanted = {"landweber_debias", "showalter_debias", "hbf_debias", "ARk_debias", "second_debias"}
        for name in ("normal", "Sparse", "sparse15"):
            for cell_index, tree in enumerate(trees(name)):
                definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
                replays = [node.name for node in definitions if node.name in wanted]
                for replay in replays:
                    env = {"np": np, "delta": -1.0, "kappa": 1.25, "max_iter": 1}
                    execute(definitions, env)
                    for steps in (0, 1, 3, 5):
                        with self.subTest(notebook=name, cell=cell_index, replay=replay, steps=steps):
                            with contextlib.redirect_stdout(io.StringIO()):
                                expected, count = env[replay.removesuffix("_debias")](X, Y, 0.1, steps)
                            self.assertEqual(count, steps)
                            actual = env[replay](X, Y, 0.1, count)
                            np.testing.assert_allclose(actual, expected, atol=1e-12)
                            np.testing.assert_array_equal(env[replay](X, Y, 0.1, 0), np.zeros(3))

    def test_landweber_dp_returns_the_iterate_satisfying_discrepancy(self):
        env = functions("NonsparseDP", {"np": np, "delta": 0.91})
        X, Y = np.eye(1), np.ones(1)
        result, count = env["landweber"](X, Y, 0.1, 5)
        self.assertEqual(count, 1)
        np.testing.assert_allclose(result, [0.1])
        self.assertLessEqual(np.linalg.norm(Y - X @ result), env["delta"])
        env["delta"] = -1.0
        for cap in (0, 1, 4):
            result, count = env["landweber"](X, Y, 0.1, cap)
            self.assertEqual(count, cap)
            np.testing.assert_allclose(result, [1 - 0.9**cap])

    def test_fractional_initial_steps_match_independent_scalar_formula(self):
        for name in ("Nonsparse", "NonsparseDP", "normal", "Sparse"):
            for theta in (0.5, 1.25):
                env = {"np": np, "gamma": math.gamma, "t": 0.025, "theta": theta,
                       "max_iter": 4, "delta": -1.0, "beta": np.ones(1)}
                fractional_tree = next(tree for tree in trees(name) if any(
                    isinstance(node, ast.FunctionDef) and node.name == "Frac" for node in tree.body))
                prefix = []
                for node in fractional_tree.body:
                    if isinstance(node, ast.FunctionDef):
                        break
                    if isinstance(node, (ast.Import, ast.ImportFrom)):
                        continue
                    if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name)
                            and target.id in {"t", "theta", "max_iter"} for target in node.targets):
                        continue
                    prefix.append(node)
                execute(prefix, env)
                execute([node for node in fractional_tree.body if isinstance(node, ast.FunctionDef)], env)
                ca = env["t"]**theta / (theta * (theta + 1))
                cb = env["t"]**theta / theta
                g = 1 / math.gamma(theta)
                first = g * (ca * (1 - g * cb) + ca * theta)
                predictor_second = g * (cb * (2**theta - 1) + cb * (1 - first))
                a0 = ca * (1 - (1 - theta) * 2**theta)
                a1 = ca * (2**(theta + 1) - 2)
                second = g * (ca * (1 - predictor_second) + a0 + a1 * (1 - first))
                for steps, expected in ((0, 0.0), (1, first), (2, second)):
                    with self.subTest(notebook=name, theta=theta, steps=steps):
                        fitted, count = env["Frac"](np.eye(1), np.ones(1), steps)
                        self.assertEqual(count, steps)
                        np.testing.assert_allclose(fitted, [expected], atol=1e-13)
                        replay = env["Frac_debias"](np.eye(1), np.ones(1), steps)
                        np.testing.assert_allclose(replay, [expected], atol=1e-13)
                        second_fit = env["Frac_debias"](np.eye(1), fitted, count)
                        np.testing.assert_allclose(2 * fitted - second_fit,
                                                   [2 * expected - expected**2], atol=1e-13)
                capped, count = env["Frac"](np.eye(1), np.ones(1), env["max_iter"])
                self.assertEqual(count, env["max_iter"])
                self.assertTrue(np.all(np.isfinite(capped)))
                for node in ast.walk(fractional_tree):
                    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "Frac_debias":
                        self.assertIsInstance(node.args[-1], ast.Name)
                        self.assertEqual(node.args[-1].id, "k_Frac")

    def test_all_bootstraps_use_fixed_operator_and_refit_support(self):
        for name, config in BOOTSTRAPS.items():
            label, _, _, _, _, selector = config
            for noise_scale in (0.01, 4.0):
                with self.subTest(notebook=name, noise_scale=noise_scale):
                    env, L = bootstrap_environment(name)
                    G = 2 * L - L @ env["X"] @ L
                    noises = [noise_scale * np.linspace(-0.2, 0.3, 7),
                              noise_scale * np.linspace(0.4, -0.1, 7)]
                    expected = []
                    for suffix, operator, noise in zip(("", "_debias"), (L, G), noises):
                        center = env["beta_" + label + suffix]
                        fitted = operator @ (env["X"] @ center + noise)
                        threshold = float(np.median(np.abs(fitted)))
                        env["bn_" + label + suffix] = threshold
                        fitted[np.abs(fitted) <= threshold] = 0
                        expected.append(np.linalg.norm(fitted - center, np.inf))
                    snapshots = {key: value.copy() for key, value in env.items()
                                 if isinstance(value, np.ndarray)}

                    def forbidden(*args, **kwargs):
                        raise AssertionError("The original-sample parameter must stay fixed")

                    env[selector] = forbidden
                    # These would stop/restrict a fresh selection call immediately.
                    env.update(delta=np.inf, max_iter=1, k_max=1)
                    with patch.object(np.random, "default_rng") as generator:
                        generator.return_value.normal.side_effect = noises
                        actual = env["compute_norminf_" + label](19)
                    np.testing.assert_allclose(actual, expected, atol=1e-11)
                    for key, before in snapshots.items():
                        np.testing.assert_array_equal(env[key], before)

    def test_sigma_interface_returns_sd_with_sample_size_denominator(self):
        X = np.zeros((7, 4))
        Y = np.arange(1.0, 8.0)
        for name in BOOTSTRAPS:
            with self.subTest(notebook=name):
                env = functions(name, {"np": np, "p": 99})
                actual = env["compute_sigma"](Y, X, np.zeros(4))
                self.assertAlmostEqual(actual, np.sqrt(np.mean(Y**2)))

    def test_nonbootstrap_sd_statements_use_the_observation_count(self):
        found = 0
        for name in ("Nonsparse", "NonsparseDP", "Sparse", "normal", "sparse15"):
            for tree in trees(name):
                for node in ast.walk(tree):
                    if not isinstance(node, ast.Assign) or not isinstance(node.value, ast.BinOp):
                        continue
                    if not isinstance(node.value.op, ast.Div):
                        continue
                    if not any(isinstance(target, ast.Name) and target.id.startswith("sigma_")
                               for target in node.targets):
                        continue
                    if not isinstance(node.value.left, ast.Name) or not node.value.left.id.startswith("sigma_"):
                        continue
                    with self.subTest(notebook=name, variable=node.value.left.id):
                        env = {"np": np, "X": np.zeros((9, 4)), "s2": 4,
                               node.value.left.id: 6.0}
                        execute([node], env)
                        self.assertEqual(env[node.targets[0].id], 2.0)
                        found += 1
        self.assertGreater(found, 30)

    def test_thresholded_sd_is_used_for_bootstrap_and_reporting(self):
        for name, config in BOOTSTRAPS.items():
            label, sigma_label = config[:2]
            env, _ = bootstrap_environment(name)
            env.update(Y=np.linspace(-0.4, 0.8, 7), std_dev=2.0,
                       variable_names=[label, label + "_debias"])
            for suffix in ("", "_debias"):
                env["sigma_" + sigma_label + suffix] = 999.0
            loops = [node for tree in trees(name) for node in ast.walk(tree)
                     if isinstance(node, ast.For) and isinstance(node.target, ast.Name)
                     and node.target.id == "name" and any(
                         isinstance(item, ast.Call) and isinstance(item.func, ast.Name)
                         and item.func.id == "compute_sigma" for item in ast.walk(node))]
            self.assertEqual(len(loops), 1, name)
            execute(loops, env)
            expected = [np.linalg.norm(env["Y"] - env["X"] @ env["beta_" + label + suffix]) / np.sqrt(7)
                        for suffix in ("", "_debias")]
            with patch.object(np.random, "default_rng") as generator:
                generator.return_value.normal.side_effect = [np.zeros(7), np.zeros(7)]
                env["compute_norminf_" + label]((1, 2, 3))
                calls = generator.return_value.normal.call_args_list
                np.testing.assert_allclose([call.args[1] for call in calls], expected)
            statements = [node for tree in trees(name) for node in ast.walk(tree)
                          if isinstance(node, ast.AugAssign) and isinstance(node.target, ast.Name)
                          and node.target.id.startswith("Sigma_")]
            for node in statements:
                env[node.target.id] = 0.0
            execute(statements, env)
            np.testing.assert_allclose([env[node.target.id] for node in statements],
                                       np.abs(np.array(expected)**2 - 4))

    def test_all_variance_summaries_accumulate_absolute_variance_error(self):
        for name in BOOTSTRAPS:
            statements = [node for tree in trees(name) for node in ast.walk(tree)
                          if isinstance(node, ast.AugAssign) and isinstance(node.target, ast.Name)
                          and node.target.id.startswith("Sigma_")]
            self.assertEqual(len(statements), 2, name)
            env = {"np": np, "std_dev": 1.7}
            for node in statements:
                env[node.target.id] = 0.0
                for item in ast.walk(node.value):
                    if isinstance(item, ast.Name) and item.id.startswith("sigma_"):
                        env[item.id] = 2.1
            execute(statements, env)
            for node in statements:
                self.assertAlmostEqual(env[node.target.id], abs(2.1**2 - 1.7**2), msg=name)

    def test_threshold_selection_uses_values_not_zero_based_row_numbers(self):
        class Column:
            def to_numpy(self):
                return np.array([0.001, 0.002, 0.003])

        class Frame:
            def __len__(self):
                return 3

            def __getitem__(self, key):
                if key != "bn":
                    raise KeyError(key)
                return Column()

        for name in BOOTSTRAPS:
            selected = [node for tree in trees(name) for node in ast.walk(tree)
                        if isinstance(node, ast.Assign) and any(
                            isinstance(target, ast.Name) and target.id == "selected_bn"
                            for target in node.targets)]
            self.assertEqual(len(selected), 1, name)
            for index, expected in ((0, 0.001), (1.5, 0.0025)):
                env = {"np": np, "avg_index": index, "df": Frame()}
                execute(selected, env)
                self.assertAlmostEqual(env["selected_bn"], expected, msg=name)

    def test_bootstrap_random_streams_vary_by_outer_replication(self):
        for name, config in BOOTSTRAPS.items():
            label = config[0]
            with self.subTest(notebook=name):
                env, _ = bootstrap_environment(name)
                sample = env["compute_norminf_" + label]
                np.random.seed(782)
                before = np.random.get_state()
                first = sample((20261006, 1, 8))
                repeated = sample((20261006, 1, 8))
                other_outer = sample((20261006, 2, 8))
                other_inner = sample((20261006, 1, 9))
                after = np.random.get_state()
                np.testing.assert_array_equal(first, repeated)
                self.assertGreater(np.linalg.norm(np.array(first) - other_outer), 1e-12)
                self.assertGreater(np.linalg.norm(np.array(first) - other_inner), 1e-12)
                self.assertEqual(before[0], after[0])
                np.testing.assert_array_equal(before[1], after[1])
                self.assertEqual(before[2:], after[2:])

    def test_coverage_includes_the_boundary_of_the_confidence_region(self):
        for name in BOOTSTRAPS:
            checks = [node for tree in trees(name) for node in ast.walk(tree)
                      if isinstance(node, ast.If) and isinstance(node.test, ast.Compare)
                      and isinstance(node.test.left, ast.Name)
                      and node.test.left.id.startswith("inf_")
                      and any(isinstance(item, ast.AugAssign) for item in node.body)]
            self.assertEqual(len(checks), 2, name)
            for check in checks:
                cutoff = check.test.comparators[0]
                counter = next(item.target.id for item in check.body
                               if isinstance(item, ast.AugAssign))
                env = {"np": np, "confidence": 0.95, counter: 0, check.test.left.id: 0.0}
                if isinstance(cutoff, ast.Name):
                    env[cutoff.id] = 0.0
                else:
                    env[cutoff.args[0].id] = np.zeros(10)
                execute([check], env)
                self.assertEqual(env[counter], 1, name)
                env[counter] = 0
                env[check.test.left.id] = np.nextafter(0.0, 1.0)
                execute([check], env)
                self.assertEqual(env[counter], 0, name)

    def test_jit_functions_have_no_data_dependent_globals(self):
        for name in BOOTSTRAPS:
            for tree in trees(name):
                for node in tree.body:
                    if not isinstance(node, ast.FunctionDef) or not node.decorator_list:
                        continue
                    local = {arg.arg for arg in node.args.args}
                    local |= {item.id for item in ast.walk(node)
                              if isinstance(item, ast.Name) and isinstance(item.ctx, ast.Store)}
                    loaded = {item.id for item in ast.walk(node)
                              if isinstance(item, ast.Name) and isinstance(item.ctx, ast.Load)}
                    self.assertFalse(loaded - local - {"np", "range", "min", "max", "jit"}, (name, node.name))

    @unittest.skipUnless(importlib.util.find_spec("numba"), "Numba is not installed in this interpreter")
    def test_compiled_iterative_functions_observe_runtime_parameters(self):
        from numba import jit

        X, Y = np.eye(1), np.ones(1)
        for name, method in (("BootstrapNesterov", "Nesterov"),
                             ("Bootstrapsecond", "second"),
                             ("Bootstrapshowalter", "showalter"),
                             ("Bootstraplandweber", "landweber")):
            with self.subTest(notebook=name):
                env = functions(name, {"np": np})
                compiled = jit(nopython=True)(env[method])
                extra = [3.0] if method == "Nesterov" else []
                if method == "landweber":
                    _, first = compiled(X, Y, 0.001, 520, np.zeros(1))
                    _, second = compiled(X, Y, 0.001, 520, np.ones(1))
                    self.assertEqual((first, second), (500, 520))
                else:
                    _, first = compiled(X, Y, 0.1, 5, 2.0, *extra)
                    _, second = compiled(X, Y, 0.1, 5, -1.0, *extra)
                    self.assertEqual((first, second), (1, 5))
                replay = env[method + "_debias"]
                compiled_replay = jit(nopython=True)(replay)
                for count in (0, 1, 5):
                    expected = replay(X, Y, 0.1, count, *extra)
                    actual = compiled_replay(X, Y, 0.1, count, *extra)
                    np.testing.assert_allclose(actual, expected, atol=1e-13)

    def test_all_notebooks_parse_and_direct_call_signatures_match(self):
        for path in ROOT.glob("*.ipynb"):
            definitions = {}
            calls = []
            for tree in trees(path.stem):
                definitions.update({node.name: node for node in tree.body
                                    if isinstance(node, ast.FunctionDef)})
                calls.extend(node for node in ast.walk(tree)
                             if isinstance(node, ast.Call) and isinstance(node.func, ast.Name))
            for call in calls:
                definition = definitions.get(call.func.id)
                if definition is None or definition.args.vararg or definition.args.kwarg:
                    continue
                required = len(definition.args.args) - len(definition.args.defaults)
                supplied = len(call.args) + len(call.keywords)
                self.assertTrue(required <= supplied <= len(definition.args.args),
                                (path.name, call.func.id, supplied))


if __name__ == "__main__":
    unittest.main()
