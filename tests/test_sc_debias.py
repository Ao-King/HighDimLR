"""Exercise SC notebook statements without running the full experiments."""

import ast
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from regression_spectral import SpectralDesign

NOTEBOOKS = ("Nonsparse", "NonsparseDP", "Sparse", "normal")


def notebook_nodes(name):
    with (ROOT / (name + ".ipynb")).open(encoding="utf-8") as handle:
        notebook = json.load(handle)
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            yield from ast.walk(ast.parse("".join(cell["source"])))


def assignment(name, target):
    matches = [node for node in notebook_nodes(name)
               if isinstance(node, ast.Assign)
               and any(isinstance(item, ast.Name) and item.id == target
                       for item in node.targets)]
    if len(matches) != 1:
        raise AssertionError((name, target, len(matches)))
    return matches[0]


def execute(nodes, env):
    module = ast.Module(body=nodes, type_ignores=[])
    exec(compile(ast.fix_missing_locations(module), "<notebook SC>", "exec"), env)


def fixture(n, p, k):
    rng = np.random.default_rng(42)
    X = rng.normal(size=(n, p))
    U, S, Vt = np.linalg.svd(X, full_matrices=True)
    Sigma = np.zeros((n, p))
    np.fill_diagonal(Sigma, S)
    inverse = np.zeros(p)
    inverse[:k] = 1 / S[:k] ** 2
    S_inv = np.diag(inverse)
    Y = rng.normal(size=n)
    beta_cut = Vt.T @ S_inv @ Vt @ X.T @ Y
    return dict(np=np, X=X, Y=Y, U=U, S=S, Vt=Vt, Sigma=Sigma,
                S_inv=S_inv, I=np.eye(p), n=n, p=p, beta_cut=beta_cut)


class SCDebiasTests(unittest.TestCase):
    def test_accelerated_bootstrap_selection_keeps_exact_sc_copy(self):
        statement = next(node for node in notebook_nodes('Bootstrapcut')
                         if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Tuple)
                         and any(isinstance(item, ast.Name) and item.id == 'beta_cut_debias'
                                 for item in node.targets[0].elts))
        for n, p in ((5, 5), (7, 4), (4, 7)):
            for k in (1, min(n, p)):
                with self.subTest(shape=(n, p), rank=k):
                    env = fixture(n, p, k)
                    expected = env['beta_cut'].copy()
                    env.update(spectral_design=SpectralDesign(env['X'], svd=(env['U'], env['S'], env['Vt'])),
                               selection_settings={'rank_max': min(n, p)}, beta=expected.copy())
                    execute([statement], env)
                    self.assertEqual(env['selected_parameter'], k)
                    np.testing.assert_allclose(env['beta_cut'], expected, atol=1e-12)
                    np.testing.assert_array_equal(env['beta_cut'], env['beta_cut_debias'])
                    self.assertFalse(np.shares_memory(env['beta_cut'], env['beta_cut_debias']))

    def test_all_main_notebooks_match_spectral_projector(self):
        for name in NOTEBOOKS:
            statement = assignment(name, "beta_cut_debias")
            for n, p in ((5, 5), (7, 4), (4, 7)):
                for k in (0, 1, min(n, p)):
                    with self.subTest(notebook=name, shape=(n, p), k=k):
                        env = fixture(n, p, k)
                        self.assertGreater(np.linalg.norm(
                            env["Vt"] - np.diag(np.diag(env["Vt"]))), 0.1)
                        original = env["beta_cut"].copy()
                        projector = env["Vt"][:k].T @ env["Vt"][:k]
                        expected = (2 * env["I"] - projector) @ original
                        execute([statement], env)
                        result = env["beta_cut_debias"]
                        np.testing.assert_allclose(result, expected, atol=1e-12)
                        np.testing.assert_array_equal(result, original)
                        self.assertFalse(np.shares_memory(result, env["beta_cut"]))
                        result[:] = 99
                        np.testing.assert_array_equal(env["beta_cut"], original)

    def test_bootstrap_debias_copies_its_own_sample(self):
        env = fixture(7, 4, 2)
        env.update(beta_cutdestar=np.array([0.1, -0.4, 0.8, 2.0]),
                   beta_cutstar=np.array([3.0, 4.0, 5.0, 6.0]))
        original = env["beta_cutdestar"].copy()
        execute([assignment("Bootstrapcut", "beta_cutstar_debias")], env)
        result = env["beta_cutstar_debias"]
        np.testing.assert_array_equal(result, original)
        self.assertFalse(np.shares_memory(result, env["beta_cutdestar"]))
        result[np.abs(result) < 0.5] = 0
        np.testing.assert_array_equal(env["beta_cutdestar"], original)

    def test_bootstrap_function_same_and_different_noise(self):
        function = next(node for node in notebook_nodes("Bootstrapcut")
                        if isinstance(node, ast.FunctionDef)
                        and node.name == "compute_norminf_cut")
        for different in (False, True):
            with self.subTest(different_noise=different):
                env = fixture(7, 4, 2)
                env.update(beta_cut_debias=env["beta_cut"].copy(),
                           sigma_cut=0.2, sigma_cut_debias=0.2,
                           bn_cut=0.3, bn_cut_debias=0.3)
                snapshots = {key: value.copy() for key, value in env.items()
                             if isinstance(value, np.ndarray)}
                noise = np.linspace(-0.2, 0.3, env["n"])
                noise_debias = noise + (1.5 if different else 0)
                estimator = env["Vt"].T @ env["S_inv"] @ env["Vt"] @ env["X"].T
                expected = []
                for suffix, error in (("", noise), ("_debias", noise_debias)):
                    center = env["beta_cut" + suffix]
                    candidate = estimator @ (env["X"] @ center + error)
                    candidate[np.abs(candidate) <= env["bn_cut" + suffix]] = 0
                    expected.append(np.linalg.norm(candidate - center, ord=np.inf))
                execute([function], env)
                with patch.object(np.random, "default_rng") as generator:
                    generator.return_value.normal.side_effect = [noise, noise_debias]
                    actual = env["compute_norminf_cut"](10)
                np.testing.assert_allclose(actual, expected, atol=1e-12)
                if different:
                    self.assertNotAlmostEqual(expected[0], expected[1])
                for key, value in snapshots.items():
                    np.testing.assert_array_equal(env[key], value)

    def test_normal_threshold_grid_preserves_sc_and_records_candidates(self):
        env = dict(np=np, beta_cut=np.array([0.1, -0.4, 0.8, 2.0]),
                   beta_cut_debias=np.array([0.1, -0.4, 0.8, 2.0]),
                   beta=np.array([0.0, -0.5, 1.0, 1.5]), results=[])
        execute([assignment("normal", name) for name in
                 ("beta_cutcopy", "beta_cut_debiascopy")], env)
        snapshots = {key: value.copy() for key, value in env.items()
                     if isinstance(value, np.ndarray)}
        self.assertFalse(np.shares_memory(env["beta_cutcopy"], env["beta_cut"]))
        self.assertFalse(np.shares_memory(env["beta_cut_debiascopy"],
                                        env["beta_cut_debias"]))
        statements = [assignment("normal", name) for name in
                      ("beta_cut_bn", "beta_cut_debias_bn", "norm_cutbn",
                       "norm_cut_debiasbn")]
        append = next(node for node in notebook_nodes("normal")
                      if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                      and isinstance(node.value.func, ast.Attribute)
                      and node.value.func.attr == "append"
                      and isinstance(node.value.args[0], ast.List)
                      and len(node.value.args[0].elts) > 20)
        for element in append.value.args[0].elts:
            if isinstance(element, ast.Name):
                env.setdefault(element.id, float("nan"))
        for threshold in (0.5, 0.0, 1.0, 0.2):
            env["bn"] = threshold
            execute(statements + [append], env)
            expected = snapshots["beta_cut"].copy()
            expected[np.abs(expected) < threshold] = 0
            error = np.linalg.norm(expected - env["beta"])
            np.testing.assert_allclose(env["results"][-1][2:4], [error, error])
            for key, value in snapshots.items():
                np.testing.assert_array_equal(env[key], value)


if __name__ == "__main__":
    unittest.main()
