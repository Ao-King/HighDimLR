"""Cold-kernel, reduced-size execution of every experiment entry point."""

import ast
import contextlib
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from regression_common import discrepancy_candidate, lasso_coefficients, ridge_pair, sparse_indices


def cells(name):
    return json.loads((ROOT / (name + '.ipynb')).read_text(encoding='utf-8'))['cells']


class SmallExperiment(ast.NodeTransformer):
    def __init__(self, p, noise, output_root, resume_from):
        self.p, self.noise = p, noise
        self.output_root, self.resume_from = output_root, resume_from

    def visit_Assign(self, node):
        self.generic_visit(node)
        sizes = dict(n=24, p=self.p, max_iter=7, k_max=7, C=2, m=1, N=2, B=3,
                     cpu_workers=2, blas_threads=1, batch_size=2)
        if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            key = node.targets[0].id
            if key in sizes and isinstance(node.value, ast.Constant):
                node.value = ast.Constant(sizes[key])
            elif key == 'noise_model':
                node.value = ast.Constant(self.noise)
            elif key == 'output_root':
                node.value = ast.Constant(str(self.output_root))
            elif key == 'resume_from':
                node.value = ast.Constant(None if self.resume_from is None else str(self.resume_from))
        return node

    def visit_comprehension(self, node):
        self.generic_visit(node)
        if isinstance(node.target, ast.Name) and node.target.id == 'seed':
            node.iter = ast.Call(ast.Name('range', ast.Load()), [ast.Constant(3)], [])
        return node

    def visit_Call(self, node):
        self.generic_visit(node)
        if isinstance(node.func, ast.Name) and node.func.id == 'Parallel':
            for keyword in node.keywords:
                if keyword.arg == 'n_jobs':
                    keyword.value = ast.Constant(2)
        if isinstance(node.func, ast.Name) and node.func.id == 'range':
            if ast.dump(node) == ast.dump(ast.parse('range(1000, 0, -1)', mode='eval').body):
                node.args[0] = ast.Constant(24)
        return node


def run_small(name, noise='normal', output_root=None, resume_from=None):
    width = 36 if name == 'sparse15' else 24
    env = {'__name__': '__main__'}
    with contextlib.ExitStack() as stack:
        stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
        stack.enter_context(threadpool_limits(limits=1))
        stack.enter_context(patch.object(plt, 'savefig'))
        stack.enter_context(patch.object(plt, 'show'))
        if not name.startswith('Bootstrap'):
            stack.enter_context(patch.object(pd.DataFrame, 'to_csv'))
        if output_root is None:
            output_root = stack.enter_context(tempfile.TemporaryDirectory())
        for index, cell in enumerate(cells(name)):
            if cell['cell_type'] == 'code':
                tree = SmallExperiment(width, noise, output_root, resume_from).visit(
                    ast.parse(''.join(cell['source'])))
                exec(compile(ast.fix_missing_locations(tree), f'{name}:cell{index}', 'exec'), env)
    plt.close('all')
    return env


class SharedRegressionTests(unittest.TestCase):
    def test_ridge_matches_spectral_formula_including_wide_and_singular_designs(self):
        rng = np.random.default_rng(19)
        for shape in ((7, 4), (4, 7), (4, 4)):
            X = rng.normal(size=shape)
            X[:, -1] = X[:, 0]
            Y = rng.normal(size=shape[0])
            U, s, Vt = np.linalg.svd(X, full_matrices=False)
            for alpha in (0.0, 0.1, 10.0):
                original, debiased = ridge_pair(X, Y, alpha)
                if alpha == 0:
                    expected = np.linalg.pinv(X) @ Y
                    expected_debiased = expected
                else:
                    expected = Vt.T @ ((s / (s*s + alpha)) * (U.T @ Y))
                    expected_debiased = Vt.T @ ((s * (s*s + 2*alpha) /
                                                   (s*s + alpha)**2) * (U.T @ Y))
                np.testing.assert_allclose(original, expected, atol=1e-12)
                np.testing.assert_allclose(debiased, expected_debiased, atol=1e-12)

    def test_lasso_uses_no_intercept_and_consistent_scale(self):
        X = np.array([-2., -1., 1., 2.])[:, None]
        Y = 2 * X[:, 0]
        np.testing.assert_allclose(lasso_coefficients(X, Y, 2 / 4), [1.8])
        np.testing.assert_allclose(lasso_coefficients(X, Y, 0), [2.0])
        np.testing.assert_allclose(lasso_coefficients(np.ones((4, 1)), np.full(4, 2.), 0.5), [1.5])

    def test_sparse_support_is_distinct_and_preserves_collision_free_legacy_draws(self):
        for seed in (18, 24, 921):
            rng = np.random.RandomState(seed)
            result = sparse_indices(rng, 1000)
            self.assertEqual(len(np.unique(result)), 20)
            raw = np.random.RandomState(seed).randint(0, 1000, size=20)
            if len(np.unique(raw)) == 20:
                np.testing.assert_array_equal(result, raw)
        self.assertEqual(len(np.unique(sparse_indices(np.random.default_rng(1), 20))), 20)
        with self.assertRaises(ValueError):
            sparse_indices(np.random.default_rng(1), 19)

    def test_discrepancy_returns_actual_grid_value_and_rejects_missing_candidates(self):
        grid = np.array([21, 11, 1])
        self.assertEqual(discrepancy_candidate(grid, [3, 2, 1], 2), 11)
        with self.assertRaises(ValueError):
            discrepancy_candidate(grid, [3, 2, 1], 0.5)
        with self.assertRaises(FloatingPointError):
            discrepancy_candidate(grid, [3, np.nan, 1], 2)


class ExperimentEntryTests(unittest.TestCase):
    def test_full_size_nesterov_steps_and_actual_discrepancy_stopping(self):
        with threadpool_limits(limits=1):
            rng = np.random.RandomState(24)
            covariance = np.full((1000, 1000), 0.5)
            np.fill_diagonal(covariance, 2.0)
            X = rng.multivariate_normal(np.zeros(1000), covariance, 1000)
            singular = np.linalg.svd(X, compute_uv=False)
            sparse_beta = np.zeros(1000)
            sparse_beta[sparse_indices(rng, 1000)] = np.repeat([2., -2., 1., -1.], 5)
            dense_beta = np.random.RandomState(25).uniform(-2, 2, 1000)
            for name, beta, noise, multiplier in (
                ('NonsparseDP', dense_beta, np.random.RandomState(26).normal(0, 2, 1000), 0.6),
                ('BootstrapNesterov', sparse_beta,
                 np.random.default_rng((20261006, 0, 0, 0)).normal(0, 2, 1000), 1.01),
            ):
                env = {'np': np, 'S': singular}
                nodes = []
                step_name = 't1' if name.startswith('Bootstrap') else 't'
                for cell in cells(name):
                    if cell['cell_type'] != 'code':
                        continue
                    text = ''.join(cell['source'])
                    for node in ast.walk(ast.parse(text)):
                        if isinstance(node, ast.Assign) and len(node.targets) == 1:
                            if isinstance(node.targets[0], ast.Name):
                                key = node.targets[0].id
                                if key in ('nominal_nesterov_step', step_name) and (
                                        name.startswith('Bootstrap') or 'def Nesterov(' in text):
                                    nodes.append(node)
                        elif isinstance(node, ast.FunctionDef) and node.name in ('Nesterov', 'Nesterov_debias'):
                            node.decorator_list = []
                            nodes.append(node)
                exec(compile(ast.Module(body=nodes, type_ignores=[]), '<default Nesterov>', 'exec'), env)
                step = env[step_name]
                self.assertLessEqual(step * singular[0]**2, 1)
                Y = X @ beta + noise
                delta = multiplier * np.linalg.norm(noise)
                estimate, count = env['Nesterov'](X, Y, step, 10000, delta, 3)
                self.assertLess(count, 10000, name)
                self.assertLessEqual(np.linalg.norm(Y - X @ estimate), delta, name)
                refit = env['Nesterov_debias'](X, X @ estimate, step, count, 3)
                self.assertTrue(np.all(np.isfinite(2*estimate - refit)), name)

    def check_lasso(self, env):
        actual = env['beta_lasso']
        expected = lasso_coefficients(env['X'], env['Y'], env['best_c'] / env['n'])
        np.testing.assert_allclose(actual, expected, atol=1e-12)
        self.assertEqual(env['lasso_alpha'], env['best_c'] / env['n'])
        original, debiased = ridge_pair(env['X'], env['Y'], env['ridge_alpha'])
        np.testing.assert_allclose(env['beta_Ridge'], original, atol=1e-12)
        np.testing.assert_allclose(env['beta_Ridge_debias'], debiased, atol=1e-12)

    def test_complete_sparse_entries_both_noise_models_and_rectangular_case(self):
        for name in ('normal', 'Sparse', 'sparse15'):
            reference_design = None
            for noise in ('normal', 'laplace'):
                with self.subTest(notebook=name, noise=noise):
                    env = run_small(name, noise)
                    self.assertEqual(env['X'].shape, (24, 36 if name == 'sparse15' else 24))
                    self.assertEqual(np.count_nonzero(env['beta']), 20)
                    for value in (2., -2., 1., -1.):
                        self.assertEqual(np.count_nonzero(env['beta'] == value), 5)
                    self.check_lasso(env)
                    self.assertEqual(env['df'].shape, (1000, 21))
                    labels = ('cut', 'cut_debias', 'Lasso', 'Ridge', 'Ridge_debias',
                              'Landweber', 'Landweber_debias', 'showalter', 'showalter_debias',
                              'hbf', 'hbf_debias', 'ARk', 'ARk_debias', 'second', 'second_debias',
                              'Nesterov', 'Nesterov_debias', 'Frac', 'Frac_debias')
                    for row, label in enumerate(('LS',) + labels):
                        error, threshold = env['results_array'][row]
                        restored = env['beta_' + label + 'copy']
                        np.testing.assert_allclose(np.linalg.norm(restored - env['beta']), error,
                                                   atol=1e-12, err_msg=label)
                        self.assertGreaterEqual(threshold, 0.001)
                        self.assertLessEqual(threshold, 1)
                    for label in ('LS', 'Ridge', 'hbf', 'ARk', 'second', 'Nesterov', 'Frac'):
                        self.assertTrue(np.isfinite(env['norm_' + label]))
                    np.testing.assert_allclose(env['Y'] - env['X'] @ env['beta'], env['e'])
                    if reference_design is None:
                        reference_design = (env['X'].copy(), env['beta'].copy())
                    else:
                        np.testing.assert_array_equal(env['X'], reference_design[0])
                        np.testing.assert_array_equal(env['beta'], reference_design[1])

    def test_complete_nonsparse_entries(self):
        results = []
        for name in ('Nonsparse', 'NonsparseDP'):
            with self.subTest(notebook=name):
                env = run_small(name)
                self.check_lasso(env)
                if name == 'NonsparseDP':
                    for key in ('beta_Ridge', 'beta_lasso'):
                        self.assertLessEqual(np.linalg.norm(env['Y'] - env['X'] @ env[key]),
                                             env['delta'] + 1e-10)
                results.append(env)
        for key in ('X', 'Y', 'beta'):
            np.testing.assert_array_equal(results[0][key], results[1][key])

    def test_all_bootstrap_entries_with_real_parallel_workers(self):
        from regression_results import export_regression_table
        from regression_execution import threshold_pair

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            completed = []
            expected = {}
            for path in ROOT.glob('Bootstrap*.ipynb'):
                with self.subTest(notebook=path.stem):
                    env = run_small(path.stem, output_root=root / 'runs')
                    self.assertEqual(env['x'], 1)
                    self.assertEqual(len(env['results']), 3)
                    self.assertEqual(np.count_nonzero(env['beta']), 20)
                    self.assertTrue(np.all(np.isfinite(np.asarray(env['results']))))
                    ordinary, debiased = env['spectral_design'].pair(
                        env['Y'], env['run_config']['method'], env['selected_parameter'],
                        env['selection_settings'])
                    centers, thresholds, scan = threshold_pair(
                        np.column_stack([ordinary, debiased]), env['beta'],
                        env['run_config']['threshold_percentile'])
                    np.testing.assert_array_equal(scan, env['df'].to_numpy())
                    for branch, label in enumerate(env['variable_names']):
                        np.testing.assert_array_equal(centers[:, branch], env['beta_' + label])
                        self.assertEqual(thresholds[branch], env['bn_' + label])
                        self.assertEqual(env['sigma_' + label], env['compute_sigma'](
                            env['Y'], env['X'], centers[:, branch]))
                    self.assertNotIn('D_array', env)
                    completed.append(env['run'].directory)
                    expected[env['run_config']['method']] = env['run'].summary
            with self.assertRaises(ValueError):
                export_regression_table(completed, root / 'formal')
            exported = export_regression_table(completed, root / 'tables', allow_nonproduction=True)
            self.assertEqual(exported['csv'].name, 'regression_table.csv')
            self.assertEqual(exported['tex'].name, 'regression_table.tex')
            self.assertNotIn('NON_MANUSCRIPT', exported['tex'].read_text(encoding='utf-8'))
            manifest = json.loads(exported['manifest'].read_text(encoding='utf-8'))
            for run in manifest['runs'].values():
                self.assertEqual(run['config']['N'], 2)
                self.assertEqual(run['config']['B'], 3)
            with exported['csv'].open(newline='', encoding='utf-8') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 9)
            for row in rows[1:]:
                for method, summary in expected.items():
                    self.assertEqual(float(row[method]), summary[row['metric']])

    def test_bootstrap_checkpoint_resume_matches_uninterrupted_run(self):
        from regression_results import BootstrapRun

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            baseline = run_small('Bootstrapshowalter', output_root=root / 'baseline')
            committed = []
            original_record = BootstrapRun.record

            def interrupt_after_commit(run, index, **kwargs):
                original_record(run, index, **kwargs)
                committed.append(run.directory)
                raise KeyboardInterrupt('simulated interruption after a complete replication')

            with patch.object(BootstrapRun, 'record', interrupt_after_commit):
                with self.assertRaises(KeyboardInterrupt):
                    run_small('Bootstrapshowalter', output_root=root / 'interrupted')
            self.assertEqual(len(committed), 1)
            resumed = run_small('Bootstrapshowalter', output_root=root / 'interrupted',
                                resume_from=committed[0])
            self.assertEqual(resumed['run'].summary, baseline['run'].summary)
            self.assertEqual(resumed['i_showalter'], baseline['i_showalter'])
            for key in ('Sigma_showalter', 'Norm_showalter', 'Inf_norm_showalter'):
                self.assertEqual(resumed[key], baseline[key])
            self.assertTrue((Path(resumed['run'].directory) / 'summary.csv').exists())

    def test_bootstrap_threshold_percentiles_match_current_notebook_settings(self):
        from regression_results import PERCENTILES

        expected = {
            'Bootstrapcut': ('SC', 50), 'Bootstrapridge': ('Ridge', 0),
            'Bootstraplandweber': ('Landweber', 72.5),
            'Bootstrapshowalter': ('Showalter', 90), 'BootstrapARk': ('ARk', 0),
            'Bootstrapsecond': ('SOAR', 60), 'BootstrapNesterov': ('Nesterov', 95),
        }
        self.assertEqual(PERCENTILES, dict(expected.values()))
        for name, (method, percentile) in expected.items():
            with self.subTest(notebook=name):
                nodes = [node for cell in cells(name) if cell['cell_type'] == 'code'
                         for node in ast.walk(ast.parse(''.join(cell['source'])))]
                configuration = next(node.value for node in nodes
                    if isinstance(node, ast.Assign) and len(node.targets) == 1
                    and isinstance(node.targets[0], ast.Name) and node.targets[0].id == 'run_config')
                stored = {key.value: value for key, value in zip(configuration.keys, configuration.values)}
                self.assertEqual(ast.literal_eval(stored['method']), method)
                self.assertEqual(ast.literal_eval(stored['threshold_percentile']), percentile)
                selections = [node for node in nodes if isinstance(node, ast.Assign)
                    and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)
                    and node.targets[0].id == 'avg_index'
                    and isinstance(node.value, ast.Call)]
                if method == 'Ridge':
                    self.assertEqual(selections, [])
                    selections = [node for node in nodes if isinstance(node, ast.Assign)
                        and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)
                        and node.targets[0].id == 'avg_index'][:1]
                self.assertEqual(len(selections), 1)
                env = {'np': np, 'min_indices': np.arange(3, 14)}
                exec(compile(ast.Module(body=selections, type_ignores=[]), '<threshold choice>', 'exec'), env)
                self.assertEqual(env['avg_index'], np.percentile(env['min_indices'], percentile))

    def test_ark_kappa_retains_each_experiment_code_setting(self):
        expected = {'normal': 1.5, 'Sparse': 1.5, 'sparse15': 1.5,
                    'Nonsparse': 0.5, 'NonsparseDP': 0.5, 'BootstrapARk': 0.5}
        for name, kappa in expected.items():
            with self.subTest(notebook=name):
                values = [ast.literal_eval(node.value)
                          for cell in cells(name) if cell['cell_type'] == 'code'
                          for node in ast.walk(ast.parse(''.join(cell['source'])))
                          if isinstance(node, ast.Assign) and len(node.targets) == 1
                          and isinstance(node.targets[0], ast.Name) and node.targets[0].id == 'kappa']
                self.assertEqual(values, [kappa])

    def test_notebook_configuration_contract(self):
        for name in ('normal', 'Sparse', 'sparse15', 'Nonsparse', 'NonsparseDP'):
            nonsparse = name.startswith('Nonsparse')
            constants = {}
            for cell in cells(name):
                if cell['cell_type'] != 'code':
                    continue
                tree = ast.parse(''.join(cell['source']))
                constants.update({node.targets[0].id: node.value.value for node in tree.body
                             if isinstance(node, ast.Assign) and len(node.targets) == 1
                             and isinstance(node.targets[0], ast.Name)
                             and isinstance(node.value, ast.Constant)})
                methods = {node.name for node in tree.body if isinstance(node, ast.FunctionDef)}
                if 'p' in constants:
                    self.assertEqual(constants['p'], 1500 if name == 'sparse15' else 1000)
                if 'Nesterov' in methods:
                    self.assertEqual(constants['omega'], 3 if nonsparse else 5)
                if 'Frac' in methods:
                    self.assertEqual(constants['t'], 5e-4 if nonsparse else 5e-5)
                    self.assertEqual(constants['theta'], 1.8 if nonsparse else 1.5)
                    self.assertEqual(constants['max_iter'], 2000 if nonsparse else 500)
                if 'ARk' in methods:
                    self.assertEqual(constants['max_iter'], 2000 if nonsparse else 1000)


if __name__ == '__main__':
    unittest.main()
