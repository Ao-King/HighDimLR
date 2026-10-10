"""Section 6.2 spectral selection against both original notebook recurrences."""

import ast
import copy
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nonsparse_experiment import METHODS, ITERATIVE_METHODS, formal_settings, select_grid, select_iterative
from regression_common import ridge_pair
from regression_spectral import SpectralDesign


NAMES = {"Landweber": "landweber", "Showalter": "showalter", "HBF": "hbf",
         "ARk": "ARk", "SOAR": "second", "Nesterov": "Nesterov", "FAR": "Frac"}
NOTEBOOKS = {"AOSR": "Nonsparse.ipynb", "TDP": "NonsparseDP.ipynb"}


def cells(rule):
    notebook = json.loads((ROOT / NOTEBOOKS[rule]).read_text(encoding="utf-8"))
    return [ast.parse("".join(cell["source"])) for cell in notebook["cells"] if cell["cell_type"] == "code"]


def direct(rule, method, X, Y, truth, count, settings, fixed=False):
    namespace = {"np": np, "X": X, "beta": truth, "delta": settings.get("delta", 0),
                 "t": settings["step"], "kappa": settings.get("kappa", .5),
                 "s2": X.shape[1], "I": np.eye(X.shape[1])}
    function_name = NAMES[method]
    for tree in cells(rule):
        definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
        if method == "ARk" and rule == "AOSR":
            for node in definitions:
                if node.name == "D":
                    exec(compile(ast.Module(body=[copy.deepcopy(node)], type_ignores=[]), "<notebook D>", "exec"), namespace)
                    build_d = namespace["D"]
                    namespace["D"] = np.asarray([build_d(k) for k in range(count+1)])
        if not any(node.name == function_name for node in definitions):
            continue
        if method == "FAR":
            prefix = []
            for node in tree.body:
                if isinstance(node, ast.FunctionDef):
                    break
                node = copy.deepcopy(node)
                if isinstance(node, ast.Assign) and len(node.targets) == 1:
                    target = node.targets[0]
                    replacements = {"t": settings["step"], "max_iter": count, "theta": settings.get("theta", 1.8)}
                    if isinstance(target, ast.Name) and target.id in replacements:
                        node.value = ast.copy_location(ast.Constant(replacements[target.id]), node.value)
                prefix.append(node)
            exec(compile(ast.fix_missing_locations(ast.Module(body=prefix, type_ignores=[])), "<notebook FAR weights>", "exec"), namespace)
        for node in definitions:
            if node.name in (function_name, function_name+"_debias"):
                exec(compile(ast.Module(body=[copy.deepcopy(node)], type_ignores=[]), "<notebook recurrence>", "exec"), namespace)
    arguments = [X, Y, count] if method == "FAR" else [X, Y, settings["step"], count]
    if method == "Nesterov":
        if not fixed:
            arguments.append(truth if rule == "AOSR" else settings["delta"])
        arguments.append(settings.get("omega", 3))
    return namespace[function_name+("_debias" if fixed else "")](*arguments)


def options(method, **changes):
    result = {"step": .05, "max_steps": 7, "minimum_steps": 500, "kappa": .5,
              "s_star": .5, "omega": 3, "eta": 5., "theta": 1.8, "delta": 0.}
    result.update(changes)
    return result


class NonsparseSettingsTests(unittest.TestCase):
    def test_formal_settings_match_both_notebooks_ast_constants(self):
        design = SpectralDesign(np.diag(np.sqrt([537050.2834, 1000., 1.])))
        for rule in NOTEBOOKS:
            actual = formal_settings(design, rule)
            self.assertEqual(tuple(actual), METHODS)
            constants = {}
            seen = set()
            for tree in cells(rule):
                for node in tree.body:
                    if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
                        try:
                            constants[node.targets[0].id] = ast.literal_eval(node.value)
                        except (ValueError, TypeError):
                            pass
                for method, name in NAMES.items():
                    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == name]
                    if not definitions:
                        continue
                    seen.add(method)
                    settings = actual[method]
                    self.assertEqual(settings["max_steps"], constants["max_iter"])
                    if method == "Nesterov":
                        self.assertEqual(settings["nominal_step"], constants["nominal_nesterov_step"])
                        self.assertEqual(settings["step"], min(constants["nominal_nesterov_step"], .99/design.s[0]**2))
                        self.assertEqual(settings["omega"], constants["omega"])
                        assignment = next(node for node in tree.body if isinstance(node, ast.Assign)
                                          and any(isinstance(target, ast.Name) and target.id == "t" for target in node.targets))
                        self.assertEqual(ast.dump(assignment.value), ast.dump(ast.parse("min(nominal_nesterov_step, 0.99 / S[0]**2)", mode="eval").body))
                    else:
                        self.assertEqual(settings["step"], constants["t"])
                    for setting in ("theta", "kappa"):
                        if setting in settings:
                            self.assertEqual(settings[setting], constants[setting])
                    for name, setting in (("eta", "eta"), ("s", "s_star")):
                        if setting in settings:
                            local = next(node for node in definitions[0].body if isinstance(node, ast.Assign)
                                         and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)
                                         and node.targets[0].id == name)
                            self.assertEqual(settings[setting], ast.literal_eval(local.value))
                    if rule == "AOSR":
                        self.assertEqual(settings["minimum_steps"], 500)
                        literals = [node.value for node in ast.walk(definitions[0]) if isinstance(node, ast.Constant)]
                        self.assertIn(500, literals)
                    else:
                        self.assertEqual(settings["tau"], constants["tau"])
                calls = [node.func.id for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)]
                for method, function_name in (("Ridge", "ridge_pair"), ("Lasso", "lasso_coefficients")):
                    if function_name in calls and any(isinstance(node, ast.FunctionDef) and node.name == "compute_mse" for node in tree.body):
                        self.assertEqual(actual[method]["C"], constants["C"])
                        self.assertEqual(actual[method]["m"], constants["m"])
                        self.assertEqual(actual[method]["alpha_divisor"], design.n)
            self.assertEqual(seen, set(ITERATIVE_METHODS))
            self.assertEqual(actual["SC"]["selection"], "oracle_rank")
            self.assertEqual(actual["LS"]["rcond"], None)


class NonsparseIterativeTests(unittest.TestCase):
    def test_formal_steps_and_orders_match_direct_large_mode_recursions(self):
        X = np.diag(np.sqrt([537050.2834, 1000., 1.]))
        truth, noise = np.array([.4, -.2, 1.]), np.array([2., -2., 2.])
        Y, design = X@truth+noise, SpectralDesign(X)
        for rule in NOTEBOOKS:
            for method, setting in formal_settings(design, rule).items():
                if method not in ITERATIVE_METHODS:
                    continue
                settings = dict(setting, max_steps=80, delta=.6*np.linalg.norm(noise))
                expected, count = direct(rule, method, X, Y, truth, 80, settings)
                refit = direct(rule, method, X, X@expected, truth, count, settings, fixed=True)
                ordinary, debiased, selected, _, _ = select_iterative(design, Y, truth, method, settings, rule)
                self.assertEqual(selected, count, (method, rule))
                np.testing.assert_allclose(ordinary, expected, rtol=5e-10, atol=2e-12)
                np.testing.assert_allclose(debiased, 2*expected-refit, rtol=5e-10, atol=2e-12)

    def test_both_rules_and_both_fits_match_original_notebook_functions(self):
        rng = np.random.default_rng(913)
        for X in (rng.normal(size=(6, 4))*.2, rng.normal(size=(4, 6))*.2):
            Y, truth = rng.normal(size=X.shape[0]), rng.normal(size=X.shape[1])
            design = SpectralDesign(X)
            for rule in NOTEBOOKS:
                for method in ITERATIVE_METHODS:
                    settings = options(method)
                    if rule == "TDP":
                        b2 = direct(rule, method, X, Y, truth, 2, settings, fixed=True)
                        b3 = direct(rule, method, X, Y, truth, 3, settings, fixed=True)
                        settings["delta"] = (np.linalg.norm(Y-X@b2)+np.linalg.norm(Y-X@b3))/2
                    with self.subTest(shape=X.shape, method=method, rule=rule):
                        expected, count = direct(rule, method, X, Y, truth, settings["max_steps"], settings)
                        refit = direct(rule, method, X, X@expected, truth, count, settings, fixed=True)
                        ordinary, debiased, selected, diagnostic, trace = select_iterative(design, Y, truth, method, settings, rule)
                        self.assertEqual(selected, count)
                        self.assertEqual(trace.shape[1], 3)
                        self.assertEqual(diagnostic["recorded_steps"], len(trace))
                        np.testing.assert_allclose(ordinary, expected, rtol=2e-11, atol=2e-13)
                        np.testing.assert_allclose(debiased, 2*expected-refit, rtol=2e-11, atol=2e-13)
                        np.testing.assert_allclose(trace[count-1, 1:], [np.linalg.norm(expected-truth), np.linalg.norm(Y-X@expected)], rtol=2e-11, atol=2e-13)

    def test_actual_minimum_500_and_first_turn_match_notebooks(self):
        X, Y = np.array([[.5]]), np.array([2.])
        design = SpectralDesign(X)
        for method in ITERATIVE_METHODS:
            settings = options(method, step=.01, max_steps=600)
            truth = direct("AOSR", method, X, Y, np.zeros(1), 10, settings, fixed=True)
            expected, count = direct("AOSR", method, X, Y, truth, 600, settings)
            ordinary, _, selected, diagnostic, trace = select_iterative(design, Y, truth, method, settings, "AOSR")
            self.assertEqual(count, 500, method)
            self.assertEqual(selected, count)
            self.assertEqual(diagnostic["first_turn_count"], 10)
            self.assertEqual(len(trace), 500)
            np.testing.assert_allclose(ordinary, expected, rtol=2e-11, atol=2e-13)

    def test_zero_cap_nonhit_cap_and_no_stop_before_first_update(self):
        X, Y, truth = np.eye(2), np.ones(2), np.zeros(2)
        design = SpectralDesign(X)
        for rule in NOTEBOOKS:
            for method in ITERATIVE_METHODS:
                for cap in (0, 1, 7):
                    settings = options(method, max_steps=cap, delta=0.)
                    expected, count = direct(rule, method, X, Y, truth, cap, settings)
                    ordinary, _, selected, diagnostic, trace = select_iterative(design, Y, truth, method, settings, rule)
                    self.assertEqual(selected, count)
                    self.assertEqual(selected, cap)
                    self.assertTrue(diagnostic["at_iteration_cap"])
                    np.testing.assert_allclose(ordinary, expected, rtol=2e-11, atol=2e-13)
                    self.assertEqual(trace.shape, (cap, 3))
                result = select_iterative(design, Y, truth, method, options(method, delta=100.), "TDP")
                self.assertEqual(result[2], 1)

    def test_trigger_iterate_and_direct_boundary_values_are_recorded(self):
        X, Y = np.eye(2), np.ones(2)
        design = SpectralDesign(X)
        truth = direct("AOSR", "Landweber", X, Y, np.zeros(2), 4, options("Landweber"), fixed=True)
        result = select_iterative(design, Y, truth, "Landweber", options("Landweber", minimum_steps=0), "AOSR")
        self.assertEqual(result[2], 4)
        self.assertEqual(result[4][-1, 0], 5)
        self.assertTrue(result[3]["trace_includes_trigger_after_selected"])
        delta = np.linalg.norm(Y-X@(0.05*Y))
        result = select_iterative(design, Y, truth, "Landweber", options("Landweber", delta=delta), "TDP")
        self.assertEqual(result[2], 1)
        self.assertEqual(result[3]["boundary_checks"][0]["direct_values"], [delta])
        result = select_iterative(design, np.zeros(2), np.ones(2), "Landweber", options("Landweber"), "AOSR")
        self.assertEqual(result[3]["first_turn_count"], 1)
        self.assertEqual(result[3]["boundary_checks"][0]["steps"], [1, 2])


class NonsparseGridTests(unittest.TestCase):
    def test_grid_choices_and_complete_curves_match_direct_fits(self):
        rng = np.random.default_rng(574)
        X, Y, truth = rng.normal(size=(4, 6)), rng.normal(size=4), rng.normal(size=6)
        design = SpectralDesign(X)
        for rule in NOTEBOOKS:
            settings = {"C": 8, "m": 2, "alpha_divisor": 4}
            numerators = list(range(0, 9, 2)) if rule == "AOSR" else list(range(8, -1, -2))
            fits = [ridge_pair(X, Y, c/4) for c in numerators]
            errors = np.array([np.linalg.norm(fit[0]-truth) for fit in fits])
            residuals = np.array([np.linalg.norm(Y-X@fit[0]) for fit in fits])
            if rule == "AOSR":
                selected = int(np.argmin(errors))
            else:
                settings["delta"] = (residuals[1]+residuals[2])/2
                selected = int(np.flatnonzero(residuals <= settings["delta"])[0])
            ordinary, debiased, parameter, diagnostic, trace = select_grid(design, Y, truth, "Ridge", settings, rule)
            self.assertEqual(parameter, numerators[selected]/4)
            self.assertEqual(diagnostic["selected_grid_numerator"], numerators[selected])
            np.testing.assert_allclose(ordinary, fits[selected][0], rtol=2e-11, atol=2e-13)
            np.testing.assert_allclose(debiased, fits[selected][1], rtol=2e-11, atol=2e-13)
            np.testing.assert_allclose(trace, np.column_stack((np.array(numerators)/4, errors, residuals)), rtol=2e-11, atol=2e-13)
            sc = select_grid(design, Y, truth, "SC", {}, rule)
            expected_rank = min(range(4, 0, -1), key=lambda k: np.linalg.norm(design.pair(Y, "SC", k, {})[0]-truth))
            self.assertEqual(sc[2], expected_rank)
            np.testing.assert_array_equal(sc[0], sc[1])
            self.assertFalse(np.shares_memory(sc[0], sc[1]))
            self.assertEqual(sc[3]["selection_exception"], "oracle_rank_under_both_rules")

    def test_ties_singular_zero_alpha_and_missing_tdp_candidate(self):
        X = np.diag([1., np.finfo(float).eps/8, 0.])
        design = SpectralDesign(X)
        for rule in NOTEBOOKS:
            self.assertEqual(select_grid(design, np.zeros(3), np.zeros(3), "SC", {}, rule)[2], 3)
        ridge = select_grid(design, np.zeros(3), np.zeros(3), "Ridge", {"C": 4, "m": 2}, "AOSR")
        self.assertEqual(ridge[2], 0)
        Y = np.ones(3)
        ridge = select_grid(design, Y, np.zeros(3), "Ridge", {"C": 0, "m": 1}, "AOSR")
        np.testing.assert_array_equal(ridge[0], np.linalg.lstsq(X, Y, rcond=None)[0])
        with self.assertRaisesRegex(ValueError, "No regularization-grid"):
            select_grid(design, Y, np.zeros(3), "Ridge", {"C": 4, "m": 2, "delta": 0.}, "TDP")

    def test_invalid_inputs_and_rule_mismatches_fail(self):
        design = SpectralDesign(np.eye(2))
        with self.assertRaises(ValueError):
            formal_settings(design, "wrong")
        for changes in ({"max_steps": -1}, {"max_steps": True}, {"step": 0}, {"minimum_steps": -1}):
            with self.assertRaises(ValueError):
                select_iterative(design, np.ones(2), np.zeros(2), "Landweber", options("Landweber", **changes), "AOSR")
        with self.assertRaises(ValueError):
            select_iterative(design, np.ones(2), np.zeros(2), "FAR", {"step": .1, "max_steps": 1}, "TDP")
        with self.assertRaises(ValueError):
            select_iterative(design, np.ones(2), np.zeros(2), "HBF", options("HBF", selection="AOSR"), "TDP")
        with self.assertRaises(FloatingPointError):
            select_iterative(design, np.ones(2), np.array([np.nan, 0]), "Landweber", options("Landweber"), "AOSR")
        for method, settings in (("SC", {"rank_max": 0}), ("SC", {"rank_max": 3}),
                                  ("Ridge", {"C": 4, "m": 0}), ("Ridge", {"C": 4, "m": 1, "alpha_divisor": 0})):
            with self.assertRaises(ValueError):
                select_grid(design, np.ones(2), np.zeros(2), method, settings, "AOSR")


if __name__ == "__main__":
    unittest.main()
