"""Exercise notebook AST functions against the manuscript's AOSR index rule.

Only function definitions and FAR coefficient construction are executed. Full
Monte Carlo cells are not run. Fixed-iteration trajectories retain the actual
notebook recurrences, with an AST-only hook recording every completed update.
"""

import ast
import copy
import json
import math
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
METHOD_STEPS = {
    "landweber": 0.001,
    "showalter": 0.001,
    "hbf": 0.001,
    "ARk": 0.00005,
    "second": 0.0005,
    "Nesterov": 0.000003,
    "Frac": 0.0005,
}


def notebook_code(path):
    notebook = json.loads(path.read_text(encoding="utf-8"))
    return [ast.parse("".join(cell["source"])) for cell in notebook["cells"]
            if cell["cell_type"] == "code"]


def define(node, namespace):
    node = copy.deepcopy(node)
    node.decorator_list = []
    module = ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[]))
    exec(compile(module, "<notebook function>", "exec"), namespace)


def define_trajectory(node, namespace):
    node = copy.deepcopy(node)
    node.name = "record_trajectory"
    node.decorator_list = []
    node.body.insert(0, ast.parse("trajectory = [np.zeros(X.shape[1])]").body[0])
    for statement in node.body:
        if isinstance(statement, ast.For):
            statement.body.extend(ast.parse("trajectory.append(beta_cur.copy())").body)
        if isinstance(statement, ast.Return):
            statement.value = ast.Name(id="trajectory", ctx=ast.Load())
    define(node, namespace)


def make_case(name, notebook="Nonsparse.ipynb", cap=700, step=None):
    cells = notebook_code(ROOT / notebook)
    definitions = {node.name: node for cell in cells for node in cell.body
                   if isinstance(node, ast.FunctionDef)}
    step = METHOD_STEPS[name] if step is None else step
    namespace = {"np": np, "gamma": math.gamma, "X": np.ones((1, 1)),
                 "Y": np.array([2.0]), "beta": np.zeros(1), "s2": 1,
                 "t": step, "kappa": 0.5, "I": np.eye(1), "theta": 1.8,
                 "max_iter": cap, "omega": 3}
    define(definitions[name], namespace)
    define(definitions[name + "_debias"], namespace)
    define_trajectory(definitions[name + "_debias"], namespace)
    if name == "ARk":
        define(definitions["D"], namespace)
        coefficient = namespace["D"]
        namespace["D"] = np.array([coefficient(k) for k in range(cap + 1)])
    if name == "Frac":
        coefficient_cell = next(cell for cell in cells if any(
            isinstance(node, ast.FunctionDef) and node.name == "Frac" for node in cell.body))
        construction = []
        for node in coefficient_cell.body:
            if isinstance(node, ast.FunctionDef):
                break
            if isinstance(node, ast.For):
                construction.append(node)
            elif isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id in {"A", "B", "g", "k1", "k2"}
                for target in node.targets
            ):
                construction.append(node)
        exec(compile(ast.Module(body=construction, type_ignores=[]),
                     "<notebook FAR coefficients>", "exec"), namespace)
    return Case(name, namespace, step, notebook.startswith("Bootstrap"))


class Case:
    def __init__(self, name, namespace, step, bootstrap):
        self.name = name
        self.ns = namespace
        self.step = step
        self.bootstrap = bootstrap

    def trajectory(self, cap, response=None):
        response = self.ns["Y"] if response is None else response
        function = self.ns["record_trajectory"]
        if self.name == "Frac":
            values = function(self.ns["X"], response, cap)
        elif self.name == "Nesterov":
            values = function(self.ns["X"], response, self.step, cap, 3)
        else:
            values = function(self.ns["X"], response, self.step, cap)
        return np.asarray(values)

    def fit(self, cap, reference, response=None):
        self.ns["beta"] = reference
        response = self.ns["Y"] if response is None else response
        function = self.ns[self.name]
        if self.name == "Frac":
            return function(self.ns["X"], response, cap)
        if self.name == "Nesterov":
            return function(self.ns["X"], response, self.step, cap, reference, 3)
        if self.bootstrap:
            return function(self.ns["X"], response, self.step, cap, reference)
        return function(self.ns["X"], response, self.step, cap)


def aosr_index(trajectory, reference, minimum=500):
    cap = len(trajectory) - 1
    errors = np.linalg.norm(trajectory - reference, axis=1)
    # beta_0 is excluded from the article's comparison set, which starts at 1.
    first_turn = next((step - 1 for step in range(2, cap + 1)
                       if errors[step] >= errors[step - 1]), cap)
    return min(cap, max(first_turn, minimum))


class AOSRSelectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = {name: make_case(name) for name in METHOD_STEPS}
        cls.cases["bootstrap_landweber"] = make_case("landweber", "Bootstraplandweber.ipynb")
        cls.paths = {name: case.trajectory(700) for name, case in cls.cases.items()}

    def assert_selected(self, case, trajectory, reference, expected=None, response=None):
        expected = aosr_index(trajectory, reference) if expected is None else expected
        estimate, steps = case.fit(len(trajectory) - 1, reference, response)
        self.assertEqual(steps, expected)
        np.testing.assert_allclose(estimate, trajectory[expected], rtol=1e-12, atol=1e-14)

    def test_zero_cap_returns_initial_estimate_and_zero(self):
        for name, case in self.cases.items():
            with self.subTest(method=name):
                self.assert_selected(case, self.paths[name][:1], np.array([1.0]), 0)

    def test_cap_before_minimum_returns_actual_cap_iterate(self):
        for name, case in self.cases.items():
            with self.subTest(method=name):
                self.assert_selected(case, self.paths[name][:8], self.paths[name][600], 7)

    def test_first_optimum_before_minimum_is_clipped_to_500(self):
        for name, case in self.cases.items():
            with self.subTest(method=name):
                reference = self.paths[name][10]
                self.assertEqual(aosr_index(self.paths[name], reference), 500)
                self.assert_selected(case, self.paths[name], reference, 500)

    def test_late_optimum_returns_best_not_first_worse_iterate(self):
        for name, case in self.cases.items():
            with self.subTest(method=name):
                reference = self.paths[name][600]
                self.assertEqual(aosr_index(self.paths[name], reference), 600)
                self.assert_selected(case, self.paths[name], reference, 600)

    def test_flat_trajectory_does_not_select_zero(self):
        for name, case in self.cases.items():
            with self.subTest(method=name):
                trajectory = np.zeros((701, 1))
                self.assert_selected(case, trajectory, np.zeros(1), 500, np.zeros(1))

    def test_first_early_turn_is_not_reselected_after_500(self):
        # An oscillatory path can be improving again at the lower cutoff.
        for step in (0.0001, 0.0003, 0.001, 0.003, 0.01):
            case = make_case("Nesterov", step=step)
            path = case.trajectory(700)
            reference = path[10]
            errors = np.linalg.norm(path - reference, axis=1)
            if errors[500] < errors[499] and aosr_index(path, reference) == 500:
                self.assert_selected(case, path, reference, 500)
                return
        self.fail("The test fixture must improve at k_min after an earlier turning point")


if __name__ == "__main__":
    unittest.main()
