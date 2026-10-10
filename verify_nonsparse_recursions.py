"""Verify archived 6.2 recurrences in full SVD coordinates and short dense runs.

The full-coordinate comparison preserves the original notebook functions and
stopping tests; it is not a complete long-trajectory dense-BLAS replay.
"""

import argparse
import ast
import copy
from datetime import datetime, timezone
import hashlib
import inspect
import json
from pathlib import Path
import sys
import time

import numpy as np
from threadpoolctl import threadpool_limits


METHODS = {"Landweber": "landweber", "Showalter": "showalter", "HBF": "hbf",
           "ARk": "ARk", "SOAR": "second", "Nesterov": "Nesterov", "FAR": "Frac"}
NOTEBOOKS = {"AOSR": "Nonsparse.ipynb", "TDP": "NonsparseDP.ipynb"}
SETTING_NAMES = {"t", "max_iter", "kappa", "omega", "theta", "nominal_nesterov_step", "effective_nesterov_step"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def compare(actual, expected, rtol=1e-8, atol=1e-9):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if actual.shape != expected.shape or not np.isfinite(actual).all() or not np.isfinite(expected).all():
        return {"close": False, "shape_actual": list(actual.shape), "shape_expected": list(expected.shape)}
    error = actual - expected
    return {"close": bool(np.allclose(actual, expected, rtol=rtol, atol=atol)),
            "max_absolute_difference": float(np.max(np.abs(error), initial=0)),
            "l2_difference": float(np.linalg.norm(error))}


class Diagonal:
    """Exact diagonal matrix algebra, without allocating dense zero entries."""
    __array_priority__ = 10000

    def __init__(self, values):
        self.values = np.asarray(values, dtype=float)
        require(self.values.ndim == 1, "A diagonal must be a vector")
        self.shape = (len(self.values), len(self.values))

    @property
    def T(self):
        return self

    def __matmul__(self, other):
        if isinstance(other, Diagonal):
            return Diagonal(self.values * other.values)
        other = np.asarray(other)
        return self.values * other if other.ndim == 1 else self.values[:, None] * other

    def __mul__(self, scalar):
        require(np.ndim(scalar) == 0, "Only scalar matrix scaling is expected")
        return Diagonal(self.values * scalar)

    __rmul__ = __mul__

    def __add__(self, other):
        require(isinstance(other, Diagonal), "Only diagonal matrix addition is expected")
        return Diagonal(self.values + other.values)

    def __neg__(self):
        return Diagonal(-self.values)


class Horizontal:
    def __init__(self, blocks):
        self.blocks = blocks


class BlockMatrix:
    def __init__(self, rows):
        self.rows = rows

    def __matmul__(self, vector):
        require(len(self.rows) == 2 and all(len(row.blocks) == 2 for row in self.rows), "Expected the HBF 2x2 block matrix")
        n = self.rows[0].blocks[0].shape[0]
        first, second = vector[:n], vector[n:]
        return np.concatenate([row.blocks[0] @ first + row.blocks[1] @ second for row in self.rows])


class LinalgProxy:
    def __getattr__(self, name):
        return getattr(np.linalg, name)

    @staticmethod
    def inv(value):
        require(isinstance(value, Diagonal), "Coordinate inverse must be diagonal")
        with np.errstate(divide="raise", invalid="raise"):
            return Diagonal(1.0 / value.values)


class NumpyProxy:
    """Only matrix algebra is specialized; vector arithmetic and norms stay NumPy."""
    linalg = LinalgProxy()

    def __getattr__(self, name):
        return getattr(np, name)

    @staticmethod
    def dot(first, second):
        return first @ second if isinstance(first, (Diagonal, BlockMatrix)) else np.dot(first, second)

    @staticmethod
    def eye(size):
        return Diagonal(np.ones(size))

    @staticmethod
    def hstack(items):
        return Horizontal(items) if all(isinstance(item, Diagonal) for item in items) else np.hstack(items)

    @staticmethod
    def vstack(items):
        return BlockMatrix(items) if all(isinstance(item, Horizontal) for item in items) else np.vstack(items)


class MatrixZeros(ast.NodeTransformer):
    """Replace only the HBF square zero block, not FAR's history arrays."""
    def visit_Call(self, node):
        node = self.generic_visit(node)
        if (isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "np" and node.func.attr == "zeros" and len(node.args) == 1
                and isinstance(node.args[0], ast.Tuple) and len(node.args[0].elts) == 2
                and all(isinstance(item, ast.Name) and item.id == "m" for item in node.args[0].elts)):
            return ast.copy_location(ast.Call(func=ast.Name(id="_zero_matrix", ctx=ast.Load()),
                                             args=[ast.Name(id="m", ctx=ast.Load())], keywords=[]), node)
        return node


class SmokeMinimum(ast.NodeTransformer):
    """Change only the explicit AOSR floor for labelled reduced test fixtures."""
    def __init__(self, minimum):
        self.minimum = minimum

    def visit_Call(self, node):
        node = self.generic_visit(node)
        if (isinstance(node.func, ast.Name) and node.func.id == "max" and len(node.args) == 2
                and isinstance(node.args[0], ast.BinOp) and isinstance(node.args[0].op, ast.Sub)
                and isinstance(node.args[0].left, ast.Name) and node.args[0].left.id == "completed_steps"
                and isinstance(node.args[1], ast.Constant) and node.args[1].value == 500):
            node.args[1] = ast.copy_location(ast.Constant(self.minimum), node.args[1])
        return node


class LazyMatrices:
    def __init__(self, function):
        self.function, self.cache = function, {}

    def __getitem__(self, index):
        if index not in self.cache:
            self.cache[index] = self.function(index)
        return self.cache[index]


def load_kernels(notebook_path, method, X, truth, singular_values, delta, *, diagonal=False, cap=None,
                 minimum_steps=None):
    """Load original functions; only replace square-zero construction in diagonal mode."""
    name = METHODS[method]
    notebook = json.loads(Path(notebook_path).read_text(encoding="utf-8"))
    namespace = {"np": NumpyProxy() if diagonal else np, "X": X, "beta": truth,
                 "S": singular_values, "delta": delta, "s2": X.shape[1],
                 "_zero_matrix": lambda size: Diagonal(np.zeros(size))}
    namespace["I"] = namespace["np"].eye(X.shape[1])
    inverse_definition = None
    for index, cell in enumerate(notebook["cells"]):
        if cell.get("cell_type") != "code":
            continue
        tree = ast.parse("".join(cell["source"]))
        definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
        for definition in definitions:
            if definition.name == "D":
                inverse_definition = copy.deepcopy(definition)
        for statement in tree.body:
            if (isinstance(statement, ast.Assign) and len(statement.targets) == 1
                    and isinstance(statement.targets[0], ast.Name)
                    and statement.targets[0].id in SETTING_NAMES):
                exec(compile(ast.Module(body=[copy.deepcopy(statement)], type_ignores=[]), str(notebook_path), "exec"), namespace)
        if not any(definition.name == name for definition in definitions):
            continue
        options = {key: namespace[key] for key in SETTING_NAMES if key in namespace}
        if cap is not None:
            options["max_iter"] = int(cap)
        if minimum_steps is not None:
            options["minimum_steps"] = int(minimum_steps)
        if method == "FAR":
            prefix = []
            for node in tree.body:
                if isinstance(node, ast.FunctionDef):
                    break
                node = copy.deepcopy(node)
                if (isinstance(node, ast.Assign) and len(node.targets) == 1
                        and isinstance(node.targets[0], ast.Name) and node.targets[0].id == "max_iter"):
                    node.value = ast.copy_location(ast.Constant(options["max_iter"]), node.value)
                prefix.append(node)
            exec(compile(ast.fix_missing_locations(ast.Module(body=prefix, type_ignores=[])),
                         f"{notebook_path}:FAR-weights", "exec"), namespace)
        if method == "ARk" and inverse_definition is not None:
            exec(compile(ast.Module(body=[inverse_definition], type_ignores=[]), str(notebook_path), "exec"), namespace)
            namespace["D"] = LazyMatrices(namespace["D"])
        for definition in definitions:
            if definition.name not in (name, name + "_debias"):
                continue
            definition = copy.deepcopy(definition)
            if minimum_steps is not None:
                definition = SmokeMinimum(minimum_steps).visit(definition)
            if diagonal:
                definition = MatrixZeros().visit(definition)
            exec(compile(ast.fix_missing_locations(ast.Module(body=[definition], type_ignores=[])),
                         f"{notebook_path}:cell-{index}:{definition.name}", "exec"), namespace)
        return namespace[name], namespace[name + "_debias"], options
    raise ValueError(f"Missing {method} functions in {notebook_path}")


def invoke(function, X, Y, options, count, truth, delta):
    arguments = {"X": X, "Y": Y, "t": options["t"], "max_iter": count, "k_max": count,
                 "beta_reference": truth, "delta": delta, "omega": options.get("omega", 3)}
    if "minimum_steps" in options:
        arguments["k_min"] = options["minimum_steps"]
    return function(**{name: arguments[name] for name in inspect.signature(function).parameters
                       if name in arguments})


def verify_method(notebook, rule, method, X, Y, truth, U, s, Vt, saved_ordinary, saved_debiased,
                  expected_count, delta, *, dense_steps=(1, 2, 5), cap=None, minimum_steps=None,
                  rtol=1e-8, atol=1e-9, full_dense_methods=()):
    """Full original selection in coordinates, fixed-k debias, and short dense checks."""
    coordinates = U.T @ Y
    truth_coordinates = Vt @ truth
    diagonal = Diagonal(s)
    selected, fixed, options = load_kernels(notebook, method, diagonal, truth_coordinates, s, delta,
                                            diagonal=True, cap=cap, minimum_steps=minimum_steps)
    timings = {}
    start = time.perf_counter()
    ordinary, count = invoke(selected, diagonal, coordinates, options, options["max_iter"], truth_coordinates, delta)
    timings["full_coordinate_selection"] = time.perf_counter() - start
    count = int(count)
    start = time.perf_counter()
    if method == "Landweber":
        debiased = invoke(fixed, diagonal, coordinates, options, 2 * count, truth_coordinates, delta)
    else:
        refit = invoke(fixed, diagonal, diagonal @ ordinary, options, count, truth_coordinates, delta)
        debiased = 2 * ordinary - refit
    timings["fixed_parameter_coordinate_debias"] = time.perf_counter() - start
    checks = {"selected_count": {"close": count == expected_count, "reference": count, "archived": expected_count},
              "ordinary": compare(Vt.T @ ordinary, saved_ordinary, rtol, atol),
              "debiased": compare(Vt.T @ debiased, saved_debiased, rtol, atol)}
    # Short dense checks run the untouched fixed-iteration functions, not selectors.
    short_cap = max(dense_steps)
    _, dense_fixed, dense_options = load_kernels(notebook, method, X, truth, s, delta, cap=short_cap)
    _, diagonal_fixed, diagonal_options = load_kernels(notebook, method, diagonal, truth_coordinates, s, delta,
                                                       diagonal=True, cap=short_cap)
    start = time.perf_counter()
    dense_checks = {}
    for steps in dense_steps:
        direct = invoke(dense_fixed, X, Y, dense_options, steps, truth, delta)
        coordinate = invoke(diagonal_fixed, diagonal, coordinates, diagonal_options, steps, truth_coordinates, delta)
        dense_checks[str(steps)] = compare(direct, Vt.T @ coordinate, rtol, atol)
    timings["short_dense_recurrences"] = time.perf_counter() - start
    full_dense = {"enabled": method in full_dense_methods}
    if full_dense["enabled"]:
        # A separate namespace retains the full cap, including FAR's full weights.
        start = time.perf_counter()
        dense_selected, dense_refit, full_options = load_kernels(
            notebook, method, X, truth, s, delta, cap=cap, minimum_steps=minimum_steps)
        timings["full_dense_setup"] = time.perf_counter() - start
        start = time.perf_counter()
        dense_ordinary, dense_count = invoke(dense_selected, X, Y, full_options,
                                            full_options["max_iter"], truth, delta)
        dense_count = int(dense_count)
        timings["full_dense_selection"] = time.perf_counter() - start
        start = time.perf_counter()
        if method == "Landweber":
            dense_debiased = invoke(dense_refit, X, Y, full_options, 2 * dense_count, truth, delta)
        else:
            dense_second = invoke(dense_refit, X, X @ dense_ordinary, full_options, dense_count, truth, delta)
            dense_debiased = 2 * dense_ordinary - dense_second
        timings["full_dense_fixed_parameter_debias"] = time.perf_counter() - start
        full_dense.update(selected_count=dense_count, iteration_cap=full_options["max_iter"], checks={
            "selected_count_archived": {"close": dense_count == expected_count,
                                        "reference": dense_count, "archived": expected_count},
            "selected_count_coordinate": {"close": dense_count == count,
                                          "dense": dense_count, "coordinate": count},
            "ordinary_archived": compare(dense_ordinary, saved_ordinary, rtol, atol),
            "debiased_archived": compare(dense_debiased, saved_debiased, rtol, atol),
            "ordinary_coordinate": compare(dense_ordinary, Vt.T @ ordinary, rtol, atol),
            "debiased_coordinate": compare(dense_debiased, Vt.T @ debiased, rtol, atol)})
        full_dense["passed"] = all(value["close"] for value in full_dense["checks"].values())
    residual = float(np.linalg.norm(Y - X @ (Vt.T @ ordinary)))
    result = {"rule": rule, "method": method, "selected_count": count, "archived_count": expected_count,
              "iteration_cap": options["max_iter"], "at_cap": count == options["max_iter"],
              "residual_norm": residual, "delta": delta if rule in ("TDP", "DP") else None,
              "tdp_satisfied": bool(residual <= delta) if rule in ("TDP", "DP") else None,
              "checks": checks, "dense_prefix_checks": dense_checks, "full_dense_reference": full_dense,
              "timings_seconds": timings}
    result["passed"] = (all(value["close"] for value in checks.values())
                        and all(value["close"] for value in dense_checks.values())
                        and full_dense.get("passed", True))
    return result


def verify_run(run_directory, *, dense_steps=(1, 2, 5), rtol=1e-8, atol=1e-9, allow_nonformal=False, cap=None,
               minimum_steps=None, full_dense_methods=()):
    directory = Path(run_directory).resolve()
    summary = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
    require(summary.get("status") == "complete", "Verification requires a completed run")
    config = summary.get("config", summary.get("configuration", {}))
    require(allow_nonformal or (config.get("n"), config.get("p")) == (1000, 1000), "Nonformal dimensions require explicit permission")
    require(allow_nonformal or not config.get("smoke_test", False), "Smoke results cannot be verified as formal")
    require(cap is None or allow_nonformal, "Changing the reference cap is only allowed for smoke fixtures")
    require(cap is None or (isinstance(cap, int) and cap > 0), "Invalid reference iteration cap")
    require(minimum_steps is None or (allow_nonformal and config.get("smoke_test") is True),
            "Changing the reference AOSR floor is only allowed for smoke fixtures")
    require(minimum_steps is None or (isinstance(minimum_steps, int) and minimum_steps > 0),
            "Invalid reference AOSR floor")
    require(dense_steps and all(isinstance(value, int) and value > 0 for value in dense_steps), "Invalid dense prefix lengths")
    require(all(method in METHODS for method in full_dense_methods), "Unknown full dense reference method")
    require(len(set(full_dense_methods)) == len(full_dense_methods), "Duplicate full dense reference methods")
    for field in ("sources", "artifact_files"):
        for relative, digest in summary.get(field, {}).items():
            path = (directory / relative).resolve()
            require(path.is_relative_to(directory), f"Archived path escapes run directory: {relative}")
            require(hashlib.sha256(path.read_bytes()).hexdigest() == digest,
                    f"Archived file hash mismatch: {relative}")
    notebooks = {}
    for rule, filename in NOTEBOOKS.items():
        path = directory / "sources" / filename
        expected = summary["sources"].get("sources/" + filename, summary["sources"].get(filename))
        require(expected is not None and hashlib.sha256(path.read_bytes()).hexdigest() == expected,
                f"Archived notebook hash mismatch: {filename}")
        notebooks[rule] = path
    with np.load(directory / "data.npz", allow_pickle=False) as data:
        X, Y, truth, e, U, s, Vt = (data[name].copy() for name in ("X", "Y", "beta", "e", "U", "S", "Vt"))
    n, p = X.shape
    require((n, p) == (config.get("n"), config.get("p")), "Configuration/data dimensions differ")
    require(n == p == len(s) and U.shape == (n, n) and Vt.shape == (p, p), "Full square SVD required for this coordinate verification")
    require(Y.shape == e.shape == (n,) and truth.shape == (p,), "Invalid archived vector shape")
    require(all(np.isfinite(value).all() for value in (X, Y, truth, e, U, s, Vt)), "Nonfinite archived data")
    close_data = {"design_reconstruction": compare((U * s) @ Vt, X, rtol, atol),
                  "response_reconstruction": compare(X @ truth + e, Y, rtol, atol),
                  "left_orthogonality": compare(U.T @ U, np.eye(n), rtol, atol),
                  "right_orthogonality": compare(Vt @ Vt.T, np.eye(p), rtol, atol)}
    rows = {("TDP" if row["rule"] == "DP" else row["rule"], row["method"]): row for row in summary["rows"]}
    require(len(rows) == len(summary["rows"]), "Duplicate archived method/rule rows")
    require(all((rule, method) in rows for rule in NOTEBOOKS for method in METHODS), "Missing iterative rows")
    reports = []
    with np.load(directory / "estimates.npz", allow_pickle=False) as estimates:
        for rule, notebook in notebooks.items():
            for method in METHODS:
                row = rows[rule, method]
                parameter = row.get("k0") if row.get("k0") is not None else row.get("parameter")
                require(parameter is not None and int(parameter) == parameter and parameter >= 1, "Invalid archived iteration count")
                print(f"Verifying {rule} {method}, archived k={parameter}", file=sys.stderr, flush=True)
                saved_rule = row["rule"]
                report = verify_method(notebook, rule, method, X, Y, truth, U, s, Vt,
                                       estimates[f"{saved_rule}__{method}__ordinary"], estimates[f"{saved_rule}__{method}__debiased"],
                                       int(parameter), .6 * float(np.linalg.norm(e)), dense_steps=dense_steps,
                                       cap=cap, minimum_steps=minimum_steps, rtol=rtol, atol=atol,
                                       full_dense_methods=full_dense_methods)
                reports.append(report)
    return {"version": "nonsparse-recursion-verification-v1", "run_directory": str(directory),
            "n": n, "p": p, "rtol": rtol, "atol": atol, "data_checks": close_data, "methods": reports,
            "nonformal_reference_overrides": {"iteration_cap": cap, "minimum_steps": minimum_steps},
            "full_dense_methods": list(full_dense_methods),
            "reference_notebook_sha256": {rule: hashlib.sha256(path.read_bytes()).hexdigest() for rule, path in notebooks.items()},
            "scope": "Original notebook functions: full selection and fixed-count debias in complete SVD coordinates; untouched dense fixed-iteration prefixes for all methods, plus complete original dense selection and fixed-count debias for the explicitly requested full_dense_methods.",
            "matrix_specializations": ["Diagonal products/inverses", "HBF 2x2 diagonal-block action", "Lazy evaluation of original ARk D(k)"],
            "limitations": "Does not replay every long-trajectory step with original dense BLAS for methods outside full_dense_methods; coordinate rounding can expose stopping-boundary differences, which are reported as failures.",
            "passed": all(item["close"] for item in close_data.values()) and all(item["passed"] for item in reports)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dense-steps", nargs="+", type=int, default=[1, 2, 5])
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--allow-nonformal", action="store_true")
    parser.add_argument("--max-steps", type=int, help="Reference cap override for explicitly permitted smoke fixtures")
    parser.add_argument("--minimum-steps", type=int, help="AOSR floor override for explicitly permitted smoke fixtures")
    parser.add_argument("--full-dense-methods", nargs="*", choices=tuple(METHODS), default=[],
                        help="Also replay complete original dense selection and fixed-count debias for these methods")
    args = parser.parse_args(argv)
    output = args.output or args.run_directory / "recursion-verification.json"
    if output.exists():
        parser.error("Verification output exists; choose a new path")
    if not np.isfinite([args.rtol, args.atol]).all() or min(args.rtol, args.atol) < 0:
        parser.error("Invalid comparison tolerance")
    started = time.monotonic()
    try:
        with threadpool_limits(limits=1):
            result = verify_run(args.run_directory, dense_steps=tuple(args.dense_steps), rtol=args.rtol, atol=args.atol,
                                allow_nonformal=args.allow_nonformal, cap=args.max_steps, minimum_steps=args.minimum_steps,
                                full_dense_methods=tuple(args.full_dense_methods))
        result["status"] = "verified" if result["passed"] else "failed"
    except Exception as error:
        result = {"status": "failed", "passed": False, "error": f"{type(error).__name__}: {error}"}
    result.update(created_at=datetime.now(timezone.utc).isoformat(), elapsed_seconds=time.monotonic() - started,
                  verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"Verification {result['status']}: {output}")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
