"""Versioned, resumable persistence for the regression bootstrap experiments."""

import ast
import csv
import hashlib
from importlib.metadata import PackageNotFoundError, version as package_version
import io
import json
import os
from pathlib import Path
import platform
import tempfile
from datetime import datetime, timezone

import numpy as np

from regression_execution import method_settings


VERSION = "regression-bootstrap-v2"
METHODS = ("SC", "Ridge", "Landweber", "Showalter", "ARk", "SOAR", "Nesterov")
METRICS = ("coverage_original", "coverage_debiased", "l2_original", "l2_debiased",
           "linf_original", "linf_debiased", "variance_absolute_error_original",
           "variance_absolute_error_debiased")
PERCENTILES = dict(zip(METHODS, (50, 0, 72.5, 90, 0, 60, 95)))


def _json(value):
    def convert(item):
        if isinstance(item, np.generic):
            return item.item()
        if isinstance(item, np.ndarray):
            return item.tolist()
        if isinstance(item, Path):
            return str(item)
        raise TypeError(f"Cannot serialize {type(item).__name__}")
    return json.dumps(value, sort_keys=True, allow_nan=False, default=convert)


def _digest(value):
    return hashlib.sha256(_json(value).encode("utf-8")).hexdigest()


def _array_hash(value):
    value = np.ascontiguousarray(value)
    digest = hashlib.sha256(_json([value.dtype.str, value.shape]).encode("ascii"))
    digest.update(memoryview(value).cast("B"))
    return digest.hexdigest()


def _now():
    return datetime.now(timezone.utc).isoformat()


def _runtime():
    result = {"python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform()}
    for package in ("scipy", "scikit-learn", "numba", "joblib", "threadpoolctl"):
        try:
            result[package] = package_version(package)
        except PackageNotFoundError:
            result[package] = None
    return result


def _atomic(path, writer):
    path = Path(path)
    handle, temporary = tempfile.mkstemp(prefix="." + path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            writer(stream)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _write_json(path, value):
    _atomic(path, lambda stream: stream.write((_json(value) + "\n").encode("utf-8")))


def _write_csv(path, columns, rows):
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer)
    writer.writerow(columns)
    writer.writerows(rows)
    _atomic(path, lambda stream: stream.write(buffer.getvalue().encode("utf-8")))


def _load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _source(notebook_path):
    notebook_path = Path(notebook_path)
    notebook = _load_json(notebook_path)
    class NormalizeIO(ast.NodeTransformer):
        def visit_Assign(self, node):
            if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id in {"resume_from", "output_root"}:
                node.value = ast.Constant(None)
            return self.generic_visit(node)

        def visit_AnnAssign(self, node):
            if isinstance(node.target, ast.Name) and node.target.id in {"resume_from", "output_root"}:
                node.value = ast.Constant(None)
            return self.generic_visit(node)

    code = []
    for cell in notebook["cells"]:
        if cell.get("cell_type") == "code":
            source = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
            tree = NormalizeIO().visit(ast.parse(source))
            code.append(ast.unparse(ast.fix_missing_locations(tree)))
    helpers = {}
    helper_paths = [notebook_path.parent / name for name in (
        "regression_common.py", "regression_execution.py", "regression_spectral.py",
        "regression_bootstrap.py")]
    for path in (*helper_paths, Path(__file__)):
        if path.exists():
            helpers[path.name] = path.read_text(encoding="utf-8")
    return {"notebook_code": code, "helpers": helpers}


def _config(config, X, beta):
    result = json.loads(_json(config))
    result.pop("is_manuscript_configuration", None)
    aliases = {"cut": "SC", "second": "SOAR"}
    result["method"] = aliases.get(result.get("method"), result.get("method"))
    if result["method"] not in METHODS:
        raise ValueError("Unknown regression method")
    for key in ("n", "p", "N", "B"):
        if isinstance(result.get(key), bool) or not isinstance(result.get(key), int) or result[key] < 1:
            raise ValueError(f"{key} must be a positive integer")
    if X.shape != (result["n"], result["p"]) or beta.shape != (result["p"],):
        raise ValueError("Configuration dimensions do not match X and beta")
    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(beta)):
        raise ValueError("X and beta must be finite")
    if not np.isfinite(result.get("sigma", np.nan)) or result["sigma"] < 0:
        raise ValueError("sigma must be a finite nonnegative standard deviation")
    if not 0 < result.get("confidence", 0) < 1:
        raise ValueError("confidence must be strictly between zero and one")
    if not 0 <= result.get("threshold_percentile", -1) <= 100:
        raise ValueError("threshold_percentile must lie in [0, 100]")
    if not isinstance(result.get("method_settings"), dict) or "seed_root" not in result:
        raise ValueError("method_settings and seed_root are required")
    return result


def is_production_config(config, X):
    """Check the manuscript defaults independently of any caller-provided flag."""
    if np.shape(X) != (1000, 1000):
        return False
    required = {"n": 1000, "p": 1000, "N": 1000, "B": 500, "sigma": 2,
                "confidence": 0.95, "seed_root": 20261006, "design_seed": 24,
                "noise_model": "normal", "threshold_grid": [0.001, 0.001, 1.0]}
    if any(config.get(key) != value for key, value in required.items()):
        return False
    method = config.get("method")
    if method not in METHODS or config.get("threshold_percentile") != PERCENTILES[method]:
        return False
    largest = config.get("largest_singular_value")
    if method == "Nesterov":
        if largest is None:
            largest = np.linalg.norm(X, ord=2)
        if not np.isfinite(largest) or largest <= 0:
            return False
    expected = method_settings(method, 1000, largest)
    actual = config.get("method_settings", {})
    return set(actual) == set(expected) and all(
        np.isclose(actual[key], value, rtol=1e-12, atol=0) if isinstance(value, float)
        else actual[key] == value for key, value in expected.items())


def _statistics(config, beta, estimates, sigma, bootstrap_errors):
    error = estimates - beta[:, None]
    l2 = np.linalg.norm(error, axis=0)
    linf = np.max(np.abs(error), axis=0)
    quantiles = np.quantile(bootstrap_errors, config["confidence"], axis=0)
    variance_error = np.abs(sigma**2 - config["sigma"]**2)
    metrics = np.concatenate(((linf <= quantiles).astype(float), l2, linf, variance_error))
    return metrics, quantiles


def _metric_mean(records):
    # Scaling prevents the sum of individually finite, nonnegative errors from overflowing.
    values = np.asarray([record[0] for record in records])
    scales = np.max(values, axis=0)
    scales[scales == 0] = 1
    return np.mean(values / scales, axis=0) * scales


def _record_arrays(config, beta, values):
    shapes = {"estimates": (config["p"], 2), "sigma": (2,), "thresholds": (2,),
              "bootstrap_errors": (config["B"], 2)}
    arrays = {}
    for name, shape in shapes.items():
        arrays[name] = np.asarray(values[name], dtype=float)
        if arrays[name].shape != shape or not np.all(np.isfinite(arrays[name])):
            raise ValueError(f"Invalid {name}: expected finite values with shape {shape}")
        if name != "estimates" and np.any(arrays[name] < 0):
            raise ValueError(f"{name} must be nonnegative")
    parameter = np.asarray(values["selected_parameter"], dtype=float)
    if parameter.ndim != 0 or not np.isfinite(parameter) or parameter < 0:
        raise ValueError("selected_parameter must be a finite nonnegative scalar")
    arrays["selected_parameter"] = parameter
    scan = np.asarray(values["threshold_scan"], dtype=float)
    columns = list(values["threshold_columns"])
    if (scan.ndim != 2 or not scan.shape[0] or scan.shape[1] != len(columns)
            or not all(isinstance(column, str) for column in columns)
            or len(set(columns)) != len(columns) or not np.all(np.isfinite(scan))):
        raise ValueError("Invalid threshold scan or column names")
    arrays["threshold_scan"] = scan
    arrays["threshold_columns"] = np.asarray(columns, dtype=str)
    try:
        with np.errstate(over="raise", invalid="raise"):
            arrays["metrics"], arrays["quantiles"] = _statistics(
                config, beta, arrays["estimates"], arrays["sigma"], arrays["bootstrap_errors"])
            arrays["variance"] = arrays["sigma"]**2
    except (FloatingPointError, OverflowError) as error:
        raise ValueError("Derived replication statistics must be finite") from error
    if any(not np.all(np.isfinite(arrays[name])) for name in ("metrics", "quantiles", "variance")):
        raise ValueError("Derived replication statistics must be finite")
    return arrays


def _read_records(directory, config, beta, identity):
    records = []
    for expected, path in enumerate(sorted((directory / "replications").glob("rep_*.npz"))):
        if path.name != f"rep_{expected:06d}.npz" or expected >= config["N"]:
            raise ValueError("Replication files are not a contiguous valid prefix")
        with np.load(path, allow_pickle=False) as saved:
            if str(saved["identity"].item()) != identity or int(saved["index"].item()) != expected:
                raise ValueError("Replication identity/index mismatch")
            arrays = _record_arrays(config, beta, saved)
            diagnostics = json.loads(str(saved["diagnostics"].item()))
            if not isinstance(diagnostics, dict):
                raise ValueError("Invalid replication diagnostics")
            _json(diagnostics)
            records.append((arrays["metrics"], float(arrays["selected_parameter"]),
                            arrays["thresholds"], arrays["quantiles"], diagnostics))
    return records


class BootstrapRun:
    """Commit each completed outer replication before updating derived summaries."""

    def __init__(self, config, X, beta, notebook_path,
                 output_root="outputs/regression-bootstrap-v2", resume_from=None):
        self.X = np.array(X, dtype=float, copy=True)
        self.beta = np.array(beta, dtype=float, copy=True)
        self.config = _config(config, self.X, self.beta)
        source = _source(notebook_path)
        fingerprints = {"config": _digest(self.config), "source": _digest(source),
                        "X": _array_hash(self.X), "beta": _array_hash(self.beta)}
        self.identity = _digest(fingerprints)
        self._records = []
        self._closed = False
        self._entered = False
        self._lock = None
        if resume_from is None:
            root = Path(output_root)
            root.mkdir(parents=True, exist_ok=True)
            prefix = self.config["method"] + "-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-"
            self.directory = Path(tempfile.mkdtemp(prefix=prefix, dir=root)).resolve()
            (self.directory / "replications").mkdir()
            self.metadata = {"version": VERSION, "config": self.config, "fingerprints": fingerprints,
                             "identity": self.identity, "created_at": _now(), "status": "created",
                             "is_manuscript_configuration": is_production_config(self.config, self.X),
                             "runtime": _runtime(),
                             "notebook_path": str(Path(notebook_path).resolve())}
            _write_json(self.directory / "source.json", source)
            _atomic(self.directory / "data.npz", lambda stream: np.savez_compressed(stream, X=self.X, beta=self.beta))
            self._persist("created")
        else:
            self.directory = Path(resume_from).resolve()
            self.metadata = _load_json(self.directory / "metadata.json")
            if self.metadata.get("version") != VERSION or self.metadata.get("fingerprints") != fingerprints:
                raise ValueError("Resume requires identical version, configuration, source code, and data")
            if self.metadata.get("identity") != self.identity or self.metadata.get("runtime") != _runtime():
                raise ValueError("Resume requires identical run identity and runtime versions")
            self._check_archived_data()
            if _digest(_load_json(self.directory / "source.json")) != fingerprints["source"]:
                raise ValueError("Archived source fingerprint mismatch")
            self._records = _read_records(self.directory, self.config, self.beta, self.identity)

    def _check_archived_data(self):
        with np.load(self.directory / "data.npz", allow_pickle=False) as data:
            if any(_array_hash(data[name]) != self.metadata["fingerprints"][name] for name in ("X", "beta")):
                raise ValueError("Archived data fingerprint mismatch")

    @property
    def start_index(self):
        return len(self._records)

    @property
    def summary(self):
        if not self._records:
            return dict.fromkeys(METRICS, None)
        means = _metric_mean(self._records)
        return dict(zip(METRICS, map(float, means)))

    def _persist(self, status, error=None):
        self.metadata.update(status=status, completed_replications=self.start_index,
                             updated_at=_now(), summary=self.summary)
        if error is not None:
            self.metadata["error"] = error
        elif status in ("running", "complete"):
            self.metadata.pop("error", None)
        checkpoint = {"version": VERSION, "identity": self.identity, "status": status,
                      "completed_replications": self.start_index, "target_replications": self.config["N"],
                      "summary": self.summary}
        _write_json(self.directory / "summary.json", checkpoint)
        _write_csv(self.directory / "summary.csv", ["metric", "value", "status", "completed_replications", "target_replications"],
                   [[name, value, status, self.start_index, self.config["N"]] for name, value in self.summary.items()])
        columns = ["index", "selected_parameter", "threshold_original", "threshold_debiased",
                   "quantile_original", "quantile_debiased", *METRICS, "diagnostics"]
        rows = [[i, record[1], *record[2], *record[3], *record[0], _json(record[4])]
                for i, record in enumerate(self._records)]
        _write_csv(self.directory / "replications.csv", columns, rows)
        _write_json(self.directory / "checkpoint.json", checkpoint)
        _write_json(self.directory / "metadata.json", self.metadata)

    def __enter__(self):
        if self._entered or self._closed:
            raise RuntimeError("A BootstrapRun context cannot be entered twice")
        lock = self.directory / ".active.lock"
        try:
            descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError as error:
            raise RuntimeError("Run is locked; another process may be using it. A stale lock after a hard crash must be removed explicitly.") from error
        self._lock = lock
        try:
            with os.fdopen(descriptor, "w") as stream:
                stream.write(_json({"pid": os.getpid(), "created_at": _now()}))
            self._records = _read_records(self.directory, self.config, self.beta, self.identity)
            self._persist("running")
            self._entered = True
            return self
        except BaseException:
            lock.unlink(missing_ok=True)
            self._lock = None
            raise

    def record(self, x, *, selected_parameter, thresholds, estimates, sigma,
               bootstrap_errors, threshold_scan, threshold_columns, diagnostics):
        if not self._entered or self._closed:
            raise RuntimeError("record must be called inside the active context")
        if isinstance(x, bool) or not isinstance(x, (int, np.integer)) or x != self.start_index or x >= self.config["N"]:
            raise ValueError("record index must be the next uncommitted replication")
        if not isinstance(diagnostics, dict):
            raise ValueError("diagnostics must be a dictionary")
        arrays = _record_arrays(self.config, self.beta, locals())
        arrays.update(index=np.asarray(x), identity=np.asarray(self.identity),
                      diagnostics=np.asarray(_json(diagnostics)), metric_names=np.asarray(METRICS))
        path = self.directory / "replications" / f"rep_{x:06d}.npz"
        if path.exists():
            raise FileExistsError("Never overwrite a committed replication")
        _atomic(path, lambda stream: np.savez_compressed(stream, **arrays))
        self._records.append((arrays["metrics"], float(arrays["selected_parameter"]), arrays["thresholds"],
                              arrays["quantiles"], diagnostics))
        self._persist("running")

    def __exit__(self, exc_type, exc, traceback):
        try:
            self._records = _read_records(self.directory, self.config, self.beta, self.identity)
            status = "failed" if exc_type else "complete" if self.start_index == self.config["N"] else "incomplete"
            self._persist(status, None if exc_type is None else f"{exc_type.__name__}: {exc}")
        finally:
            self._closed = True
            if self._lock is not None:
                self._lock.unlink(missing_ok=True)
                self._lock = None
        return False


def export_regression_table(run_dirs, destination, allow_nonproduction=False):
    """Export all seven completed methods, recomputing statistics from raw records."""
    loaded = {}
    provenance = {}
    shared = None
    for directory in map(Path, run_dirs):
        metadata = _load_json(directory / "metadata.json")
        if metadata.get("version") != VERSION or metadata.get("status") != "complete":
            raise ValueError("Only complete current-version runs can be exported")
        if (directory / ".active.lock").exists():
            raise ValueError("Cannot export an active run")
        with np.load(directory / "data.npz", allow_pickle=False) as data:
            X, beta = data["X"], data["beta"]
        config = _config(metadata["config"], X, beta)
        fingerprints = metadata["fingerprints"]
        source = _load_json(directory / "source.json")
        if (fingerprints.get("config") != _digest(config) or fingerprints.get("X") != _array_hash(X)
                or fingerprints.get("beta") != _array_hash(beta)
                or fingerprints.get("source") != _digest(source)
                or metadata.get("identity") != _digest(fingerprints)):
            raise ValueError("Run fingerprints do not match the saved artifacts")
        method = config["method"]
        if method in loaded:
            raise ValueError("Duplicate method in table export")
        common = {key: config.get(key) for key in ("n", "p", "N", "B", "sigma", "confidence",
                  "seed_root", "design_seed", "noise_model", "threshold_grid", "execution", "backend_info",
                  "preprocess_blas_threads")}
        common.update(X=fingerprints["X"], beta=fingerprints["beta"], helpers=_digest(source["helpers"]),
                      runtime=metadata["runtime"])
        if shared is not None and common != shared:
            raise ValueError("Methods must share the design, truth, noise, simulation settings, helper code, and runtime")
        shared = common
        production_config = config.copy()
        if not allow_nonproduction and method == "Nesterov":
            production_config["largest_singular_value"] = float(np.linalg.norm(X, ord=2))
        if not allow_nonproduction and not is_production_config(production_config, X):
            raise ValueError("Nonproduction configurations require allow_nonproduction=True")
        records = _read_records(directory, config, beta, metadata["identity"])
        if len(records) != config["N"] or metadata.get("completed_replications") != config["N"]:
            raise ValueError("Incomplete replication set")
        loaded[method] = (config, _metric_mean(records))
        provenance[method] = {"directory": str(directory.resolve()), "identity": metadata["identity"],
                              "config": config, "runtime": metadata["runtime"], "fingerprints": fingerprints}
    if set(loaded) != set(METHODS):
        raise ValueError("Exactly one completed run for each of the seven methods is required")
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    csv_path = destination / "regression_table.csv"
    tex_path = destination / "regression_table.tex"
    manifest_path = destination / "regression_table_manifest.json"
    if any(path.exists() for path in (csv_path, tex_path, manifest_path)):
        raise FileExistsError("Table export never overwrites existing files")
    rows = [["threshold_percentile", *[loaded[method][0]["threshold_percentile"] for method in METHODS]]]
    rows.extend([[name, *[float(loaded[method][1][j]) for method in METHODS]] for j, name in enumerate(METRICS)])
    _write_csv(csv_path, ["metric", *METHODS], rows)
    labels = {
        "threshold_percentile": r"$b_n^\circ$",
        "coverage_original": "Coverage I", "coverage_debiased": "Coverage II",
        "l2_original": r"$\overline{\|\hat{\boldsymbol{\theta}}-\boldsymbol{\beta}\|_2}$",
        "l2_debiased": r"$\overline{\|\tilde{\boldsymbol{\theta}}-\boldsymbol{\beta}\|_2}$",
        "linf_original": r"$\overline{\|\hat{\boldsymbol{\theta}}-\boldsymbol{\beta}\|_\infty}$",
        "linf_debiased": r"$\overline{\|\tilde{\boldsymbol{\theta}}-\boldsymbol{\beta}\|_\infty}$",
        "variance_absolute_error_original": r"$\overline{|\widehat{\sigma}^2-\sigma^2|}$",
        "variance_absolute_error_debiased": r"$\overline{|\widetilde{\sigma}^2-\sigma^2|}$",
    }
    headers = [r"$\mathrm{AR}^{\kappa}$" if method == "ARk" else method for method in METHODS]
    lines = ["% Generated from complete versioned bootstrap records; CSV coverage is a probability, TeX coverage is a percentage.",
             "\\begin{tabular}{lrrrrrrr}", "\\hline", "Metric & " + " & ".join(headers) + r" \\"]
    for row in rows:
        values = [f"{100*value:.1f}\\%" if row[0].startswith("coverage_") else f"{value:.6g}" for value in row[1:]]
        lines.append(labels[row[0]] + " & " + " & ".join(values) + r" \\")
    lines.extend(["\\hline", "\\end{tabular}"])
    _atomic(tex_path, lambda stream: stream.write(("\n".join(lines) + "\n").encode("utf-8")))
    manifest = {"version": VERSION, "created_at": _now(), "non_manuscript_export": bool(allow_nonproduction),
                "common_settings": shared, "outputs": {"csv": csv_path.name, "tex": tex_path.name},
                "runs": {method: provenance[method] for method in METHODS}}
    _write_json(manifest_path, manifest)
    return {"csv": csv_path, "tex": tex_path, "manifest": manifest_path}
