"""Checkpointed, independent threshold calibration for nine regression methods."""

import argparse
import ast
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import traceback

import numpy as np
from threadpoolctl import threadpool_limits

from regression_results import METRICS, _array_hash, _atomic, _digest, _json, _runtime, _statistics, _write_json


ROOT = Path(__file__).resolve().parent
VERSION = "regression-threshold-scan-v1"
METHODS = ("Ridge", "SC", "Landweber", "Showalter", "ARk", "SOAR", "Nesterov", "HBF", "FAR")
SEEDS = {"baseline": 20261006, "calibration": 20261007, "refinement": 20261007, "validation": 20261008}


def now():
    return datetime.now(timezone.utc).isoformat()


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def source_snapshot():
    """Archive the actual local Python import closure, including conditional imports."""
    pending, sources = [Path(__file__).resolve()], {}
    while pending:
        path = pending.pop()
        if path.name in sources:
            continue
        text = path.read_text(encoding="utf-8")
        sources[path.name] = text
        for node in ast.walk(ast.parse(text)):
            names = ([node.module] if isinstance(node, ast.ImportFrom) and node.module
                     else [item.name for item in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                candidate = ROOT / (name.split(".")[0] + ".py")
                if candidate.is_file() and candidate.name not in sources:
                    pending.append(candidate)
    return {"files": dict(sorted(sources.items()))}


def validate_plan(plan, method, percentiles, seed_root):
    if plan.get("version") != "regression-threshold-selection-v1" or plan.get("status") != "ready_for_validation":
        raise ValueError("Validation requires a frozen, ready selection plan")
    if seed_root != plan.get("validation_seed_root") or seed_root == plan.get("calibration_seed_root"):
        raise ValueError("Validation must use the plan's independent validation seed")
    if list(percentiles) != plan.get("validation_percentiles", {}).get(method):
        raise ValueError("Validation percentiles must exactly match the frozen method plan")


def _check_config(config):
    if config.get("stage") not in SEEDS or config.get("method") not in METHODS:
        raise ValueError("Unknown threshold-scan stage or method")
    for key in ("n", "p", "N", "B"):
        if isinstance(config.get(key), bool) or not isinstance(config.get(key), int) or config[key] < 1:
            raise ValueError(f"Invalid positive integer: {key}")
    percentiles = config.get("percentiles", [])
    if (not percentiles or percentiles != sorted(set(percentiles))
            or not all(np.isfinite(q) and 0 <= q <= 100 for q in percentiles)):
        raise ValueError("Percentiles must be a sorted, unique, nonempty list in [0,100]")
    if config.get("sigma") != 2.0 or config.get("confidence") != 0.95 or config.get("noise_model") != "normal":
        raise ValueError("Unexpected statistical model")
    if config.get("threshold_grid") != [0.001, 0.001, 1.0] or config.get("preprocess_blas_threads") != 1:
        raise ValueError("Unexpected threshold grid or preprocessing settings")
    if (isinstance(config.get("seed_root"), bool) or not isinstance(config.get("seed_root"), int)
            or config["seed_root"] < 0 or not isinstance(config.get("method_settings"), dict)):
        raise ValueError("Invalid seed or method settings")


def _record(config, beta, values):
    q, p, b = len(config["percentiles"]), config["p"], config["B"]
    shapes = {"estimates": (q, p, 2), "sigma": (q, 2), "thresholds": (q, 2),
              "bootstrap_errors": (q, b, 2), "threshold_scan": (1000, 3), "percentiles": (q,)}
    result = {}
    for name, shape in shapes.items():
        result[name] = np.asarray(values[name], dtype=float)
        if result[name].shape != shape or not np.isfinite(result[name]).all():
            raise ValueError(f"Invalid {name}: expected finite {shape}")
        if name != "estimates" and np.any(result[name] < 0):
            raise ValueError(f"Negative {name}")
    if not np.array_equal(result["percentiles"], config["percentiles"]):
        raise ValueError("Replication percentile order does not match the configuration")
    if not np.array_equal(result["threshold_scan"][:, 0], 0.001 * np.arange(1, 1001)):
        raise ValueError("Threshold scan grid changed")
    parameter = np.asarray(values["selected_parameter"], dtype=float)
    if parameter.ndim != 0 or not np.isfinite(parameter) or parameter < 0:
        raise ValueError("Invalid selected parameter")
    result["selected_parameter"] = parameter
    diagnostics = values["diagnostics"]
    if not isinstance(diagnostics, dict):
        diagnostics = json.loads(str(np.asarray(diagnostics).item()))
    if not isinstance(diagnostics, dict):
        raise ValueError("Diagnostics must be an object")
    result["diagnostics"] = np.asarray(_json(diagnostics))
    with np.errstate(over="raise", invalid="raise"):
        pairs = [_statistics(config, beta, result["estimates"][j], result["sigma"][j],
                             result["bootstrap_errors"][j]) for j in range(q)]
        result["metrics"] = np.stack([pair[0] for pair in pairs])
        result["quantiles"] = np.stack([pair[1] for pair in pairs])
        result["variance"] = result["sigma"] ** 2
    if any(not np.isfinite(result[name]).all() for name in ("metrics", "quantiles", "variance")):
        raise ValueError("Nonfinite derived statistics")
    return result


def validate_run(directory, require_complete=True):
    """Validate immutable run metadata, source, data and any frozen selection plan."""
    directory = Path(directory)
    meta = load_json(directory / "metadata.json")
    if meta.get("version") != VERSION:
        raise ValueError("Unsupported threshold-scan version")
    if require_complete and (meta.get("status") != "complete" or (directory / ".active.lock").exists()):
        raise ValueError("Only complete, inactive scan runs may be exported")
    config = meta["config"]
    _check_config(config)
    source = load_json(directory / "source.json")
    fingerprints = meta["fingerprints"]
    if (_digest(config) != fingerprints.get("config") or _digest(source) != fingerprints.get("source")
            or _digest(fingerprints) != meta.get("identity")):
        raise ValueError("Configuration/source fingerprint mismatch")
    with np.load(directory / "data.npz", allow_pickle=False) as data:
        for name, shape in (("X", (config["n"], config["p"])), ("beta", (config["p"],))):
            if (data[name].shape != shape or not np.isfinite(data[name]).all()
                    or _array_hash(data[name]) != fingerprints.get(name)):
                raise ValueError(f"Saved {name} fingerprint/shape mismatch")
    if config["stage"] == "validation":
        plan = load_json(directory / "selection_plan.json")
        if _digest(plan) != config.get("selection_plan_sha256"):
            raise ValueError("Selection plan fingerprint mismatch")
        validate_plan(plan, config["method"], config["percentiles"], config["seed_root"])
    formal = (config["n"], config["p"], config["N"], config["B"]) == (1000, 1000, 1000, 500)
    if meta.get("is_manuscript_configuration") != formal:
        raise ValueError("Incorrect formal-size marker")
    return meta


def read_scan_records(directory, metadata=None):
    """Stream validated raw records without keeping all bootstrap draws in memory."""
    directory = Path(directory)
    meta = validate_run(directory, require_complete=False) if metadata is None else metadata
    config = meta["config"]
    with np.load(directory / "data.npz", allow_pickle=False) as data:
        beta = data["beta"].copy()
    files = sorted((directory / "replications").glob("*.npz"))
    if len(files) > config["N"]:
        raise ValueError("Too many replication files")
    if meta["status"] == "complete" and (len(files) != config["N"] or meta.get("completed_replications") != len(files)):
        raise ValueError("Complete run has missing records")
    for index, path in enumerate(files):
        if path.name != f"rep_{index:06d}.npz":
            raise ValueError("Replication files must be a contiguous prefix")
        with np.load(path, allow_pickle=False) as saved:
            if int(saved["index"].item()) != index or str(saved["identity"].item()) != meta["identity"]:
                raise ValueError("Replication identity/index mismatch")
            result = _record(config, beta, saved)
            for name in ("metrics", "quantiles", "variance"):
                if saved[name].shape != result[name].shape or not np.allclose(saved[name], result[name], rtol=1e-12, atol=1e-12):
                    raise ValueError(f"Corrupt derived {name}")
            result.update(index=np.asarray(index), identity=np.asarray(meta["identity"]))
            yield result


class ScanRun:
    def __init__(self, config, X, beta, source, output_root, resume_from=None, selection_plan=None, data_origin=None):
        _check_config(config)
        if (np.shape(X) != (config["n"], config["p"]) or np.shape(beta) != (config["p"],)
                or not np.isfinite(X).all() or not np.isfinite(beta).all()):
            raise ValueError("Input data does not match the configuration")
        self.config, self.X, self.beta = config, X, beta
        fingerprints = {"config": _digest(config), "source": _digest(source),
                        "X": _array_hash(X), "beta": _array_hash(beta)}
        self.identity = _digest(fingerprints)
        self.records, self.lock = [], None
        if resume_from is None:
            root = Path(output_root)
            root.mkdir(parents=True, exist_ok=True)
            prefix = f"{config['stage']}-{config['method']}-"
            self.directory = Path(tempfile.mkdtemp(prefix=prefix, dir=root)).resolve()
            (self.directory / "replications").mkdir()
            self.metadata = {"version": VERSION, "config": config, "fingerprints": fingerprints,
                             "identity": self.identity, "created_at": now(), "runtime": _runtime(),
                             "data_origin": data_origin,
                             "is_manuscript_configuration": (config["n"], config["p"], config["N"], config["B"]) == (1000, 1000, 1000, 500)}
            _write_json(self.directory / "source.json", source)
            _atomic(self.directory / "data.npz", lambda stream: np.savez_compressed(stream, X=X, beta=beta))
            if selection_plan is not None:
                _write_json(self.directory / "selection_plan.json", selection_plan)
            self.persist("created")
        else:
            self.directory = Path(resume_from).resolve()
            self.metadata = validate_run(self.directory, require_complete=False)
            if self.metadata["fingerprints"] != fingerprints or self.metadata["runtime"] != _runtime():
                raise ValueError("Resume requires identical configuration, data, source, and runtime")

    @property
    def start_index(self):
        return len(self.records)

    def summary(self):
        rows = []
        if self.records:
            metrics = np.stack([record[0] for record in self.records])
            radii = np.stack([record[1] for record in self.records])
            means, radius_means = metrics.mean(axis=0), radii.mean(axis=0)
            for j, percentile in enumerate(self.config["percentiles"]):
                row = dict(zip(METRICS, map(float, means[j])))
                row.update(percentile=percentile, mean_radius_original=float(radius_means[j, 0]),
                           mean_radius_debiased=float(radius_means[j, 1]))
                for branch, index in (("original", 0), ("debiased", 1)):
                    rate = means[j, index]
                    row["mcse_coverage_" + branch] = float(np.sqrt(rate * (1 - rate) / self.start_index))
                rows.append(row)
        return rows

    def persist(self, status, error=None):
        self.metadata.update(status=status, completed_replications=self.start_index,
                             target_replications=self.config["N"], updated_at=now())
        if error is None:
            self.metadata.pop("error", None)
        else:
            self.metadata["error"] = error
        summary = {"version": VERSION, "status": status, "config": self.config, "identity": self.identity,
                   "completed_replications": self.start_index, "target_replications": self.config["N"],
                   "results": self.summary()}
        _write_json(self.directory / "summary.json", summary)
        _write_json(self.directory / "checkpoint.json", {key: value for key, value in summary.items() if key != "results"})
        _write_json(self.directory / "metadata.json", self.metadata)

    def __enter__(self):
        lock = self.directory / ".active.lock"
        try:
            descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError as error:
            raise RuntimeError("Run is active or has a stale crash lock; never remove a live process lock") from error
        self.lock = lock
        try:
            with os.fdopen(descriptor, "w") as stream:
                stream.write(_json({"pid": os.getpid(), "created_at": now()}))
            self.records = [(record["metrics"], record["quantiles"])
                            for record in read_scan_records(self.directory, self.metadata)]
            self.persist("running")
            return self
        except BaseException:
            lock.unlink(missing_ok=True)
            self.lock = None
            raise

    def record(self, index, **values):
        if self.lock is None or index != self.start_index or index >= self.config["N"]:
            raise ValueError("Record must be the next index inside the active context")
        arrays = _record(self.config, self.beta, values)
        arrays.update(index=np.asarray(index), identity=np.asarray(self.identity))
        path = self.directory / "replications" / f"rep_{index:06d}.npz"
        if path.exists():
            raise FileExistsError("Never overwrite committed records")
        _atomic(path, lambda stream: np.savez_compressed(stream, **arrays))
        self.records.append((arrays["metrics"], arrays["quantiles"]))
        self.persist("running")

    def __exit__(self, exc_type, exc, tb):
        try:
            self.records = [(record["metrics"], record["quantiles"])
                            for record in read_scan_records(self.directory, self.metadata)]
            self.persist("failed" if exc_type else "complete" if self.start_index == self.config["N"] else "incomplete",
                         None if exc_type is None else f"{exc_type.__name__}: {exc}")
        finally:
            if self.lock is not None:
                self.lock.unlink(missing_ok=True)
                self.lock = None
        return False


def positive(value):
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return result


def ridge_grid_stride(alpha_step, settings):
    """Map the requested alpha spacing to an integer grid with unchanged endpoints."""
    scaled = alpha_step * settings["alpha_divisor"]
    if not np.isfinite(scaled) or scaled <= 0:
        raise ValueError("Ridge alpha step must map to a finite positive stride")
    stride = int(round(scaled))
    tolerance = 8 * np.finfo(float).eps * max(1.0, abs(scaled))
    if stride < 1 or abs(scaled - stride) > tolerance:
        raise ValueError("Ridge alpha step times alpha_divisor must be an integer")
    if settings["C"] % stride:
        raise ValueError("Ridge alpha step must exactly divide the original grid maximum")
    return stride


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=SEEDS, required=True)
    parser.add_argument("--methods", nargs="+", choices=METHODS)
    parser.add_argument("--percentiles", nargs="+", type=float)
    parser.add_argument("--N", type=positive, default=1000)
    parser.add_argument("--B", type=positive, default=500)
    parser.add_argument("--seed-root", type=int)
    parser.add_argument("--data-from", type=Path, required=True)
    parser.add_argument("--selection-plan", type=Path)
    parser.add_argument("--output-root", type=Path, default=ROOT / "outputs" / VERSION)
    parser.add_argument("--resume-from", type=Path)
    parser.add_argument("--smoke-test", action="store_true", help="Allow nonproduction input dimensions")
    parser.add_argument("--backend", choices=("cpu", "gpu"), default="gpu")
    parser.add_argument("--gpu-devices", default="0")
    parser.add_argument("--cpu-workers", type=positive, default=4)
    parser.add_argument("--blas-threads", type=positive, default=8)
    parser.add_argument("--batch-size", type=positive, default=64)
    parser.add_argument("--ridge-alpha-step", type=float,
                        help="Override the alpha spacing, retaining grid endpoints, for a Ridge-only sensitivity run")
    parser.add_argument("--nesterov-tau", type=float,
                        help="Override the original-sample TDP multiplier for a Nesterov-only sensitivity run")
    parser.add_argument("--hbf-tau", type=float,
                        help="Override the original-sample TDP multiplier for an HBF-only sensitivity run")
    parser.add_argument("--hbf-eta", type=float,
                        help="Override the positive damping coefficient for an HBF-only sensitivity run")
    parser.add_argument("--far-theta", type=float,
                        help="Override the order in (1,2) for an FAR-only sensitivity run")
    parser.add_argument("--far-tau", type=float,
                        help="Override the original-sample TDP multiplier for an FAR-only sensitivity run")
    parser.add_argument("--soar-tau", type=float,
                        help="Override the original-sample TDP multiplier for a SOAR-only sensitivity run")
    parser.add_argument("--landweber-selection", choices=("AOSR", "TDP"),
                        help="Override the original-sample stopping rule for a Landweber-only sensitivity run")
    parser.add_argument("--landweber-tau", type=float,
                        help="Landweber TDP multiplier; defaults to 0.9 when TDP is selected")
    parser.add_argument("--_worker-method", choices=METHODS, help=argparse.SUPPRESS)
    parser.add_argument("--_batch-directory", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    args.methods = args.methods or (["HBF", "FAR"] if args.stage == "baseline" else list(METHODS))
    if len(set(args.methods)) != len(args.methods):
        parser.error("Methods must be unique")
    if args.ridge_alpha_step is not None:
        if not np.isfinite(args.ridge_alpha_step) or args.ridge_alpha_step <= 0:
            parser.error("Ridge alpha step must be finite and positive")
        if args.methods != ["Ridge"] or args._worker_method not in (None, "Ridge"):
            parser.error("--ridge-alpha-step requires a Ridge-only run")
    if args.nesterov_tau is not None:
        if not np.isfinite(args.nesterov_tau) or args.nesterov_tau <= 0:
            parser.error("Nesterov tau must be finite and positive")
        if args.methods != ["Nesterov"] or args._worker_method not in (None, "Nesterov"):
            parser.error("--nesterov-tau requires a Nesterov-only run")
    for method, name in (("HBF", "hbf_tau"), ("SOAR", "soar_tau")):
        value = getattr(args, name)
        if value is not None:
            if not np.isfinite(value) or value <= 0:
                parser.error(f"{method} tau must be finite and positive")
            if args.methods != [method] or args._worker_method not in (None, method):
                parser.error(f"--{name.replace('_', '-')} requires a {method}-only run")
    if args.hbf_eta is not None:
        if not np.isfinite(args.hbf_eta) or args.hbf_eta <= 0:
            parser.error("HBF eta must be finite and positive")
        if args.methods != ["HBF"] or args._worker_method not in (None, "HBF"):
            parser.error("--hbf-eta requires an HBF-only run")
    if args.far_theta is not None or args.far_tau is not None:
        if args.methods != ["FAR"] or args._worker_method not in (None, "FAR"):
            parser.error("FAR overrides require an FAR-only run")
        if args.far_theta is not None and (not np.isfinite(args.far_theta) or not 1 < args.far_theta < 2):
            parser.error("FAR theta must be finite and lie in (1,2)")
        if args.far_tau is not None and (not np.isfinite(args.far_tau) or args.far_tau <= 0):
            parser.error("FAR tau must be finite and positive")
    if args.landweber_selection is not None or args.landweber_tau is not None:
        if args.methods != ["Landweber"] or args._worker_method not in (None, "Landweber"):
            parser.error("Landweber overrides require a Landweber-only run")
        if args.landweber_tau is not None:
            if not np.isfinite(args.landweber_tau) or args.landweber_tau <= 0:
                parser.error("Landweber tau must be finite and positive")
            if args.landweber_selection != "TDP":
                parser.error("--landweber-tau requires --landweber-selection TDP")
    if args.resume_from is not None and len(args.methods) != 1:
        parser.error("Resume requires exactly one method")
    if args.stage == "validation" and (args.percentiles is None or args.selection_plan is None):
        parser.error("Validation requires explicit percentiles and a frozen selection plan")
    if args.stage != "validation" and args.selection_plan is not None:
        parser.error("Selection plans are only used by the independent validation stage")
    if args.percentiles is None:
        args.percentiles = ([50.0] if args.stage == "baseline" else list(map(float, range(5, 100, 10)))
                            if args.stage == "refinement" else list(map(float, range(0, 101, 10))))
    if (not args.percentiles or len(args.percentiles) != len(set(args.percentiles))
            or not all(np.isfinite(q) and 0 <= q <= 100 for q in args.percentiles)):
        parser.error("Percentiles must be unique values in [0,100]")
    args.percentiles = sorted(args.percentiles)
    args.seed_root = SEEDS[args.stage] if args.seed_root is None else args.seed_root
    if args.seed_root < 0:
        parser.error("Seed must be nonnegative")
    try:
        args.gpu_devices = tuple(map(int, args.gpu_devices.split(",")))
        if not args.gpu_devices or min(args.gpu_devices) < 0 or len(set(args.gpu_devices)) != len(args.gpu_devices):
            raise ValueError
    except ValueError:
        parser.error("GPU devices must be unique nonnegative integers")
    for key in ("data_from", "output_root", "resume_from", "selection_plan", "_batch_directory"):
        value = getattr(args, key)
        if value is not None:
            setattr(args, key, value.resolve())
    if bool(args._worker_method) != bool(args._batch_directory):
        parser.error("Internal worker arguments must be provided together")
    return args


def run_worker(args):
    from regression_bootstrap import validate_backend
    from regression_spectral import SpectralDesign
    from regression_threshold_scan import scan_method_settings, threshold_grid_pair, bootstrap_errors_grid

    method = args._worker_method
    method_report = args._batch_directory / f"{method}.json"
    report = {"method": method, "stage": args.stage, "status": "starting", "pid": os.getpid(), "started_at": now()}
    _write_json(method_report, report)
    try:
        parent = load_json(args.data_from / "metadata.json")
        if parent.get("status") != "complete":
            raise ValueError("Input design must come from a completed archived run")
        with np.load(args.data_from / "data.npz", allow_pickle=False) as data:
            X, beta = data["X"].copy(), data["beta"].copy()
        if _array_hash(X) != parent["fingerprints"]["X"] or _array_hash(beta) != parent["fingerprints"]["beta"]:
            raise ValueError("Input design/truth archive fingerprint mismatch")
        if not args.smoke_test and X.shape != (1000, 1000):
            raise ValueError("Nonproduction dimensions require --smoke-test")
        backend_info = validate_backend(args.backend, gpu_devices=args.gpu_devices)
        with threadpool_limits(limits=1):
            design = SpectralDesign(X)
        settings = scan_method_settings(method, X.shape[0], float(design.s[0]))
        if method == "Ridge" and args.ridge_alpha_step is not None:
            settings["m"] = ridge_grid_stride(args.ridge_alpha_step, settings)
        if args.nesterov_tau is not None:
            settings["tau"] = args.nesterov_tau
        for selected_method, name in (("HBF", "hbf_tau"), ("SOAR", "soar_tau")):
            if method == selected_method and getattr(args, name) is not None:
                settings["tau"] = getattr(args, name)
        if method == "HBF" and args.hbf_eta is not None:
            settings["eta"] = args.hbf_eta
        if method == "FAR":
            if args.far_theta is not None:
                settings["theta"] = args.far_theta
            if args.far_tau is not None:
                settings["tau"] = args.far_tau
        if args.landweber_selection is not None:
            settings["selection"] = args.landweber_selection
            if args.landweber_selection == "TDP":
                settings.pop("minimum_steps", None)
                settings["tau"] = 0.9 if args.landweber_tau is None else args.landweber_tau
        config = {"stage": args.stage, "method": method, "percentiles": args.percentiles,
                  "n": X.shape[0], "p": X.shape[1], "N": args.N, "B": args.B,
                  "sigma": 2.0, "confidence": 0.95, "seed_root": args.seed_root,
                  "design_seed": parent["config"]["design_seed"], "noise_model": "normal",
                  "threshold_grid": [0.001, 0.001, 1.0], "preprocess_blas_threads": 1,
                  "method_settings": settings, "largest_singular_value": float(design.s[0]),
                  "execution": {"backend": args.backend, "gpu_devices": list(args.gpu_devices),
                                "cpu_workers": args.cpu_workers, "blas_threads": args.blas_threads,
                                "batch_size": args.batch_size}, "backend_info": backend_info}
        plan = None
        if args.selection_plan is not None:
            plan = load_json(args.selection_plan)
            validate_plan(plan, method, args.percentiles, args.seed_root)
            config["selection_plan_sha256"] = _digest(plan)
        source = source_snapshot()
        origin = {"directory": str(args.data_from), "identity": parent.get("identity"),
                  "metadata_sha256": _digest(parent)}
        run = ScanRun(config, X, beta, source, args._batch_directory / "runs", args.resume_from, plan, origin)
        report.update(run_directory=str(run.directory), status="running", target_replications=args.N)
        _write_json(method_report, report)
        started = time.monotonic()
        with threadpool_limits(limits=args.blas_threads), run:
            print(f"{method} {args.stage}: {run.start_index}/{args.N}; {run.directory}", flush=True)
            for index in range(run.start_index, args.N):
                rng = np.random.default_rng((args.seed_root, index, 0, 0))
                noise = rng.normal(0, 2.0, X.shape[0])
                response = X @ beta + noise
                selection_settings = settings.copy()
                if "tau" in selection_settings:
                    selection_settings["delta"] = selection_settings["tau"] * float(np.linalg.norm(noise))
                ordinary, debiased, parameter, diagnostics = design.select(response, method, selection_settings, beta)
                centers, thresholds, scan = threshold_grid_pair(np.column_stack((ordinary, debiased)), beta, args.percentiles)
                sigma_cache = {}
                sigma_rows = []
                for center in centers:
                    key = center.tobytes()
                    if key not in sigma_cache:
                        sigma_cache[key] = np.linalg.norm(response[:, None] - X @ center, axis=0) / np.sqrt(X.shape[0])
                    sigma_rows.append(sigma_cache[key])
                sigma = np.stack(sigma_rows)
                errors = bootstrap_errors_grid(design, method, parameter, selection_settings, centers, sigma, thresholds,
                                               args.B, args.seed_root, index, **config["execution"])
                run.record(index, estimates=centers, sigma=sigma, thresholds=thresholds,
                           bootstrap_errors=errors, selected_parameter=parameter, diagnostics=diagnostics,
                           threshold_scan=scan, percentiles=args.percentiles)
                report.update(completed_replications=run.start_index, elapsed_seconds=time.monotonic() - started,
                              updated_at=now())
                _write_json(method_report, report)
                if run.start_index == 1 or run.start_index % 10 == 0 or run.start_index == args.N:
                    print(f"[{now()}] {method}: checkpoint {run.start_index}/{args.N}", flush=True)
        report.update(status="complete", completed_replications=run.start_index, finished_at=now())
        _write_json(method_report, report)
        return 0
    except BaseException as error:
        report.update(status="failed", error=f"{type(error).__name__}: {error}", finished_at=now())
        _write_json(method_report, report)
        raise


def child_arguments(args, method, batch):
    result = [sys.executable, "-u", str(Path(__file__).resolve()), "--stage", args.stage,
              "--methods", method, "--percentiles", *map(str, args.percentiles),
              "--N", str(args.N), "--B", str(args.B), "--seed-root", str(args.seed_root),
              "--data-from", str(args.data_from), "--backend", args.backend,
              "--gpu-devices", ",".join(map(str, args.gpu_devices)), "--cpu-workers", str(args.cpu_workers),
              "--blas-threads", str(args.blas_threads), "--batch-size", str(args.batch_size),
              "--_worker-method", method, "--_batch-directory", str(batch)]
    if args.smoke_test:
        result.append("--smoke-test")
    if args.nesterov_tau is not None:
        result.extend(["--nesterov-tau", str(args.nesterov_tau)])
    for name in ("ridge_alpha_step", "hbf_tau", "hbf_eta", "far_theta", "far_tau", "soar_tau"):
        if getattr(args, name) is not None:
            result.extend(["--" + name.replace("_", "-"), str(getattr(args, name))])
    if args.landweber_selection is not None:
        result.extend(["--landweber-selection", args.landweber_selection])
    if args.landweber_tau is not None:
        result.extend(["--landweber-tau", str(args.landweber_tau)])
    for name in ("resume_from", "selection_plan"):
        if getattr(args, name) is not None:
            result.extend(["--" + name.replace("_", "-"), str(getattr(args, name))])
    return result


def run_batch(args):
    args.output_root.mkdir(parents=True, exist_ok=True)
    batch = Path(tempfile.mkdtemp(prefix=args.stage + "-", dir=args.output_root)).resolve()
    configuration = {key: getattr(args, key) for key in ("stage", "methods", "percentiles", "N", "B", "seed_root",
                                                       "backend", "gpu_devices", "batch_size", "cpu_workers", "blas_threads",
                                                       "ridge_alpha_step", "nesterov_tau", "hbf_tau", "hbf_eta", "far_theta", "far_tau", "soar_tau",
                                                       "landweber_selection", "landweber_tau")}
    manifest = {"version": VERSION, "status": "running", "configuration": configuration,
                "started_at": now(), "pid": os.getpid(), "runs": {}}
    path = batch / "manifest.json"
    _write_json(path, manifest)
    print(f"Threshold batch directory: {batch}", flush=True)
    environment = dict(os.environ, PYTHONUNBUFFERED="1", MPLBACKEND="Agg")
    for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[name] = str(args.blas_threads)
    try:
        for method in args.methods:
            with (batch / f"{method}.log").open("x", encoding="utf-8") as log:
                process = subprocess.Popen(child_arguments(args, method, batch), cwd=ROOT, env=environment,
                                           stdout=log, stderr=subprocess.STDOUT)
                manifest["runs"][method] = {"status": "running", "pid": process.pid, "log": str(batch / f"{method}.log")}
                _write_json(path, manifest)
                try:
                    code = process.wait()
                except BaseException:
                    process.terminate()
                    process.wait()
                    raise
            status_path = batch / f"{method}.json"
            status = load_json(status_path) if status_path.exists() else {}
            manifest["runs"][method].update(status, returncode=code)
            _write_json(path, manifest)
            if code != 0 or status.get("status") != "complete":
                raise RuntimeError(f"{method} failed; see {batch / (method + '.log')}")
            print(f"[{now()}] Complete {method}: {status['run_directory']}", flush=True)
        manifest.update(status="complete", finished_at=now())
        _write_json(path, manifest)
        print(f"Threshold batch complete: {path}", flush=True)
        return 0
    except BaseException as error:
        manifest.update(status="failed", finished_at=now(), error=f"{type(error).__name__}: {error}")
        _write_json(path, manifest)
        raise


def main(argv=None):
    args = parse_args(argv)
    try:
        return run_worker(args) if args._worker_method else run_batch(args)
    except Exception:
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
