"""Run the existing bootstrap notebooks as isolated, checkpointed terminal jobs."""

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import traceback


ROOT = Path(__file__).resolve().parent
VERSION = "regression-terminal-runner-v1"
NOTEBOOKS = {
    "Ridge": "Bootstrapridge.ipynb",
    "SC": "Bootstrapcut.ipynb",
    "Landweber": "Bootstraplandweber.ipynb",
    "Showalter": "Bootstrapshowalter.ipynb",
    "ARk": "BootstrapARk.ipynb",
    "SOAR": "Bootstrapsecond.ipynb",
    "Nesterov": "BootstrapNesterov.ipynb",
}
CONTROLS = frozenset({"n", "p", "N", "B", "backend", "gpu_devices",
                      "cpu_workers", "blas_threads", "batch_size",
                      "output_root", "resume_from"})


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def write_json(path, value):
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix="." + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def prepare_cells(notebook, overrides):
    """Change only simple module-level configuration assignments, exactly once."""
    if not set(overrides) <= CONTROLS:
        raise ValueError("Only execution/output and explicit sample-size controls may be overridden")
    counts = dict.fromkeys(overrides, 0)
    cells = []
    for index, cell in enumerate(notebook["cells"]):
        if cell.get("cell_type") != "code":
            continue
        source = cell["source"]
        tree = ast.parse("".join(source) if isinstance(source, list) else source)
        for node in tree.body:
            if (isinstance(node, ast.Assign) and len(node.targets) == 1
                    and isinstance(node.targets[0], ast.Name)
                    and node.targets[0].id in overrides):
                name = node.targets[0].id
                node.value = ast.parse(repr(overrides[name]), mode="eval").body
                counts[name] += 1
        cells.append((index, ast.fix_missing_locations(tree)))
    if any(count != 1 for count in counts.values()):
        raise ValueError(f"Configuration assignment count changed: {counts}")
    return cells


def cells_source(cells):
    return "\n\n".join(f"# Cell {index}\n{ast.unparse(tree)}" for index, tree in cells) + "\n"


def positive(value):
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return result


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methods", nargs="+", choices=NOTEBOOKS, default=list(NOTEBOOKS))
    parser.add_argument("--backend", choices=("cpu", "gpu"), default="gpu")
    parser.add_argument("--gpu-devices", default="0")
    parser.add_argument("--cpu-workers", type=positive, default=4)
    parser.add_argument("--blas-threads", type=positive, default=8)
    parser.add_argument("--batch-size", type=positive, default=64)
    parser.add_argument("--N", type=positive)
    parser.add_argument("--B", type=positive)
    parser.add_argument("--n", type=positive)
    parser.add_argument("--p", type=positive)
    parser.add_argument("--smoke", action="store_true",
                        help="Explicit nonproduction run; defaults n=p=40, N=2, B=4")
    parser.add_argument("--output-root", type=Path, default=ROOT / "outputs" / "section61-bootstrap")
    parser.add_argument("--resume-from", type=Path, help="Existing BootstrapRun directory; one method only")
    parser.add_argument("--_worker-method", choices=NOTEBOOKS, help=argparse.SUPPRESS)
    parser.add_argument("--_batch-directory", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if len(args.methods) != len(set(args.methods)):
        parser.error("Duplicate methods are not allowed")
    if args.resume_from is not None and len(args.methods) != 1:
        parser.error("--resume-from requires exactly one method")
    if not args.smoke and (args.n is not None or args.p is not None):
        parser.error("Dimension overrides require --smoke")
    for name, formal, small in (("n", 1000, 40), ("p", 1000, 40), ("N", 1000, 2), ("B", 500, 4)):
        if getattr(args, name) is None:
            setattr(args, name, small if args.smoke else formal)
    if args.p < 20:
        parser.error("p must be at least 20 for the notebook's sparse truth")
    if args.smoke and (args.n, args.p, args.N, args.B) == (1000, 1000, 1000, 500):
        parser.error("--smoke must use a nonproduction sample size")
    try:
        args.gpu_devices = tuple(int(value) for value in args.gpu_devices.split(","))
        if not args.gpu_devices or min(args.gpu_devices) < 0 or len(set(args.gpu_devices)) != len(args.gpu_devices):
            raise ValueError
    except ValueError:
        parser.error("--gpu-devices must list unique nonnegative device IDs, for example 0,1")
    args.output_root = args.output_root.resolve()
    if args.resume_from is not None:
        args.resume_from = args.resume_from.resolve()
    if bool(args._worker_method) != bool(args._batch_directory):
        parser.error("Internal worker arguments must be supplied together")
    return args


def run_worker(args):
    from regression_results import BootstrapRun

    method = args._worker_method
    batch = args._batch_directory.resolve()
    status_path = batch / f"{method}.json"
    status = {"method": method, "status": "starting", "pid": os.getpid(), "started_at": now()}
    write_json(status_path, status)
    notebook_path = ROOT / NOTEBOOKS[method]
    notebook_text = notebook_path.read_text(encoding="utf-8")
    notebook = json.loads(notebook_text)
    overrides = {key: getattr(args, key) for key in CONTROLS - {"output_root", "resume_from"}}
    overrides.update(output_root=str(batch / "runs"),
                     resume_from=None if args.resume_from is None else str(args.resume_from))
    cells = prepare_cells(notebook, overrides)
    effective_source = cells_source(cells)
    normalized = dict(overrides, output_root=None, resume_from=None)
    normalized_source = cells_source(prepare_cells(notebook, normalized))
    runner_source = Path(__file__).read_text(encoding="utf-8")
    provenance = {"version": VERSION, "runner_sha256": digest(runner_source),
                  "effective_code_sha256": digest(normalized_source), "smoke": args.smoke}
    started = time.monotonic()

    class LoggedRun(BootstrapRun):
        def __init__(self, config, *positional, **keywords):
            super().__init__(dict(config, terminal_runner=provenance), *positional, **keywords)
            status.update(run_directory=str(self.directory), completed_replications=self.start_index,
                          target_replications=self.config["N"], status="prepared",
                          is_manuscript_configuration=self.metadata["is_manuscript_configuration"])
            write_json(status_path, status)
            archive = self.directory / "terminal_runner"
            if args.resume_from is None:
                archive.mkdir()
                (archive / "runner.py").write_text(runner_source, encoding="utf-8")
                (archive / notebook_path.name).write_text(notebook_text, encoding="utf-8")
                (archive / "effective_code.py").write_text(effective_source, encoding="utf-8")
                (archive / "normalized_code.py").write_text(normalized_source, encoding="utf-8")
                write_json(archive / "invocation.json", {"overrides": overrides, "provenance": provenance,
                                                        "python": sys.executable, "created_at": now()})
            else:
                if (digest((archive / "runner.py").read_text(encoding="utf-8")) != provenance["runner_sha256"]
                        or digest((archive / "normalized_code.py").read_text(encoding="utf-8"))
                        != provenance["effective_code_sha256"]):
                    raise ValueError("Archived terminal runner source does not match")
            print(f"[{now()}] {method}: {self.start_index}/{self.config['N']} committed; {self.directory}", flush=True)

        def record(self, index, **values):
            super().record(index, **values)
            status.update(status="running", completed_replications=self.start_index,
                          updated_at=now(), elapsed_seconds=time.monotonic() - started)
            write_json(status_path, status)
            print(f"[{now()}] {method}: checkpoint {self.start_index}/{self.config['N']} committed", flush=True)

    namespace = {"__name__": "__main__", "__file__": str(notebook_path)}
    previous_directory = Path.cwd()
    try:
        os.chdir(ROOT)
        for index, tree in cells:
            exec(compile(tree, f"{notebook_path}:cell-{index}", "exec"), namespace)
            # Imports in the first cell use the original class; only persistence/logging is wrapped.
            if namespace.get("BootstrapRun") is BootstrapRun:
                namespace["BootstrapRun"] = LoggedRun
        run = namespace.get("run")
        if run is None or run.metadata["status"] != "complete":
            raise RuntimeError("Notebook did not produce a complete BootstrapRun")
        status.update(status="complete", completed_replications=run.start_index,
                      finished_at=now(), elapsed_seconds=time.monotonic() - started,
                      summary=run.summary)
        write_json(status_path, status)
        return 0
    except BaseException as error:
        status.update(status="failed", error=f"{type(error).__name__}: {error}", finished_at=now())
        write_json(status_path, status)
        raise
    finally:
        os.chdir(previous_directory)


def child_arguments(args, method, batch):
    command = [sys.executable, "-u", str(Path(__file__).resolve()), "--methods", method,
               "--_worker-method", method, "--_batch-directory", str(batch),
               "--backend", args.backend, "--gpu-devices", ",".join(map(str, args.gpu_devices)),
               "--N", str(args.N), "--B", str(args.B), "--cpu-workers", str(args.cpu_workers),
               "--blas-threads", str(args.blas_threads), "--batch-size", str(args.batch_size)]
    if args.smoke:
        command += ["--smoke", "--n", str(args.n), "--p", str(args.p)]
    if args.resume_from is not None:
        command += ["--resume-from", str(args.resume_from)]
    return command


def run_batch(args):
    args.output_root.mkdir(parents=True, exist_ok=True)
    prefix = "bootstrap-"
    batch = Path(tempfile.mkdtemp(prefix=prefix + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ-") ,
                                  dir=args.output_root)).resolve()
    settings = {key: getattr(args, key) for key in ("methods", "backend", "gpu_devices", "N", "B", "n", "p",
                                                   "cpu_workers", "blas_threads", "batch_size", "smoke")}
    manifest = {"version": VERSION, "status": "running", "started_at": now(), "pid": os.getpid(),
                "configuration": settings, "runs": {}, "batch_directory": str(batch)}
    write_json(batch / "manifest.json", manifest)
    print(f"Batch directory: {batch}", flush=True)
    environment = dict(os.environ, MPLBACKEND="Agg", PYTHONUNBUFFERED="1")
    for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[variable] = str(args.blas_threads)
    try:
        for method in args.methods:
            print(f"[{now()}] Starting {method}", flush=True)
            command = child_arguments(args, method, batch)
            with (batch / f"{method}.log").open("w", encoding="utf-8") as log:
                process = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=log, stderr=subprocess.STDOUT)
                manifest["runs"][method] = {"status": "running", "pid": process.pid,
                                            "log": str(batch / f"{method}.log")}
                write_json(batch / "manifest.json", manifest)
                try:
                    result = process.wait()
                except BaseException:
                    process.terminate()
                    process.wait()
                    raise
            report_path = batch / f"{method}.json"
            report = json.loads(report_path.read_text(encoding="utf-8")) if report_path.exists() else {}
            manifest["runs"][method].update(report, returncode=result)
            write_json(batch / "manifest.json", manifest)
            if result != 0 or report.get("status") != "complete":
                raise RuntimeError(f"{method} failed; see {batch / (method + '.log')}")
            print(f"[{now()}] Finished {method}: {report['run_directory']}", flush=True)
        manifest.update(status="complete", finished_at=now())
        write_json(batch / "manifest.json", manifest)
        print(f"Batch complete: {batch / 'manifest.json'}", flush=True)
        return 0
    except BaseException as error:
        manifest.update(status="failed", error=f"{type(error).__name__}: {error}", finished_at=now())
        write_json(batch / "manifest.json", manifest)
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
