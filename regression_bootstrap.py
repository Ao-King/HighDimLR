"""Deterministic, float64 batched bootstrap execution on CPU or optional CUDA."""

from concurrent.futures import ThreadPoolExecutor, as_completed
import importlib
import threading
import weakref

import numpy as np
from threadpoolctl import threadpool_limits


_GPU_CACHE = weakref.WeakKeyDictionary()
_GPU_CACHE_LOCK = threading.Lock()


def _positive_integer(value, name, minimum=1):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _devices(gpu_devices):
    try:
        devices = tuple(_positive_integer(device, "GPU device ID", 0) for device in gpu_devices)
    except TypeError as error:
        raise ValueError("gpu_devices must be a nonempty sequence of device IDs") from error
    if not devices or len(set(devices)) != len(devices):
        raise ValueError("gpu_devices must be nonempty and contain no duplicate IDs")
    return devices


def _cupy():
    try:
        return importlib.import_module("cupy")
    except (ImportError, OSError) as error:
        raise RuntimeError("GPU backend requires CuPy and a working CUDA installation; no CPU fallback is performed") from error


def validate_backend(backend="cpu", gpu_devices=(0,)):
    """Fail before fitting when the requested backend/device is unavailable."""
    if isinstance(backend, dict):
        gpu_devices = backend.get("gpu_devices", gpu_devices)
        backend = backend.get("backend", "cpu")
    if backend == "cpu":
        return {"backend": "cpu", "numpy": np.__version__, "dtype": "float64"}
    if backend != "gpu":
        raise ValueError("backend must be 'cpu' or 'gpu'")
    devices = _devices(gpu_devices)
    cp = _cupy()
    try:
        count = int(cp.cuda.runtime.getDeviceCount())
        if count == 0:
            raise RuntimeError("No CUDA GPU is available")
        if any(device >= count for device in devices):
            raise RuntimeError(f"Requested GPU device IDs {devices}, but only {count} device(s) are available")
        details = []
        for device in devices:
            with cp.cuda.Device(device):
                properties = cp.cuda.runtime.getDeviceProperties(device)
                name = properties["name"]
                if isinstance(name, bytes):
                    name = name.decode("utf-8", errors="replace")
                probe = cp.ones(1, dtype=cp.float64)
                if float(cp.asnumpy(probe @ probe)) != 1.0:
                    raise RuntimeError(f"GPU {device} failed the float64 initialization check")
                cp.cuda.get_current_stream().synchronize()
                details.append({"id": device, "name": name})
        return {"backend": "gpu", "numpy": np.__version__, "cupy": cp.__version__,
                "cuda_runtime": int(cp.cuda.runtime.runtimeGetVersion()), "devices": details,
                "dtype": "float64"}
    except Exception as error:
        raise RuntimeError(f"GPU backend initialization failed: {error}; no CPU fallback is performed") from error


def _noise(n, sigma, B, seed_root, outer_index):
    noise = np.empty((2, n, B), dtype=np.float64)
    for seed in range(B):
        rng = np.random.default_rng((seed_root, outer_index, 1, seed))
        noise[0, :, seed] = rng.normal(0, sigma[0], n)
        noise[1, :, seed] = rng.normal(0, sigma[1], n)
    if not np.all(np.isfinite(noise)):
        raise FloatingPointError("Nonfinite bootstrap noise")
    return noise


def _arguments(design, centers, sigma, thresholds):
    X, U, s, Vt = (np.asarray(getattr(design, name), dtype=np.float64) for name in ("X", "U", "s", "Vt"))
    if X.ndim != 2 or not min(X.shape):
        raise ValueError("design.X must be a nonempty matrix")
    n, p = X.shape
    if s.ndim != 1 or U.shape != (n, len(s)) or Vt.shape != (len(s), p):
        raise ValueError("The spectral design must contain compatible thin SVD arrays")
    if any(not np.all(np.isfinite(value)) for value in (X, U, s, Vt)) or np.any(s < 0):
        raise ValueError("Spectral design values must be finite with nonnegative singular values")
    centers, sigma, thresholds = (np.asarray(value, dtype=np.float64) for value in (centers, sigma, thresholds))
    if centers.shape != (p, 2) or sigma.shape != (2,) or thresholds.shape != (2,):
        raise ValueError("centers, sigma, thresholds must have shapes (p,2), (2,), (2,)")
    if any(not np.all(np.isfinite(value)) for value in (centers, sigma, thresholds)):
        raise ValueError("Bootstrap inputs must be finite")
    if np.any(sigma < 0) or np.any(thresholds < 0):
        raise ValueError("Noise standard deviations and thresholds must be nonnegative")
    return X, U, s, Vt, centers, sigma, thresholds


def _cpu_chunk(U, Vt, a, debiased, means, centers, thresholds, noise, start, stop):
    count = stop - start
    with np.errstate(over="raise", invalid="raise"):
        responses = np.concatenate((means[:, :1] + noise[0, :, start:stop],
                                    means[:, 1:] + noise[1, :, start:stop]), axis=1)
        coordinates = U.T @ responses
        coordinates[:, :count] *= a[:, None]
        coordinates[:, count:] *= debiased[:, None]
        fitted = Vt.T @ coordinates
        # Threshold only after returning to the original coefficient coordinates.
        fitted[:, :count][np.abs(fitted[:, :count]) <= thresholds[0]] = 0
        fitted[:, count:][np.abs(fitted[:, count:]) <= thresholds[1]] = 0
        errors = np.column_stack((np.max(np.abs(fitted[:, :count] - centers[:, :1]), axis=0),
                                  np.max(np.abs(fitted[:, count:] - centers[:, 1:]), axis=0)))
    if not np.all(np.isfinite(errors)):
        raise FloatingPointError("Nonfinite CPU bootstrap result")
    return start, stop, errors


def _gpu_design(cp, design, device, U, Vt):
    with _GPU_CACHE_LOCK:
        cache = _GPU_CACHE.setdefault(design, {})
        if device not in cache:
            # The creator's stream may differ from the next call's worker stream.
            cache[device] = (cp.asarray(U), cp.asarray(Vt.T))
            cp.cuda.get_current_stream().synchronize()
        return cache[device]


def _gpu_worker(cp, design, device, U, Vt, a, debiased, means, centers, thresholds, noise, chunks, stop_event):
    results = []
    with cp.cuda.Device(device):
        device_u, device_v = _gpu_design(cp, design, device, U, Vt)
        device_a, device_debiased = cp.asarray(a), cp.asarray(debiased)
        device_means, device_centers = cp.asarray(means), cp.asarray(centers)
        for start, stop in chunks:
            if stop_event.is_set():
                break
            count = stop - start
            responses = cp.asarray(np.concatenate((noise[0, :, start:stop], noise[1, :, start:stop]), axis=1))
            responses[:, :count] += device_means[:, :1]
            responses[:, count:] += device_means[:, 1:]
            coordinates = device_u.T @ responses
            coordinates[:, :count] *= device_a[:, None]
            coordinates[:, count:] *= device_debiased[:, None]
            fitted = device_v @ coordinates
            fitted[:, :count] = cp.where(cp.abs(fitted[:, :count]) <= thresholds[0], 0., fitted[:, :count])
            fitted[:, count:] = cp.where(cp.abs(fitted[:, count:]) <= thresholds[1], 0., fitted[:, count:])
            errors = cp.column_stack((cp.max(cp.abs(fitted[:, :count] - device_centers[:, :1]), axis=0),
                                      cp.max(cp.abs(fitted[:, count:] - device_centers[:, 1:]), axis=0)))
            errors = cp.asnumpy(errors)
            if not np.all(np.isfinite(errors)):
                raise FloatingPointError(f"Nonfinite bootstrap result on GPU {device}")
            results.append((start, stop, errors))
        cp.cuda.get_current_stream().synchronize()
    return results


def bootstrap_errors(design, method, parameter, settings, centers, sigma, thresholds,
                     B, seed_root, outer_index, *, backend="cpu", batch_size=64,
                     cpu_workers=4, blas_threads=8, gpu_devices=(0,)):
    """Return B paired maximum errors, with the legacy per-replicate RNG stream."""
    B = _positive_integer(B, "B")
    batch_size = _positive_integer(batch_size, "batch_size")
    cpu_workers = _positive_integer(cpu_workers, "cpu_workers")
    blas_threads = _positive_integer(blas_threads, "blas_threads")
    seed_root = _positive_integer(seed_root, "seed_root", 0)
    outer_index = _positive_integer(outer_index, "outer_index", 0)
    validate_backend(backend, gpu_devices)
    X, U, s, Vt, centers, sigma, thresholds = _arguments(design, centers, sigma, thresholds)
    with np.errstate(over="raise", invalid="raise"):
        a = np.asarray(design.filter(method, parameter, settings), dtype=np.float64)
        if a.shape != s.shape or not np.all(np.isfinite(a)):
            raise ValueError("The spectral filter must be a finite vector matching the singular values")
        # SC is exactly idempotent: retain identical floating-point gains as well.
        debiased = a.copy() if str(method).lower() in ("sc", "cut") else 2 * a - s * a**2
    if not np.all(np.isfinite(debiased)):
        raise FloatingPointError("Nonfinite debiased spectral filter")
    noise = _noise(X.shape[0], sigma, B, seed_root, outer_index)
    chunks = [(start, min(B, start + batch_size)) for start in range(0, B, batch_size)]
    result = np.empty((B, 2), dtype=np.float64)
    # One outer limit context avoids racing process-global BLAS settings in workers.
    with threadpool_limits(limits=blas_threads, user_api="blas"):
        with np.errstate(over="raise", invalid="raise"):
            means = X @ centers
        if not np.all(np.isfinite(means)):
            raise FloatingPointError("Nonfinite bootstrap response center")
        if backend == "cpu":
            workers = min(cpu_workers, len(chunks))
            if workers == 1:
                for start, stop in chunks:
                    _, _, result[start:stop] = _cpu_chunk(U, Vt, a, debiased, means, centers, thresholds, noise, start, stop)
            else:
                with ThreadPoolExecutor(max_workers=workers) as executor:
                    futures = [executor.submit(_cpu_chunk, U, Vt, a, debiased, means, centers,
                                               thresholds, noise, start, stop) for start, stop in chunks]
                    try:
                        for future in as_completed(futures):
                            start, stop, errors = future.result()
                            result[start:stop] = errors
                    except BaseException:
                        for future in futures:
                            future.cancel()
                        raise
        else:
            cp = _cupy()
            devices = _devices(gpu_devices)
            stop_event = threading.Event()
            try:
                with ThreadPoolExecutor(max_workers=len(devices)) as executor:
                    futures = [executor.submit(_gpu_worker, cp, design, device, U, Vt, a, debiased,
                                               means, centers, thresholds, noise, chunks[index::len(devices)], stop_event)
                               for index, device in enumerate(devices)]
                    try:
                        for future in as_completed(futures):
                            for start, stop, errors in future.result():
                                result[start:stop] = errors
                    except BaseException:
                        stop_event.set()
                        for future in futures:
                            future.cancel()
                        raise
            except Exception as error:
                raise RuntimeError(f"GPU bootstrap failed: {error}; no CPU fallback is performed") from error
    return result
