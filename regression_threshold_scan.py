"""Paired threshold calibration with the original discrete regression filters."""

from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

import numpy as np
from threadpoolctl import threadpool_limits

from regression_bootstrap import (_arguments, _cupy, _devices, _gpu_design,
                                  _noise, _positive_integer, validate_backend)
from regression_execution import method_settings


METHODS = ("SC", "Ridge", "Landweber", "Showalter", "ARk", "SOAR", "Nesterov", "HBF", "FAR")
BASELINE_PERCENTILES = dict(zip(METHODS, (50, 0, 72.5, 90, 0, 5, 95, 50, 50)))
VERSION = "regression-threshold-grid-v1"


def scan_method_settings(method, n, largest_singular_value):
    if method == "HBF":
        return dict(selection="TDP", step=5e-4, max_steps=10000, tau=1.0, eta=5.0)
    if method == "FAR":
        return dict(selection="TDP", step=5e-5, max_steps=500, tau=1.0, theta=1.5)
    return method_settings(method, n, largest_singular_value)


def threshold_grid_pair(estimates, truth, percentiles):
    """Compute one oracle error curve, then apply each requested tie rule."""
    if any(np.iscomplexobj(value) for value in (estimates, truth, percentiles)):
        raise ValueError("Threshold inputs must be real")
    estimates = np.asarray(estimates, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    percentiles = np.asarray(percentiles, dtype=np.float64)
    if truth.ndim != 1 or estimates.shape != (truth.size, 2):
        raise ValueError("Expected estimates (p,2) and truth (p,)")
    if (percentiles.ndim != 1 or not percentiles.size
            or np.any((percentiles < 0) | (percentiles > 100))
            or not np.all(np.isfinite(percentiles))
            or len(np.unique(percentiles)) != len(percentiles)):
        raise ValueError("Percentiles must be distinct finite values in [0,100]")
    if not np.all(np.isfinite(estimates)) or not np.all(np.isfinite(truth)):
        raise ValueError("Estimates and truth must be finite")
    grid = 0.001 * np.arange(1, 1001)
    scan = np.empty((grid.size, 3))
    scan[:, 0] = grid
    centers = np.broadcast_to(estimates, (len(percentiles), *estimates.shape)).copy()
    thresholds = np.empty((len(percentiles), 2))
    for branch in range(2):
        coefficients = estimates[:, branch]
        for index, threshold in enumerate(grid):
            candidate = np.where(np.abs(coefficients) <= threshold, 0, coefficients)
            scan[index, branch + 1] = np.linalg.norm(candidate - truth)
        errors = scan[:, branch + 1]
        if not np.all(np.isfinite(errors)):
            raise FloatingPointError("Nonfinite oracle threshold errors")
        indices = np.flatnonzero(errors == errors.min())
        selected = (np.percentile(indices, percentiles) if len(indices) > 1
                    and np.all(np.diff(indices) == 1)
                    else np.full(len(percentiles), indices[0]))
        thresholds[:, branch] = np.interp(selected, np.arange(len(grid)), grid)
        centers[:, :, branch] = np.where(np.abs(coefficients)[None, :] <= thresholds[:, branch, None],
                                         0, coefficients[None, :])
    return centers, thresholds, scan


def _groups(centers, sigma):
    groups = {}
    for index in range(len(centers)):
        key = (centers[index].tobytes(), sigma[index].tobytes())
        groups.setdefault(key, []).append(index)
    return list(groups.values())


def _cpu_grid_chunk(U, Vt, a, debiased, means, center, thresholds, noise, start, stop):
    count = stop - start
    with np.errstate(over="raise", invalid="raise"):
        responses = np.concatenate((means[:, :1] + noise[0, :, start:stop],
                                    means[:, 1:] + noise[1, :, start:stop]), axis=1)
        coordinates = U.T @ responses
        coordinates[:, :count] *= a[:, None]
        coordinates[:, count:] *= debiased[:, None]
        fitted = Vt.T @ coordinates
        errors = np.empty((len(thresholds), count, 2))
        for branch, values in enumerate((fitted[:, :count], fitted[:, count:])):
            candidates = np.where(np.abs(values)[None, :, :] <= thresholds[:, branch, None, None],
                                  0, values[None, :, :])
            errors[:, :, branch] = np.max(np.abs(candidates-center[None, :, branch, None]), axis=1)
    if not np.all(np.isfinite(errors)):
        raise FloatingPointError("Nonfinite threshold-grid CPU bootstrap result")
    return start, stop, errors


def _gpu_grid_worker(cp, design, device, U, Vt, a, debiased, means, center,
                     thresholds, noise, chunks, stop_event):
    results = []
    with cp.cuda.Device(device):
        device_u, device_v = _gpu_design(cp, design, device, U, Vt)
        device_a, device_debiased = cp.asarray(a), cp.asarray(debiased)
        device_means, device_center = cp.asarray(means), cp.asarray(center)
        device_thresholds = cp.asarray(thresholds)
        for start, stop in chunks:
            if stop_event.is_set():
                break
            count = stop-start
            responses = cp.asarray(np.concatenate((noise[0, :, start:stop],
                                                     noise[1, :, start:stop]), axis=1))
            responses[:, :count] += device_means[:, :1]
            responses[:, count:] += device_means[:, 1:]
            coordinates = device_u.T @ responses
            coordinates[:, :count] *= device_a[:, None]
            coordinates[:, count:] *= device_debiased[:, None]
            fitted = device_v @ coordinates
            errors = cp.empty((len(thresholds), count, 2), dtype=cp.float64)
            for branch, values in enumerate((fitted[:, :count], fitted[:, count:])):
                candidates = cp.where(cp.abs(values)[None, :, :] <= device_thresholds[:, branch, None, None],
                                      0, values[None, :, :])
                errors[:, :, branch] = cp.max(cp.abs(candidates-device_center[None, :, branch, None]), axis=1)
            errors = cp.asnumpy(errors)
            if not np.all(np.isfinite(errors)):
                raise FloatingPointError(f"Nonfinite threshold-grid bootstrap on GPU {device}")
            results.append((start, stop, errors))
        cp.cuda.get_current_stream().synchronize()
    return results


def bootstrap_errors_grid(design, method, parameter, settings, centers, sigma, thresholds,
                          B, seed_root, outer_index, *, backend="cpu", batch_size=64,
                          cpu_workers=4, blas_threads=8, gpu_devices=(0,)):
    """Return (Q,B,2) errors; reuse fits only for exactly equal center/noise pairs."""
    B = _positive_integer(B, "B")
    batch_size = _positive_integer(batch_size, "batch_size")
    cpu_workers = _positive_integer(cpu_workers, "cpu_workers")
    blas_threads = _positive_integer(blas_threads, "blas_threads")
    seed_root = _positive_integer(seed_root, "seed_root", 0)
    outer_index = _positive_integer(outer_index, "outer_index", 0)
    validate_backend(backend, gpu_devices)
    if any(np.iscomplexobj(value) for value in (centers, sigma, thresholds)):
        raise ValueError("Bootstrap inputs must be real")
    centers, sigma, thresholds = (np.asarray(value, dtype=np.float64) for value in (centers, sigma, thresholds))
    if (centers.ndim != 3 or not len(centers) or centers.shape[2] != 2
            or sigma.shape != (len(centers), 2) or thresholds.shape != (len(centers), 2)):
        raise ValueError("Expected centers (Q,p,2), sigma (Q,2), thresholds (Q,2)")
    if (any(not np.all(np.isfinite(value)) for value in (centers, sigma, thresholds))
            or np.any(sigma < 0) or np.any(thresholds < 0)):
        raise ValueError("Bootstrap inputs must be finite; sigma/thresholds nonnegative")
    X, U, s, Vt, _, _, _ = _arguments(design, centers[0], sigma[0], thresholds[0])
    if centers.shape[1] != X.shape[1]:
        raise ValueError("Centers do not match the design")
    with np.errstate(over="raise", invalid="raise"):
        a = np.asarray(design.filter(method, parameter, settings), dtype=np.float64)
        if a.shape != s.shape or not np.all(np.isfinite(a)):
            raise ValueError("Invalid spectral filter")
        debiased = a.copy() if str(method).lower() in ("sc", "cut") else 2*a-s*a**2
    if not np.all(np.isfinite(debiased)):
        raise FloatingPointError("Nonfinite debiased filter")
    chunks = [(start, min(B, start+batch_size)) for start in range(0, B, batch_size)]
    result = np.empty((len(centers), B, 2), dtype=np.float64)
    with threadpool_limits(limits=blas_threads, user_api="blas"):
        for indices in _groups(centers, sigma):
            index = indices[0]
            center, selected_thresholds = centers[index], thresholds[indices]
            noise = _noise(X.shape[0], sigma[index], B, seed_root, outer_index)
            with np.errstate(over="raise", invalid="raise"):
                means = X @ center
            if not np.all(np.isfinite(means)):
                raise FloatingPointError("Nonfinite response center")
            if backend == "cpu":
                workers = min(cpu_workers, len(chunks))
                if workers == 1:
                    values = [_cpu_grid_chunk(U,Vt,a,debiased,means,center,selected_thresholds,noise,start,stop)
                              for start,stop in chunks]
                else:
                    with ThreadPoolExecutor(max_workers=workers) as executor:
                        futures = [executor.submit(_cpu_grid_chunk,U,Vt,a,debiased,means,center,
                                                   selected_thresholds,noise,start,stop) for start,stop in chunks]
                        values = [future.result() for future in as_completed(futures)]
            else:
                cp, devices = _cupy(), _devices(gpu_devices)
                stop_event = threading.Event()
                try:
                    with ThreadPoolExecutor(max_workers=len(devices)) as executor:
                        futures = [executor.submit(_gpu_grid_worker,cp,design,device,U,Vt,a,debiased,
                                                   means,center,selected_thresholds,noise,
                                                   chunks[position::len(devices)],stop_event)
                                   for position,device in enumerate(devices)]
                        values = []
                        try:
                            for future in as_completed(futures):
                                values.extend(future.result())
                        except BaseException:
                            stop_event.set()
                            for future in futures:
                                future.cancel()
                            raise
                except Exception as error:
                    raise RuntimeError(f"GPU threshold-grid bootstrap failed: {error}") from error
            for start, stop, errors in values:
                result[np.asarray(indices), start:stop, :] = errors
    return result
