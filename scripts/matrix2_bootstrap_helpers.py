"""Exact sparse-center acceleration for the matrix-2 percentile bootstrap.

This module changes only the maximum-error reduction after the same float64
spectral fits. It preserves the canonical independent branch draws generated
by regression_bootstrap._noise. Fits are discarded after each batch.
"""

import numpy as np
from threadpoolctl import threadpool_limits

from regression_bootstrap import _arguments, _noise, _positive_integer


DEFAULT_NOISE_ROOT = 2026100802


def max_errors_from_fits(fitted, center, thresholds, *, max_active_cells=4_000_000):
    """Return exact (Q,B) hard-thresholded max errors, using sparse centers.

    At coordinates whose center is zero, the maximum error is either the
    largest absolute fitted coefficient or zero, according as that largest
    coefficient exceeds the threshold. Only nonzero-center coordinates need
    the complete Q-by-coordinate-by-bootstrap calculation. Dense centers are
    supported by processing a bounded number of thresholds at a time.
    """
    if any(np.iscomplexobj(value) for value in (fitted, center, thresholds)):
        raise ValueError("Inputs must be real")
    fitted, center, thresholds = (np.asarray(value, dtype=np.float64)
                                  for value in (fitted, center, thresholds))
    if (fitted.ndim != 2 or not min(fitted.shape)
            or center.shape != (fitted.shape[0],)
            or thresholds.ndim != 1 or not thresholds.size):
        raise ValueError("Expected fitted (p,B), center (p,), thresholds (Q,)")
    if any(not np.all(np.isfinite(value)) for value in (fitted, center, thresholds)):
        raise ValueError("Inputs must be finite")
    if np.any(thresholds < 0):
        raise ValueError("Thresholds must be nonnegative")
    max_active_cells = _positive_integer(max_active_cells, "max_active_cells")
    active = center != 0
    B = fitted.shape[1]
    if np.any(~active):
        inactive_max = np.max(np.abs(fitted[~active]), axis=0)
    else:
        inactive_max = np.zeros(B, dtype=np.float64)
    result = np.where(inactive_max[None, :] <= thresholds[:, None],
                      0.0, inactive_max[None, :])
    if not np.any(active):
        return result
    values = fitted[active]
    absolute = np.abs(values)
    active_center = center[active]
    block = max(1, max_active_cells // (len(active_center) * B))
    with np.errstate(over="raise", invalid="raise"):
        for first in range(0, len(thresholds), block):
            stop = min(len(thresholds), first + block)
            thresholded = np.where(absolute[None, :, :] <= thresholds[first:stop, None, None],
                                   0.0, values[None, :, :])
            active_error = np.max(np.abs(thresholded - active_center[None, :, None]), axis=1)
            result[first:stop] = np.maximum(result[first:stop], active_error)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("Nonfinite sparse-center maximum error")
    return result


def bootstrap_errors_grid_fast(design, method, parameter, settings, centers,
                               sigma, thresholds, B, seed_root=DEFAULT_NOISE_ROOT,
                               outer_index=0, *, batch_size=500, blas_threads=1,
                               max_active_cells=4_000_000):
    """Return (Q,B,2), reusing fits only for byte-identical center/sigma pairs.

    This is a CPU implementation. The caller can parallelize outer replications
    in separate processes. One process-global BLAS limit covers the complete
    call, avoiding competing thread limits within a process.
    """
    B = _positive_integer(B, "B")
    batch_size = _positive_integer(batch_size, "batch_size")
    blas_threads = _positive_integer(blas_threads, "blas_threads")
    seed_root = _positive_integer(seed_root, "seed_root", 0)
    outer_index = _positive_integer(outer_index, "outer_index", 0)
    max_active_cells = _positive_integer(max_active_cells, "max_active_cells")
    if any(np.iscomplexobj(value) for value in (centers, sigma, thresholds)):
        raise ValueError("Bootstrap inputs must be real")
    centers, sigma, thresholds = (np.asarray(value, dtype=np.float64)
                                  for value in (centers, sigma, thresholds))
    if (centers.ndim != 3 or not len(centers) or centers.shape[2] != 2
            or sigma.shape != (len(centers), 2)
            or thresholds.shape != (len(centers), 2)):
        raise ValueError("Expected centers (Q,p,2), sigma (Q,2), thresholds (Q,2)")
    if (any(not np.all(np.isfinite(value)) for value in (centers, sigma, thresholds))
            or np.any(sigma < 0) or np.any(thresholds < 0)):
        raise ValueError("Bootstrap inputs must be finite with nonnegative sigma/thresholds")
    X, U, s, Vt, _, _, _ = _arguments(design, centers[0], sigma[0], thresholds[0])
    if centers.shape[1] != X.shape[1]:
        raise ValueError("Centers do not match the design")
    with np.errstate(over="raise", invalid="raise"):
        ordinary_gain = np.asarray(design.filter(method, parameter, settings), dtype=np.float64)
        if ordinary_gain.shape != s.shape or not np.all(np.isfinite(ordinary_gain)):
            raise ValueError("Invalid spectral filter")
        debiased_gain = (ordinary_gain.copy() if str(method).lower() in ("sc", "cut")
                         else 2 * ordinary_gain - s * ordinary_gain**2)
    if not np.all(np.isfinite(debiased_gain)):
        raise FloatingPointError("Nonfinite debiased spectral filter")
    groups = {}
    for index in range(len(centers)):
        groups.setdefault((centers[index].tobytes(), sigma[index].tobytes()), []).append(index)
    result = np.empty((len(centers), B, 2), dtype=np.float64)
    with threadpool_limits(limits=blas_threads, user_api="blas"):
        for indices in groups.values():
            first = indices[0]
            center = centers[first]
            selected_thresholds = thresholds[indices]
            noise = _noise(X.shape[0], sigma[first], B, seed_root, outer_index)
            with np.errstate(over="raise", invalid="raise"):
                means = X @ center
            if not np.all(np.isfinite(means)):
                raise FloatingPointError("Nonfinite response center")
            for start in range(0, B, batch_size):
                stop = min(B, start + batch_size)
                count = stop - start
                with np.errstate(over="raise", invalid="raise"):
                    responses = np.concatenate((means[:, :1] + noise[0, :, start:stop],
                                                means[:, 1:] + noise[1, :, start:stop]), axis=1)
                    coordinates = U.T @ responses
                    coordinates[:, :count] *= ordinary_gain[:, None]
                    coordinates[:, count:] *= debiased_gain[:, None]
                    fitted = Vt.T @ coordinates
                for branch, values in enumerate((fitted[:, :count], fitted[:, count:])):
                    errors = max_errors_from_fits(values, center[:, branch],
                                                  selected_thresholds[:, branch],
                                                  max_active_cells=max_active_cells)
                    result[np.asarray(indices), start:stop, branch] = errors
    return result


bootstrap_errors_grid_sparse = bootstrap_errors_grid_fast
