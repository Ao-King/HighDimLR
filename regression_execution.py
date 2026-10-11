"""Shared defaults and threshold selection for accelerated regression runs."""

import numpy as np


def method_settings(method, n, largest_singular_value):
    """Return the existing bootstrap experiment's numerical settings."""
    defaults = {
        'SC': dict(selection='oracle_rank', rank_max=n),
        'Ridge': dict(selection='oracle_grid', C=200000, m=100, alpha_divisor=n),
        'Landweber': dict(selection='AOSR', step=2.5e-6, max_steps=10000, minimum_steps=500),
        'Showalter': dict(selection='TDP', step=2.5e-6, max_steps=10000, tau=0.9),
        'ARk': dict(selection='TDP', step=5e-5, max_steps=1000, tau=0.9, kappa=0.5),
        'SOAR': dict(selection='TDP', step=2.5e-4, max_steps=10000, tau=1.01, s_star=1.5),
        'Nesterov': dict(selection='TDP', step=2.5e-6, max_steps=10000, tau=1.01,
                         omega=5, nominal_step=2.5e-6),
    }
    if method not in defaults:
        raise ValueError(f'Unknown method: {method}')
    settings = defaults[method].copy()
    if method == 'Nesterov':
        if not np.isfinite(largest_singular_value) or largest_singular_value <= 0:
            raise ValueError('The largest singular value must be finite and positive')
        settings['step'] = min(settings['nominal_step'], 0.99 / largest_singular_value**2)
    return settings


def threshold_pair(estimates, truth, percentile):
    """Apply the notebooks' original-coordinate threshold grid to both branches."""
    estimates = np.asarray(estimates, dtype=float)
    truth = np.asarray(truth, dtype=float)
    if truth.ndim != 1 or estimates.shape != (truth.size, 2):
        raise ValueError('Expected a coefficient pair of shape (p, 2) and truth of shape (p,)')
    if not np.all(np.isfinite(estimates)) or not np.all(np.isfinite(truth)):
        raise ValueError('Coefficients and truth must be finite')
    if not np.isfinite(percentile) or not 0 <= percentile <= 100:
        raise ValueError('The threshold percentile must lie in [0, 100]')
    grid = 0.001 * np.arange(1, 1001)
    scan = np.empty((grid.size, 3))
    scan[:, 0] = grid
    centers = estimates.copy()
    thresholds = np.empty(2)
    for branch in range(2):
        coefficients = estimates[:, branch]
        for index, threshold in enumerate(grid):
            candidate = np.where(np.abs(coefficients) <= threshold, 0, coefficients)
            scan[index, branch + 1] = np.linalg.norm(candidate - truth)
        errors = scan[:, branch + 1]
        indices = np.flatnonzero(errors == errors.min())
        selected = (np.percentile(indices, percentile) if len(indices) > 1
                    and np.all(np.diff(indices) == 1) else indices[0])
        thresholds[branch] = np.interp(selected, np.arange(len(grid)), grid)
        centers[np.abs(centers[:, branch]) <= thresholds[branch], branch] = 0
    return centers, thresholds, scan
