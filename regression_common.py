"""Shared numerical and sampling conventions for the regression experiments."""

import numpy as np


def ridge_pair(X, Y, alpha):
    """Original and debiased ridge, including the minimum-norm alpha=0 limit."""
    if not np.isfinite(alpha) or alpha < 0:
        raise ValueError("alpha must be finite and nonnegative")
    if alpha == 0:
        beta = np.linalg.lstsq(X, Y, rcond=None)[0]
        return beta, beta.copy()
    n, p = X.shape
    if p <= n:
        system = X.T @ X + alpha * np.eye(p)
        beta = np.linalg.solve(system, X.T @ Y)
        debiased = beta + alpha * np.linalg.solve(system, beta)
    else:
        system = X @ X.T + alpha * np.eye(n)
        dual = np.linalg.solve(system, Y)
        beta = X.T @ dual
        debiased = beta + alpha * X.T @ np.linalg.solve(system, dual)
    return beta, debiased


def lasso_coefficients(X, Y, alpha):
    """Return a converged no-intercept fit, or fail after a bounded retry budget.

    Each retry warm-starts the previous coefficients. The three iteration
    budgets total at most 140000 coordinate-descent sweeps, at tol=1e-6.
    """
    import warnings

    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import Lasso

    if not np.isfinite(alpha) or alpha < 0:
        raise ValueError("alpha must be finite and nonnegative")
    if alpha == 0:
        coefficients = np.linalg.lstsq(X, Y, rcond=None)[0]
        if not np.all(np.isfinite(coefficients)):
            raise FloatingPointError("least-squares coefficients are nonfinite")
        return coefficients

    tolerance = 1e-6
    iteration_budgets = (10000, 30000, 100000)
    model = Lasso(alpha=alpha, fit_intercept=False, tol=tolerance,
                  warm_start=True, selection="cyclic")
    completed_iterations = 0
    for budget in iteration_budgets:
        model.set_params(max_iter=budget)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", ConvergenceWarning)
            model.fit(X, Y)
        convergence_warning = False
        for warning in caught:
            if issubclass(warning.category, ConvergenceWarning):
                convergence_warning = True
            else:
                warnings.warn_explicit(warning.message, warning.category,
                                       warning.filename, warning.lineno)

        coefficients = np.asarray(model.coef_)
        dual_gap = np.asarray(model.dual_gap_, dtype=float)
        if not np.all(np.isfinite(coefficients)) or not np.all(np.isfinite(dual_gap)):
            raise FloatingPointError("Lasso returned nonfinite coefficients or dual gap")
        # sklearn reports dual_gap_ per sample; its stopping threshold is
        # tol * ||Y||^2 / n when no intercept or sample weights are fitted.
        response_energy = np.sum(np.asarray(Y, dtype=float)**2, axis=0) / X.shape[0]
        dual_tolerance = tolerance * response_energy
        roundoff = 100 * np.finfo(coefficients.dtype).eps * np.maximum(1.0, response_energy)
        if not np.all(np.isfinite(dual_tolerance)):
            raise FloatingPointError("Lasso convergence threshold is nonfinite")
        if np.any(dual_gap < -roundoff):
            raise FloatingPointError("Lasso returned an invalid negative dual gap")
        # A zero response and zero coefficients attain the exact minimum.
        # sklearn's strict gap < tolerance test can still warn when both are 0.
        if np.all(np.asarray(Y) == 0) and np.all(coefficients == 0):
            return coefficients.copy()
        completed_iterations += int(np.max(model.n_iter_))
        if not convergence_warning and np.all(dual_gap <= dual_tolerance + roundoff):
            return coefficients.copy()

    raise RuntimeError(
        f"Lasso did not converge for alpha={alpha} after {completed_iterations} "
        f"iterations (budgets={iteration_budgets}, tol={tolerance}); "
        f"dual_gap={dual_gap}, allowed_gap={dual_tolerance}. "
        "No unconverged coefficients were returned."
    )


def sparse_indices(rng, p, size=20):
    """Sample distinct positions, retaining legacy draws when already distinct."""
    if size < 0 or size > p:
        raise ValueError("the sparse support size must be between zero and p")
    draw = rng.integers if hasattr(rng, "integers") else rng.randint
    indices = np.asarray(draw(0, p, size=size), dtype=int)
    seen = set()
    for j in range(size):
        while int(indices[j]) in seen:
            indices[j] = draw(0, p)
        seen.add(int(indices[j]))
    return indices


def discrepancy_candidate(candidates, residuals, delta):
    """Return the actual first grid value satisfying the closed TDP inequality."""
    candidates = np.asarray(candidates)
    residuals = np.asarray(residuals, dtype=float)
    if candidates.ndim != 1 or residuals.shape != candidates.shape:
        raise ValueError("candidates and residuals must be equally sized vectors")
    if not np.all(np.isfinite(residuals)):
        raise FloatingPointError("nonfinite residual in the parameter grid")
    acceptable = np.flatnonzero(residuals <= delta)
    if not len(acceptable):
        raise ValueError("no regularization-grid candidate satisfies the discrepancy rule")
    return candidates[acceptable[0]]
