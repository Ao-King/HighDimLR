"""Reproduce the non-sparse matrix-02 Lasso TDP calculation.

LARS only supplies initial estimates; every examined grid candidate is then
accepted only after converged cyclic coordinate descent with tol=1e-6.
No true coefficient values are used for selecting the parameter.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import warnings

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

import numpy as np
import sklearn
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import Lasso, lars_path
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "outputs/regression/nonsparse/Lasso"
DATA = ROOT / "data/regression/nonsparse/data.npz"
TOL = 1e-6


def emit(**message):
    line = json.dumps(message, ensure_ascii=True)
    print(line, flush=True)
    with (DEST / "events.jsonl").open("a", encoding="utf-8") as stream:
        stream.write(line + "\n")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    started = time.perf_counter()
    DEST.mkdir(parents=True, exist_ok=True)
    (DEST / "events.jsonl").write_text("", encoding="utf-8")
    data_sha = sha(DATA)
    with np.load(DATA) as data:
        X = np.array(data["X"], dtype=np.float64, order="F")
        y = np.asarray(data["Y"], dtype=np.float64).reshape(-1)
        truth = np.asarray(data["beta"], dtype=np.float64).reshape(-1)
        noise = np.asarray(data["e"], dtype=np.float64).reshape(-1)
    delta = float(.6 * np.linalg.norm(noise))
    assert np.allclose(y, X @ truth + noise, rtol=0, atol=1e-10)
    energy = float(y @ y / X.shape[0])
    allowed = TOL * energy
    roundoff = 100 * np.finfo(float).eps * max(1., energy)
    trace, gaps, iterations, retries, initial_norm_changes = [], [], [], [], []
    with threadpool_limits(limits=1):
        gram = np.ascontiguousarray(X.T @ X)
        xy = X.T @ y
        emit(event="initializing_lars", data_sha256=data_sha,
             delta=delta, allowed_dual_gap=allowed)
        lars_started = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            knots, active, knot_coefs = lars_path(
                X, y, Xy=xy, Gram=gram, alpha_min=.001, method="lasso",
                max_iter=20000, return_path=True)
        path_warnings = [str(w.message) for w in caught]
        np.savez_compressed(DEST / "initialization.npz", alpha=knots,
                            coefficients=knot_coefs)
        emit(event="initialized_lars", knots=len(knots),
             minimum_alpha=float(knots[-1]), warnings=path_warnings,
             seconds=time.perf_counter()-lars_started)
        if knots[-1] > .001 + 1e-12:
            raise RuntimeError("LARS initialization did not reach alpha=.001")
        asc_alpha, asc_coef = knots[::-1], knot_coefs[:, ::-1]

        def initial(alpha):
            position = int(np.searchsorted(asc_alpha, alpha, side="right"))
            if position == 0:
                return asc_coef[:, 0].copy()
            if position == len(asc_alpha):
                return asc_coef[:, -1].copy()
            low, high = asc_alpha[position-1], asc_alpha[position]
            fraction = (alpha - low) / (high - low)
            return (asc_coef[:, position-1] + fraction * (
                asc_coef[:, position]-asc_coef[:, position-1])).copy()

        model = Lasso(alpha=1., fit_intercept=False, tol=TOL,
                      warm_start=True, selection="cyclic", precompute=gram)
        checkpoint = time.perf_counter()
        selected = None
        for numerator in range(1000, -1, -1):
            alpha = numerator / 1000.
            if numerator == 0:
                coef = np.linalg.lstsq(X, y, rcond=None)[0]
                gap, count, retry, change = 0., 0, 0, 0.
            else:
                start_coef = initial(alpha)
                model.set_params(alpha=alpha)
                model.coef_ = start_coef.copy()
                count = 0
                for retry, budget in enumerate((10000, 30000, 100000)):
                    model.set_params(max_iter=budget)
                    with warnings.catch_warnings(record=True) as fit_warnings:
                        warnings.simplefilter("always", ConvergenceWarning)
                        model.fit(X, y)
                    unexpected = [str(w.message) for w in fit_warnings
                                  if not issubclass(w.category, ConvergenceWarning)]
                    if unexpected:
                        emit(event="fit_nonconvergence_warnings", alpha=alpha,
                             warnings=unexpected)
                    failed = any(issubclass(w.category, ConvergenceWarning)
                                 for w in fit_warnings)
                    gap = float(model.dual_gap_)
                    count += int(model.n_iter_)
                    coef = model.coef_.copy()
                    if (not np.isfinite(gap) or not np.isfinite(coef).all()
                            or gap < -roundoff):
                        raise FloatingPointError("Invalid Lasso coefficients or gap")
                    if not failed and gap <= allowed + roundoff:
                        break
                else:
                    raise RuntimeError(f"CD did not converge at alpha={alpha}")
                change = float(np.linalg.norm(coef - start_coef))
            residual = float(np.linalg.norm(y-X@coef))
            error = float(np.linalg.norm(coef-truth))
            trace.append([alpha, error, residual])
            gaps.append(gap)
            iterations.append(count)
            retries.append(retry)
            initial_norm_changes.append(change)
            if len(trace) % 100 == 0 or time.perf_counter()-checkpoint > 40:
                emit(event="progress", examined=len(trace), alpha=alpha,
                     residual=residual, dual_gap=gap, iterations=count,
                     elapsed_seconds=time.perf_counter()-started)
                checkpoint = time.perf_counter()
            if residual <= delta:
                selected = coef.copy()
                emit(event="tdp_selected", alpha=alpha, residual=residual,
                     l2=error, dual_gap=gap, iterations=count,
                     candidates_examined=len(trace))
                break
        if selected is None:
            raise RuntimeError("No grid candidate met TDP")

        # Independent calculations from the saved estimate, including KKT and
        # dual-gap checks with direct matrix products rather than Gram products.
        alpha, saved_error, saved_residual = trace[-1]
        direct_residual = y-X@selected
        recalculated_error = float(np.linalg.norm(selected-truth))
        recalculated_residual = float(np.linalg.norm(direct_residual))
        gradient = X.T @ direct_residual / X.shape[0]
        active = selected != 0.
        kkt_active = float(np.max(np.abs(gradient[active]
                    - alpha*np.sign(selected[active])))) if active.any() else 0.
        kkt_inactive = float(np.max(np.maximum(np.abs(gradient[~active])
                      - alpha, 0.))) if (~active).any() else 0.
        dual_scale = min(1., alpha / float(np.max(np.abs(gradient))))
        dual = direct_residual * dual_scale
        primal_objective = float(direct_residual@direct_residual/(2*X.shape[0])
                                 + alpha*np.sum(np.abs(selected)))
        dual_objective = float((y@dual - .5*dual@dual)/X.shape[0])
        direct_gap = primal_objective - dual_objective
        assert direct_gap >= -roundoff
        assert direct_gap <= allowed+roundoff
        assert recalculated_residual <= delta
        assert all(row[2] > delta for row in trace[:-1])
        assert sha(DATA) == data_sha

        np.savez_compressed(DEST / "fit.npz", estimate=selected,
                            trace=np.asarray(trace), dual_gap=np.asarray(gaps),
                            iterations=np.asarray(iterations),
                            retries=np.asarray(retries),
                            initialization_norm_change=np.asarray(initial_norm_changes))
        result = dict(status="complete", method="Lasso", rule="TDP",
            data_file_sha256=data_sha, script_sha256=sha(__file__), parameter=alpha,
            ordinary_l2=recalculated_error,
            ordinary_linf=float(np.linalg.norm(selected-truth, ord=np.inf)),
            ordinary_residual=recalculated_residual,
            ordinary_sigma=recalculated_residual/np.sqrt(X.shape[0]),
            elapsed_seconds=time.perf_counter()-started,
            diagnostics=dict(selection="discrepancy_grid", fit_intercept=False,
                parameter_kind="alpha", selected_alpha=alpha,
                selected_grid_numerator=int(round(1000*alpha)),
                selected_trace_index=len(trace)-1, recorded_candidates=len(trace),
                trace_columns=["alpha", "coefficient_l2_error", "residual_norm"],
                descending_grid=dict(start=1., stop=0., step=-.001),
                tdp_satisfied=True, stopping_threshold=delta,
                residual_norm=recalculated_residual,
                coefficient_l2_error=recalculated_error,
                previous_candidate_alpha=trace[-2][0] if len(trace)>1 else None,
                previous_candidate_residual=trace[-2][2] if len(trace)>1 else None,
                all_examined_candidates_converged=True, selection_rule_uses_truth=False,
                solver="sklearn Lasso cyclic coordinate descent",
                sklearn_version=sklearn.__version__, convergence_tolerance=TOL,
                allowed_dual_gap=allowed, dual_gap=gaps[-1],
                direct_dual_gap=direct_gap,
                max_dual_gap=float(max(gaps)), iterations=iterations[-1],
                total_iterations=int(sum(iterations)),
                max_iterations=int(max(iterations)),
                kkt_active_max_error=kkt_active, kkt_inactive_max_violation=kkt_inactive,
                initialization="LARS path interpolation; every examined point validated by converged cyclic CD",
                initialization_path_knots=len(knots),
                initialization_path_warnings=path_warnings,
                initialization_max_coefficient_norm_change=float(max(initial_norm_changes)),
                blas_threads=1,
                boundary_checks=[dict(check="previous_candidate_rejects_tdp", passed=True),
                                 dict(check="all_prior_candidates_reject_tdp", passed=True),
                                 dict(check="selected_candidate_satisfies_tdp", passed=True),
                                 dict(check="direct_dual_gap_within_tolerance", passed=True),
                                 dict(check="immutable_data_hash", passed=True)]))
        (DEST / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        emit(event="complete", **result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=DATA)
    parser.add_argument('--output', type=Path, default=DEST)
    args = parser.parse_args()
    DATA, DEST = args.data, args.output
    main()
