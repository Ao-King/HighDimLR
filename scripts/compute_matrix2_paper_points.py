"""Complete fixed-index Normal/Laplace point tables on immutable matrix 02."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time
import warnings

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from regression_common import lasso_coefficients
from regression_results import _array_hash
from regression_spectral import SpectralDesign
from regression_threshold_scan import threshold_grid_pair
from paper_regression_config import SETTINGS, PERCENTILES, load_design

CASE = ROOT / 'data/regression/caseI'
DEST = ROOT / 'outputs/regression/caseI/point'
DATA = CASE / 'data.npz'
SVD = CASE / 'svd.npz'
METHODS = ('LS', 'SC', 'Lasso', 'Ridge', 'Landweber', 'Showalter', 'HBF', 'ARk', 'SOAR', 'Nesterov', 'FAR')
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lasso_grid(X, y, truth, label, least_squares):
    """LARS is initialization only; every grid point must pass cyclic-CD checks."""
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import Lasso, lars_path
    Xf = np.asfortranarray(X)
    gram = np.ascontiguousarray(X.T @ X)
    started = time.monotonic()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        knots, active, knot_coefficients = lars_path(Xf, y, Xy=X.T@y, Gram=gram,
            alpha_min=.001, method='lasso', max_iter=10000, return_path=True)
    path_warnings = [str(w.message) for w in caught]
    if knots[-1] > .001 + 1e-12:
        raise RuntimeError(f'LARS initialization did not reach grid lower bound: {knots[-1]}')
    print(json.dumps(dict(event='lasso_initialization_path', distribution=label,
                          knots=len(knots), seconds=time.monotonic()-started,
                          path_warnings=path_warnings)), flush=True)
    grid = .001*np.arange(1001)
    coefficients = np.empty((1001, X.shape[1]))
    coefficients[0] = least_squares
    errors = np.empty(1001)
    errors[0] = np.linalg.norm(least_squares-truth)
    dual_gaps = np.zeros(1001)
    iterations = np.zeros(1001, dtype=int)
    retry_counts = np.zeros(1001, dtype=int)
    energy = float(y@y / X.shape[0])
    allowed = 1e-6*energy
    roundoff = 100*np.finfo(float).eps*max(1., energy)
    model = Lasso(alpha=1., fit_intercept=False, tol=1e-6, warm_start=True,
                  selection='cyclic', precompute=gram)
    reversed_knots = knots[::-1]
    reversed_coefs = knot_coefficients[:, ::-1]
    checkpoint = time.monotonic()
    for j in range(1000, 0, -1):
        alpha = grid[j]
        pos = int(np.searchsorted(reversed_knots, alpha, side='right'))
        if pos == 0:
            initial = reversed_coefs[:, 0]
        elif pos == len(reversed_knots):
            initial = reversed_coefs[:, -1]
        else:
            fraction = (alpha-reversed_knots[pos-1])/(reversed_knots[pos]-reversed_knots[pos-1])
            initial = reversed_coefs[:, pos-1] + fraction*(reversed_coefs[:, pos]-reversed_coefs[:, pos-1])
        model.set_params(alpha=alpha)
        model.coef_ = initial.copy()
        converged = False
        for ri, budget in enumerate((10000, 30000, 100000)):
            model.set_params(max_iter=budget)
            with warnings.catch_warnings(record=True) as fit_warnings:
                warnings.simplefilter('always', ConvergenceWarning)
                model.fit(Xf, y)
            warning = any(issubclass(w.category, ConvergenceWarning) for w in fit_warnings)
            gap = float(model.dual_gap_)
            iterations[j] += int(model.n_iter_)
            if not np.isfinite(gap) or not np.isfinite(model.coef_).all() or gap < -roundoff:
                raise FloatingPointError('Invalid Lasso coefficients or dual gap')
            if not warning and gap <= allowed+roundoff:
                converged = True
                retry_counts[j] = ri
                break
        if not converged:
            raise RuntimeError(f'Lasso did not converge at alpha={alpha}')
        coefficients[j], dual_gaps[j] = model.coef_, gap
        errors[j] = np.linalg.norm(coefficients[j]-truth)
        if j % 100 == 0 or time.monotonic()-checkpoint >= 50:
            print(json.dumps(dict(event='lasso_grid_progress', distribution=label,
                                  completed=1001-j, total=1001,
                                  seconds=time.monotonic()-started)), flush=True)
            checkpoint = time.monotonic()
    warm_errors = errors.copy()
    warm_index = int(errors.argmin())
    verified = {}
    radius = 4
    while True:
        target = range(max(1, warm_index-radius), min(1000, warm_index+radius)+1)
        for j in target:
            if j in verified:
                continue
            cold = lasso_coefficients(X, y, grid[j])
            verified[j] = dict(alpha=float(grid[j]), cold_l2=float(np.linalg.norm(cold-truth)),
                               warm_l2=float(warm_errors[j]),
                               max_coefficient_difference=float(np.max(np.abs(cold-coefficients[j]))))
            coefficients[j] = cold
            errors[j] = verified[j]['cold_l2']
        best_index = int(errors.argmin())
        if best_index in verified and abs(best_index-warm_index) < radius:
            break
        radius *= 2
        if radius > 1000:
            raise RuntimeError('Lasso optimal-index neighborhood failed to stabilize')
    outside = np.ones(1001, dtype=bool)
    outside[list(verified)] = False
    outside_min = float(warm_errors[outside].min())
    observed_norm_discrepancy = max(abs(v['cold_l2']-v['warm_l2']) for v in verified.values())
    separation = outside_min-float(errors[best_index])
    if separation <= 10*observed_norm_discrepancy:
        # Expand cold checking to every numerically plausible candidate.
        plausible = np.flatnonzero(warm_errors <= errors[best_index]+10*observed_norm_discrepancy)
        for j in plausible:
            if j == 0 or j in verified:
                continue
            cold = lasso_coefficients(X, y, grid[j])
            verified[int(j)] = dict(alpha=float(grid[j]), cold_l2=float(np.linalg.norm(cold-truth)),
                                   warm_l2=float(warm_errors[j]),
                                   max_coefficient_difference=float(np.max(np.abs(cold-coefficients[j]))))
            coefficients[j], errors[j] = cold, verified[int(j)]['cold_l2']
        best_index = int(errors.argmin())
        if best_index not in verified:
            raise AssertionError('Final Lasso winner not cold-start validated')
    diagnostic = dict(grid=[0.,.001,1.], grid_size=1001, all_grid_candidates_converged=True,
                      convergence_tolerance=1e-6, allowed_dual_gap=allowed,
                      max_dual_gap=float(dual_gaps.max()), max_iterations=int(iterations.max()),
                      initialization='LARS exact-target path interpolated at each alpha; then converged cyclic CD',
                      precompute_gram=True, warm_best_alpha=float(grid[warm_index]),
                      selected_alpha=float(grid[best_index]), cold_start_checks=list(verified.values()),
                      outside_checked_window_warm_error_min=outside_min,
                      selected_l2=float(errors[best_index]), outside_window_error_separation=separation,
                      maximum_observed_cold_warm_l2_difference=observed_norm_discrepancy,
                      path_warning_messages=path_warnings, elapsed_seconds=time.monotonic()-started,
                      numerical_scope='Finite-grid numerical selection with complete converged path and cold validation near the optimum; not a symbolic certificate of identical floating-point CD histories.')
    np.savez(DEST / f'{label.lower()}_lasso_grid.npz', alpha=grid, coefficients=coefficients,
             l2=errors, warm_l2=warm_errors, dual_gap=dual_gaps, iterations=iterations,
             retry_counts=retry_counts, lars_knots=knots, lars_coefficients=knot_coefficients)
    (DEST / f'{label.lower()}_lasso_grid.json').write_text(json.dumps(diagnostic, indent=2), encoding='utf-8')
    return coefficients[best_index], float(grid[best_index]), diagnostic


def main(methods=METHODS):
    DEST.mkdir(parents=True, exist_ok=True)
    with np.load(DATA, allow_pickle=False) as archive:
        X, truth = archive['X'], archive['beta']
    hashes = dict(X_hash=_array_hash(X), beta_hash=_array_hash(truth))
    source_files = [Path(__file__), ROOT/'regression_common.py', ROOT/'regression_spectral.py',
                    ROOT/'regression_threshold_scan.py', ROOT/'regression_results.py', ROOT/'paper_regression_config.py']
    source_hashes = {str(path):sha(path) for path in source_files}
    source_directory=DEST/'sources'
    source_directory.mkdir(exist_ok=True)
    for source in source_files:
        shutil.copy2(source,source_directory/source.name)
    metadata = dict(version='matrix2-complete-paper-point-v1', status='running', **hashes,
                    noise_root=2026100802, index=0, n=1000, p=1000, rows=[],
                    noise_variance=4., fixed_before_evaluation=True, methods=list(methods),
                    noise_rule='default_rng((2026100802,0,0,0)) freshly initialized per distribution',
                    normal_noise=dict(distribution='Normal', mean=0., standard_deviation=2.),
                    laplace_noise=dict(distribution='Laplace', location=0., scale=float(np.sqrt(2))),
                    source_sha256=source_hashes,
                    limitations=['Fixed index 0, not a sample chosen for favorable outcomes.',
                                 'True beta is used in oracle rank/alpha/threshold selection.',
                                 'Iterative settings and threshold percentiles are frozen from earlier Normal calibration; no Laplace retuning.',
                                 'LS and Lasso have no debiased branch in this table, so II is not applicable.',
                                 'Point estimates contain no bootstrap coverage calculation.'])
    checks = []
    with threadpool_limits(limits=2):
        design = load_design(X, SVD)
        for distribution in ('Normal','Laplace'):
            rng = np.random.default_rng((2026100802,0,0,0))
            noise = rng.normal(0.,2.,1000) if distribution=='Normal' else rng.laplace(0.,np.sqrt(2),1000)
            response = X@truth+noise
            np.savez(DEST/f'{distribution.lower()}_data.npz', X=X, beta=truth, noise=noise, response=response)
            least_squares = np.linalg.lstsq(X,response,rcond=None)[0] if any(m in ('LS','Lasso') for m in methods) else None
            for method in methods:
                started = time.monotonic()
                applicable = method not in ('LS','Lasso')
                diagnostics = {}
                lasso_diagnostic = None
                if method == 'LS':
                    settings, percentile, parameter, parameter_kind = dict(selection='least_squares', rcond=None),50.,None,'none'
                    raw = np.column_stack((least_squares,least_squares))
                elif method == 'Lasso':
                    fit, parameter, lasso_diagnostic = lasso_grid(X,response,truth,distribution,least_squares)
                    settings, percentile, parameter_kind = dict(selection='oracle_grid',alpha_grid=[0.,.001,1.],fit_intercept=False,tol=1e-6),50.,'alpha'
                    raw = np.column_stack((fit,fit))
                else:
                    if method == 'SC':
                        settings, percentile, parameter_kind = dict(selection='oracle_rank',rank_max=1000),0.,'rank'
                    elif method == 'Ridge':
                        settings, percentile, parameter_kind = dict(selection='oracle_grid',C=200000,m=100,alpha_divisor=1000),20.,'alpha'
                    else:
                        settings, percentile, parameter_kind = SETTINGS[method].copy(),PERCENTILES[method],'steps'
                        settings['delta'] = float(settings['tau']*np.linalg.norm(noise))
                    ordinary,debiased,parameter,diagnostics = design.select(response,method,settings,truth)
                    raw = np.column_stack((ordinary,debiased))
                centers,thresholds,scan = threshold_grid_pair(raw,truth,[percentile])
                center,bn = centers[0],thresholds[0]
                reuse = dict(performed=False, reason='Recomputed with frozen paper settings.')
                raw_errors,errors = raw-truth[:,None],center-truth[:,None]
                raw_l2,raw_linf = np.linalg.norm(raw_errors,axis=0),np.max(np.abs(raw_errors),axis=0)
                l2,linf = np.linalg.norm(errors,axis=0),np.max(np.abs(errors),axis=0)
                sigma = np.linalg.norm(response[:,None]-X@center,axis=0)/np.sqrt(1000)
                row = dict(distribution=distribution,method=method,raw_l2_I=float(raw_l2[0]),
                    raw_l2_II=float(raw_l2[1]) if applicable else None,
                    l2_I=float(l2[0]),l2_II=float(l2[1]) if applicable else None,
                    threshold_I=float(bn[0]),threshold_II=float(bn[1]) if applicable else None,
                    parameter=int(parameter) if parameter_kind in ('rank','steps') else parameter,
                    parameter_kind=parameter_kind, steps=int(parameter) if parameter_kind=='steps' else None,
                    percentile=percentile,tau=settings.get('tau'),settings={k:v for k,v in settings.items() if k!='delta'},
                    raw_linf_I=float(raw_linf[0]),raw_linf_II=float(raw_linf[1]) if applicable else None,
                    linf_I=float(linf[0]),linf_II=float(linf[1]) if applicable else None,
                    sigma_I=float(sigma[0]),sigma_II=float(sigma[1]) if applicable else None,
                    variance_error_I=float(abs(sigma[0]**2-4)),variance_error_II=float(abs(sigma[1]**2-4)) if applicable else None,
                    debiased_applicable=applicable, tdp_satisfied=diagnostics.get('tdp_satisfied'),
                    at_iteration_cap=diagnostics.get('at_iteration_cap'), stopping_threshold=settings.get('delta'))
                filename=f'{distribution.lower()}_{method.lower()}'
                branches=2 if applicable else 1
                np.savez(DEST/(filename+'.npz'), raw_estimates=raw[:,:branches], estimates=center[:,:branches],
                         thresholds=bn[:branches], sigma=sigma[:branches], threshold_scan=scan,
                         noise=noise, response=response, beta=truth, parameter=np.asarray(parameter if parameter is not None else np.nan))
                check=dict(distribution=distribution,method=method,status='passed',
                           canonical_parameter_and_pair=True,canonical_threshold_grid=True,
                           normal_archive_reuse=reuse,lasso_full_grid=lasso_diagnostic,
                           tdp_diagnostics=diagnostics)
                checks.append(check)
                metadata['rows'].append(row)
                (DEST/(filename+'.json')).write_text(json.dumps(dict(row=row,validation=check),indent=2,allow_nan=False),encoding='utf-8')
                (DEST/'summary.json').write_text(json.dumps(metadata,indent=2,allow_nan=False),encoding='utf-8')
                print(json.dumps(dict(event='point_row_complete',distribution=distribution,method=method,
                                      parameter=row['parameter'],raw_l2=[row['raw_l2_I'],row['raw_l2_II']],
                                      l2=[row['l2_I'],row['l2_II']],seconds=time.monotonic()-started)),flush=True)
    expected_rows = 2*len(methods)
    if len(metadata['rows'])!=expected_rows:
        raise AssertionError('Incomplete point table')
    metadata['status']='complete'
    metadata['completed_rows']=expected_rows
    (DEST/'summary.json').write_text(json.dumps(metadata,indent=2,allow_nan=False),encoding='utf-8')
    with (DEST/'statistics_table.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(metadata['rows'][0]))
        writer.writeheader()
        for row in metadata['rows']:
            writer.writerow(dict(row,settings=json.dumps(row['settings'],sort_keys=True)))
    audit=dict(status='passed',rows=expected_rows,**hashes,noise_root=2026100802,index=0,noise_variance=4.,
               checks=checks,source_sha256=source_hashes,validation_tolerance=1e-10,
               no_existing_archives_modified=True)
    (DEST/'audit.json').write_text(json.dumps(audit,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(dict(event='point_table_complete',rows=expected_rows,path=str(DEST/'summary.json'))),flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=DATA)
    parser.add_argument('--svd', type=Path, default=SVD)
    parser.add_argument('--output', type=Path, default=DEST)
    parser.add_argument('--methods', nargs='+', choices=METHODS, default=list(METHODS))
    args = parser.parse_args()
    if len(set(args.methods)) != len(args.methods):
        parser.error('Methods must not repeat')
    DATA, SVD, DEST = args.data, args.svd, args.output
    main(args.methods)
