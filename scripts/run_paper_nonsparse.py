"""Recompute the non-sparse point table with tau=0.6 and HBF eta=10."""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nonsparse_experiment import METHODS, select_grid, select_iterative
from paper_regression_config import NONSPARSE_SETTINGS
from regression_results import _runtime, _write_json
from regression_spectral import SpectralDesign
import run_nonsparse_lasso as lasso


def result_row(method, ordinary, debiased, parameter, diagnostics, data, elapsed):
    X, response, truth = (data[k] for k in ('X', 'Y', 'beta'))
    row = dict(rule='TDP', method=method, parameter=parameter,
               k0=int(parameter) if method not in ('LS', 'SC', 'Lasso', 'Ridge') else None,
               diagnostics=diagnostics, elapsed_seconds=elapsed)
    for branch, estimate in (('ordinary', ordinary), ('debiased', debiased)):
        if estimate is None:
            for field in ('l2', 'linf', 'residual', 'sigma'):
                row[branch + '_' + field] = None
            continue
        if estimate.shape != truth.shape or not np.all(np.isfinite(estimate)):
            raise FloatingPointError(f'Invalid estimate: {method}/{branch}')
        row[branch + '_l2'] = float(np.linalg.norm(estimate-truth))
        row[branch + '_linf'] = float(np.max(np.abs(estimate-truth)))
        row[branch + '_residual'] = float(np.linalg.norm(response-X @ estimate))
        row[branch + '_sigma'] = row[branch + '_residual'] / np.sqrt(X.shape[0])
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=ROOT / 'data/regression/nonsparse/data.npz')
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/regression/nonsparse')
    parser.add_argument('--methods', nargs='+', choices=METHODS, default=list(METHODS))
    args = parser.parse_args()
    if len(set(args.methods)) != len(args.methods):
        parser.error('Methods must not repeat')
    args.output.mkdir(parents=True, exist_ok=True)
    with np.load(args.data, allow_pickle=False) as saved:
        data = {name: saved[name].copy() for name in saved.files}
    rows, estimates, traces = [], {}, {}
    with threadpool_limits(limits=1, user_api='blas'):
        design = SpectralDesign(data['X'], svd=(data['U'], data['S'], data['Vt']))
        for method in args.methods:
            started = time.perf_counter()
            settings = NONSPARSE_SETTINGS[method].copy()
            if 'tau' in settings:
                settings['delta'] = float(settings['tau'] * np.linalg.norm(data['e']))
            if method == 'Lasso':
                lasso.DATA, lasso.DEST = args.data, args.output / 'Lasso'
                lasso.main()
                result = json.loads((lasso.DEST / 'result.json').read_text(encoding='utf-8'))
                row = dict(result, k0=None, debiased_l2=None, debiased_linf=None)
            else:
                if method == 'LS':
                    ordinary = np.linalg.lstsq(data['X'], data['Y'], rcond=None)[0]
                    debiased, parameter, diagnostics, trace = None, None, dict(selection='minimum_norm'), np.empty((0, 3))
                elif method in ('SC', 'Ridge'):
                    ordinary, debiased, parameter, diagnostics, trace = select_grid(
                        design, data['Y'], data['beta'], method, settings, 'TDP')
                else:
                    ordinary, debiased, parameter, diagnostics, trace = select_iterative(
                        design, data['Y'], data['beta'], method, settings, 'TDP')
                row = result_row(method, ordinary, debiased, parameter, diagnostics, data, time.perf_counter()-started)
                estimates[method + '__ordinary'] = ordinary
                if debiased is not None:
                    estimates[method + '__debiased'] = debiased
                traces[method] = trace
            row['settings'] = settings
            rows.append(row)
            print(json.dumps(row, allow_nan=False), flush=True)
    np.savez_compressed(args.output / 'estimates.npz', **estimates)
    np.savez_compressed(args.output / 'traces.npz', **traces)
    _write_json(args.output / 'summary.json', dict(status='complete', tau=.6, HBF_eta=10.,
        thresholding=False, bootstrap=False, backend='cpu', runtime=_runtime(), rows=rows))


if __name__ == '__main__':
    main()
