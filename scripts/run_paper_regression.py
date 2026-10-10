"""Recompute Table 3 using the fixed published regression settings."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from paper_regression_config import SETTINGS, PERCENTILES, NOISE_ROOT, load_design
from regression_results import METRICS, _array_hash, _atomic, _runtime, _statistics, _write_json
from regression_threshold_scan import threshold_grid_pair
from matrix2_bootstrap_helpers import bootstrap_errors_grid_fast
from run_matrix02_sc_paper_statistics_fast import rank_screen

METHODS = ('SC', 'Ridge', 'Landweber', 'Showalter', 'ARk', 'SOAR', 'Nesterov', 'HBF', 'FAR')


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=ROOT / 'data/regression/caseI/data.npz')
    parser.add_argument('--svd', type=Path, default=ROOT / 'data/regression/caseI/svd.npz')
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/regression/table3')
    parser.add_argument('--methods', nargs='+', choices=METHODS, default=list(METHODS))
    parser.add_argument('--N', type=int, default=1000)
    parser.add_argument('--B', type=int, default=500)
    parser.add_argument('--blas-threads', type=int, default=1)
    args = parser.parse_args(argv)
    if min(args.N, args.B, args.blas_threads) < 1:
        parser.error('N, B and blas-threads must be positive')
    if len(set(args.methods)) != len(args.methods):
        parser.error('Methods must not repeat')
    return args


def summarize(method, metadata, records):
    means = np.mean([r['metrics'] for r in records], axis=0)
    result = dict(metadata, status='complete', completed_replications=len(records))
    result.update(dict(zip(METRICS, map(float, means))))
    table = dict(method=method, percentile=PERCENTILES[method],
                 parameter_kind='rank' if method == 'SC' else 'alpha' if method == 'Ridge' else 'steps',
                 tau=SETTINGS[method].get('tau'),
                 mean_parameter=float(np.mean([r['parameter'] for r in records])))
    for branch, roman in enumerate(('I', 'II')):
        table['coverage_' + roman] = float(means[branch])
        table['l2_' + roman] = float(means[2 + branch])
        table['linf_' + roman] = float(means[4 + branch])
        table['variance_error_' + roman] = float(means[6 + branch])
        result['coverage_count_' + roman] = int(sum(r['coverage'][branch] for r in records))
        for field, source in (('threshold', 'thresholds'), ('sigma', 'sigma'), ('radius', 'quantiles'),
                              ('raw_l2', 'raw_l2'), ('raw_linf', 'raw_linf')):
            table[field + '_' + roman] = float(np.mean([r[source][branch] for r in records]))
    result['table_row'] = table
    result['iteration_cap_hits'] = int(sum(bool(r.get('at_iteration_cap', False)) for r in records))
    return result


def run(args):
    with np.load(args.data, allow_pickle=False) as saved:
        X, truth = saved['X'].copy(), saved['beta'].copy()
    if X.shape != (1000, 1000) or truth.shape != (1000,):
        raise ValueError('Table 3 requires the archived 1000 by 1000 design and truth')
    input_hashes = dict(X_hash=_array_hash(X), beta_hash=_array_hash(truth))
    if input_hashes['X_hash'] != 'eabdee9dcb96a8312ca7ea0897a56e0a5af8392f7e48eb2e2cdc4a304261e9ef':
        raise ValueError('The design does not match the published experiment')
    paths = [Path(__file__), ROOT / 'paper_regression_config.py',
             ROOT / 'scripts/matrix2_bootstrap_helpers.py',
             ROOT / 'scripts/run_matrix02_sc_paper_statistics_fast.py',
             ROOT / 'scripts/run_matrix02_sc_paper_statistics.py']
    paths += [ROOT / name for name in ('regression_spectral.py', 'regression_results.py',
              'regression_common.py', 'regression_execution.py', 'regression_bootstrap.py',
              'regression_threshold_scan.py')]
    source_hashes = {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in paths}
    config = dict(confidence=.95, sigma=2.)
    tables = []
    with threadpool_limits(limits=args.blas_threads, user_api='blas'):
        design = load_design(X, args.svd)
        signal = X @ truth
        ranks = None
        if 'SC' in args.methods:
            responses = np.column_stack([signal + np.random.default_rng((NOISE_ROOT, i, 0, 0)).normal(0, 2., X.shape[0])
                                         for i in range(args.N)])
            ranks, _ = rank_screen(design, responses, truth)
        for method in args.methods:
            started = time.perf_counter()
            output = args.output / method
            (output / 'replications').mkdir(parents=True, exist_ok=True)
            metadata = dict(version='paper-regression-fixed-v1', method=method,
                            N=args.N, B=args.B, noise_root=NOISE_ROOT, percentile=PERCENTILES[method],
                            settings=SETTINGS[method], backend='cpu', dtype='float64',
                            confidence=.95, sigma=2., **input_hashes,
                            source_sha256=source_hashes, runtime=_runtime())
            if (output / 'metadata.json').exists():
                previous = json.loads((output / 'metadata.json').read_text(encoding='utf-8'))
                for key in metadata.keys() - {'runtime'}:
                    if metadata[key] != previous[key]:
                        raise ValueError(f'Resume configuration changed: {method}/{key}')
            else:
                _write_json(output / 'metadata.json', metadata)
            records = []
            for index in range(args.N):
                destination = output / 'replications' / f'rep_{index:06d}.npz'
                if destination.exists():
                    with np.load(destination, allow_pickle=False) as saved:
                        row = {name: saved[name].copy() for name in saved.files}
                    if int(row['index']) != index or row['bootstrap_errors'].shape != (args.B, 2):
                        raise ValueError(f'Invalid saved replication: {destination}')
                else:
                    noise = np.random.default_rng((NOISE_ROOT, index, 0, 0)).normal(0, 2., X.shape[0])
                    response = signal + noise
                    settings = SETTINGS[method].copy()
                    if settings.get('selection') == 'TDP':
                        settings['delta'] = float(settings['tau'] * np.linalg.norm(noise))
                    if method == 'SC':
                        parameter = int(ranks[index])
                        ordinary, debiased = design.pair(response, method, parameter, settings)
                        diagnostics = {}
                    else:
                        ordinary, debiased, parameter, diagnostics = design.select(response, method, settings, truth)
                    raw = np.column_stack((ordinary, debiased))
                    centers, thresholds, _ = threshold_grid_pair(raw, truth, [PERCENTILES[method]])
                    center, threshold = centers[0], thresholds[0]
                    sigma = np.linalg.norm(response[:, None] - X @ center, axis=0) / np.sqrt(X.shape[0])
                    errors = bootstrap_errors_grid_fast(design, method, parameter, settings, centers,
                        sigma[None], thresholds, args.B, NOISE_ROOT, index,
                        batch_size=args.B, blas_threads=args.blas_threads)[0]
                    metrics, quantiles = _statistics(config, truth, center, sigma, errors)
                    row = dict(index=np.asarray(index), parameter=np.asarray(parameter),
                        raw_estimates=raw, estimates=center, thresholds=threshold, sigma=sigma,
                        bootstrap_errors=errors, metrics=metrics, quantiles=quantiles,
                        coverage=metrics[:2].astype(bool), raw_l2=np.linalg.norm(raw-truth[:, None], axis=0),
                        raw_linf=np.max(np.abs(raw-truth[:, None]), axis=0),
                        at_iteration_cap=np.asarray(diagnostics.get('at_iteration_cap', False)))
                    _atomic(destination, lambda stream: np.savez_compressed(stream, **row))
                metrics, quantiles = _statistics(config, truth, row['estimates'], row['sigma'], row['bootstrap_errors'])
                if not np.array_equal(metrics, row['metrics']) or not np.array_equal(quantiles, row['quantiles']):
                    raise ValueError(f'Saved metrics failed verification: {destination}')
                records.append(row)
                if (index + 1) % 100 == 0 or index + 1 == args.N:
                    print(f'{method}: {index+1}/{args.N}, {time.perf_counter()-started:.1f}s', flush=True)
            summary = summarize(method, metadata, records)
            summary['elapsed_seconds'] = time.perf_counter() - started
            _write_json(output / 'summary.json', summary)
            tables.append(summary['table_row'])
    _write_json(args.output / 'summary.json', dict(status='complete', N=args.N, B=args.B,
        backend='cpu', noise_root=NOISE_ROOT, **input_hashes, rows=tables))


if __name__ == '__main__':
    run(parse_args())
