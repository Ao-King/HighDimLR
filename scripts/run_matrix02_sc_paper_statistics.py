"""Complete canonical SC paper statistics on additional Gaussian matrix 02."""

from __future__ import annotations

from datetime import datetime, timezone
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from threadpoolctl import threadpool_limits

from regression_bootstrap import bootstrap_errors
from regression_results import METRICS, _array_hash, _atomic, _statistics, _write_json
from regression_spectral import SpectralDesign
from regression_threshold_scan import threshold_grid_pair
from matrix2_bootstrap_helpers import bootstrap_errors_grid_fast


DATA = ROOT / "data/regression/caseI/data.npz"
OUTPUT = ROOT / "outputs/regression/table3/SC"
EXPECTED_X_HASH = "eabdee9dcb96a8312ca7ea0897a56e0a5af8392f7e48eb2e2cdc4a304261e9ef"
NOISE_ROOT, N, B = 2026100802, 1000, 500


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def prepare(X, beta):
    assert _array_hash(X) == EXPECTED_X_HASH
    assert X.shape == (1000, 1000) and beta.shape == (1000,)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "replications").mkdir(exist_ok=True)
    paths = [Path(__file__), ROOT / "scripts/matrix2_bootstrap_helpers.py"]
    paths += [ROOT / name for name in ("regression_results.py", "regression_spectral.py",
        "regression_threshold_scan.py", "regression_bootstrap.py", "regression_execution.py", "regression_common.py")]
    hashes = {str(path.relative_to(ROOT)): sha(path) for path in paths}
    for path in paths:
        destination = OUTPUT / "sources" / path.relative_to(ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            assert sha(destination) == sha(path), f"Frozen source mismatch: {path}"
        else:
            shutil.copy2(path, destination)
    metadata = dict(version="matrix02-SC-paper-statistics-v1", method="SC", N=N, B=B,
        n=1000, p=1000, noise_root=NOISE_ROOT, confidence=.95, sigma=2., percentile=0.,
        threshold_percentile=0., X_hash=_array_hash(X), beta_hash=_array_hash(beta),
        data_origin=str(DATA), source_sha256=hashes,
        settings=dict(selection="oracle_rank", rank_max=1000), parameter_kind="rank",
        rank_rule="canonical ordinary raw oracle L2; ranks 1000 down to 1; first strict minimum",
        threshold_rule="canonical oracle L2 on .001:.001:1; q=0 first minimum",
        bootstrap_rule="500 independent original/debiased Gaussian draws with canonical (root,outer,1,b) RNG",
        debias_rule="SC idempotence: point estimators identical; bootstrap branch draws independent",
        status="running", created_at=now())
    metadata_path = OUTPUT / "metadata.json"
    if metadata_path.exists():
        old = read(metadata_path)
        for key in ("version", "N", "B", "noise_root", "X_hash", "beta_hash", "settings", "percentile", "source_sha256"):
            assert old[key] == metadata[key], f"Resume metadata mismatch: {key}"
        metadata = old
    else:
        _write_json(metadata_path, metadata)
        shutil.copy2(DATA, OUTPUT / "data.npz")
    with np.load(OUTPUT / "data.npz", allow_pickle=False) as saved:
        assert _array_hash(saved["X"]) == metadata["X_hash"]
        assert _array_hash(saved["beta"]) == metadata["beta_hash"]
    print(json.dumps(dict(initial_identity_check="passed", X_hash=metadata["X_hash"],
        beta_hash=metadata["beta_hash"], noise_root=NOISE_ROOT, N=N, B=B, percentile=0.,
        output=str(OUTPUT))), flush=True)
    return metadata


def summarize(metadata, rows, elapsed):
    means = np.mean([row["metrics"] for row in rows], axis=0)
    summary = dict(metadata)
    summary.update(status="complete", completed_replications=len(rows), completed_at=now(),
                   elapsed_seconds_this_session=elapsed, **dict(zip(METRICS, map(float, means))))
    for branch, label in enumerate(("original", "debiased")):
        summary["coverage_count_" + label] = int(sum(row["metrics"][branch] for row in rows))
        coverage = summary["coverage_" + label]
        summary["mcse_coverage_" + label] = float(np.sqrt(coverage * (1 - coverage) / N))
        for name, prefix in (("thresholds", "mean_threshold_"), ("sigma", "mean_sigma_"), ("quantiles", "mean_radius_")):
            summary[prefix + label] = float(np.mean([row[name][branch] for row in rows]))
        summary["raw_l2_" + label] = float(np.mean([row["raw_l2"][branch] for row in rows]))
        summary["raw_linf_" + label] = float(np.mean([row["raw_linf"][branch] for row in rows]))
    parameters = [int(row["parameter"]) for row in rows]
    summary["mean_parameter"] = float(np.mean(parameters))
    summary["parameter_range"] = [min(parameters), max(parameters)]
    table = dict(method="SC", percentile=0., mean_parameter=summary["mean_parameter"], tau=None, parameter_kind="rank")
    mapping = dict(coverage="coverage_", raw_l2="raw_l2_", raw_linf="raw_linf_",
        l2="l2_", linf="linf_", variance_error="variance_absolute_error_",
        sigma="mean_sigma_", threshold="mean_threshold_", radius="mean_radius_")
    for name, prefix in mapping.items():
        for roman, label in (("I", "original"), ("II", "debiased")):
            table[name + "_" + roman] = summary[prefix + label]
    summary["table_row"] = table
    _write_json(OUTPUT / "summary.json", summary)
    print(json.dumps(table), flush=True)
    return summary


def audit(metadata, summary, design, X, beta, signal):
    started = time.perf_counter()
    config = dict(confidence=.95, sigma=2.)
    maxima = {}
    checks = []
    rows = []
    for name, expected in metadata["source_sha256"].items():
        assert sha(OUTPUT / "sources" / name) == expected
        assert sha(ROOT / name) == expected
    with np.load(OUTPUT / "data.npz", allow_pickle=False) as data:
        assert _array_hash(data["X"]) == EXPECTED_X_HASH
        assert _array_hash(data["beta"]) == metadata["beta_hash"]
    for index in range(N):
        with np.load(OUTPUT / "replications" / f"rep_{index:06d}.npz", allow_pickle=False) as saved:
            row = {name: saved[name].copy() for name in saved.files}
        assert int(row["index"]) == index and str(row["X_hash"]) == EXPECTED_X_HASH
        assert row["bootstrap_errors"].shape == (B, 2)
        assert np.isfinite(row["bootstrap_errors"]).all()
        assert np.array_equal(row["estimates"][:, 0], row["estimates"][:, 1])
        assert np.array_equal(row["raw_estimates"][:, 0], row["raw_estimates"][:, 1])
        assert np.array_equal(row["sigma"][:1], row["sigma"][1:])
        assert np.array_equal(row["thresholds"][:1], row["thresholds"][1:])
        response = signal + np.random.default_rng((NOISE_ROOT, index, 0, 0)).normal(0, 2., X.shape[0])
        expected_centers, expected_thresholds, _ = threshold_grid_pair(row["raw_estimates"], beta, [0.])
        assert np.array_equal(row["estimates"], expected_centers[0])
        assert np.array_equal(row["thresholds"], expected_thresholds[0])
        sigma = np.linalg.norm(response[:, None] - X @ row["estimates"], axis=0) / np.sqrt(X.shape[0])
        assert np.array_equal(row["sigma"], sigma)
        metrics, quantiles = _statistics(config, beta, row["estimates"], sigma, row["bootstrap_errors"])
        assert np.array_equal(row["metrics"], metrics)
        assert np.array_equal(row["quantiles"], quantiles)
        assert np.array_equal(row["coverage"], metrics[:2].astype(bool))
        assert np.array_equal(row["raw_l2"], np.linalg.norm(row["raw_estimates"] - beta[:, None], axis=0))
        assert np.array_equal(row["raw_linf"], np.max(np.abs(row["raw_estimates"] - beta[:, None]), axis=0))
        ordinary, debiased = design.pair(response, "SC", int(row["parameter"]), metadata["settings"])
        raw_difference = float(np.max(np.abs(np.column_stack((ordinary, debiased)) - row["raw_estimates"])))
        assert raw_difference <= 1e-11
        maxima["all_raw_fit"] = max(maxima.get("all_raw_fit", 0.), raw_difference)
        if index in (0, 333, 999):
            ordinary, debiased, rank, _ = design.select(response, "SC", metadata["settings"], beta)
            assert int(rank) == int(row["parameter"])
            assert np.array_equal(row["raw_estimates"], np.column_stack((ordinary, debiased)))
            expected_errors = bootstrap_errors(design, "SC", rank, metadata["settings"], row["estimates"],
                sigma, row["thresholds"], B, NOISE_ROOT, index, backend="cpu", batch_size=B,
                cpu_workers=1, blas_threads=1)
            error = float(np.max(np.abs(expected_errors - row["bootstrap_errors"])))
            assert np.allclose(row["bootstrap_errors"], expected_errors, atol=1e-11, rtol=1e-11)
            maxima["canonical_bootstrap"] = max(maxima.get("canonical_bootstrap", 0.), error)
            expected_radius = np.quantile(expected_errors, .95, axis=0)
            assert np.array_equal(metrics[:2].astype(bool), np.max(np.abs(row["estimates"] - beta[:, None]), axis=0) <= expected_radius)
            checks.append(dict(index=index, rank=int(rank), raw_fit_exact=True, full_bootstrap_B=B,
                               max_bootstrap_difference=error, bootstrap_branches_independent=True))
        rows.append(row)
        if (index + 1) % 200 == 0:
            print(f"SC audit: {index+1}/{N}; {time.perf_counter()-started:.1f}s", flush=True)
    recomputed = dict(zip(METRICS, np.mean([row["metrics"] for row in rows], axis=0)))
    for name, value in recomputed.items():
        assert summary[name] == float(value)
    for branch, label in enumerate(("original", "debiased")):
        assert summary["coverage_count_" + label] == int(sum(row["metrics"][branch] for row in rows))
        for name, prefix in (("thresholds", "mean_threshold_"), ("sigma", "mean_sigma_"), ("quantiles", "mean_radius_")):
            assert summary[prefix + label] == float(np.mean([row[name][branch] for row in rows]))
        assert summary["raw_l2_" + label] == float(np.mean([row["raw_l2"][branch] for row in rows]))
        assert summary["raw_linf_" + label] == float(np.mean([row["raw_linf"][branch] for row in rows]))
    assert summary["mean_parameter"] == float(np.mean([row["parameter"] for row in rows]))
    report = dict(status="passed", method="SC", N=N, B=B, noise_root=NOISE_ROOT, X_hash=EXPECTED_X_HASH,
        beta_hash=metadata["beta_hash"], percentile=0., completed_replications=N,
        all_metrics_hash_threshold_sigma_rawfit_verified=True, frozen_sources_verified=True,
        canonical_checks=checks, maximum_absolute_differences=maxima,
        summary_sha256=sha(OUTPUT / "summary.json"), completed_at=now(), elapsed_seconds=time.perf_counter()-started)
    _write_json(OUTPUT / "audit.json", report)
    print(json.dumps(report), flush=True)


def run():
    started = time.perf_counter()
    with np.load(DATA, allow_pickle=False) as saved:
        X, beta = saved["X"].copy(), saved["beta"].copy()
    metadata = prepare(X, beta)
    rows = []
    config = dict(confidence=.95, sigma=2.)
    with threadpool_limits(limits=1, user_api="blas"):
        design = SpectralDesign(X)
        signal = X @ beta
        for index in range(N):
            path = OUTPUT / "replications" / f"rep_{index:06d}.npz"
            if path.exists():
                with np.load(path, allow_pickle=False) as saved:
                    row = {name: saved[name].copy() for name in saved.files}
                assert int(row["index"]) == index and str(row["X_hash"]) == EXPECTED_X_HASH
                assert row["bootstrap_errors"].shape == (B, 2)
            else:
                response = signal + np.random.default_rng((NOISE_ROOT, index, 0, 0)).normal(0, 2., X.shape[0])
                ordinary, debiased, rank, _ = design.select(response, "SC", metadata["settings"], beta)
                raw = np.column_stack((ordinary, debiased))
                centers, thresholds, _ = threshold_grid_pair(raw, beta, [0.])
                center, threshold = centers[0], thresholds[0]
                sigma = np.linalg.norm(response[:, None] - X @ center, axis=0) / np.sqrt(X.shape[0])
                errors = bootstrap_errors_grid_fast(design, "SC", rank, metadata["settings"], centers,
                    sigma[None], thresholds, B, NOISE_ROOT, index, batch_size=B, blas_threads=1)[0]
                metrics, quantiles = _statistics(config, beta, center, sigma, errors)
                row = dict(index=np.asarray(index), X_hash=np.asarray(EXPECTED_X_HASH), parameter=np.asarray(rank),
                    selected_parameter=np.asarray(rank), raw_estimates=raw, estimates=center, thresholds=threshold,
                    sigma=sigma, bootstrap_errors=errors, metrics=metrics, quantiles=quantiles,
                    coverage=metrics[:2].astype(bool), l2=metrics[2:4], linf=metrics[4:6],
                    raw_l2=np.linalg.norm(raw - beta[:, None], axis=0),
                    raw_linf=np.max(np.abs(raw - beta[:, None]), axis=0))
                _atomic(path, lambda stream: np.savez_compressed(stream, **row))
            rows.append(row)
            if (index + 1) % 100 == 0:
                print(f"SC complete: {index+1}/{N}; {time.perf_counter()-started:.1f}s", flush=True)
        summary = summarize(metadata, rows, time.perf_counter()-started)
        audit(metadata, summary, design, X, beta, signal)


def configure(argv=None):
    global DATA, OUTPUT, N, B
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=DATA)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--N', type=int, default=N)
    parser.add_argument('--B', type=int, default=B)
    args = parser.parse_args(argv)
    if args.N < 1 or args.B < 1:
        parser.error('N and B must be positive')
    DATA, OUTPUT, N, B = args.data, args.output, args.N, args.B


if __name__ == "__main__":
    configure()
    run()
