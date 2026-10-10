"""Resume SC with equivalent batched spectral oracle-rank screening.

The original canonical runner/source and completed records remain archived.
Only the rank search is accelerated. Close rank candidates are reevaluated
using the identical scalar norm expression and descending tie order.
"""

from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from threadpoolctl import threadpool_limits

import run_matrix02_sc_paper_statistics as original
from regression_results import _atomic, _statistics, _write_json
from regression_spectral import SpectralDesign
from regression_threshold_scan import threshold_grid_pair
from matrix2_bootstrap_helpers import bootstrap_errors_grid_fast


def rank_screen(design, responses, truth):
    projected = design.U.T @ responses
    truth_coordinates = design.Vt @ truth
    truth_perp = np.linalg.norm(truth - design.Vt.T @ truth_coordinates)
    coordinates = np.zeros_like(projected)
    np.divide(projected, design.s[:, None], out=coordinates, where=design.s[:, None] != 0)
    # Every summed term is nonnegative, avoiding subtractive cancellation.
    selected = np.cumsum((coordinates - truth_coordinates[:, None])**2, axis=0)
    suffix = np.r_[np.cumsum(truth_coordinates[::-1]**2)[::-1], 0.]
    costs = selected + suffix[1:, None] + truth_perp**2
    minima = costs.min(axis=0)
    # Float64 summation and GEMM-vs-GEMV roundoff are dominated by this guard.
    guard = 128 * np.finfo(np.float64).eps * design.p * np.maximum(1., minima)
    ranks = np.empty(responses.shape[1], dtype=int)
    candidate_counts = []
    for index in range(responses.shape[1]):
        candidates = np.flatnonzero(costs[:, index] <= minima[index] + guard[index])[::-1] + 1
        canonical_projected = design.U.T @ responses[:, index]
        minimum = np.inf
        for rank in candidates:
            gain = design.filter("SC", int(rank), {})
            error = float(np.hypot(np.linalg.norm(gain * canonical_projected - truth_coordinates), truth_perp))
            if error < minimum:
                minimum, chosen = error, int(rank)
        ranks[index] = chosen
        candidate_counts.append(len(candidates))
    return ranks, candidate_counts


def run():
    started = time.perf_counter()
    with np.load(original.DATA, allow_pickle=False) as saved:
        X, beta = saved["X"].copy(), saved["beta"].copy()
    metadata = original.prepare(X, beta)
    source = Path(__file__).resolve()
    destination = original.OUTPUT / "sources/scripts" / source.name
    source_hash = original.sha(source)
    if destination.exists():
        assert original.sha(destination) == source_hash
    else:
        destination.write_bytes(source.read_bytes())
    metadata["execution_source_sha256"] = {str(source.relative_to(ROOT)): source_hash}
    metadata["rank_acceleration"] = "nonnegative prefix/suffix spectral oracle squared errors; close candidates reevaluated by canonical scalar norms, descending ties"
    _write_json(original.OUTPUT / "metadata.json", metadata)
    rows = []
    with threadpool_limits(limits=1, user_api="blas"):
        design = SpectralDesign(X)
        signal = X @ beta
        noise = np.column_stack([np.random.default_rng((original.NOISE_ROOT, index, 0, 0)).normal(0, 2., X.shape[0])
                                 for index in range(original.N)])
        responses = signal[:, None] + noise
        # Explicit per-sample addition agrees exactly with the original response.
        assert np.array_equal(responses[:, 0], signal + noise[:, 0])
        ranks, candidate_counts = rank_screen(design, responses, beta)
        previous_count = len(list((original.OUTPUT / "replications").glob("rep_*.npz")))
        history = dict(status="running", source_sha256=metadata["execution_source_sha256"],
            prior_canonical_replications=previous_count, rank_screen_candidate_count_range=[min(candidate_counts),max(candidate_counts)],
            rank_screen_candidate_count_mean=float(np.mean(candidate_counts)), created_at=original.now())
        _write_json(original.OUTPUT / "execution_history.json", history)
        print(json.dumps(dict(resuming=previous_count, fast_rank_screen="complete", **history)), flush=True)
        config = dict(confidence=.95, sigma=2.)
        for index in range(original.N):
            path = original.OUTPUT / "replications" / f"rep_{index:06d}.npz"
            if path.exists():
                with np.load(path, allow_pickle=False) as saved:
                    row = {name: saved[name].copy() for name in saved.files}
                assert int(row["index"]) == index and str(row["X_hash"]) == original.EXPECTED_X_HASH
                assert int(row["parameter"]) == int(ranks[index]), f"Canonical/fast rank mismatch at {index}"
                assert row["bootstrap_errors"].shape == (original.B, 2)
            else:
                response, rank = responses[:, index], int(ranks[index])
                ordinary, debiased = design.pair(response, "SC", rank, metadata["settings"])
                raw = np.column_stack((ordinary, debiased))
                centers, thresholds, _ = threshold_grid_pair(raw, beta, [0.])
                center, threshold = centers[0], thresholds[0]
                sigma = np.linalg.norm(response[:, None] - X @ center, axis=0) / np.sqrt(X.shape[0])
                errors = bootstrap_errors_grid_fast(design, "SC", rank, metadata["settings"], centers,
                    sigma[None], thresholds, original.B, original.NOISE_ROOT, index,
                    batch_size=original.B, blas_threads=1)[0]
                metrics, quantiles = _statistics(config, beta, center, sigma, errors)
                row = dict(index=np.asarray(index), X_hash=np.asarray(original.EXPECTED_X_HASH), parameter=np.asarray(rank),
                    selected_parameter=np.asarray(rank), raw_estimates=raw, estimates=center, thresholds=threshold,
                    sigma=sigma, bootstrap_errors=errors, metrics=metrics, quantiles=quantiles,
                    coverage=metrics[:2].astype(bool), l2=metrics[2:4], linf=metrics[4:6],
                    raw_l2=np.linalg.norm(raw - beta[:, None], axis=0),
                    raw_linf=np.max(np.abs(raw - beta[:, None]), axis=0))
                _atomic(path, lambda stream: np.savez_compressed(stream, **row))
            rows.append(row)
            if (index + 1) % 100 == 0:
                print(f"Fast SC: {index+1}/{original.N}; {time.perf_counter()-started:.1f}s", flush=True)
        summary = original.summarize(metadata, rows, time.perf_counter()-started)
        original.audit(metadata, summary, design, X, beta, signal)
        report = original.read(original.OUTPUT / "audit.json")
        assert original.sha(destination) == source_hash
        report["additional_execution_source_verified"] = metadata["execution_source_sha256"]
        report["prior_canonical_ranks_verified"] = previous_count
        _write_json(original.OUTPUT / "audit.json", report)
        history.update(status="complete", completed_at=original.now(), final_replications=original.N,
                       prior_canonical_ranks_verified=previous_count)
        _write_json(original.OUTPUT / "execution_history.json", history)


if __name__ == "__main__":
    original.configure()
    run()
