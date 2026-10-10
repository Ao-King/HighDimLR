"""Audited exports and prespecified selection for the nine-method threshold study."""

import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np

import regression_results as rr


VERSION = "regression-threshold-export-v1"
PLAN_VERSION = "regression-threshold-selection-v1"
METHODS = (*rr.METHODS, "HBF", "FAR")
BASELINE_PERCENTILES = {**rr.PERCENTILES, "HBF": 50, "FAR": 50}
SEEDS = {"baseline": 20261006, "calibration": 20261007,
         "refinement": 20261007, "validation": 20261008}
COARSE = tuple(range(0, 101, 10))
MIDPOINTS = tuple(range(5, 100, 10))
TARGET, TOLERANCE = 0.95, 0.015
STATISTICS = (*rr.METRICS, "mean_radius_original", "mean_radius_debiased")
MCSE_NAMES = tuple("mcse_" + name.removeprefix("mean_") for name in STATISTICS)
COLUMNS = ("stage", "seed_root", "method", "percentile", "N", "B", *STATISTICS,
           *MCSE_NAMES, "coverage_count_original", "coverage_count_debiased")


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def _csv(path, columns, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _require_formal(config, stage):
    expected = {"n": 1000, "p": 1000, "N": 1000, "B": 500, "sigma": 2,
                "confidence": TARGET, "design_seed": 24, "noise_model": "normal",
                "threshold_grid": [0.001, 0.001, 1.0], "preprocess_blas_threads": 1,
                "seed_root": SEEDS[stage]}
    if any(config.get(key) != value for key, value in expected.items()):
        raise ValueError(f"Nonformal or wrong-seed configuration in {stage}")
    if config.get("method") not in METHODS:
        raise ValueError("Unknown threshold-study method")


def summarize_samples(config, percentiles, values):
    """Summarize N x Q x 10 statistics without retaining bootstrap coefficient arrays."""
    values = np.asarray(values, dtype=float)
    n = config["N"]
    if values.shape != (n, len(percentiles), len(STATISTICS)) or not np.all(np.isfinite(values)):
        raise ValueError("Invalid outer-sample statistics")
    if np.any(values < 0) or np.any((values[:, :, :2] != 0) & (values[:, :, :2] != 1)):
        raise ValueError("Coverage must be binary and errors/radii nonnegative")
    rows = []
    for index, percentile in enumerate(percentiles):
        samples = values[:, index]
        counts = samples[:, :2].sum(axis=0).astype(int)
        scale = np.max(samples, axis=0)
        scale[scale == 0] = 1
        normalized = samples / scale
        means = normalized.mean(axis=0) * scale
        means[:2] = counts / n
        errors = np.std(normalized, axis=0, ddof=1) / np.sqrt(n) * scale if n > 1 else np.full(len(STATISTICS), np.nan)
        errors[:2] = np.sqrt(means[:2] * (1 - means[:2]) / n)
        row = {"stage": config["stage"], "seed_root": config["seed_root"],
               "method": config["method"], "percentile": float(percentile), "N": n, "B": config["B"]}
        row.update({name: float(value) for name, value in zip(STATISTICS, means)})
        row.update({name: None if not np.isfinite(value) else float(value)
                    for name, value in zip(MCSE_NAMES, errors)})
        row.update(coverage_count_original=int(counts[0]), coverage_count_debiased=int(counts[1]))
        rows.append(row)
    return rows


def _load_legacy_baseline(directory):
    metadata = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
    if metadata.get("version") != rr.VERSION or metadata.get("status") != "complete" or (directory / ".active.lock").exists():
        raise ValueError("Legacy baseline must be an inactive complete v2 run")
    config = metadata["config"]
    _require_formal(config, "baseline")
    method = config["method"]
    if method not in rr.METHODS or config.get("threshold_percentile") != BASELINE_PERCENTILES[method]:
        raise ValueError("Legacy run is not a historical baseline percentile")
    source = json.loads((directory / "source.json").read_text(encoding="utf-8"))
    with np.load(directory / "data.npz", allow_pickle=False) as data:
        X, beta = data["X"], data["beta"]
        fingerprints = {"config": rr._digest(config), "source": rr._digest(source),
                        "X": rr._array_hash(X), "beta": rr._array_hash(beta)}
        if not rr.is_production_config(config, X):
            raise ValueError("Legacy baseline does not match historical method settings")
    if metadata["fingerprints"] != fingerprints or metadata["identity"] != rr._digest(fingerprints):
        raise ValueError("Legacy baseline fingerprint mismatch")
    records = rr._read_records(directory, config, beta, metadata["identity"])
    if len(records) != config["N"] or metadata.get("completed_replications") != config["N"]:
        raise ValueError("Incomplete legacy baseline")
    statistics = np.asarray([np.concatenate((record[0], record[3])) for record in records])[:, None, :]
    export_config = {**config, "stage": "baseline"}
    return metadata, summarize_samples(export_config, [config["threshold_percentile"]], statistics)


def _load_new_run(directory):
    from run_regression_threshold_scan import read_scan_records, validate_run
    from regression_threshold_scan import scan_method_settings

    metadata = validate_run(directory, require_complete=True)
    config = metadata["config"]
    stage = config.get("stage")
    if stage not in SEEDS:
        raise ValueError("Unrecognized threshold-study stage")
    _require_formal(config, stage)
    expected_settings = scan_method_settings(config["method"], config["n"], config["largest_singular_value"])
    if rr._digest(config["method_settings"]) != rr._digest(expected_settings):
        raise ValueError("Method settings differ from the prespecified scan configuration")
    percentiles = list(config["percentiles"])
    if len(set(percentiles)) != len(percentiles) or any(not 0 <= q <= 100 for q in percentiles):
        raise ValueError("Invalid or repeated percentiles")
    if stage == "calibration" and set(percentiles) != set(COARSE):
        raise ValueError("Calibration run must contain the entire coarse grid")
    if stage == "refinement" and not set(percentiles).issubset(MIDPOINTS):
        raise ValueError("Refinement can only add previously untested five-percent midpoints")
    if stage == "baseline" and (config["method"] not in ("HBF", "FAR") or percentiles != [50]):
        raise ValueError("New baseline runs are HBF/FAR at percentile 50")
    values = np.empty((config["N"], len(percentiles), len(STATISTICS)))
    completed = 0
    for completed, arrays in enumerate(read_scan_records(directory, metadata), start=1):
        if completed > config["N"]:
            raise ValueError("Unexpected extra outer replication")
        values[completed - 1] = np.concatenate((arrays["metrics"], arrays["quantiles"]), axis=1)
    if completed != config["N"]:
        raise ValueError("Incomplete threshold-study run")
    return metadata, summarize_samples(config, percentiles, values)


def _load_runs(directories, mode):
    loaded, rows, common, source_fingerprints = [], [], None, set()
    allowed = {"scan": {"calibration", "refinement"}, "baseline": {"baseline"},
               "validation": {"validation"}}[mode]
    for directory in map(lambda path: Path(path).resolve(), directories):
        metadata = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
        legacy = metadata.get("version") == rr.VERSION
        if legacy and mode != "baseline":
            raise ValueError("Historical baseline cannot be mixed into scan or independent validation")
        metadata, run_rows = _load_legacy_baseline(directory) if legacy else _load_new_run(directory)
        config = metadata["config"]
        stage = "baseline" if legacy else config["stage"]
        if stage not in allowed:
            raise ValueError("Do not combine baseline, calibration, and independent validation results")
        stat_keys = ("n", "p", "N", "B", "sigma", "confidence", "seed_root", "design_seed",
                     "noise_model", "threshold_grid", "preprocess_blas_threads")
        current = {key: config[key] for key in stat_keys}
        current.update({key: metadata["fingerprints"][key] for key in ("X", "beta")})
        if common is not None and current != common:
            raise ValueError("Runs must share statistical configuration, design, truth and seed")
        common = current
        source_fingerprints.add(metadata["fingerprints"]["source"])
        loaded.append({"directory": str(directory), "identity": metadata["identity"],
                       "metadata_sha256": _hash(directory / "metadata.json"),
                       "version": metadata["version"], "stage": stage, "method": config["method"],
                       "config": config, "fingerprints": metadata["fingerprints"]})
        rows.extend(run_rows)
    if not loaded or {row["method"] for row in rows} != set(METHODS):
        raise ValueError("Exactly the nine study methods are required")
    if mode != "baseline" and len(source_fingerprints) != 1:
        raise ValueError("Calibration/refinement or validation runs have mixed source versions")
    keys = [(row["method"], row["percentile"]) for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("Repeated method/percentile; never average duplicate scans")
    if mode == "baseline" and (len(rows) != 9 or any(row["percentile"] != BASELINE_PERCENTILES[row["method"]] for row in rows)):
        raise ValueError("Baseline must contain one original percentile per method")
    rows.sort(key=lambda row: (METHODS.index(row["method"]), row["percentile"]))
    return rows, loaded, common


def _gaps(row):
    n = Decimal(row["N"])
    target_count = Decimal(str(TARGET)) * n
    return tuple(abs(Decimal(row[name]) - target_count) / n
                 for name in ("coverage_count_original", "coverage_count_debiased"))


def _score(row, rule):
    first, second = _gaps(row)
    radius = row["mean_radius_debiased"] if rule == "bestII" else row["mean_radius_original"] + row["mean_radius_debiased"]
    primary, secondary = (second, max(first, second)) if rule == "bestII" else (max(first, second), second)
    return primary, secondary, radius, row["percentile"]


def choose_percentiles(rows):
    """Apply the fixed two selection rules; do not interpolate or assume monotonicity."""
    selections, refinement, monotonicity = {}, {}, {}
    for method in METHODS:
        candidates = sorted([row for row in rows if row["method"] == method], key=lambda row: row["percentile"])
        coarse = [row for row in candidates if row["stage"] == "calibration"]
        if {row["percentile"] for row in coarse} != set(COARSE):
            raise ValueError(f"Incomplete coarse grid for {method}")
        if any(row["seed_root"] != SEEDS["calibration"] for row in candidates):
            raise ValueError("Only the calibration seed may select percentiles")
        selections[method] = {}
        needs_refinement = False
        for rule in ("bestII", "bestBoth"):
            best_coarse = min(coarse, key=lambda row: _score(row, rule))
            needs_refinement |= _score(best_coarse, rule)[0] > Decimal(str(TOLERANCE))
            best = min(candidates, key=lambda row: _score(row, rule))
            first, second = _gaps(best)
            selections[method][rule] = {"percentile": best["percentile"],
                "coverage_original": best["coverage_original"], "coverage_debiased": best["coverage_debiased"],
                "gap_original": float(first), "gap_debiased": float(second),
                "mean_radius_original": best["mean_radius_original"], "mean_radius_debiased": best["mean_radius_debiased"],
                "meets_target": _score(best, rule)[0] <= Decimal(str(TOLERANCE)),
                "N": best["N"], "source_stage": best["stage"]}
        missing = sorted(set(MIDPOINTS) - {row["percentile"] for row in candidates})
        if needs_refinement and missing:
            refinement[method] = missing
        monotonicity[method] = {}
        for branch in ("original", "debiased"):
            differences = np.diff([row[f"coverage_count_{branch}"] for row in candidates])
            positive, negative = bool(np.any(differences > 0)), bool(np.any(differences < 0))
            label = "nonmonotone" if positive and negative else "nondecreasing" if positive else "nonincreasing" if negative else "constant"
            monotonicity[method][branch] = label
    return selections, refinement, monotonicity


def _selection_plan(rows, loaded, common):
    selections, refinement, monotonicity = choose_percentiles(rows)
    plan = {"version": PLAN_VERSION, "status": "needs_refinement" if refinement else "ready_for_validation",
            "target": TARGET, "tolerance": TOLERANCE,
            "calibration_seed_root": SEEDS["calibration"], "validation_seed_root": SEEDS["validation"],
            "calibration_configuration": common, "calibration_inputs": loaded,
            "selection_rules": {"bestII": ["absolute_CII_gap", "maximum_CI_CII_gap", "mean_radius_II", "smaller_q"],
                                "bestBoth": ["maximum_CI_CII_gap", "absolute_CII_gap", "sum_of_mean_radii", "smaller_q"]},
            "selections": selections, "refinement_percentiles": refinement,
            "observed_monotonicity": monotonicity, "interpolation_used": False,
            "validation_percentiles": {} if refinement else {method: sorted({choice["percentile"] for choice in rules.values()})
                                                             for method, rules in selections.items()}}
    return plan


def _check_validation_plan(plan, rows, loaded, common):
    if (plan.get("version") != PLAN_VERSION or plan.get("status") != "ready_for_validation"
            or plan.get("calibration_seed_root") != SEEDS["calibration"]
            or plan.get("validation_seed_root") != SEEDS["validation"]
            or common.get("seed_root") != SEEDS["validation"]
            or plan.get("target") != TARGET or plan.get("tolerance") != TOLERANCE
            or set(plan.get("validation_percentiles", {})) != set(METHODS)):
        raise ValueError("An independently seeded ready selection plan is required")
    calibration = plan["calibration_configuration"]
    if any(common[key] != value for key, value in calibration.items() if key != "seed_root"):
        raise ValueError("Validation must retain the calibration design and statistical settings")
    digest = rr._digest(plan)
    calibration_sources = {item["fingerprints"]["source"] for item in plan["calibration_inputs"]}
    if {item["fingerprints"]["source"] for item in loaded} != calibration_sources:
        raise ValueError("Validation must retain the calibration source version")
    for item in loaded:
        archived = Path(item["directory"]) / "selection_plan.json"
        if item["config"].get("selection_plan_sha256") != digest or not archived.is_file():
            raise ValueError("Validation selection-plan provenance mismatch")
        if rr._digest(json.loads(archived.read_text(encoding="utf-8"))) != digest:
            raise ValueError("Validation was run with a different selection plan")
    comparisons = []
    for method in METHODS:
        actual = {row["percentile"]: row for row in rows if row["method"] == method}
        expected = set(plan["validation_percentiles"][method])
        if (set(plan.get("selections", {}).get(method, {})) != {"bestII", "bestBoth"}
                or {item["percentile"] for item in plan["selections"][method].values()} != expected):
            raise ValueError("Selected-rule details do not match the frozen percentiles")
        if set(actual) != expected:
            raise ValueError("Validation percentiles do not match the prespecified choices")
        for rule in ("bestII", "bestBoth"):
            chosen = plan["selections"][method][rule]
            row = actual[chosen["percentile"]]
            first, second = _gaps(row)
            comparisons.append({"method": method, "selection_rule": rule, "percentile": chosen["percentile"],
                "calibration_coverage_original": chosen["coverage_original"],
                "calibration_coverage_debiased": chosen["coverage_debiased"],
                "calibration_meets_target": chosen["meets_target"],
                "validation_coverage_original": row["coverage_original"],
                "validation_coverage_debiased": row["coverage_debiased"],
                "validation_mcse_original": row["mcse_coverage_original"],
                "validation_mcse_debiased": row["mcse_coverage_debiased"],
                "validation_gap_original": float(first), "validation_gap_debiased": float(second),
                "validation_meets_target": _score(row, rule)[0] <= Decimal(str(TOLERANCE))})
    return comparisons


def _matrix_tex(rows, mode):
    headers = ["Method", "$q$", "Coverage I", "Coverage II", "$L_2$ I", "$L_2$ II",
               "$L_\\infty$ I", "$L_\\infty$ II", "$V$ I", "$V$ II", "$R$ I", "$R$ II"]
    lines = ["% Requires longtable. Coverage entries and MCSE are percentages; all other entries are means.",
             "% V = mean absolute variance error relative to sigma^2=4; R = mean bootstrap 95% Linf radius.",
             "% q is a percentile, not the coordinate threshold. HBF/FAR Gaussian bootstrap is a new empirical extension.",
             "% FAR order 1.5 is outside the manuscript's proved range 0 < order < 1.",
             r"\begin{longtable}{lr" + "r" * 10 + "}",
             " & ".join(headers) + r" \\ \hline", r"\endfirsthead",
             " & ".join(headers) + r" \\ \hline", r"\endhead"]
    for row in rows:
        fields = [r"$\mathrm{AR}^{\kappa}$" if row["method"] == "ARk" else row["method"], f"{row['percentile']:g}"]
        fields.extend(f"${100*row[name]:.1f} \\pm {100*row['mcse_'+name]:.2f}$" for name in rr.METRICS[:2])
        fields.extend(f"{row[name]:.5g}" for name in STATISTICS[2:])
        lines.append(" & ".join(fields) + r" \\")
    lines.extend([r"\hline", r"\end{longtable}", ""])
    return "\n".join(lines)


def _plot(rows, output, mode):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True, sharey=True, layout="constrained")
    for method, axis in zip(METHODS, axes.flat):
        items = [row for row in rows if row["method"] == method]
        q = [row["percentile"] for row in items]
        for branch, color, marker, label in (("original", "#0072b2", "o", "Coverage I"),
                                             ("debiased", "#d55e00", "s", "Coverage II")):
            axis.errorbar(q, [row[f"coverage_{branch}"] for row in items],
                          yerr=[1.96 * row[f"mcse_coverage_{branch}"] for row in items],
                          color=color, marker=marker, markersize=3, linewidth=1, capsize=2, label=label)
        axis.axhline(TARGET, color="black", linestyle="--", linewidth=0.8, label="Nominal 0.95")
        axis.set(title=r"AR$^\kappa$" if method == "ARk" else method, xlim=(-2, 102), ylim=(-0.04, 1.04))
        axis.set_xticks([0, 20, 40, 60, 80, 100])
        axis.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8, loc="lower left")
    figure.supxlabel("Threshold percentile q")
    figure.supylabel("Empirical coverage probability")
    figure.savefig(output / f"{mode}-coverage.pdf", bbox_inches="tight")
    figure.savefig(output / f"{mode}-coverage.png", dpi=160, bbox_inches="tight")
    plt.close(figure)


def export_threshold_scan(run_directories, output_root, *, mode="scan", plots=True, selection_plan=None):
    """Export a baseline, full scan, or independent validation without mixing stages."""
    if mode not in ("scan", "baseline", "validation"):
        raise ValueError("mode must be scan, baseline, or validation")
    rows, loaded, common = _load_runs(run_directories, mode)
    plan, comparisons = None, None
    if mode == "scan":
        plan = _selection_plan(rows, loaded, common)
    elif mode == "validation":
        if selection_plan is None:
            raise ValueError("Validation export requires its original selection plan")
        plan = json.loads(Path(selection_plan).read_text(encoding="utf-8")) if isinstance(selection_plan, (str, Path)) else selection_plan
        comparisons = _check_validation_plan(plan, rows, loaded, common)
    elif selection_plan is not None:
        raise ValueError("Baseline does not use a threshold-selection plan")
    root = Path(output_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix=f"{mode}-threshold-export-", dir=root))
    _csv(output / "aggregate.csv", COLUMNS, rows)
    (output / "aggregate.tex").write_text(_matrix_tex(rows, mode), encoding="utf-8")
    if mode == "scan":
        _json(output / "selection_report.json", plan)
        if plan["status"] == "ready_for_validation":
            _json(output / "selection_plan.json", plan)
        choices = [{"method": method, "selection_rule": rule, **choice}
                   for method, rules in plan["selections"].items() for rule, choice in rules.items()]
        _csv(output / "selected_percentiles.csv", list(choices[0]), choices)
    elif comparisons is not None:
        _csv(output / "validation_comparison.csv", list(comparisons[0]), comparisons)
        _json(output / "selection_plan.json", plan)
    if plots:
        _plot(rows, output, mode)
    manifest = {"version": VERSION, "mode": mode, "created_utc": datetime.now(timezone.utc).isoformat(),
                "status": "complete", "row_count": len(rows), "configuration": common,
                "inputs": loaded, "exporter_sha256": _hash(__file__),
                "mixed_historical_and_extended_sources": mode == "baseline",
                "notes": ["All rows and both branches are reported, including nonmonotone and poor coverage results.",
                          "HBF/FAR Gaussian regression bootstrap is a new empirical extension, not an original seven-method run.",
                          "FAR order 1.5 is outside the manuscript's proved 0 < order < 1 range.",
                          "Coverage MCSE is sqrt(p*(1-p)/N); zero MCSE at observed 0 or 1 does not imply zero uncertainty.",
                          "Other MCSE values use sample standard deviation divided by sqrt(N).",
                          "Figure error bars are 1.96 times the plug-in Monte Carlo SE, not simultaneous confidence bands.",
                          "Selecting q on calibration data does not establish nominal coverage; validation uses a separate seed."],
                "outputs": {path.name: _hash(path) for path in output.iterdir()}}
    if mode == "scan":
        manifest.update(refinement_percentiles=plan["refinement_percentiles"],
                        observed_monotonicity=plan["observed_monotonicity"])
    _json(output / "manifest.json", manifest)
    return {"directory": output, "aggregate_csv": output / "aggregate.csv", "aggregate_tex": output / "aggregate.tex",
            "manifest": output / "manifest.json",
            "selection_plan": output / "selection_plan.json" if (output / "selection_plan.json").exists() else None,
            "refinement_percentiles": plan["refinement_percentiles"] if mode == "scan" else {}}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+", type=Path)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--mode", choices=("scan", "baseline", "validation"), default="scan")
    parser.add_argument("--selection-plan", type=Path)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)
    print(json.dumps(export_threshold_scan(args.runs, args.output_root, mode=args.mode,
                                          plots=not args.no_plots, selection_plan=args.selection_plan), default=str))


if __name__ == "__main__":
    main()
