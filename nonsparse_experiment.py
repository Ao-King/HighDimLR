"""Section 6.2 selection rules over the unchanged discrete spectral iterates."""

import numpy as np

from regression_spectral import _count, _finite, _real


METHODS = ("LS", "SC", "Lasso", "Ridge", "Landweber", "Showalter", "HBF",
           "ARk", "SOAR", "Nesterov", "FAR")
ITERATIVE_METHODS = METHODS[4:]
TRACE_COLUMNS = ("parameter", "coefficient_l2_error", "residual_norm")


def _rule(rule):
    if rule not in ("AOSR", "TDP"):
        raise ValueError("rule must be AOSR or TDP")
    return rule


def formal_settings(design, rule):
    """Return unscaled-X/Y settings; TDP callers supply the sample's delta.

    SC deliberately retains oracle rank tuning in the TDP comparison. Ridge
    helper alpha is c/n in X.T@X + alpha*I, not the normalized-theory alpha.
    """
    _rule(rule)
    maximum = float(design.s[0]**2)
    nesterov_step = min(3e-6, .99/maximum) if maximum > 0 else 3e-6
    grid_selection = "oracle_grid" if rule == "AOSR" else "discrepancy_grid"
    result = {
        "LS": {"selection": "minimum_norm", "rcond": None},
        "SC": {"selection": "oracle_rank", "rank_max": design.rank_size,
               "selection_exception": "oracle_rank_under_both_rules"},
        "Lasso": {"selection": grid_selection, "C": 1000, "m": 1,
                  "alpha_divisor": design.n, "fit_intercept": False},
        "Ridge": {"selection": grid_selection, "C": 50000, "m": 10,
                  "alpha_divisor": design.n},
        "Landweber": {"step": 3e-6, "max_steps": 100000},
        "Showalter": {"step": 3e-6, "max_steps": 100000},
        "HBF": {"step": 5e-4, "max_steps": 100000, "eta": 5.0},
        "ARk": {"step": 5e-5, "max_steps": 2000, "kappa": .5},
        "SOAR": {"step": 5e-4, "max_steps": 100000, "s_star": .5},
        "Nesterov": {"step": nesterov_step, "nominal_step": 3e-6,
                     "max_steps": 100000, "omega": 3},
        "FAR": {"step": 5e-4, "max_steps": 2000, "theta": 1.8},
    }
    for method in ITERATIVE_METHODS:
        result[method]["selection"] = rule
        if rule == "AOSR":
            result[method]["minimum_steps"] = 500
        else:
            result[method]["tau"] = .6
    if rule == "TDP":
        for method in ("Ridge", "Lasso"):
            result[method]["tau"] = .6
    return result


class _Geometry:
    def __init__(self, design, Y, truth):
        self.design = design
        self.Y = design._response(Y)
        if self.Y.ndim != 1:
            raise ValueError("Selection requires a single response vector")
        if np.iscomplexobj(truth):
            raise ValueError("truth must be real")
        self.truth = np.asarray(truth, dtype=np.float64)
        if self.truth.shape != (design.p,):
            raise ValueError("truth must have one value per design column")
        _finite(self.truth, "truth")
        with np.errstate(over="raise", invalid="raise"):
            self.projected = _finite(design.U.T @ self.Y, "response coordinates")
            self.truth_coordinates = _finite(design.Vt @ self.truth, "truth coordinates")
            self.truth_perpendicular = _finite(np.linalg.norm(self.truth-design.Vt.T@self.truth_coordinates), "truth null component")
            self.response_perpendicular = _finite(np.linalg.norm(self.Y-design.U@self.projected), "response null component")
        self.response_norm = float(_finite(np.linalg.norm(self.Y), "response norm"))

    def estimate(self, gain):
        with np.errstate(over="raise", invalid="raise"):
            return _finite(self.design.Vt.T @ (gain*self.projected), "coefficient estimate")

    def errors(self, gain):
        with np.errstate(over="raise", invalid="raise"):
            error = np.hypot(np.linalg.norm(gain*self.projected-self.truth_coordinates), self.truth_perpendicular)
            residual = np.hypot(np.linalg.norm((1-self.design.s*gain)*self.projected), self.response_perpendicular)
        return [float(_finite(error, "coefficient error")), float(_finite(residual, "residual"))]

    def direct_error(self, gain):
        return float(_finite(np.linalg.norm(self.estimate(gain)-self.truth), "direct coefficient error"))

    def direct_residual(self, gain):
        with np.errstate(over="raise", invalid="raise"):
            return float(_finite(np.linalg.norm(self.Y-self.design.X@self.estimate(gain)), "direct residual"))


def _near(a, b, *scales):
    return abs(a-b) <= 1e-12*max(1.0, abs(a), abs(b), *scales)


def _delta(settings):
    if "delta" not in settings:
        raise ValueError("TDP selection requires the sample-specific delta")
    value = _real(settings["delta"], "delta")
    if value < 0:
        raise ValueError("delta must be nonnegative")
    return value


def select_iterative(design, Y, truth, method, settings, rule):
    """Return two fits, selected count, diagnostics, and an N x 3 trace.

    AOSR records the triggering next iterate even if the previous one is
    returned. TDP checks each completed update from 1, never the zero iterate.
    """
    _rule(rule)
    if method not in ITERATIVE_METHODS:
        raise ValueError("Expected one of the seven iterative methods")
    if settings.get("selection", rule) != rule:
        raise ValueError("Method settings and requested selection rule disagree")
    if "max_steps" not in settings:
        raise ValueError("Missing max_steps")
    cap = _count(settings["max_steps"], "max_steps")
    if "step" not in settings:
        raise ValueError("Missing step")
    design._key(method, cap, settings)
    minimum = _count(settings.get("minimum_steps", 500), "minimum_steps") if rule == "AOSR" else None
    delta = _delta(settings) if rule == "TDP" else None
    geometry = _Geometry(design, Y, truth)
    gain = np.zeros_like(design.s)
    previous_gain, previous_error = gain, None
    selected = 0
    target, first_turn = cap, None
    trace, boundary_checks = [], []
    for count, candidate in enumerate(design._trajectory(method, cap, settings), 1):
        error, residual = geometry.errors(candidate)
        trace.append([float(count), error, residual])
        if rule == "AOSR" and first_turn is None and count >= 2:
            if _near(error, previous_error):
                error = geometry.direct_error(candidate)
                previous_error = geometry.direct_error(previous_gain)
                trace[-1][1], trace[-2][1] = error, previous_error
                boundary_checks.append({"kind": "coefficient_error", "steps": [count-1, count],
                                        "direct_values": [previous_error, error]})
            if error >= previous_error:
                first_turn = count-1
                target = min(cap, max(first_turn, minimum))
                if target < count:
                    selected, gain = count-1, previous_gain
                    break
        selected, gain = count, candidate
        if rule == "TDP":
            if _near(residual, delta, geometry.response_norm):
                residual = geometry.direct_residual(candidate)
                trace[-1][2] = residual
                boundary_checks.append({"kind": "residual_norm", "steps": [count],
                                        "direct_values": [residual]})
            if residual <= delta:
                break
        elif count >= target:
            break
        previous_gain, previous_error = candidate, error
    ordinary, debiased = design._pair_from_gain(geometry.Y, gain)
    replay_difference = None
    if method == "Landweber":
        twice_gain = np.zeros_like(gain)
        for twice_gain in design._trajectory(method, 2*selected, settings):
            pass
        twice_fit = geometry.estimate(twice_gain)
        replay_difference = float(np.max(np.abs(twice_fit-debiased), initial=0))
        if not np.allclose(twice_fit, debiased, rtol=1e-9, atol=1e-10):
            raise FloatingPointError("Landweber 2k replay disagrees with 2L-LXL")
        debiased = twice_fit
    residual = geometry.direct_residual(gain)
    trace_array = np.asarray(trace, dtype=np.float64).reshape(-1, 3)
    diagnostics = {"selection": rule, "parameter_kind": "steps", "selected_count": selected,
                   "first_turn_count": first_turn, "minimum_steps": minimum, "max_steps": cap,
                   "recorded_steps": len(trace), "trace_columns": list(TRACE_COLUMNS),
                   "trace_includes_trigger_after_selected": bool(len(trace) and trace[-1][0] > selected),
                   "at_iteration_cap": selected == cap, "residual_norm": residual,
                   "coefficient_l2_error": geometry.direct_error(gain),
                   "stopping_threshold": delta,
                   "tdp_satisfied": None if delta is None else bool(residual <= delta),
                   "boundary_checks": boundary_checks,
                   "landweber_2k_max_difference": replay_difference}
    return ordinary, debiased, selected, diagnostics, trace_array


def select_grid(design, Y, truth, method, settings, rule):
    """SC uses descending oracle ranks in both groups; Ridge uses its own grid."""
    _rule(rule)
    if method not in ("SC", "Ridge"):
        raise ValueError("select_grid handles SC and Ridge only")
    geometry = _Geometry(design, Y, truth)
    oracle = rule == "AOSR" or method == "SC"
    expected = "oracle_rank" if method == "SC" else "oracle_grid" if oracle else "discrepancy_grid"
    if settings.get("selection", expected) != expected:
        raise ValueError("Grid settings and requested selection rule disagree")
    delta = _delta(settings) if not oracle else None
    if method == "SC":
        rank_max = _count(settings.get("rank_max", design.rank_size), "rank_max")
        if not 1 <= rank_max <= design.rank_size:
            raise ValueError("rank_max must lie in [1,min(n,p)]")
        candidates = [(rank, None) for rank in range(rank_max, 0, -1)]
    else:
        if "C" not in settings or "m" not in settings:
            raise ValueError("Ridge requires C and m")
        maximum, stride = _count(settings["C"], "C"), _count(settings["m"], "m")
        if not stride:
            raise ValueError("m must be positive")
        divisor = _real(settings.get("alpha_divisor", design.n), "alpha_divisor", positive=True)
        numerators = range(0, maximum+1, stride) if oracle else range(maximum, -1, -stride)
        candidates = [(c/divisor, c) for c in numerators]
    trace, boundary_checks = [], []
    best_error, selected_gain, selected, selected_index, selected_c = np.inf, None, None, None, None
    for parameter, numerator in candidates:
        gain = design.filter(method, parameter, settings)
        error, residual = geometry.errors(gain)
        trace.append([float(parameter), error, residual])
        if oracle:
            if selected_gain is not None and _near(error, best_error):
                error, best_error = geometry.direct_error(gain), geometry.direct_error(selected_gain)
                trace[-1][1], trace[selected_index][1] = error, best_error
                boundary_checks.append({"kind": "coefficient_error", "parameters": [selected, parameter],
                                        "direct_values": [best_error, error]})
            choose = error < best_error
        else:
            if _near(residual, delta, geometry.response_norm):
                residual = geometry.direct_residual(gain)
                trace[-1][2] = residual
                boundary_checks.append({"kind": "residual_norm", "parameters": [parameter],
                                        "direct_values": [residual]})
            choose = selected is None and residual <= delta
        if choose:
            selected, selected_c, selected_gain = parameter, numerator, gain
            selected_index, best_error = len(trace)-1, error
    if selected is None:
        raise ValueError("No regularization-grid candidate satisfies the discrepancy rule")
    ordinary, debiased = design._pair_from_gain(geometry.Y, selected_gain, identity_debias=method == "SC")
    residual = geometry.direct_residual(selected_gain)
    diagnostics = {"selection": expected, "comparison_rule": rule,
                   "parameter_kind": "rank" if method == "SC" else "alpha",
                   "selected_count": int(selected) if method == "SC" else None,
                   "selected_grid_numerator": selected_c, "selected_trace_index": selected_index,
                   "recorded_candidates": len(trace), "trace_columns": list(TRACE_COLUMNS),
                   "residual_norm": residual, "coefficient_l2_error": geometry.direct_error(selected_gain),
                   "stopping_threshold": delta, "tdp_satisfied": None if oracle else bool(residual <= delta),
                   "selection_exception": "oracle_rank_under_both_rules" if method == "SC" else None,
                   "boundary_checks": boundary_checks}
    return ordinary, debiased, selected, diagnostics, np.asarray(trace, dtype=np.float64).reshape(-1, 3)
