"""Discrete regression filters in the SVD basis, without changing estimators.

The default formulas are algebraically equivalent to the notebook recurrences;
HBF also offers the paper's damped Stormer--Verlet discretization. Neither is
a continuous-time substitute. Reordered float64 arithmetic cannot guarantee
bitwise agreement at stopping, oracle-tie, or threshold boundaries.
"""

import numbers

import numpy as np
from scipy.special import gamma


_ALIASES = {"sc": "SC", "cut": "SC", "ridge": "Ridge",
            "landweber": "Landweber", "showalter": "Showalter",
            "ark": "ARk", "soar": "SOAR", "second": "SOAR",
            "nesterov": "Nesterov", "hbf": "HBF", "far": "FAR", "frac": "FAR"}


def _method(method):
    try:
        return _ALIASES[str(method).lower()]
    except KeyError as error:
        raise ValueError(f"Unknown spectral regression method: {method}") from error


def _count(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise ValueError(f"{name} must be a nonnegative integer")
    if not np.isfinite(value) or value < 0 or int(value) != value:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _finite(value, name):
    if not np.all(np.isfinite(value)):
        raise FloatingPointError(f"Nonfinite {name} in spectral regression")
    return value


def _real(value, name, positive=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise ValueError(f"{name} must be a finite real number")
    value = float(value)
    if not np.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"Invalid {name}")
    return value


class SpectralDesign:
    """Immutable float64 design and thin SVD; arithmetic is not bitwise legacy.

    A supplied ``svd=(U, s, Vt)`` must be an SVD of X. Full U/Vt matrices are
    accepted and trimmed to min(n, p); no singular directions are discarded.
    """

    def __init__(self, X, svd=None):
        if np.iscomplexobj(X):
            raise ValueError("X must be real")
        self.X = np.array(X, dtype=np.float64, copy=True)
        if self.X.ndim != 2 or min(self.X.shape) < 1:
            raise ValueError("X must be a nonempty matrix")
        _finite(self.X, "design")
        self.n, self.p = self.X.shape
        self.rank_size = min(self.n, self.p)
        if svd is None:
            U, s, Vt = np.linalg.svd(self.X, full_matrices=False)
        else:
            if len(svd) != 3 or any(np.iscomplexobj(part) for part in svd):
                raise ValueError("svd must contain real (U, s, Vt)")
            U, s, Vt = (np.asarray(part, dtype=np.float64) for part in svd)
        r = self.rank_size
        if (U.ndim != 2 or U.shape[0] != self.n or U.shape[1] < r
                or s.shape != (r,) or Vt.ndim != 2 or Vt.shape[0] < r
                or Vt.shape[1] != self.p or np.any(s < 0)
                or np.any(s[:-1] < s[1:])):
            raise ValueError("SVD dimensions/order do not match X")
        self.U = _finite(np.array(U[:, :r], dtype=np.float64, copy=True), "left singular vectors")
        self.s = _finite(np.array(s, dtype=np.float64, copy=True), "singular values")
        self.Vt = _finite(np.array(Vt[:r], dtype=np.float64, copy=True), "right singular vectors")
        with np.errstate(over="raise", invalid="raise"):
            self._eigenvalues = self.s * self.s
        for array in (self.X, self.U, self.s, self.Vt, self._eigenvalues):
            array.setflags(write=False)
        self._filters = {}
        self._trajectories = {}

    def _response(self, Y):
        if np.iscomplexobj(Y):
            raise ValueError("Y must be real")
        Y = np.asarray(Y, dtype=np.float64)
        if Y.ndim not in (1, 2) or Y.shape[0] != self.n:
            raise ValueError("Y must have one row per observation")
        return _finite(Y, "response")

    def _key(self, method, parameter, settings):
        if method == "SC":
            rank = _count(parameter, "SC rank")
            if rank > self.rank_size:
                raise ValueError("SC rank exceeds the thin SVD size")
            return method, rank
        if method == "Ridge":
            alpha = _real(parameter, "alpha")
            if alpha < 0:
                raise ValueError("alpha must be nonnegative")
            return method, alpha
        count = _count(parameter, "iteration count")
        step = _real(settings["step"], "step", positive=True)
        extra = ()
        if method == "ARk":
            extra = (_real(settings.get("kappa", .5), "kappa"),)
        elif method == "SOAR":
            value = _real(settings.get("s_star", .5), "s_star")
            if value <= -.5:
                raise ValueError("s_star must exceed -0.5")
            extra = (value,)
        elif method == "Nesterov":
            value = _real(settings.get("omega", 3), "omega")
            if value <= -1:
                raise ValueError("omega must exceed -1")
            extra = (value,)
        elif method == "HBF":
            eta = _real(settings.get("eta", 5), "eta")
            if eta < 0:
                raise ValueError("eta must be nonnegative")
            integrator = settings.get("integrator", "rk4")
            if integrator not in ("rk4", "stormer-verlet"):
                raise ValueError("HBF integrator must be rk4 or stormer-verlet")
            extra = (eta,)
            if integrator == "stormer-verlet":
                # Zhang--Hofmann (2020), Proposition 6.1, equation (46):
                # use the sufficient bound for their damped scheme (45).
                if eta <= 0:
                    raise ValueError("HBF stormer-verlet requires positive eta")
                slack = 1 + 8*np.finfo(np.float64).eps
                if (step*float(self.s[0]) > np.sqrt(2)*slack
                        or step*eta > 2*slack):
                    raise ValueError("HBF stormer-verlet step must satisfy "
                                     "step <= min(sqrt(2)/s_max, 2/eta)")
                # Preserve the legacy RK4 key, while isolating both caches.
                extra += (integrator,)
        elif method == "FAR":
            extra = (_real(settings.get("theta", 1.5), "theta", positive=True),)
        return method, count, step, *extra

    def _advance_hbf(self, state, h, eta):
        gain, velocity = state["gains"][-1], state["velocity"]
        eigenvalues = self._eigenvalues
        k1_b = velocity
        k1_v = self.s - eigenvalues*gain - eta*velocity
        k2_b = velocity + .5*h*k1_v
        k2_v = self.s - eigenvalues*(gain+.5*h*k1_b) - eta*k2_b
        k3_b = velocity + .5*h*k2_v
        k3_v = self.s - eigenvalues*(gain+.5*h*k2_b) - eta*k3_b
        k4_b = velocity + h*k3_v
        k4_v = self.s - eigenvalues*(gain+h*k3_b) - eta*k4_b
        next_gain = _finite(gain + h/6*(k1_b+2*k2_b+2*k3_b+k4_b), "HBF filter")
        next_velocity = _finite(velocity + h/6*(k1_v+2*k2_v+2*k3_v+k4_v), "HBF velocity")
        state["velocity"] = next_velocity
        return next_gain

    def _advance_hbf_verlet(self, state, h, eta):
        """Damped Stormer--Verlet, Zhang--Hofmann (2020), equation (45)."""
        gain, velocity = state["gains"][-1], state["velocity"]
        half_velocity = (velocity + .5*h*(self.s-self._eigenvalues*gain))/(1+.5*eta*h)
        next_gain = _finite(gain + h*half_velocity, "HBF filter")
        next_velocity = _finite(half_velocity + .5*h*(
            self.s-self._eigenvalues*next_gain-eta*half_velocity), "HBF velocity")
        state["velocity"] = next_velocity
        return next_gain

    def _advance_far(self, state, h, theta):
        # The notebook stores beta_0,...,beta_k and uses the residual at every
        # completed iterate in both the predictor and corrector history sums.
        k = len(state["gains"]) - 1
        residuals = state["residuals"]
        k1 = h**theta / (theta*(theta+1))
        k2 = h**theta / theta
        reciprocal_gamma = 1 / _finite(gamma(theta), "FAR gamma")
        predictor_sum = np.zeros_like(self.s)
        corrector_sum = np.zeros_like(self.s)
        for j, residual in enumerate(residuals):
            b_weight = k2*((k-j+1)**theta-(k-j)**theta)
            if j == 0:
                a_weight = k1*(k**(theta+1)-(k-theta)*(k+1)**theta)
            else:
                a_weight = k1*((k-j+2)**(theta+1)+(k-j)**(theta+1)
                              - 2*(k-j+1)**(theta+1))
            predictor_sum = predictor_sum + b_weight*residual
            corrector_sum = corrector_sum + a_weight*residual
        predictor = reciprocal_gamma*predictor_sum
        next_gain = _finite(reciprocal_gamma*(
            k1*(self.s-self._eigenvalues*predictor) + corrector_sum), "FAR filter")
        next_residual = _finite(self.s-self._eigenvalues*next_gain, "FAR residual")
        state["residuals"].append(next_residual)
        return next_gain

    def _cached_trajectory(self, method, count, settings):
        key = self._key(method, count, settings)
        trajectory_key = (method, *key[2:])
        if trajectory_key not in self._trajectories:
            zero = np.zeros_like(self.s)
            zero.setflags(write=False)
            self._trajectories[trajectory_key] = {
                "gains": [zero], "velocity": zero, "residuals": [self.s]}
        state = self._trajectories[trajectory_key]
        # Extend lazily: an early TDP hit must not evaluate a later unstable
        # iterate, and subsequent responses reuse all already completed steps.
        for k in range(1, count + 1):
            if len(state["gains"]) <= k:
                with np.errstate(over="raise", invalid="raise", divide="raise"):
                    if method == "HBF":
                        advance = (self._advance_hbf_verlet if len(key) == 5
                                   else self._advance_hbf)
                        gain = advance(state, key[2], key[3])
                    else:
                        gain = self._advance_far(state, key[2], key[3])
                gain.setflags(write=False)
                state["gains"].append(gain)
            yield state["gains"][k]

    def _trajectory(self, method, count, settings):
        if method in ("HBF", "FAR"):
            yield from self._cached_trajectory(method, count, settings)
            return
        key = self._key(method, count, settings)
        h = key[2]
        gain = np.zeros_like(self.s)
        previous = gain
        velocity = gain
        q = gain
        eigenvalues = self._eigenvalues
        for k in range(1, count + 1):
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                if method == "Landweber":
                    gain = gain + h * (self.s - eigenvalues * gain)
                elif method == "Showalter":
                    k1 = self.s - eigenvalues * gain
                    k2 = self.s - eigenvalues * (gain + .5*h*k1)
                    k3 = self.s - eigenvalues * (gain + .5*h*k2)
                    k4 = self.s - eigenvalues * (gain + h*k3)
                    gain = gain + h/6 * (k1 + 2*k2 + 2*k3 + k4)
                elif method == "ARk":
                    kappa = key[3]
                    previous = gain
                    gain = previous + k*h*velocity
                    # Preserve reference (51), including k*h and the old state.
                    inverse = 1 / (h*(k*h)**kappa * eigenvalues
                                   + 1 + ((k*h)**(-kappa)-kappa)/k)
                    velocity = inverse*velocity + h*inverse*(self.s-eigenvalues*previous)
                elif method == "SOAR":
                    s_star = key[3]
                    coefficient = 2*k/(2*k+2*s_star+1)
                    half_q = coefficient*q + coefficient*(h/2)*(self.s-eigenvalues*gain)
                    gain = gain + h*half_q
                    acceleration = (2*k-2*s_star+1)*(k+1)/(2*k+2*s_star+3)/k
                    extrapolated = gain + 2*h*acceleration*half_q
                    q = (2*k+1-2*s_star)/(2*k+2)*half_q + (h/2)*(self.s-eigenvalues*extrapolated)
                else:
                    omega = key[3]
                    momentum = 0.0 if k == 1 else (k-2)/(k-1+omega)
                    extrapolated = gain + momentum*(gain-previous)
                    previous = gain
                    gain = extrapolated + h*(self.s-eigenvalues*extrapolated)
            _finite(gain, "filter")
            _finite(velocity, "velocity")
            _finite(q, "SOAR momentum")
            yield gain

    def filter(self, method, parameter, settings):
        """Return discrete gain a; stopping/threshold boundaries are not bitwise.

        SC takes an integer rank, Ridge its actual alpha (not c), and iterative
        methods their completed update count. No response-dependent selection
        occurs here. Zero-alpha Ridge uses NumPy's default lstsq rank cutoff.
        HBF defaults to ``integrator="rk4"``; ``"stormer-verlet"`` selects
        Zhang--Hofmann's damped scheme (45), with positive eta and the
        sufficient step bound h <= min(sqrt(2)/s_max, 2/eta).
        """
        method = _method(method)
        key = self._key(method, parameter, settings)
        if key not in self._filters:
            if method == "SC":
                gain = np.zeros_like(self.s)
                rank = key[1]
                with np.errstate(over="raise", invalid="raise", divide="raise"):
                    np.divide(1.0, self.s[:rank], out=gain[:rank], where=self.s[:rank] != 0)
            elif method == "Ridge":
                alpha = key[1]
                if alpha == 0:
                    cutoff = np.finfo(np.float64).eps * max(self.n, self.p) * self.s[0]
                    gain = np.zeros_like(self.s)
                    np.divide(1.0, self.s, out=gain, where=self.s > cutoff)
                else:
                    with np.errstate(over="raise", invalid="raise", divide="raise"):
                        gain = self.s / (self._eigenvalues + alpha)
            else:
                gain = np.zeros_like(self.s)
                for gain in self._trajectory(method, key[1], settings):
                    pass
            self._filters[key] = _finite(gain.copy(), "filter")
        return self._filters[key].copy()

    def _pair_from_gain(self, Y, gain, identity_debias=False):
        projected = self.U.T @ Y
        with np.errstate(over="raise", invalid="raise"):
            ordinary_gain = gain[:, None] if Y.ndim == 2 else gain
            ordinary = _finite(self.Vt.T @ (ordinary_gain * projected), "estimate")
            # SC has G=L; avoid cancellation roundoff and keep branch storage separate.
            if identity_debias:
                return ordinary, ordinary.copy()
            debiased_gain = 2*gain - self.s*gain*gain
            if Y.ndim == 2:
                debiased_gain = debiased_gain[:, None]
            debiased = self.Vt.T @ (debiased_gain * projected)
        return ordinary, _finite(debiased, "debiased estimate")

    def pair(self, Y, method, parameter, settings):
        """Apply L and 2L-LXL to vector/batch Y; not bitwise at numeric boundaries."""
        method = _method(method)
        return self._pair_from_gain(self._response(Y), self.filter(method, parameter, settings),
                                    identity_debias=method == "SC")

    def select(self, Y, method, settings, truth):
        """Select the legacy parameter, then return its two fixed linear fits.

        Oracle SC scans descending ranks; Ridge scans c=0,m,...,C and returns
        alpha=c/alpha_divisor. Landweber defaults to first-turn AOSR and also
        accepts TDP. TDP uses the first residual hit starting at update 1,
        with explicit cap-zero support.
        Near oracle ties and stopping boundaries float64 order can change a
        discrete decision; algebraic equivalence is not bitwise equivalence.
        """
        Y = self._response(Y)
        if Y.ndim != 1:
            raise ValueError("Parameter selection requires one response vector")
        method = _method(method)
        expected_selection = {"SC": "oracle_rank", "Ridge": "oracle_grid",
                              "Landweber": "AOSR"}.get(method, "TDP")
        selection = settings.get("selection", expected_selection)
        if method == "Landweber":
            if selection not in ("AOSR", "TDP"):
                raise ValueError("Landweber selection must be AOSR or TDP")
            if selection == "TDP" and "delta" not in settings:
                raise ValueError("Landweber TDP selection requires delta")
        elif selection != expected_selection:
            raise ValueError(f"{method} selection must be {expected_selection}")
        if np.iscomplexobj(truth):
            raise ValueError("truth must be real")
        truth = np.asarray(truth, dtype=np.float64)
        if truth.shape != (self.p,):
            raise ValueError("truth must have one coefficient per column")
        _finite(truth, "truth")
        projected = self.U.T @ Y
        truth_coordinates = self.Vt @ truth
        truth_perpendicular = np.linalg.norm(truth-self.Vt.T@truth_coordinates)
        response_perpendicular = np.linalg.norm(Y-self.U@projected)

        def error_norm(gain):
            return float(_finite(np.hypot(np.linalg.norm(gain*projected-truth_coordinates),
                                          truth_perpendicular), "oracle error"))

        def residual_norm(gain):
            value = np.hypot(np.linalg.norm((1-self.s*gain)*projected), response_perpendicular)
            return float(_finite(value, "selection residual"))

        count = None
        cap = None
        delta = None
        gain = np.zeros_like(self.s)
        if method in ("SC", "Ridge"):
            if method == "SC":
                rank_max = _count(settings.get("rank_max", self.rank_size), "rank_max")
                if not 1 <= rank_max <= self.rank_size:
                    raise ValueError("rank_max must lie in [1, min(n,p)]")
                candidates = range(rank_max, 0, -1)
            else:
                C = _count(settings["C"], "C")
                stride = _count(settings["m"], "m")
                if stride == 0:
                    raise ValueError("m must be positive")
                divisor = _real(settings.get("alpha_divisor", self.n), "alpha_divisor", positive=True)
                candidates = (c/divisor for c in range(0, C+1, stride))
            minimum = np.inf
            for candidate in candidates:
                candidate_gain = self.filter(method, candidate, settings)
                error = error_norm(candidate_gain)
                if error < minimum:
                    minimum, parameter, gain = error, candidate, candidate_gain
            count = int(parameter) if method == "SC" else None
        else:
            cap = _count(settings["max_steps"], "max_steps")
            self._key(method, cap, settings)
            count = 0
            if method == "Landweber" and selection == "AOSR":
                minimum_steps = _count(settings.get("minimum_steps", 500), "minimum_steps")
                turning_found = False
                target = cap
                previous_error = None
                previous_gain = gain
                for count, candidate_gain in enumerate(self._trajectory(method, cap, settings), 1):
                    error = error_norm(candidate_gain)
                    if not turning_found and count >= 2 and error >= previous_error:
                        target = min(cap, max(count-1, minimum_steps))
                        turning_found = True
                        if target < count:
                            count, gain = count-1, previous_gain
                            break
                    gain = candidate_gain
                    if count >= target:
                        break
                    previous_error, previous_gain = error, gain
            else:
                delta = _real(settings["delta"], "delta")
                if delta < 0:
                    raise ValueError("delta must be nonnegative")
                for count, gain in enumerate(self._trajectory(method, cap, settings), 1):
                    residual = residual_norm(gain)
                    # At a close boundary, check the returned coefficient vector.
                    if abs(residual-delta) <= 1e-12*max(1.0, np.linalg.norm(Y), delta):
                        estimate = self.Vt.T@(gain*projected)
                        residual = float(_finite(np.linalg.norm(Y-self.X@estimate), "selection residual"))
                    if residual <= delta:
                        break
            parameter = count
            self._filters[self._key(method, parameter, settings)] = gain.copy()
        ordinary, debiased = self._pair_from_gain(Y, gain, identity_debias=method == "SC")
        residual = float(_finite(np.linalg.norm(Y-self.X@ordinary), "selected residual"))
        diagnostics = {
            "selected_count": count,
            "at_iteration_cap": cap is not None and count == cap,
            "residual_norm": residual,
            "tdp_satisfied": None if delta is None else bool(residual <= delta),
            "parameter_kind": "alpha" if method == "Ridge" else "rank" if method == "SC" else "steps",
        }
        if delta is not None:
            diagnostics["stopping_threshold"] = delta
        return ordinary, debiased, parameter, diagnostics
