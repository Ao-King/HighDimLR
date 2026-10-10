"""Frozen numerical settings for the published regression experiments."""

from pathlib import Path
import numpy as np
from regression_spectral import SpectralDesign

ROOT = Path(__file__).resolve().parent
NOISE_ROOT = 2026100802
N = 1000
B = 500
SETTINGS = {'SC': {'selection': 'oracle_rank', 'rank_max': 1000}, 'Ridge': {'selection': 'oracle_grid', 'C': 200000, 'm': 100, 'alpha_divisor': 1000}, 'Landweber': {'max_steps': 10000, 'selection': 'TDP', 'step': 2.5e-06, 'tau': 0.9}, 'Showalter': {'max_steps': 10000, 'selection': 'TDP', 'step': 2.5e-06, 'tau': 0.9}, 'ARk': {'kappa': 0.5, 'max_steps': 1000, 'selection': 'TDP', 'step': 5e-05, 'tau': 0.9}, 'SOAR': {'max_steps': 10000, 's_star': 0.5, 'selection': 'TDP', 'step': 0.00025, 'tau': 1.1}, 'Nesterov': {'max_steps': 10000, 'nominal_step': 2.5e-06, 'omega': 5.0, 'selection': 'TDP', 'step': 2.0574059723851243e-06, 'tau': 0.9}, 'HBF': {'eta': 50.0, 'integrator': 'rk4', 'max_steps': 10000, 'selection': 'TDP', 'step': 0.0005, 'tau': 0.9}, 'FAR': {'max_steps': 500, 'selection': 'TDP', 'step': 5e-05, 'tau': 0.9, 'theta': 1.3}}
PERCENTILES = {'SC': 0.0, 'Ridge': 20.0, 'Landweber': 34.0, 'Showalter': 34.0, 'ARk': 35.0, 'SOAR': 53.0, 'Nesterov': 35.0, 'HBF': 32.0, 'FAR': 30.0}
NONSPARSE_SETTINGS = {'LS': {'selection': 'minimum_norm', 'rcond': None}, 'SC': {'selection': 'oracle_rank', 'rank_max': 1000, 'selection_exception': 'oracle_rank_reference'}, 'Lasso': {'selection': 'discrepancy_grid', 'C': 1000, 'm': 1, 'alpha_divisor': 1000, 'fit_intercept': False, 'tau': 0.6, 'delta': 37.398009755964544}, 'Ridge': {'selection': 'discrepancy_grid', 'C': 50000, 'm': 10, 'alpha_divisor': 1000, 'tau': 0.6, 'delta': 37.398009755964544}, 'Landweber': {'step': 3e-06, 'max_steps': 100000, 'selection': 'TDP', 'tau': 0.6, 'delta': 37.398009755964544}, 'Showalter': {'step': 3e-06, 'max_steps': 100000, 'selection': 'TDP', 'tau': 0.6, 'delta': 37.398009755964544}, 'HBF': {'step': 0.0005, 'max_steps': 100000, 'eta': 10.0, 'selection': 'TDP', 'tau': 0.6, 'delta': 37.398009755964544}, 'ARk': {'step': 5e-05, 'max_steps': 2000, 'kappa': 0.5, 'selection': 'TDP', 'tau': 0.6, 'delta': 37.398009755964544}, 'SOAR': {'step': 0.0005, 'max_steps': 100000, 's_star': 0.5, 'selection': 'TDP', 'tau': 0.6, 'delta': 37.398009755964544}, 'Nesterov': {'step': 2.0574059723851243e-06, 'nominal_step': 3e-06, 'max_steps': 100000, 'omega': 3, 'selection': 'TDP', 'tau': 0.6, 'delta': 37.398009755964544}, 'FAR': {'step': 0.0005, 'max_steps': 2000, 'theta': 1.8, 'selection': 'TDP', 'tau': 0.6, 'delta': 37.398009755964544}}

def load_design(X, svd_path=None):
    """Use an archived SVD when provided, or compute the same reduced SVD."""
    if svd_path is not None and Path(svd_path).is_file():
        with np.load(svd_path, allow_pickle=False) as saved:
            singular_values = saved['s'] if 's' in saved else saved['S']
            return SpectralDesign(X, svd=(saved['U'], singular_values, saved['Vt']))
    return SpectralDesign(X)
