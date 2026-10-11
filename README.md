# HighDimLR

Code and data for *On a class of high dimensional linear regression methods
with debiasing and thresholding*.

The commands below use the settings and saved inputs of the submitted tables
and figures. Run them from the repository root.

## Installation

Python 3.12 is supported. Install the CPU dependencies in an isolated environment:

```text
python -m pip install -r requirements.txt
```

`requirements-reference.txt` records the available package versions of the
archived non-sparse Python experiment. GPU execution of the shared regression bootstrap
is optional; install `requirements-gpu.txt` only when using a compatible CUDA
environment. The publication entry points use CPU double precision by default.

The PDE calculations require MATLAB and Partial Differential Equation Toolbox.
The current formal PDE experiments used MATLAB R2025a Update 1. The public entry points
and figure renderer were also checked with MATLAB R2022b. The saved PDE figure
data can be rendered without regenerating the finite element mesh.

## Reproduce the regression experiments

```text
python scripts/compute_matrix2_paper_points.py
python scripts/compute_matrix2_caseII_points.py
python scripts/run_paper_nonsparse.py
python scripts/run_paper_regression.py
```

The first two commands compute the Case I and Case II sparse point experiments.
The non-sparse experiment uses tau=0.6, SOAR s=1.5 and HBF eta=10. The Table 3 bootstrap entry
uses the nine published methods, N=1000 and B=500. The design, singular factors,
noise seeds, method settings and threshold percentiles are fixed by the saved
inputs and `paper_regression_config.py`.

The regression SOAR settings are:

| Experiment | s | tau | h | Iteration cap | Threshold percentile |
|---|---:|---:|---:|---:|---:|
| Sparse Cases I-II and Table 3 | 1.5 | 1.01 | 0.00025 | 10000 | 60 |
| Non-sparse point experiment | 1.5 | 0.6 | 0.0005 | 100000 | None |

The supplementary independent validation uses the same fixed configuration
and fresh outer and bootstrap noise. Run or export it separately:

```text
python scripts/run_paper_regression.py --experiment independent
python scripts/export_paper_regression.py --experiment independent
```

Its saved results use N=1000 and B=500. The noise seed is 2026101101 for SOAR
and 2026101802 for the other eight methods. Coverage Monte Carlo standard
errors are included in the independent-validation export.

For a short bootstrap calculation:

```text
python scripts/run_paper_regression.py --N 1 --B 8
```

The notebooks provide interactive versions of the methods. Their saved cell
outputs are cleared; use the commands above for the submitted experiment settings.

## Reproduce PDE Table 6

In MATLAB, set the current folder to this repository:

```matlab
results = run_paper_table6('all');
```

To run one method or a short calculation:

```matlab
result = run_paper_table6('SOAR');
result = run_paper_table6('SOAR',struct('N',1,'B',8));
```

The common settings are rt=2, N=1000, B=500, sigma=0.002, confidence=0.95 and
max_steps=100000. The seven `Bootstrap_*.m` entries use the same configuration.

| Method | tau | h | Other settings | Seed |
|---|---:|---:|---|---:|
| Landweber | 3.25 | 2 | | 20261011 |
| Showalter | 3.25 | 1 | | 20261012 |
| HBF | 3.2 | 0.5 | eta=5 | 20261013 |
| ARk | 4.1 | 0.725 | kappa=1.5 | 20261015 |
| SOAR | 3.4 | 0.25 | s=1.5, initial_time=1 | 20261009 |
| Nesterov | 3.3 | 1 | omega=5 | 20261010 |
| FAR | 3.35 | 0.5 | theta=1.25 | 20261014 |

The PDE bootstrap reselects the discrepancy stopping step in each resample.
Debiasing then replays that step count from zero. Regression bootstrap fits keep
the original-sample selected regression parameters and thresholds fixed.
FAR uses the validated full spectral implementation of its discrete recurrence.

## Tables and figures from the archived results

```text
python scripts/export_paper_regression.py
python scripts/export_table6.py
python plotting/plot_regression.py
python plotting/plot_pde.py
```

The PDE plotting command locates MATLAB automatically; `--matlab` can supply its
executable. The seven saved plot inputs reproduce the current PDE figures.
The regression figures read the saved point estimates and threshold scans.
These commands do not rerun the Monte Carlo calculations.

## Data

- `data/regression/caseI/` and `caseII/`: designs, singular factors and sparse point results.
- `data/regression/nonsparse/`: non-sparse input and current results with SOAR s=1.5 and HBF eta=10.
- `data/regression/table3/`: final method configurations and summaries.
- `data/regression/independent_validation/`: supplementary validation summaries, configurations and Monte Carlo standard errors.
- `data/pde/`: the common finite element operator, Table 6 statistics and original run-source hashes.
- `data/plotting/pde/`: the seven figure inputs.
- `data/input_checksums.json`: hashes of copied numerical inputs.

The repository includes the selected final statistics and figure inputs.
The full simulation records can be regenerated using the fixed configurations,
saved inputs and reproducible random streams in the publication entry points.

## Verified AMAX environment

The AMAX machine used for the latest formal PDE experiments was checked on
2026-10-10:

| Component | Configuration |
|---|---|
| Operating system | Ubuntu 20.04.6 LTS |
| CPU | Two Intel Xeon Platinum 8358P processors at 2.60 GHz |
| Cores | 64 physical cores, 128 logical threads |
| OS-visible memory | 503.5 GiB |
| GPUs | Two NVIDIA GeForce RTX 4090, 24564 MiB each |
| NVIDIA driver | 560.35.03 |
| MATLAB | R2025a Update 1 |
| Scientific Python | 3.12.2 |
| System orchestration Python | 3.8.10 |

The selected regression and PDE runs used double-precision CPU arithmetic.
The non-sparse regression and PDE runs used one computational thread; the
sparse and regression-bootstrap runs used no more than four BLAS threads per
process. The two GPUs were available but were not used by those runs.
The verified scientific Python
environment contains NumPy 1.26.4, SciPy 1.14.1, scikit-learn 1.4.2, matplotlib
3.8.4, Numba 0.59.1 and threadpoolctl 2.2.0. Its optional GPU environment adds
CuPy 13.6.0.

The original non-sparse Python records specify NumPy 1.26.4, SciPy 1.13.1,
scikit-learn 1.5.1 and threadpoolctl 3.5.0. These recorded versions are distinct
from the current AMAX installation. The saved designs and singular factors are
included to support comparisons across environments.
The updated non-sparse SOAR result was computed on AMAX with Python 3.12.2,
NumPy 1.26.4 and SciPy 1.14.1, using one BLAS thread.

## Validation

```text
python -m unittest discover -s tests -p "test_*.py"
```

```matlab
results = runtests('tests');
assertSuccess(results);
```

The release preparation checked all seven PDE methods against the first formal
AMAX record and its first four bootstrap draws. Stopping steps and random
substreams matched, and the maximum transformed-source difference was
2.36e-16. The PDE figure inputs and public renderer were checked against the
submitted figures. The table exporters and regression figure data were checked
against the current manuscript values.

Author: Wang Ying-Ao. Portions of the MATLAB finite element implementation
were adapted from work by Gong Rongfang.
