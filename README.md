# brainshuttle_QSP

A quantitative systems pharmacology (QSP) model for anti-amyloid immunotherapy.

The model includes amyloid production, aggregation, dynamics, antibody pharmacokinetics across brain, CSF, and plasma compartments, and antibody-dependent clearance of amyloid aggregates.  

This repository also documents data used for model calibration, in particular agains SILK and PK data, and reproduces all figures contained within the manuscript: ''

## Models

The installable package is `QSP_models`:

| Class | Role |
|-------|------|
| `NoAbModel` | Amyloid-only model (allowing for leucine-labeled SILK variants); disease or healthy parameter sets |
| `OneAbModel` | Single-antibody model on top of the amyloid backbone; dosing, pharmacokinetics, and antibody-dependent amyloid clearance|
| `Solution` | ODE integration wrapper (`scipy.integrate.solve_ivp`, LSODA) |
| `NoAbParameters`, `Lecanemab` | Default parameters for both models, including lecanemab-specific affinities |

Typical usage:

```python
from QSP_models import OneAbModel, Solution

model = OneAbModel(dose=10, bodyweight=83, gpm=147181.62,
                   final_dose=24 * 364, dose_interval=14 * 24)
solver = Solution(model, 0, int(24 * 364 * 1.5), 1)
solutions = solver.solve()
```

## Repository layout

```
brainshuttle_QSP/
├── QSP_models/          # Core ODE models and parameters
├── Calibration/         # SILK, steady-state, and PK calibration scripts
│   ├── Data/            # SILK and PK input data
│   └── Plots/           # Saved calibration figures
├── Examples/            # Example model simulation script and paper figures
│   ├── Data/            # PK, clearance, model inputs
│   └── Plots/           # Saved example outputs
├── Parameters/          # Fitted parameters for healthy base model and OneAbModel
└── Sensitivity/         # eFAST sensitivity, profile-likelihood, and clearance parameter sweeps
    ├── Results/         # sensitivity indices, profile-likelihood results, sweep results
    └── Plots/           # Saved analysis outputs
```

## Requirements

- Python ≥ 3.10
- Dependencies listed in `requirements.txt` (NumPy, SciPy, pandas, matplotlib, Plotly, Jupyter, etc.)

## Installation

From the repository root:

```bash
python -m venv QSP_env
source QSP_env/bin/activate
pip install -r requirements.txt
pip install -e .
```

This installs `QSP_models` in editable mode so local changes are picked up immediately.

## Running analyses

Run scripts from the directory noted in each file header (or from the repo root with paths adjusted). Jupyter notebooks should use that folder as the working directory.

### Calibration

| Script / notebook | Description |
|-------------------|-------------|
| `Calibration/SILK_simulations.py` | Disease and healthy CSF/plasma SILK fits (2×2 figure) |
| `Calibration/steady_state_simulations.py` | One-year steady-state amyloid trajectories, disease vs healthy (3×2 figure) |
| `Calibration/PK_profiles.ipynb` | Comparison of PK profiles for different anti-amyloid mAbs |
| `Calibration/Data/CSF_SILK_data.ipynb` | Process raw SILK cohort data |

### Examples

| Script / notebook | Description |
|-------------------|-------------|
| `Examples/paper_simulations.ipynb` | key paper figures: PK, brain:plasma ratios, CBL clearance vs trials |
| `Examples/lecanemab_clarity_simulation.py` | Example lecanemab CLARITY AD simulation |

### Sensitivity

| Script / notebook | Description |
|-------------------|-------------|
| `Sensitivity/sweep_clearance_params.py` | Grid search over `k_ADCP` and coclearance `scale` |
| `Sensitivity/plotting_clearance_sweep.py` | Heatmaps of sweep outputs |
| `Sensitivity/plotting_sensitivity.ipynb` | eFAST bar charts and profile-likelihood scans |

## Data notes

Summaries of the contents of Data/ are provided within each directory.
