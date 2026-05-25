#!/usr/bin/env python3
"""
Lecanemab CLARITY-dose simulation example
=========================================

This script reproduces a single treatment scenario from the paper
simulations notebook (``Examples/paper_simulations.ipynb``): lecanemab
10 mg/kg every two weeks (Q2W), the regimen used in the CLARITY-AD trial
(model key ``lec_CLARITY_10bw`` in ``drug_data.csv``).

It demonstrates the standard workflow for running the one-antibody QSP model:

1. Configure dose, body weight, molecular weight, schedule, and binding affinities.
2. Instantiate :class:`QSP_models.one_ab_model.OneAbModel`.
3. Integrate the ODEs with :class:`QSP_models.solution.Solution`.
4. Extract the brain plaque readout and plot simulation output.

Requirements
------------
Install the package from the repository root (recommended)::

    pip install -e .

Or ensure the repo root is on ``PYTHONPATH`` when running this file.

Usage
-----
From the repository root::

    python Examples/lecanemab_clarity_simulation.py

Optional: save the figure to a custom path::

    python Examples/lecanemab_clarity_simulation.py --output path/to/plot.png

Time units
----------
All model times are in **hours**. Dose intervals in ``drug_data.csv`` are stored
in **days** and are converted to hours (× 24) before building the model,
matching ``paper_simulations.ipynb``.

Plaque readout
--------------
Total brain plaque is the sum of free plaque and antibody-bound plaque
(ODE state variables 3 and 12). For plotting, this is expressed as percent
change from the model baseline burden (1300 nM), i.e. % fall in plaque (CBL),
as in the paper figures.
"""

import matplotlib.pyplot as plt
import numpy as np

from QSP_models.one_ab_model import OneAbModel
from QSP_models.solution import Solution

# ---------------------------------------------------------------------------
# CLARITY regimen (lec_CLARITY_10bw row in Examples/drug_data.csv)
# ---------------------------------------------------------------------------
DOSE_MG_PER_KG = 10.0
BODYWEIGHT_KG = 83.0
GPM = 147181.62  # grams per mole (lecanemab)
FINAL_DOSE_H = 13104  # last dose time (h); ~18 months of Q2W dosing
DOSE_INTERVAL_DAYS = 14  # Q2W: 14 days between doses in drug_data.csv
DOSE_INTERVAL_H = DOSE_INTERVAL_DAYS * 24
SIMULATION_TIME_H = 13104
STEP_SIZE_H = 1

# Binding affinities (nM): monomer, oligomer, plaque — passed via z= to OneAbModel
AFFINITIES_NM = [2300.0, 67.3, 1.8]

# Baseline brain plaque used for % CBL conversion (matches paper_simulations.ipynb)
BASELINE_PLAQUE_NM = 1300.0

# Plot styling
PLOT_COLOUR = "#C0251A"
HOURS_PER_MONTH = 728  # 2184 h = 3 months in paper tick marks


# Construct the lecanemab 10 mg/kg Q2W model used in paper simulations."""
model = OneAbModel(
        dose=DOSE_MG_PER_KG,
        bodyweight=BODYWEIGHT_KG,
        gpm=GPM,
        final_dose=FINAL_DOSE_H,
        dose_interval=DOSE_INTERVAL_H,
        z=AFFINITIES_NM,
    )


# Integrate the model and return time (hours) and % fall in brain plaque.
solver = Solution(model, t_0=0, t_end=SIMULATION_TIME_H, step_size=STEP_SIZE_H)
solution = solver.solve()

# Extract brain plaque result for plotting
total_brain_plaque = solution.y[3] + solution.y[12]
pct_fall_plaque = -((BASELINE_PLAQUE_NM - total_brain_plaque) / BASELINE_PLAQUE_NM) * 100
time_months = solution.t / HOURS_PER_MONTH

# plot simulation % fall in brain plaque
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(
    time_months,
    pct_fall_plaque,
    color=PLOT_COLOUR,
    label="Lecanemab 10 mg/kg Q2W (simulation)",
)

ax.set_xlabel("Time, months")
ax.set_ylabel("% Fall in brain plaque")
ax.set_title("Lecanemab CLARITY dose simulation")
ax.set_xlim(0, SIMULATION_TIME_H / HOURS_PER_MONTH)
ax.set_ylim(-110, 10)
ax.legend(fontsize=10)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
plt.show()
