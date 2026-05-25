"""
Parameter sweep over antibody-dependent clearance rate (k_ADCP) and coclearance scale factor (s).

Runs the one-antibody QSP model on a 2D grid, storing simulated change from baseline in 
plaque level (CBL) and an error score comparing CBL to published reductions in amyloid PET signal
for two lecanemab doses. Results can be visualised with plotting_clearance_sweep.py.

Run from the Examples/ directory so CSV paths resolve correctly.
"""

import random
import pandas as pd
import numpy as np
from tqdm import tqdm

from QSP_models.one_ab_model import OneAbModel
from QSP_models.solution import Solution

# --- Load baseline PK parameters from prior model fit ---
# Column order in model1_PK_0_1.csv is rearranged to match OneAbModel state vector y
params = pd.read_csv("../Parameters/one_ab.csv", header=None).values
params = [x[0] for x in params]

# indices 7 and 8 are log10(k_ADCP) and scale (s); defaults above are the fitted values

random.seed(1)  # reproducible ODE solves if the model uses randomness internally

# --- Sweep grid (must match plotting_clearance_sweep.py) ---
range_k_adcp = np.linspace(-5, 0, 25)   # log10(k_ADCP)
range_scale = np.linspace(0, 50, 25)    # coclearancescale factor s

outputs = np.zeros((len(range_k_adcp), len(range_scale)))
error_scores = np.zeros((len(range_k_adcp), len(range_scale)))

for i, adcp in enumerate(tqdm(range_k_adcp)):
    for j, scale in enumerate(range_scale):
        params_use = params.copy()
        params_use[7] = adcp
        params_use[8] = scale

        # High-dose regimen (10 mg/kg Q2W) — primary lecanemab calibration arm
        model = OneAbModel(
            dose=10,
            bodyweight=83,
            gpm=147181.62,
            final_dose=24 * 364 * 1,
            dose_interval=14 * 24,
            y=params_use,
        )

        # Low-dose regimen (2.5 mg/kg Q2W) — second arm to better bound the error score
        model_2 = OneAbModel(
            dose=2.5,
            bodyweight=82,
            gpm=147181.62,
            final_dose=24 * 364 * 1,
            dose_interval=14 * 24,
            y=params_use,
        )

        # Integrate 1.5 years in 1 h steps (LSODA in Solution.solve)
        solver = Solution(model, 0, int(24 * 364 * 1.5), 1)
        solver_2 = Solution(model_2, 0, int(24 * 364 * 1.5), 1)
        solutions_plaque = solver.solve()
        solutions_plaque_2 = solver_2.solve()

        # Total brain plaque = soluble + insoluble compartments (states 3 and 12)
        brain_plaque = [
            x + y
            for x, y in zip(solutions_plaque.y[3], solutions_plaque.y[12])
        ]
        brain_plaque_2 = [
            x + y
            for x, y in zip(solutions_plaque_2.y[3], solutions_plaque_2.y[12])
        ]

        # CBL (% reduction from baseline centiloid load 1300) at protocol time points
        CBL_3months = ((1300 - brain_plaque[13 * 7 * 24]) / 1300) * 100
        CBL_6months = ((1300 - brain_plaque[25 * 7 * 24]) / 1300) * 100
        CBL_12months = ((1300 - brain_plaque[53 * 7 * 24]) / 1300) * 100
        CBL_end = ((1300 - brain_plaque[-1]) / 1300) * 100
        CBL_12months_2 = ((1300 - brain_plaque_2[53 * 7 * 24]) / 1300) * 100
        CBL_end_2 = ((1300 - brain_plaque_2[-1]) / 1300) * 100

        # Normalised squared error vs trial CBL data
        error = (
            ((CBL_3months - 28.7) / 28.7) ** 2
            + ((CBL_6months - 50.1) / 50.1) ** 2
            + ((CBL_12months - 72.9) / 72.9) ** 2
            + ((CBL_end - 81.7) / 81.7) ** 2
            + ((CBL_12months_2 - 19.7) / 19.7) ** 2
            + ((CBL_end_2 - 29.1) / 29.1) ** 2
        )
        error_scores[i, j] = error

        # Store combined end-of-study CBL
        outputs[i, j] = CBL_end + CBL_end_2

# Write results to CSV files
df_outputs = pd.DataFrame(outputs)
df_outputs.to_csv("Results/sweep_outputs_{}.csv".format(len(range_k_adcp)), header=False, index=False)

df_error = pd.DataFrame(error_scores)
df_error.to_csv("Results/sweep_error_scores_{}.csv".format(len(range_k_adcp)), header=False, index=False)
