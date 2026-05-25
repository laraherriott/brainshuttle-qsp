# Combined SILK calibration: disease/healthy CSF and plasma in one 2x2 figure.
# Run from repo root:  python Calibration/realistic_SILK_combined.py
# Paths assume Calibration/ is the working directory.

import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from QSP_models.abeta42_and_40 import NoAbModel
from QSP_models.solution import Solution

random.seed(1)

# read in parameters for healthy simulations
parameters = pd.read_csv('../Parameters/healthy.csv',
                         header=None).values
parameters = [x[0] for x in parameters]

# =============================================================================
# CSF SILK — mean leucine input, 37 h (realistic_SILK + healthy_CSF_SILK)
# =============================================================================

# read in labeled leucine fractions for plasma and CSF in CSF SILK experiment
f_leu_plasma_csf = pd.read_csv("Data/mean_f_leu.csv", header=None).values
f_leu_plasma_csf = np.append([0], f_leu_plasma_csf)
f_leu_csf_csf = pd.read_csv("Data/mean_f_leu_csf.csv", header=None).values
f_leu_csf_csf = np.append([0], f_leu_csf_csf)

# CSF SILK simulation in disease model
CSF_model = NoAbModel(
    leucine=[f_leu_plasma_csf, f_leu_csf_csf],
    time_points=37,
)
CSF_solver = Solution(CSF_model, 0, 37, 1)
disease_csf_sol = CSF_solver.solve()
time_csf = list(disease_csf_sol.t)
# save labeled CSF concentration
disease_csf_leu_total = [
    x + y for x, y in zip(disease_csf_sol.y[11], disease_csf_sol.y[13])
]

# CSF SILK simulation in healthy model
CSF_healthy_model = NoAbModel(
    leucine=[f_leu_plasma_csf, f_leu_csf_csf],
    x=parameters,
    time_points=37,
    healthy=True,
)
CSF_healthy_solver = Solution(CSF_healthy_model, 0, 37, 1)
healthy_csf_sol = CSF_healthy_solver.solve()
# save labeled CSF concentration
healthy_csf_leu_total = [
    x + y for x, y in zip(healthy_csf_sol.y[11], healthy_csf_sol.y[13])
]

# CSF observed data
CSF_times = range(37)

# read in observed data from amyloid positive individuals
CSF_df = pd.read_csv("Data/mean_Ab42_CSF_SILK.csv")
CSF_vals_disease = CSF_df.iloc[:, 0]
error_bars = pd.read_csv("Data/csf_SILK_std.csv", header=None)
std = error_bars.iloc[:, 0]
std_adjusted = [((((x * 1000) * 1e-9) / 4514) / 1e-9) for x in std]
ci_95_disease_csf = [2 * x for x in std_adjusted]

# read in observed data from amyloid negative individuals
CSF_df_healthy = pd.read_csv("Data/mean_nM_Ab42_young_CSF_SILK.csv")
CSF_vals_healthy = CSF_df_healthy.iloc[:, 0]
error_bars_healthy = pd.read_csv("Data/SILK_sd_young.csv", header=None)
std_healthy = error_bars_healthy.iloc[:, 0]
std_adjusted_healthy = [((((x * 1000) * 1e-9) / 4514) / 1e-9) for x in std_healthy]
ci_95_healthy_csf = [2 * x for x in std_adjusted_healthy]

# =============================================================================
# Plasma SILK — bolus leucine input, 25 h (realistic_plasma_SILK + healthy_plasma_SILK)
# =============================================================================

# read in labeled leucine fractions for plasma and CSF in Plasma SILK experiment
f_leu_plasma = [
    x[0] for x in pd.read_csv("Data/bolus_800_plasma_new.csv", header=None).values
]
f_leu_csf = [
    x[0] for x in pd.read_csv("Data/bolus_800_csf_new.csv", header=None).values
]

# Plasma SILK simulation in disease model
Plasma_model = NoAbModel(leucine=[f_leu_plasma, f_leu_csf], time_points=25)
Plasma_solver = Solution(Plasma_model, 0, 25, 1)
disease_plasma_sol = Plasma_solver.solve()
time_plasma = list(disease_plasma_sol.t)
# save labeled plasma concentration
disease_plasma_leu_total = [
    x + y for x, y in zip(disease_plasma_sol.y[7], disease_plasma_sol.y[9])
]

# Plasma SILK simulation in healthy model
Plasma_healthy_model = NoAbModel(
    leucine=[f_leu_plasma, f_leu_csf],
    x=parameters,
    time_points=25,
    healthy=True,
)
Plasma_healthy_solver = Solution(Plasma_healthy_model, 0, 25, 1)
healthy_plasma_sol = Plasma_healthy_solver.solve()
# save labeled plasma concentration
healthy_plasma_leu_total = [
    x + y for x, y in zip(healthy_plasma_sol.y[7], healthy_plasma_sol.y[9])
]

# read in observed data from plasma SILK experiment
plasma_df = pd.read_csv("Data/plasma_silk_digitized.csv")
error_bars_plasma = pd.read_csv("Data/silk_ci_std.csv")
x_error = error_bars_plasma["time_plasma"]
ci_95 = error_bars_plasma[" ci_plasma"][:14]

# convert reported normalized MFLs to labeled plasma concentrations
plasma_vals_disease = [
    (x * 0.8799657293261557) * (6.936e-3 + 1.092e-3) for x in plasma_df.iloc[:, 1]
]
y_error_disease = plasma_vals_disease[1:]
ci_adjusted_disease = [
    (x * 0.8799657293261557) * (6.936e-3 + 1.092e-3) for x in ci_95
]
errors_disease = [y_error_disease[i] - ci_adjusted_disease[i] for i in range(14)]

plasma_vals_healthy = [
    (x * 0.8799657293261557) * (2.5357e-3 + 1.72238e-3) for x in plasma_df.iloc[:, 1]
]
y_error_healthy = plasma_vals_healthy[1:]
ci_adjusted_healthy = [
    (x * 0.8799657293261557) * (2.5357e-3 + 1.72238e-3) for x in ci_95
]
errors_healthy = [y_error_healthy[i] - ci_adjusted_healthy[i] for i in range(14)]

# =============================================================================
# 2x2 plot: A disease CSF, B healthy CSF, C disease plasma, D healthy plasma
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# A — disease CSF
axes[0, 0].plot(time_csf, disease_csf_leu_total, marker="x", label="Fitted (disease)")
axes[0, 0].errorbar(
    CSF_times, CSF_vals_disease, yerr=ci_95_disease_csf,
    label="Observed (disease)", fmt="o", capsize=5,
)
axes[0, 0].set_xlabel("Time, h")
axes[0, 0].set_ylabel("Labeled CSF amyloid, nM")
axes[0, 0].set_ylim(-0.0025, 0.04)
axes[0, 0].legend()
axes[0, 0].text(-0.15, 1.1, r"$\mathbf{A}$", transform=axes[0, 0].transAxes,
                fontsize=18, va="top", ha="left")

# B — healthy CSF
axes[0, 1].errorbar(
    CSF_times, CSF_vals_healthy, yerr=ci_95_healthy_csf,
    label="Observed (healthy)", fmt="o", capsize=5, color="tab:red",
)
axes[0, 1].plot(
    time_csf, healthy_csf_leu_total, marker="x",
    label="Fitted (healthy)", color="tab:green",
)
axes[0, 1].set_xlabel("Time, h")
axes[0, 1].set_ylabel("Labeled CSF amyloid, nM")
axes[0, 1].set_ylim(-0.0025, 0.04)
axes[0, 1].legend()
axes[0, 1].text(-0.15, 1.1, r"$\mathbf{B}$", transform=axes[0, 1].transAxes,
                fontsize=18, va="top", ha="left")

# C — disease plasma
axes[1, 0].plot(time_plasma, disease_plasma_leu_total, marker="x", label="Fitted (disease)")
axes[1, 0].errorbar(
    x_error[:14], y_error_disease, yerr=errors_disease,
    label="Observed (disease)", fmt="o", capsize=5, color="tab:orange",
)
axes[1, 0].set_xlabel("Time, h")
axes[1, 0].set_ylabel("Labeled plasma amyloid, nM")
axes[1, 0].set_ylim(0, 0.0008)
axes[1, 0].legend()
axes[1, 0].text(-0.15, 1.1, r"$\mathbf{C}$", transform=axes[1, 0].transAxes,
                fontsize=18, va="top", ha="left")

# D — healthy plasma
axes[1, 1].plot(
    time_plasma, healthy_plasma_leu_total, marker="x",
    label="Fitted (healthy)", color="tab:green",
)
axes[1, 1].errorbar(
    x_error[:14], y_error_healthy, yerr=errors_healthy,
    label="Observed (healthy)", fmt="o", capsize=5, color="tab:red",
)
axes[1, 1].set_xlabel("Time, h")
axes[1, 1].set_ylabel("Labeled plasma amyloid, nM")
axes[1, 1].set_ylim(0, 0.0008)
axes[1, 1].legend()
axes[1, 1].text(-0.15, 1.1, r"$\mathbf{D}$", transform=axes[1, 1].transAxes,
                fontsize=18, va="top", ha="left")

for ax in axes.flat:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)

fig.tight_layout()
fig.savefig("Plots/SILK_combined_2x2.png", dpi=300, bbox_inches="tight")
plt.show()
