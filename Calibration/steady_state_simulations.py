# Steady-state no-antibody model: disease (MCI) vs healthy, 1-year simulation.
# Run from inside Calibration/

import matplotlib.pyplot as plt
import pandas as pd

from QSP_models.abeta42_and_40 import NoAbModel
from QSP_models.solution import Solution

simulation_length = 24 * 364
save_interval = 24

# read in parameters for healthy simulations
parameters = pd.read_csv('../Parameters/healthy.csv',
                         header=None).values
parameters = [x[0] for x in parameters]

# =============================================================================
# Disease — default NoAb parameters (no_ab_example)
# =============================================================================

disease_model = NoAbModel(time_points=simulation_length)
disease_solver = Solution(disease_model, 0, simulation_length, save_interval)
disease_sol = disease_solver.solve()
time = disease_sol.t

dis_brain_monomer = disease_sol.y[0]
dis_brain_oligomer = disease_sol.y[2]
dis_brain_plaque = disease_sol.y[4]
dis_plasma_monomer = disease_sol.y[6]
dis_plasma_oligomer = disease_sol.y[8]
dis_csf_monomer = disease_sol.y[10]
dis_csf_oligomer = disease_sol.y[12]

# =============================================================================
# Healthy — fitted parameters (baseline_healthy)
# =============================================================================

healthy_model = NoAbModel(x=parameters, time_points=simulation_length, healthy=True)
healthy_solver = Solution(healthy_model, 0, simulation_length, save_interval)
healthy_sol = healthy_solver.solve()

hel_brain_monomer = healthy_sol.y[0]
hel_brain_oligomer = healthy_sol.y[2]
hel_brain_plaque = healthy_sol.y[4]
hel_plasma_monomer = healthy_sol.y[6]
hel_plasma_oligomer = healthy_sol.y[8]
hel_csf_monomer = healthy_sol.y[10]
hel_csf_oligomer = healthy_sol.y[12]

# =============================================================================
# Print steady-state (endpoint) concentrations, nM
# =============================================================================

with open("Plots/steady_state_concentrations.txt", "w") as f:
    print("Disease (MCI) — endpoint concentrations (nM)", file=f)
    print("  Brain monomer:", dis_brain_monomer[-1], file=f)
    print("  Brain oligomer:", dis_brain_oligomer[-1], file=f)
    print("  Brain plaque:", dis_brain_plaque[-1], file=f)
    print("  Plasma monomer:", dis_plasma_monomer[-1], file=f)
    print("  Plasma oligomer:", dis_plasma_oligomer[-1], file=f)
    print("  CSF monomer:", dis_csf_monomer[-1], file=f)
    print("  CSF oligomer:", dis_csf_oligomer[-1], file=f)
    print(file=f)
    print("Healthy — endpoint concentrations (nM)", file=f)
    print("  Brain monomer:", hel_brain_monomer[-1], file=f)
    print("  Brain oligomer:", hel_brain_oligomer[-1], file=f)
    print("  Brain plaque:", hel_brain_plaque[-1], file=f)
    print("  Plasma monomer:", hel_plasma_monomer[-1], file=f)
    print("  Plasma oligomer:", hel_plasma_oligomer[-1], file=f)
    print("  CSF monomer:", hel_csf_monomer[-1], file=f)
    print("  CSF oligomer:", hel_csf_oligomer[-1], file=f)

print("Wrote Plots/steady_state_concentrations.txt")

# =============================================================================
# 3x2 figure: rows = brain / plasma / CSF, cols = disease / healthy
# =============================================================================

# month ticks shared by all panels (28-day months, hourly time axis)
no_months = 12
xticks = [i * (28 * 24) for i in range(0, no_months + 1, 3)]
xtick_labels = [i for i in range(0, no_months + 1, 3)]

fig, axes = plt.subplots(3, 2, figsize=(12, 12))

# A — disease brain
axes[0, 0].plot(time, dis_brain_monomer, label="Monomer")
axes[0, 0].plot(time, dis_brain_oligomer, label="Oligomer")
axes[0, 0].plot(time, dis_brain_plaque, label="Plaque")
axes[0, 0].set_xlabel("Time, months")
axes[0, 0].set_ylabel("Species (brain), nM")
axes[0, 0].set_xticks(xticks)
axes[0, 0].set_xticklabels(xtick_labels)
axes[0, 0].legend()
axes[0, 0].text(0.01, 0.99, r"$\mathbf{A}$", transform=axes[0, 0].transAxes,
                fontsize=18, va="top", ha="left")

# B — healthy brain
axes[0, 1].plot(time, hel_brain_monomer, label="Monomer")
axes[0, 1].plot(time, hel_brain_oligomer, label="Oligomer")
axes[0, 1].plot(time, hel_brain_plaque, label="Plaque")
axes[0, 1].set_xlabel("Time, months")
axes[0, 1].set_ylabel("Species (brain), nM")
axes[0, 1].set_xticks(xticks)
axes[0, 1].set_xticklabels(xtick_labels)
axes[0, 1].legend()
axes[0, 1].text(0.01, 0.99, r"$\mathbf{B}$", transform=axes[0, 1].transAxes,
                fontsize=18, va="top", ha="left")

# C — disease plasma
axes[1, 0].plot(time, dis_plasma_monomer, label="Monomer")
axes[1, 0].plot(time, dis_plasma_oligomer, label="Oligomer")
axes[1, 0].set_xlabel("Time, months")
axes[1, 0].set_ylabel("Species (plasma), nM")
axes[1, 0].set_xticks(xticks)
axes[1, 0].set_xticklabels(xtick_labels)
axes[1, 0].legend()
axes[1, 0].text(0.01, 0.99, r"$\mathbf{C}$", transform=axes[1, 0].transAxes,
                fontsize=18, va="top", ha="left")

# D — healthy plasma
axes[1, 1].plot(time, hel_plasma_monomer, label="Monomer")
axes[1, 1].plot(time, hel_plasma_oligomer, label="Oligomer")
axes[1, 1].set_xlabel("Time, months")
axes[1, 1].set_ylabel("Species (plasma), nM")
axes[1, 1].set_xticks(xticks)
axes[1, 1].set_xticklabels(xtick_labels)
axes[1, 1].legend()
axes[1, 1].text(0.01, 0.99, r"$\mathbf{D}$", transform=axes[1, 1].transAxes,
                fontsize=18, va="top", ha="left")

# E — disease CSF
axes[2, 0].plot(time, dis_csf_monomer, label="Monomer")
axes[2, 0].plot(time, dis_csf_oligomer, label="Oligomer")
axes[2, 0].set_xlabel("Time, months")
axes[2, 0].set_ylabel("Species (CSF), nM")
axes[2, 0].set_xticks(xticks)
axes[2, 0].set_xticklabels(xtick_labels)
axes[2, 0].legend()
axes[2, 0].text(0.01, 0.99, r"$\mathbf{E}$", transform=axes[2, 0].transAxes,
                fontsize=18, va="top", ha="left")

# F — healthy CSF
axes[2, 1].plot(time, hel_csf_monomer, label="Monomer")
axes[2, 1].plot(time, hel_csf_oligomer, label="Oligomer")
axes[2, 1].set_xlabel("Time, months")
axes[2, 1].set_ylabel("Species (CSF), nM")
axes[2, 1].set_xticks(xticks)
axes[2, 1].set_xticklabels(xtick_labels)
axes[2, 1].legend()
axes[2, 1].text(0.01, 0.99, r"$\mathbf{F}$", transform=axes[2, 1].transAxes,
                fontsize=18, va="top", ha="left")

fig.tight_layout()
fig.savefig("Plots/steady_state_combined_3x2.png", dpi=300, bbox_inches="tight")
plt.show()
