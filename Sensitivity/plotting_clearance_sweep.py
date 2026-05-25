"""
Visualise clearance-parameter sweep results from sweep_clearance_params.py.

Reads sweep_outputs_<N>.csv and sweep_error_scores_<N>.csv (N = grid resolution)
and plots side-by-side heatmaps of combined end-of-study CBL and error score.

Run from Examples/ after the sweep script has finished.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Grid definition (must match sweep_clearance_params.py) ---
range_k_adcp = np.linspace(-5, 0, 25)
range_scale = np.linspace(0, 50, 25)
grid_n = len(range_k_adcp)

outputs = pd.read_csv(f"Results/sweep_outputs_{grid_n}.csv", header=None).values
error_scores = pd.read_csv(f"Results/sweep_error_scores_{grid_n}.csv", header=None).values

# Fitted parameter values to overlay on the error heatmap (axis coords: s, log10 k_ADCP)
scale_target = 10.5
k_adcp_target = -1.2

fig, axs = plt.subplots(1, 2, figsize=(18, 8))

# Left: combined end-of-study CBL (%) for 10 mg/kg + 2.5 mg/kg lecanemab
im0 = axs[0].imshow(
    outputs,
    extent=[range_scale[0], range_scale[-1], range_k_adcp[0], range_k_adcp[-1]],
    origin="lower",
    aspect="auto",
    cmap="viridis",
)
axs[0].set_xlabel("s", fontsize=12)
axs[0].set_ylabel(r"$\log_{10}$(k$_{\mathrm{ADCP}}$)", fontsize=12)
axs[0].text(-4, 0, "A", fontsize=18, fontweight="bold", color="black")
plt.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04, label="CBL")

# Right: normalised squared error vs CLARITY AD CBL targets (lower is better)
im1 = axs[1].imshow(
    error_scores,
    extent=[range_scale[0], range_scale[-1], range_k_adcp[0], range_k_adcp[-1]],
    origin="lower",
    aspect="auto",
    cmap="viridis",
)
axs[1].set_xlabel("s", fontsize=12)
axs[1].set_ylabel(r"$\log_{10}$(k$_{\mathrm{ADCP}}$)", fontsize=12)
# Red cross marks the fitted (s, log10 k_ADCP) pair in axis coordinates
axs[1].plot(scale_target, k_adcp_target, "rx", markersize=15, markeredgewidth=4)
axs[1].text(-4, 0, "B", fontsize=18, fontweight="bold", color="black")
plt.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04, label="Error score")

plt.tight_layout()
fig.savefig("Plots/param_sweep_subplots_{}.png".format(grid_n), dpi=300, bbox_inches="tight")
