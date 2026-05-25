# Example/Data README

This folder contains input and processed data used for model simulation, validation, and plotting.

---

## Folder Overview

- Observed pharmacokinetic profiles.
- Observed percentage change from baseline in amyloid PET signal from clinical trials.
- Model inputs to simulate various drugs in various clinical trial contexts.

---

## File/Folder Listings

| File/Folder Name         | Type         | Description / Source / Notes                         |
|-------------------------|--------------|------------------------------------------------------|
| `clearance_data.csv`        | Digitized  | Percentage change from baseline in amyloid PET signals       |
| `drug_data.csv`        | -  | Model inputs to simulate different drugs in different clinical trial scenarios   |
| `PK_data.csv`        | -  | Model inputs to simulate different drugs in different clinical trial scenarios        |
| `PK.csv`        | Digitized  | Plasma drug concentrations from SAD studies of aducanumab, lecanemab, and donanemab       |

---

## Data Format Description

- `clearance_data.csv` columns: time_points (in days); months; one column for the %CBL of each trial scenario
- `drug_data.csv` columns: model; dose; bodyweight; gpm; final_dose (h); dose_interval (days); monomer (affinity); oligomer (affinity); plaque (affinity); simulation_time (h); colour (for plotting)
- `PK_data.csv` columns: columns: model; dose; bodyweight; gpm; final_dose (h); dose_interval (days); monomer (affinity); oligomer (affinity); plaque (affinity); simulation_time (h)
- `PK.csv` columns: time (in h); value (in ug/ml); Dose (3 or 10 mg/kg); Drug (aducanumab, lecanemab, or donanemab)

---
