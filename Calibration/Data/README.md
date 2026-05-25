# Calibration/Data README

This folder contains input and processed data used for model calibration.

---

## Folder Overview
- CSF and plasma SILK data used for model calibration
- Pharmacokinetic profiles for three anti-amyloid mAbs used for model calibration

---

## File Listings

| File         | Type         | Description                         |
|-------------------------|--------------|------------------------------------------------------|
| `bolus_800_csf_new.csv`        | Derived  | Mean fraction CSF leucine labeled in plasma SILK experiments        | 
| `bolus_800_plasma_new.csv`        | Derived  | Mean fraction plasma leucine labeled in plasma SILK experiments, derived from a profile digitized from [4]        |
| `csf_SILK_std.csv`        | Derived  | Standard deviation of labeled CSF Ab42 concentration in amyloid positive individuals at each timepoint    (concentrations in ng/ml)    |
| `mean_Ab42_CSF_SILK.csv`        | Derived  | Mean labeled CSF Ab42 concentration in amyloid positive individuals at each timepoint     |
| `mean_f_leu_csf.csv`        | Derived  | Mean fraction CSF leucine labeled at each timepoint        |
| `mean_f_leu.csv`        | Derived  | Mean fraction plasma leucine labeled at each timepoint        |
| `mean_nM_Ab42_young_CSF_SILK.csv`        | Derived  | Mean labeled CSF Ab42 concentration in young, amyloid negative individuals at each timepoint    (concentrations in nM)    |
| `PK.csv`        | Digitized  | Plasma drug concentrations from SAD studies of aducanumab, lecanemab, and donanemab        |
| `plasma_silk_digitized.csv`        | Digitized  | Mean fraction plasma Ab42 labeled, digitized from [3]         |
| `silk_ci_std.csv`        | Digitized  | 95% CI of fraction plasma Ab42 labeled, digitized from [3]         | 
| `SILK_sd_young.csv`        | Derived  | Standard deviation of labeled CSF Ab42 concentration in young, amyloid negative individuals at each timepoint   (concentrations in ng/ml)     |
| `CSF_SILK_data.ipynb`      | Notebook     | Details processing of data to obtain derived data used for model simulations     |


## Notes on Data Formats

Describe expected columns, units, and formats for each type of data file, e.g.:

- `mean_Ab42_CSF_SILK.csv` columns: amyloid positive; amyloid negative (concentrations in nM)
- `PK.csv` columns: time (in h); value (in ug/ml); Dose (3 or 10 mg/kg); Drug (aducanumab, lecanemab, or donanemab)
- `plasma_silk_digitized.csv` columns: time, h; labeled fraction
- `silk_ci_std.csv` columns: time_plasma (in h); ci_plasma (as MFL)

---


## Raw data not provided here

| File         | Type         | Description / Source / Notes                         |
|-------------------------|--------------|------------------------------------------------------|
| `Ab42_MFL_timepoints.csv`        | Raw  | Timepoints for SILK measurements for all individuals [1]       |
| `Ab42_mol_fraction_labeled.csv`        | Raw  | MFL CSF Ab42 at each timepoint for all individuals   [1]     |
| `leucine_mol_fraction_labeled.csv`        | Raw  | MFL plasma leucine at each timepoint for all individuals   [1]     |
| `patient_metadata.csv`        | Raw  | Includes demographic information on individual participants   [1]     |
| `patient_processing_order.csv` | Raw | Patient IDs matching order of data in each raw file [1] |
| `plasma_vs_csf_leucineTTR.csv`        | Digitized  | Mean fraction leucine labeled in CSF and plasma  [2]      |
| `total_Ab42_concentration.csv`        | Raw  | CSF Ab42 concentration at each timepoint for all individuals  [1]      |

---

## Processing Notes

- Individual .mat files from Supplementary Data 1 were downloaded and processed to collate relevant data across individuals. Each .mat file corresponds to an individual participant. The following files were generated during data collation: patient_processing_order, total_Ab42_concentration (Ab42conc in .mat), Ab42_mol_fraction_labeled (MFL42 in .mat), leucine_mol_fraction_labeled (f in .mat), and Ab42_MFL_timepoints (T42 in .mat). 
- `patient_metadata.csv` corresponds to Supplementary Data 3 and was used to subset individuals based on age and amyloid status.

---

## References

- <sup>1</sup>Elbert et al. 2022. 
- <sup>2</sup>Bateman et al. 2006.
- <sup>3</sup>Ovod et al. 2017.  
- <sup>4</sup>Sato et al. 2018.

---
