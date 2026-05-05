# # A Realistic Simulation Framework for Evaluating Microbiome Normalization in Sample Stratification and Taxa-Level Analysis

This repository contains MATLAB and R scripts for evaluating normalization methods in simulated microbiome count data using alpha-diversity screening, PCA, and oMEDA.

The framework compares raw and normalized datasets against known simulated ground truth at phylum and genus levels.

---

## Overview

This repository contains scripts and data for evaluating normalization methods in simulated metagenomics count data using alpha-diversity screening, PCA, and observation-based Missing-data methods for Exploratory Data Analysis (oMEDA).

The framework compares raw and normalized microbiome datasets against known simulated ground truth at both phylum and genus levels.

## Contact

For questions or support, please contact:

**Amen A. Khabeer**  
amen.a.khabeer@uotechnology.edu.iq

Last document update: **16/05/2026**


## Requirements
- MATLAB
- R / RStudio
- MEDA Toolbox v1.10

MEDA Toolbox:

https://github.com/josecamachop/MEDA-Toolbox

---
## Installation 
addpath(genpath('path_to/MEDA-Toolbox-1.10'));
savepath; 
Replace path_to with your local folder.

---

## Repository Structure 
.
├── README.md
├── LICENSE
├── Data/
├── PCA_oMEDA/
├── Result/
└── Simulation_Normalization/
---

## Folder Description

| Folder                   | Description                            |
| ------------------------ | -------------------------------------- |
| Data                     | Simulated datasets and ground truth    |
| PCA_oMEDA                | MATLAB analysis functions              |
| Result                   | Figures and performance outputs        |
| Simulation_Normalization | R simulation and normalization scripts |

---

## General Workflow

Generate or load simulated data
        ↓
Apply normalization methods
        ↓
Run alpha-diversity screening
        ↓
Run PCA
        ↓
Run oMEDA
        ↓
Compare with ground truth
        ↓
Summarize performance

---

## Usage 

### Option 1: Use Existing Data
Use files inside Data/.

Example
phylums = readtable('your_phylum_dataset.csv');
phylumg = readmatrix('your_phylum_ground_truth.csv');

Val_PCAOmeda(phylums, phylumg);

### Option 2: Generate New Simulations

Open R scripts in:

Simulation_Normalization/

Example:

true.R

Run in RStudio, then move outputs to MATLAB working directory.

---
Additional Scripts
Alpha Diversity
alpha_code.m
oMEDA
omeda_run.m

---

## Output 

Results are stored in:

Result/
├── high/
├── low/
└── whole/

Outputs may include:

PCA score plots
PCA loading plots
oMEDA plots
Distance metrics
Performance summaries

--- 

## Supported Normalization Methods

Raw counts
Total-sum scaling (TSS)
Rarefaction
Centered log-ratio (CLR)
Other methods from included R scripts

---


## Citation

If you use this repository, please cite:

Al Khafaji A, Vallejo-España D, Gómez-Llorente C, Camacho J. *A Realistic Simulation Framework for Evaluating Microbiome Normalization in Sample Stratification and Taxa-Level Analysis*. Manuscript in preparation.
https://github.com/AmenALkhafaji/PCA_oMEDA-for-assesst_norm_metageon

---

## Acknowledgements

We thank Alejandro García Vázquez for contributions to simulation development.

For MEDA-related issues:

Dr. José Camacho
josecamacho@ugr.es

.



## Acknowledgements
## License
