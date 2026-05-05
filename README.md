#A Realistic Simulation Framework for Evaluating Microbiome Normalization in Sample Stratification and Taxa-Level Analysis
This repository contains scripts and data for evaluating normalization methods in simulated metagenomics count data using alpha-diversity screening, PCA, and observation-based Missing-data methods for Exploratory Data Analysis (oMEDA).

The framework compares raw and normalized microbiome datasets against known simulated ground truth at both phylum and genus levels.

## Contact

For questions or support, please contact:

**Amen A. Khabeer**  
amen.a.khabeer@uotechnology.edu.iq

Last document update: **15/04/2026**

---

## Requirements

The following software is required:

- MATLAB
- RStudio
- MEDA Toolbox version 1.10

Download MEDA Toolbox from the official GitHub repository:

https://github.com/josecamachop/MEDA-Toolbox

---

## Installing MEDA Toolbox in MATLAB

After downloading MEDA Toolbox version 1.10, unzip the folder and add it to the MATLAB path.

---

##In MATLAB, run:

```matlab
addpath(genpath('path_to/MEDA-Toolbox-1.10'));
savepath;

---

##Data/

Contains the simulated datasets and the corresponding ground-truth files for both phylum-level and genus-level analyses.

This folder can be used directly if you do not want to regenerate the simulations in R.

---

##PCA_oMEDA/

Contains MATLAB functions used to run PCA and oMEDA analyses and to compare the recovered structure against the simulated ground truth.

---

##Result/

Contains output figures and results from PCA-oMEDA analyses.

Results are organized by sequencing-depth subset:

Result/
├── high/
├── low/
└── whole/
Each folder contains outputs such as:

PCA score plots
PCA loading plots
oMEDA plots
distance/error results
Simulation_Normalization/

Contains the R scripts used to generate the simulated datasets and apply normalization methods.

---

##How to Use This Repository

You have two options:

Use the pre-generated datasets available in the Data/ folder.
Regenerate the simulated datasets from scratch using the R scripts.
Option 1: Use the Existing Data

This is the recommended option if you want to reproduce the analysis directly.

Step 1: Unzip the Data Files

Unzip all data files inside the Data/ folder for both:

Phylum simulations
Genus simulations
Step 2: Place the Data in the MATLAB Working Directory

Place the unzipped datasets in the same working directory where you placed or added the MEDA Toolbox.

The working directory should contain:

MEDA-Toolbox-1.10/
Data/
PCA_oMEDA/
Step 3: Load the Data in MATLAB

Load the simulated dataset and its corresponding ground-truth file into the MATLAB workspace.

For example, for the phylum-level simulation:

phylums = readtable('your_phylum_dataset.csv');
phylumg = readmatrix('your_phylum_ground_truth.csv');

The simulated dataset should contain the taxa abundance table.

The ground-truth file should be loaded as numeric data.

Step 4: Run the PCA-oMEDA Analysis

Open the PCA-oMEDA function in MATLAB. The main function is located in:

PCA_oMEDA/

Run the analysis using:

Val_PCAOmeda(phylums, phylumg);

where:

phylums is the simulated or normalized phylum-level dataset.
phylumg is the phylum-level ground truth.

Repeat the same procedure for each normalization method by replacing phylums with the corresponding normalized dataset.

Step 5: Repeat for Genus-Level Data

For genus-level simulations, repeat the same workflow using the genus dataset and genus ground-truth file.

For example:

Val_PCAOmeda(genuss, genusg);

where:

genuss is the simulated or normalized genus-level dataset.
genusg is the genus-level ground truth.
Option 2: Regenerate the Simulations from Scratch

Use this option if you want to recreate the simulated datasets using the R scripts.

Step 1: Open the Simulation Scripts in RStudio

Go to the folder:

Simulation_Normalization/

Open the R script corresponding to the scenario that you want to regenerate.

For example, to regenerate the true positive phylum simulation, open:

true.R
Step 2: Run the R Script

Run the selected R script in RStudio.

For the phylum true-positive scenario, running true.R will generate:

100 raw simulated datasets
600 normalized datasets across the selected normalization methods

The same logic applies to the other simulation scenarios and to the genus-level simulations.

Step 3: Move the Generated Data to MATLAB

After generating the datasets in R, place the generated files in the same MATLAB working directory used for the MEDA Toolbox and PCA-oMEDA functions.

The working directory should contain:

MEDA-Toolbox-1.10/
PCA_oMEDA/
Generated simulation datasets/
Ground truth files/
Step 4: Run Alpha-Diversity Analysis

After placing the generated datasets in the MATLAB working directory, run:

alpha_code.m

This script performs the alpha-diversity analysis used to screen the simulated datasets.

Step 5: Run PCA-oMEDA Analysis

For the oMEDA analysis, place the following files in the same MATLAB working directory:

Ground truth files
omeda_run.m
MEDA-Toolbox-1.10/

Then run the oMEDA analysis script:

omeda_run.m

This script performs PCA-oMEDA analysis and compares the recovered multivariate structure with the known simulated ground truth.

Step 6: Run ANOVA Analysis

After completing the alpha-diversity and PCA-oMEDA analyses, the ANOVA analysis can be run using the corresponding MATLAB scripts.

---

##General Workflow

The full workflow is:

Generate or load simulated count data
        ↓
Apply normalization methods
        ↓
Run alpha-diversity screening
        ↓
Run PCA-oMEDA analysis
        ↓
Compare recovered taxa-level structure with ground truth
        ↓
Summarize performance across normalization methods
Normalization Methods

The repository supports analysis of raw and normalized datasets, including:

Raw counts
Total-sum scaling
Rarefaction
Centered log-ratio transformation
Other normalization outputs generated by the R scripts

The exact number of datasets depends on the selected simulation scenario and taxonomic level.

---

##Notes on Data Placement

For the MATLAB scripts to run correctly, the following files should be in the same working directory or added to the MATLAB path:

MEDA-Toolbox-1.10/
PCA_oMEDA/
Ground truth files/
Simulation or normalized datasets/

If MATLAB cannot find a function or dataset, check that the correct folders have been added to the MATLAB path:

addpath(genpath(pwd));
savepath;
---

##Citation

Please acknowledge the use of MEDA Toolbox by citing:

Camacho, J., Pérez, A., Rodríguez, R., Jiménez-Mañas, E.
Multivariate Exploratory Data Analysis (MEDA) Toolbox.
Chemometrics and Intelligent Laboratory Systems, 2015, 143: 49–57.
Available at: https://github.com/josecamachop/MEDA-Toolbox

Please also check the documentation of the MEDA routines used in this repository for additional relevant references.

---

##Acknowledgements

We acknowledge the direct and indirect contributions of:

Alejandro García Vázquez, for contributions to the simulation design.

For issues related specifically to MEDA Toolbox, please contact:

Dr. José Camacho
josecamacho@ugr.es

---

##License

Copyright (C) 2026 Universidad de Granada
Copyright (C) 2026 Amen Adnan and aniel Vallejo-Espa ˜na

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or, at your option, any later version.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see:

http://www.gnu.org/licenses/



One important correction: do **not** write “normalization methods for metagenomics data” without specifying the data type. This framework is not assessing all metagenomics normalization in general; it assesses normalization of **simulated microbiome/metagenomic count tables** under known ground truth. That wording is more defensible.
