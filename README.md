# Fibroblast Mechano-Chemo Signaling Network Model Analysis

08.05.2021 Jesse Rogers  

This repository provides code for fitting a Netflux-generated ODE network of fibroblast signaling to patent-derived plasma biomarker and cell gene expression data. Included here are scripts for performing model fitting, model validation, network perturbation analyses, as well as current drug case studies comparing indiviual patient responses to simulated inhibitors.

## Required Programs/Toolboxes

- **MATLAB 2017a or newer**: needed to generate heatmap objects (MATLAB 2020b used for study)
  - *Global Optimization Toolbox*: needed for genetic algorithm-based model fitting
  - *Parallel Computing Toolbox*: needed for drug screens on high performance computing environment
- **High performance computing environment/resource**: needed for model fitting only
  - Provided bash scripts rely on Portable Batch System (PBS) job scheduler as implemented by the [Palmetto cluster at Clemson University](http://www.palmetto.clemson.edu/palmetto/), and other environments may require alternative scripts

## Included Analyses

### 1. Model Fitting

Fitting of network reaction weights to patient-derived data using a genetic algorithm (see `src\model_fitting\`).

- `ModelFittingEnsemble.m`: Runtime script for implementing genetic algorithm with the provided data.
  - `ModelFitness.m`: Dependency for calculating fitness scores based on network simulations.

### 2. Model Validation

Qualitative validation of input-output predictions with patient proteomic/transcriptomic data.

- `EnsembleAggregation_training.m`: Runtime script for aggregating reaction weights from fitting. Uses training data (patients A-H) only.
  - `EnsembleFitness_training.m`: Dependency for predicting model output expression for individual reaction weight sets.
- `EnsembleAggregation_testing.m`: Runtime script for aggregating reaction weights from fitting. Uses testing data (patients I-L) only.
  - `EnsembleFitness_testing.m`: Dependency for predicting model output expression for individual reaction weight sets.
- `EnsembleAnalysis.m`: Runtime script for comparing model-predicted changes between pre-/post-TAVR conditions for individual patients (Figure 2) and comparing model predictions to clinical parameters (Figure 3).

### 3. Network Sensitivity Analysis

Simulation of individual node knockdowns under individual patient plasma biomarker profiles, and identification of node sensitivity towards knockdown ("knockdown sensitivity") and influence on network-wide activity with knockdown ("knockdown influence")

- `SensitivityAggregation.m`: Runtime script for simulating steady-state network activation in response to knockdown of individual nodes under specified patient plasma biomarker levels
  - `SensitivitySimulations.m`: Dependency for individual model predictions per set of reaction weights
- `SensitivityVisualization.m`: Runtime script for comparing top-ranking nodes in sensitivity and influence for each patient (Figure 4)

### 4. Patient-Specific Drug Screens

Simulation and comparison of patient-specific changes to matrix-related output expression with simulated inhibitor dosing (via node knockdown).

- `StratificationAggregation.m`: Runtime script for simulating single-dose or dose-response changes in network activation in response to individual node knockdown
  - `StratificationSimulation_singledose.m`: Dependency for model simulations per set of reaction weights
  - `StratificationSimulation_doseresponse.m`: Dependency for model simulations per set of reaction weights
- `StratificationVisualization.m`: Runtime script for analyzing/plotting patient responses to node knockdowns (Figure 5)

### 5. HPC Batch Scripts

Includes bash shell scripts for submitting Model Fitting scripts to PBS scheduler. Recommended resources per script: 24 CPUs, 120GB memory, 48:00 runtime