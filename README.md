# Tumor-Growth-PDE-Model-Relevant-Code
# A Three-dimensional Tumor Growth Model and Its Boundary Instability

This repository contains the MATLAB code accompanying the paper "A Three-dimensional tumor growth model and its boundary instability". The code includes both analytical and numerical techniques used to validate and extend the results from Feng et al., 2023. The focus is on boundary instability using perturbation, asymptotic analysis with spherical harmonics, and numerical simulations with the Alternating Directional Implicit (ADI) method.

## Code Files

### 3D_tumor_evolution_visualization.ipynb
Visualization of the 3D tumor evolution model. This notebook contains scripts for generating Figure 7: "3D tumor evolution speed in vitro and in vivo (unperturbed)".

### full_coupling_adaptive_lambda_1_2.m
MATLAB script used for generating Figure 3: "2D tumor growth boundary evolution with λ = 0.5 (three time profiles evolving from left to right) for an initial perturbation with wave number m = 8".

### full_coupling_adaptive_lambda_40.m
MATLAB script for Figure 4: "2D tumor growth boundary evolution with λ = 40 (three time profiles evolving from left to right) for an initial perturbation with wave number m = 8".

### drdt.m
This script corresponds to Figures 5 and 6, representing "3D in vitro tumor evolution" and "3D in vivo tumor evolution", respectively.

### vitro_evolution_function.m
Used to simulate the in vitro tumor growth profiles, resulting in Figure 8: "3D in vitro λ = 0.8, ℓ = 3,8,12,16,18" and Figure 9: "3D in vitro λ = 5, ℓ = 3,8,12,16,18".

### vivo_evolution_function_plot.m
Produces the in vivo tumor growth simulations as shown in Figures 10 to 13, demonstrating various growth parameters and their effects on tumor stability.

### axisymm_pme_2d.m
Generates the plots for Figure B1: "2D axisymmetric PME-Comparison of numerical and analytical solution" and Figure B2: "validation of 2D axisymmetric PME".

### adi_2d_compare.m
This script is used to compare the 2D ADI method against axisymmetric results, as shown in Figure B3: "Comparison of 2D ADI and axisymmetric".

## Usage

To run the scripts, MATLAB R2021a or later is recommended. Each script is associated with a specific figure from the paper. Run the script in MATLAB to reproduce the results and figures as shown in the study.

