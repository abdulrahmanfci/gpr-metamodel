# Two-level Gaussian Process Regression – Linear Response Metamodel for epidemic Intervention Simulation

This repository contains the implementation of a two-stage metamodeling framework developed for county-level opioid epidemic simulations in Pennsylvania. The framework approximates outputs of a large-scale agent-based model (FRED) under varying intervention dispensing rates of naloxone and buprenorphine.

## Overview

The method integrates:
1. **Multi-Output Gaussian Process Regression (MO–GPR):**  
   - Incorporates spatial (county centroids) and socio-economic features (population size, overdose trends, dispensing rates).  
   - Uses composite kernels (RBF + periodic components).  
   - Produces posterior mean and variance for intervention effect coefficients.

2. **Response Function:**  
   - Approximates overdose mortality outcomes as:  
     $$ z(n,b \mid c) = \mu_0(\mathbf{x}_c) + \mu_n(\mathbf{x}_c) \cdot n + \mu_b(\mathbf{x}_c) \cdot b. $$  
   - Provides interpretable coefficients for naloxone and buprenorphine effects at the county level.

3. **Sequential Design:**  
   - Allocates simulation runs adaptively across counties and treatments.  
   - Uses a signal-to-noise acquisition rule to prioritize informative simulations.  
   - Iteratively updates the metamodel after each simulation batch.

## Notebook Structure

- **Data Preprocessing:** Load calibrated county simulations and prepare feature sets.  
- **GPR-RF:** Fit and update MO–GPR using BoTorch.  
- **Sequential Design Loop:** Select counties and treatment conditions to simulate.  
- **Response Stage:** Estimate intervention coefficients from GPR posteriors.  
- **Results & Visualization:** Generate comparison plots across counties and treatment strategies.

## Dependencies

- Python 3.10+  
- [BoTorch](https://botorch.org/)  
- PyTorch  
- NumPy, Pandas, Matplotlib  

## Citation
Please cite our work:
