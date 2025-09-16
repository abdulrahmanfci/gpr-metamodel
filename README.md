# Two-level Gaussian Process Regression – Linear Response Metamodel for Epidemic Intervention Simulation

This repository contains the implementation of a two-level metamodeling framework developed for county-level opioid epidemic simulations in Pennsylvania. The framework approximates outputs of a large-scale agent-based model (FRED) under varying intervention dispensing rates of naloxone and buprenorphine.

## Overview

## Gaussian Process Regression (Background)

Gaussian Process Regression (GPR) is a Bayesian nonparametric method for learning functions from data.  
Instead of assuming a fixed functional form, GPR places a prior distribution directly over possible functions, defined by a **mean function** and a **kernel (covariance) function**.  

For a set of training data points, GPR computes a posterior distribution over functions that balance two sources of information:
- **Prior knowledge** encoded in the kernel (e.g., smoothness, periodicity).
- **Observed simulation outputs** from county–treatment experiments.  

This results in predictions that include both a **posterior mean** (the most likely function value) and a **posterior variance** (uncertainty about that value).  

In our framework:
- The GPR is **multi-output (MO–GPR)**, jointly modeling the overdose mortality and the effects of naloxone and buprenorphine.  
- The kernel is a **composite structure** (multiple RBF components plus a periodic term), which allows the model to capture both smooth spatial variation and recurrent patterns in county-level features.  
- The posterior mean function $\mathbf{\mu(x_c)} = \[\mu_0(\mathbf{x}_c), \ \mu_n(\mathbf{x}_c), \ \mu_b(\mathbf{x}_c)\]^{\top}$ represent regression coefficients for county $c$, while the posterior variances quantify uncertainty due to limited simulation runs.  

This stage provides a flexible statistical surrogate for the simulation model, enabling efficient exploration of the treatment space before the response function stage is applied.

The method integrates:
1. **Multi-Output Gaussian Process Regression (MO–GPR):**  
   - Incorporates spatial (county centroids) and socio-economic features (population density, median income, black residents percentage).  
   - Uses composite kernels (RBF + periodic components).  
   - Produces posterior mean and variance for intervention effect coefficients.

2. **Response Function:**  
   - Approximates overdose mortality outcomes as:  
     ```math
     z(n,b \mid c) = \mu_0(\mathbf{x}_c) + \mu_n(\mathbf{x}_c) \cdot n + \mu_b(\mathbf{x}_c) \cdot b.
     ```  
   - Provides interpretable coefficients for naloxone and buprenorphine effects at the county level.

3. **Sequential Design:**  
   - Allocates simulation runs adaptively across counties and treatments.  
   - Uses a signal-to-noise ratio acquisition rule to prioritize informative simulations.  
   - Iteratively updates the metamodel after each simulation batch.
  
## Data Description

The experiments in this repository are based on simulation outputs from a calibrated agent-based model of the opioid epidemic in Pennsylvania.  

- **Geographic scope:** All 67 counties in Pennsylvania.  
- **Interventions modeled (dispensing rates):**  
  - **Naloxone (overdose reversal/harm reduction):**  
    Naloxone is a short-acting opioid antagonist that rapidly reverses overdoses. In the model, naloxone intervention levels correspond to different county-level dispensing rates of naloxone kits.  
  - **Buprenorphine (treatment for opioid use disorder):**  
    Buprenorphine is a partial opioid agonist that reduces cravings and withdrawal symptoms. It is used as a maintenance treatment to mitigate the long-term burden of opioid use disorder. In the model, buprenorphine intervention levels correspond to different county-level dispensing rates of buprenorphine prescriptions.  

- **Treatment grid:** Each intervention is varied over **five discrete dispensing rate levels**, resulting in a **5 × 5 grid** (25 total treatment conditions).  
- **Simulation outputs:** County-level overdose mortality counts and rates under each intervention combination.  
- **Calibration:** Six representative counties (Allegheny, Philadelphia, Dauphin, Erie, Columbia, and Clearfield) were first calibrated to historical overdose mortality data, and the remaining counties were assigned to one of these prototypes based on demographic and epidemic similarity.  

These calibrated simulations form the input data for the two-level metamodeling framework (MO–GPR followed by a linear response surface), which approximates outcomes for all county–treatment combinations.


## Notebook Structure

- **Data Preprocessing:** Load calibrated county simulations and prepare feature sets.  
- **GPR-RF:** Fit and update MO–GPR using BoTorch.  
- **Sequential Design Loop:** Select counties and treatment conditions to simulate.  
- **Response Level:** Estimate intervention coefficients from GPR posteriors.  
- **Results & Visualization:** Generate comparison plots across counties and treatment strategies.

## Dependencies

- Python 3.10+  
- [BoTorch](https://botorch.org/)  
- PyTorch  
- NumPy, Pandas, Matplotlib  

## Citation
Please cite our work:
