# OUD Model Calibration

This folder contains R scripts for calibrating the Opioid Use Disorder (OUD) Markov model using Incremental Mixture Importance Sampling (IMIS).

## Overview

The calibration process estimates 7 transition probability parameters by matching simulated overdose death rates to observed county-level mortality data from 2015-2019.

## Files

| File | Description |
|------|-------------|
| `calibration_pa.R` | Defines the OUD Markov model, parameter bounds, likelihood function, and PSA |
| `imis_calibration.R` | Runs IMIS algorithm to obtain posterior parameter estimates |

## Workflow
```
┌─────────────────────────────────┐
│  calibration_pa.R               │
│  ─────────────────────────────  │
│  1. Define parameter bounds     │
│  2. Define OUD Markov model     │
│  3. Define likelihood function  │
│  4. Run initial PSA             │
└───────────────┬─────────────────┘
                │
                ▼
┌─────────────────────────────────┐
│  imis_calibration.R             │
│  ─────────────────────────────  │
│  1. Define prior distribution   │
│  2. Run IMIS algorithm          │
│  3. Sample from posterior       │
│  4. Evaluate model fit          │
│  5. Export calibrated params    │
└─────────────────────────────────┘
```

## Usage

### Step 1: Configure County

In `calibration_pa.R`, uncomment the parameter bounds and target for your county:
```r
# --- Allegheny ---
v1_lower <- matrix(c(-130.40, -25.5, -9.0, 0.3, 0.0, 0.4, 0.6),
                   nrow = 1, ncol = 7, byrow = TRUE)
v1_upper <- matrix(c(0.16, 0.19, -5.0, 1.0, 2.0, 1.0, 1.1),
                   nrow = 1, ncol = 7, byrow = TRUE)
target <- c(41.76, 63.09, 77.14, 46.12, 53.31)
```

### Step 2: Add Dispensing Data

Populate the standardized dispensing rate vectors in `oud_model()` with your county's IQVIA data (z-scores for 2015-2019):
```r
opioid_rx <- c(...)
bupre_dispense_rate <- c(...)
fent_dispense_rate <- c(...)
naloxone_dispense_rate <- c(...)
```

### Step 3: Run Calibration
```r
# First, run the model script
source("calibration_pa.R")

# Then, run IMIS calibration
source("imis_calibration.R")
```

### Step 4: Check Outputs

Results are saved to `routput/`:

| Output File | Description |
|-------------|-------------|
| `imis_resample.csv` | Posterior parameter samples |
| `model_pred.csv` | Model predictions from posterior |
| `posterior_summary.csv` | Parameter estimates with 95% credible intervals |
| `model_rate_with_targets.pdf` | Visual comparison of model vs targets |
| `prior_post.pdf` | Prior vs posterior distributions |

## Model Parameters

The calibration estimates 7 parameters controlling three transition probabilities:

| Parameter | Description | Transition |
|-----------|-------------|------------|
| `NU_PU_intercept` | Intercept for non-user → prescribed user | $p_1$ |
| `NU_PU_slope` | Effect of opioid prescribing rate | $p_1$ |
| `OUD_Rx_intercept` | Intercept for OUD → treatment | $p_2$ |
| `OUD_Rx_bup` | Effect of buprenorphine availability | $p_2$ |
| `OUD_DeathOD_intercept` | Intercept for OUD → overdose death | $p_3$ |
| `OUD_DeathOD_fent` | Effect of fentanyl prevalence | $p_3$ |
| `OUD_DeathOD_nal` | Effect of naloxone availability | $p_3$ |

## Calibrated Counties

Six prototype counties were calibrated:

| County | Population | Type |
|--------|------------|------|
| Allegheny | 1,248,340 | Large urban |
| Philadelphia | 1,583,219 | Large urban |
| Dauphin | 276,645 | Medium |
| Erie | 279,609 | Medium |
| Columbia | 66,647 | Small rural |
| Clearfield | 81,708 | Small rural |

## Data Requirements

- **Target mortality rates**: Annual overdose deaths per 100,000 (2015-2019) from CDC WONDER
- **Dispensing rates**: County-level standardized rates from IQVIA (not publicly available)
  - Opioid prescriptions
  - Buprenorphine dispensing
  - Naloxone dispensing
  - Fentanyl seizure rates (from NFLIS)

## Dependencies
```r
install.packages(c("lhs", "matrixStats", "tidyverse", "reshape2", "ggplot2", "mvtnorm"))

# IMIS package (install from archive)
devtools::install_version("IMIS", version = "0.1", repos = "http://cran.us.r-project.org")
```
