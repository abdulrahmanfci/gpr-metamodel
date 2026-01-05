#===============================================================================
# Calibrating PA Counties using IMIS (Incremental Mixture Importance Sampling)
# This file contains the cohort model, PSA, and likelihood functions
# We have 7 parameters to estimate for 3 transition probabilities
#===============================================================================

set.seed(1234)

#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

logi_fun <- function(val) {
  
  return((1 / (1 + exp(-1 * val))))
}

#===============================================================================
# PARAMETER BOUNDS FOR CALIBRATION
# Each county has specific lower/upper bounds for the 7 parameters:
# [1] NU_PU_intercept, [2] OUD_Rx_intercept, [3] OUD_DeathOD_intercept,
# [4] NU_PU_slope, [5] OUD_Rx_bup, [6] OUD_DeathOD_fent, [7] OUD_DeathOD_nal
#===============================================================================

# --- Allegheny ---
# v1_lower <- matrix(c(-130.40, -25.5, -9.0, 0.3, 0.0, 0.4, 0.6),
#                    nrow = 1, ncol = 7, byrow = TRUE)
# v1_upper <- matrix(c(0.16, 0.19, -5.0, 1.0, 2.0, 1.0, 1.1),
#                    nrow = 1, ncol = 7, byrow = TRUE)
# target <- c(41.76, 63.09, 77.14, 46.12, 53.31)

# --- Philadelphia ---
# v1_lower <- matrix(c(-130.40, -15.5, -9.0, 0.3, 1.1, 0.4, 0.9),
#                    nrow = 1, ncol = 7, byrow = TRUE)
# v1_upper <- matrix(c(0.16, -0.19, -4.53, 1.0, 2.9, 1.9, 1.1),
#                    nrow = 1, ncol = 7, byrow = TRUE)
# target <- c(48.33, 61.96, 88.7, 80.17, 84.07)

# --- Dauphin ---
# v1_lower <- matrix(c(-130.40, -7.5, -8.0, 0.3, 0.0, 0.0, 0.0),
#                    nrow = 1, ncol = 7, byrow = TRUE)
# v1_upper <- matrix(c(0.16, 1.19, 0.53, 1.0, 0.5, 0.9, 1.5),
#                    nrow = 1, ncol = 7, byrow = TRUE)
# target <- c(24.93, 34.07, 38.27, 48.73, 41.77)

# --- Erie ---
# v1_lower <- matrix(c(-130.40, -5.5, -15.0, 0.3, 0.0, 0.1, 0.0),
#                    nrow = 1, ncol = 7, byrow = TRUE)
# v1_upper <- matrix(c(0.16, 5.19, 0.53, 1.0, 0.9, 0.9, 1.5),
#                    nrow = 1, ncol = 7, byrow = TRUE)
# target <- c(23.3, 34.2, 48.68, 30.8, 28.17)

# --- Columbia ---
# v1_lower <- matrix(c(-130.40, -55.5, -10.0, 0.3, 0.0, 0.1, 0.0),
#                    nrow = 1, ncol = 7, byrow = TRUE)
# v1_upper <- matrix(c(0.16, 1.19, 0.53, 1.8, 0.5, 0.9, 1.5),
#                    nrow = 1, ncol = 7, byrow = TRUE)
# target <- c(20.37, 24.29, 26.37, 34.1, 26.6)

# --- Clearfield ---
v1_lower <- matrix(c(-130.40, -55.5, -10.0, 0.3, 0.0, 0.0, 0.0),
                   nrow = 1, ncol = 7, byrow = TRUE)
v1_upper <- matrix(c(0.16, 1.19, 0.53, 1.0, 0.5, 0.9, 1.5),
                   nrow = 1, ncol = 7, byrow = TRUE)
target <- c(18.23, 18.29, 12.21, 32.02, 18.32)

#===============================================================================
# OUD MARKOV MODEL
# Simulates opioid use disorder dynamics over 5 years (60 monthly cycles)
#===============================================================================

oud_model <- function(params) {
  cycles_length <- 5
  n_cycles <- 12 * 5
  pr <- matrix(0, 1, 3)
  
  #---------------------------------------------------------------------------
  # Standardized dispensing rates (z-scores) for each county
  
  # Data source: IQVIA (not publicly available)
  # Users must obtain their own county-level dispensing data and standardize
  #---------------------------------------------------------------------------
  
  # Placeholder vectors - replace with your standardized county data
  # Format: c(year1, year2, year3, year4, year5) for 2015-2019
  opioid_rx <- c(0, 0, 0, 0, 0)
  bupre_dispense_rate <- c(0, 0, 0, 0, 0)
  fent_dispense_rate <- c(0, 0, 0, 0, 0)
  naloxone_dispense_rate <- c(0, 0, 0, 0, 0)
  
  # --- Uncomment and populate for your county ---
  # Allegheny standardized rates
  # opioid_rx <- c(...)
  # bupre_dispense_rate <- c(...)
  # fent_dispense_rate <- c(...)
  # naloxone_dispense_rate <- c(...)
  
  # Philadelphia standardized rates
  # opioid_rx <- c(...)
  # bupre_dispense_rate <- c(...)
  # fent_dispense_rate <- c(...)
  # naloxone_dispense_rate <- c(...)
  
  # Dauphin standardized rates
  # opioid_rx <- c(...)
  # bupre_dispense_rate <- c(...)
  # fent_dispense_rate <- c(...)
  # naloxone_dispense_rate <- c(...)
  
  # Erie standardized rates
  # opioid_rx <- c(...)
  # bupre_dispense_rate <- c(...)
  # fent_dispense_rate <- c(...)
  # naloxone_dispense_rate <- c(...)
  
  # Columbia standardized rates
  # opioid_rx <- c(...)
  # bupre_dispense_rate <- c(...)
  # fent_dispense_rate <- c(...)
  # naloxone_dispense_rate <- c(...)
  
  # Clearfield standardized rates
  # opioid_rx <- c(...)
  # bupre_dispense_rate <- c(...)
  # fent_dispense_rate <- c(...)
  # naloxone_dispense_rate <- c(...)
  
  #---------------------------------------------------------------------------
  # Calculate initial transition probabilities
  #---------------------------------------------------------------------------
  pr[1] <- params[1] + params[4] * opioid_rx[1]
  pr[2] <- params[2] + params[5] * bupre_dispense_rate[1]
  pr[3] <- params[3] + params[6] * fent_dispense_rate[1] - params[7] * naloxone_dispense_rate[1]
  
  # Apply logistic transformation
  pr[1] <- max(0, logi_fun(pr[1]))
  pr[2] <- max(0, logi_fun(pr[2]))
  pr[3] <- max(0, logi_fun(pr[3]))
  
  #---------------------------------------------------------------------------
  # Initialize cohort
  #---------------------------------------------------------------------------
  v_m_init <- c(NU = 0.909, PU = 0.006, MU = 0.04, OUD = 0.04, 
                OUD_treat = 0.005, DeathOd = 0, Death = 0)
  n_states <- length(v_m_init)
  v_names_states <- names(v_m_init)
  
  m_M <- matrix(NA, nrow = (n_cycles + 1), ncol = n_states,
                dimnames = list(0:n_cycles, v_names_states))
  m_M[1, ] <- v_m_init
  
  #---------------------------------------------------------------------------
  # Build transition probability matrix
  #---------------------------------------------------------------------------
  m_P <- matrix(0, nrow = n_states, ncol = n_states,
                dimnames = list(v_names_states, v_names_states))
  
  zero_p <- 0.0
  one_p <- 1.0
  
  # From NU (Non-User)
  m_P["NU", "NU"] <- zero_p
  m_P["NU", "PU"] <- pr[1]
  m_P["NU", "MU"] <- 0.001
  m_P["NU", "OUD"] <- 0.0008
  m_P["NU", "OUD_treat"] <- zero_p
  m_P["NU", "DeathOd"] <- zero_p
  m_P["NU", "Death"] <- zero_p
  m_P["NU", "NU"] <- max(0, 1.0 - sum(m_P["NU", ]))
  
  # From PU (Prescribed User)
  m_P["PU", "PU"] <- zero_p
  m_P["PU", "NU"] <- 0.373
  m_P["PU", "MU"] <- 0.05
  m_P["PU", "OUD"] <- 0.0012778
  m_P["PU", "OUD_treat"] <- zero_p
  m_P["PU", "DeathOd"] <- zero_p
  m_P["PU", "Death"] <- 0.000238
  m_P["PU", "PU"] <- max(0, 1.0 - sum(m_P["PU", ]))
  
  # From MU (Misuser)
  m_P["MU", "MU"] <- zero_p
  m_P["MU", "NU"] <- 0.039
  m_P["MU", "PU"] <- 0
  m_P["MU", "OUD"] <- 0.023
  m_P["MU", "OUD_treat"] <- 0
  m_P["MU", "DeathOd"] <- 1.3e-7
  m_P["MU", "Death"] <- 0.000238
  m_P["MU", "MU"] <- max(0, 1.0 - sum(m_P["MU", ]))
  
  # From OUD (Opioid Use Disorder)
  m_P["OUD", "OUD"] <- zero_p
  m_P["OUD", "NU"] <- 0.0162
  m_P["OUD", "PU"] <- zero_p
  m_P["OUD", "MU"] <- 0.013
  m_P["OUD", "OUD_treat"] <- pr[2]
  m_P["OUD", "DeathOd"] <- pr[3]
  m_P["OUD", "Death"] <- 0.00258
  m_P["OUD", "OUD"] <- max(0, 1.0 - sum(m_P["OUD", ]))
  
  # From OUD_treat (OUD in Treatment)
  m_P["OUD_treat", "OUD_treat"] <- zero_p
  m_P["OUD_treat", "NU"] <- 0.009
  m_P["OUD_treat", "PU"] <- zero_p
  m_P["OUD_treat", "MU"] <- zero_p
  m_P["OUD_treat", "OUD"] <- 0.0089
  m_P["OUD_treat", "DeathOd"] <- zero_p
  m_P["OUD_treat", "Death"] <- 0.000238
  m_P["OUD_treat", "OUD_treat"] <- max(0, 1.0 - sum(m_P["OUD_treat", ]))
  
  # Absorbing states
  m_P["DeathOd", "DeathOd"] <- one_p
  m_P["Death", "Death"] <- one_p
  
  #---------------------------------------------------------------------------
  # Run Markov model with time-varying transition probabilities
  #---------------------------------------------------------------------------
  year_change <- c(1, 12, 24, 36, 48)
  
  for (t in 1:n_cycles) {
    m_M[t + 1, ] <- m_M[t, ] %*% m_P
    
    if (t %in% year_change) {
      index <- which(year_change == t)
      
      # Update transition probabilities based on year-specific dispensing rates
      pr[1] <- params[1] + params[4] * opioid_rx[index]
      pr[2] <- params[2] + params[5] * bupre_dispense_rate[index]
      pr[3] <- params[3] + params[6] * fent_dispense_rate[index] - params[7] * naloxone_dispense_rate[index]
      
      pr[1] <- max(0, logi_fun(pr[1]))
      pr[2] <- max(0, logi_fun(pr[2]))
      pr[3] <- max(0, logi_fun(pr[3]))
      
      # Update NU -> PU transition
      m_P["NU", "PU"] <- pr[1]
      m_P["NU", "NU"] <- zero_p
      m_P["NU", "NU"] <- max(0, 1.0 - sum(m_P["NU", ]))
      
      # Update OUD -> OUD_treat and OUD -> DeathOd transitions
      m_P["OUD", "OUD_treat"] <- pr[2]
      m_P["OUD", "DeathOd"] <- pr[3]
      m_P["OUD", "OUD"] <- zero_p
      m_P["OUD", "OUD"] <- max(0, 1.0 - sum(m_P["OUD", ]))
    }
  }
  
  #---------------------------------------------------------------------------
  # Calculate annual overdose death rates per 100,000
  #---------------------------------------------------------------------------
  m_M <- cbind(m_M, sum_alive = rowSums(m_M[, 1:5], na.rm = TRUE))
  
  rate_years <- rep(0, cycles_length)
  od_death_prob <- 0
  
  for (i in 1:cycles_length) {
    k <- i - 1
    rate_years[i] <- (m_M[(i * 12), 6] - od_death_prob) / 
      (mean(m_M[(k * 12):(i * 12), 8])) * 1e5
    od_death_prob <- m_M[(i * 12), 6]
  }
  
  return(rate_years)
}

#===============================================================================
# LATIN HYPERCUBE SAMPLING FOR PRIOR
#===============================================================================

require(lhs)
n_calibrated_params <- length(v1_lower)

sample.prior.lhs <- function(n) {
  A_samples <- randomLHS(n = n, k = n_calibrated_params)
  draws <- data.frame(
    NU_PU_intercept = qunif(A_samples[, 1], v1_lower[1], v1_upper[1]),
    OUD_Rx_intercept = qunif(A_samples[, 2], v1_lower[2], v1_upper[2]),
    OUD_DeathOD_intercept = qunif(A_samples[, 3], v1_lower[3], v1_upper[3]),
    NU_PU_slope = qunif(A_samples[, 4], v1_lower[4], v1_upper[4]),
    OUD_Rx_bup = qunif(A_samples[, 5], v1_lower[5], v1_upper[5]),
    OUD_DeathOD_fent = qunif(A_samples[, 6], v1_lower[6], v1_upper[6]),
    OUD_DeathOD_nal = qunif(A_samples[, 7], v1_lower[7], v1_upper[7])
  )
  return(as.matrix(draws))
}

sample.prior <- sample.prior.lhs

#===============================================================================
# PROBABILISTIC SENSITIVITY ANALYSIS (PSA)
#===============================================================================

size <- 250

psa_run <- function(par_vector) {
  if (is.null(dim(par_vector))) par_vector <- t(par_vector)
  n <- 5
  psa <- matrix(0, nrow = nrow(par_vector), ncol = n)
  
  for (j in 1:nrow(par_vector)) {
    res_j <- oud_model(as.numeric(par_vector[j, ]))
    psa[j, 1:5] <- res_j
  }
  return(psa)
}

set.seed(1234)
input_params <- sample.prior.lhs(size)
psa_data <- psa_run(input_params)
write.csv(psa_data, "psa_data.csv")

#===============================================================================
# VISUALIZE PSA RESULTS
#===============================================================================

library(ggplot2)
library(reshape2)

data <- read.csv("psa_data.csv")
data <- data[, -1]
colnames(data) <- c(1, 2, 3, 4, 5)

data_long <- melt(data, variable.name = "Year", value.name = "Rate")

# Calculate 95% credible intervals
summary_stats <- aggregate(Rate ~ Year, data_long, function(x) {
  c(mean = mean(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975))
})
summary_stats <- do.call(data.frame, summary_stats)
colnames(summary_stats) <- c("Year", "Mean", "Lower", "Upper")

# Plot
p <- ggplot(summary_stats, aes(x = Year, y = Mean)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "blue") +
  geom_point(data = data.frame(Year = unique(summary_stats$Year), Target = target),
             aes(x = Year, y = target), color = "red", size = 3) +
  labs(title = "PSA Results: Model vs Target Overdose Rates",
       x = "Year",
       y = "Overdose Death Rate per 100,000")
print(p)

#===============================================================================
# LIKELIHOOD FUNCTION FOR CALIBRATION
#===============================================================================

# County populations for standard error calculation
# Uncomment the appropriate line for your county
# std_target <- sqrt(target * 1248340 / 100000) / 1248340 * 1e5  # Allegheny
# std_target <- sqrt(target * 1583219 / 100000) / 1583219 * 1e5  # Philadelphia
# std_target <- sqrt(target * 276645 / 100000) / 276645 * 1e5    # Dauphin
# std_target <- sqrt(target * 279609 / 100000) / 279609 * 1e5    # Erie
# std_target <- sqrt(target * 66647 / 100000) / 66647 * 1e5      # Columbia
std_target <- sqrt(target * 81708 / 100000) / 81708 * 1e5        # Clearfield

l_likelihood <- function(par_vector) {
  if (is.null(dim(par_vector))) par_vector <- t(par_vector)
  llik <- rep(0, nrow(par_vector))
  
  for (j in 1:nrow(par_vector)) {
    jj <- tryCatch({
      res_j <- oud_model(c(as.numeric(par_vector[j, ]), 1))
      llik[j] <- llik[j] + sum(dnorm(as.numeric(target), 
                                     mean = mean(res_j), 
                                     sd = std_target, 
                                     log = TRUE))
    }, error = function(e) NA)
    
    if (is.na(jj)) {
      llik[j] <- -Inf
    }
  }
  return(llik)
}

# Test likelihood function
par_vector <- sample.prior.lhs(250)
l_likelihood_val <- l_likelihood(par_vector)
summary(l_likelihood_val)