#===============================================================================
# IMIS Calibration for OUD Model
# This script runs Incremental Mixture Importance Sampling (IMIS) to obtain
# posterior parameter estimates, following the PSA and likelihood functions
# defined in the companion script (oud_calibration_model.R)
#
# Prerequisites: Run oud_calibration_model.R first to define:
#   - v1_lower, v1_upper (parameter bounds)
#   - target (observed overdose death rates)
#   - std_target (standard errors for likelihood)
#   - sample.prior.lhs() (prior sampling function)
#   - l_likelihood() (log-likelihood function)
#   - psa_run() (model evaluation function)
#===============================================================================

library(matrixStats)
library(tidyverse)
library(reshape2)
library(mvtnorm)
library(usethis)
library(devtools)

options(dplyr.summarise.inform = FALSE)

#===============================================================================
# PRIOR DISTRIBUTION
# Log-prior function using normal distributions centered on parameter bounds
# Parameters are assumed independent with sd proportional to range
#===============================================================================

l_prior <- function(par_vector) {
  if (is.null(dim(par_vector))) par_vector <- t(par_vector)
  lprior <- rep(0, nrow(par_vector))
  
  for (i in 1:7) {
    # Normal prior centered at midpoint, sd = range/4 (or range/20 for param 1)
    sd_scale <- ifelse(i == 1, 20, 4)
    lprior <- lprior + dnorm(par_vector[, i],
                             mean = (v1_lower[i] + v1_upper[i]) / 2,
                             sd = (v1_upper[i] - v1_lower[i]) / sd_scale,
                             log = TRUE)
  }
  return(lprior)
}

#===============================================================================
# IMIS SETUP
# IMIS requires three functions:
#   1. sample.prior - draws from prior (defined in oud_calibration_model.R)
#   2. prior - evaluates prior density
#   3. likelihood - evaluates likelihood
#===============================================================================

# Load IMIS package (install from local source if needed)
# devtools::install_version("IMIS", version = "0.1", repos = "http://cran.us.r-project.org")
load_all("path/to/IMIS")  # Update path to your IMIS installation
library(IMIS)

# Wrapper functions for IMIS (converts log-scale to regular scale)
prior <- function(par_vector) {
  exp(l_prior(par_vector))
}

likelihood <- function(par_vector) {
  exp(l_likelihood(par_vector))
}

#===============================================================================
# RUN IMIS
# B: initial sample size
# B.re: resample size
# number_k: number of mixture components
# D: number of optimizers (0 = none)
#===============================================================================

set.seed(123)
imis_res <- IMIS(B = 1500, B.re = 3e3, number_k = 150, D = 0)

# Diagnostics
cat("Unique parameter sets:", length(unique(imis_res$resample[, 1])), "\n")
cat("Effective sample size:",
    sum(table(imis_res$resample[, 1]))^2 / sum(table(imis_res$resample[, 1])^2), "\n")
cat("Max weight:",
    max(table(imis_res$resample[, 1])) / sum(table(imis_res$resample[, 1])), "\n")

# Save IMIS outputs
dir.create("routput", showWarnings = FALSE)
write.csv(imis_res$stat, "routput/imis_stat.csv")
write.csv(imis_res$resample, "routput/imis_resample.csv")
write.csv(imis_res$center, "routput/imis_center.csv")

#===============================================================================
# POSTERIOR ANALYSIS
# Draw samples from posterior and evaluate model predictions
#===============================================================================

# Load posterior samples
posterior_data <- read.csv("routput/imis_resample.csv")
posterior_data <- data.frame(posterior_data[, -1])
posterior.df <- as.matrix(posterior_data)
colnames(posterior.df) <- names(posterior_data)

# Draw random subsample from posterior
smp_size <- 250
set.seed(23)
counter <- sample(seq_len(nrow(posterior.df)), size = smp_size)
samp_posterior <- posterior.df[counter, ]

# Run model with posterior samples (uses psa_run from oud_calibration_model.R)
model_pred <- psa_run(samp_posterior)
colnames(model_pred) <- 1:5
write.csv(model_pred, "routput/model_pred.csv")

#===============================================================================
# VISUALIZE MODEL FIT
# Compare posterior predictions to calibration targets
#===============================================================================

model_pred <- read.csv("routput/model_pred.csv")
model_pred <- model_pred[, -1]
colnames(model_pred) <- 1:5

# Calculate summary statistics
summary_stats <- data.frame(
  Year = colnames(model_pred),
  Mean = sapply(model_pred, mean),
  Lower = sapply(model_pred, function(x) quantile(x, probs = 0.05)),
  Upper = sapply(model_pred, function(x) quantile(x, probs = 0.95)),
  std_lower = target - 1.96 * std_target,
  std_upper = target + 1.96 * std_target,
  Target = target
)

# Plot model fit vs targets
pd <- position_dodge(0.2)
plot_model_rate <- ggplot(summary_stats, aes(x = Year, y = Mean)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, position = pd, color = "blue") +
  geom_point(position = pd, size = 2, color = "black") +
  geom_point(aes(y = Target),
             position = pd, size = 3, shape = 18, color = "red") +
  geom_errorbar(aes(ymin = std_lower, ymax = std_upper),
                width = 0.2, position = pd, color = "cyan3") +
  ylab("Deaths per 100,000 people") +
  xlab("Year") +
  theme_bw()

ggsave("routput/model_rate_with_targets.pdf", plot = plot_model_rate,
       device = "pdf", dpi = 400, width = 8, height = 4, units = "in")
print(plot_model_rate)

#===============================================================================
# PRIOR VS POSTERIOR COMPARISON
# Visualize how calibration updated parameter distributions
#===============================================================================

prior_samples <- data.frame(sample.prior.lhs(smp_size))
prior_samples$type <- 'Prior'
posterior_plot_data <- posterior_data
posterior_plot_data$type <- 'Posterior'

priorpost <- rbind(prior_samples, posterior_plot_data)
priorpost_df <- reshape2::melt(priorpost, id.vars = "type")

median_posterior.df <- data.frame(
  variable = colnames(posterior.df),
  intercept = colMedians(posterior.df)
)

prior_post_plot <- ggplot(priorpost_df, aes(value, fill = type, colour = type)) +
  geom_density(alpha = 0.1) +
  facet_wrap(~variable, scales = "free", ncol = 4) +
  geom_vline(data = median_posterior.df,
             aes(xintercept = intercept), color = "red", linetype = "dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        strip.text.x = element_text(size = 8),
        legend.position = "bottom") +
  labs(x = "Parameter Value", y = "Density")

ggsave("routput/prior_post.pdf", plot = prior_post_plot,
       device = "pdf", dpi = 400, width = 10, height = 8, units = "in")

#===============================================================================
# POSTERIOR SUMMARY STATISTICS
# Calculate credible intervals for calibrated parameters
#===============================================================================

posterior_matrix <- as.matrix(posterior_data[, sapply(posterior_data, is.numeric)])

posterior_summary <- data.frame(
  variable = colnames(posterior_matrix),
  mean = colMeans(posterior_matrix),
  median = colMedians(posterior_matrix),
  sd = apply(posterior_matrix, 2, sd),
  q2.5 = apply(posterior_matrix, 2, quantile, probs = 0.025),
  q97.5 = apply(posterior_matrix, 2, quantile, probs = 0.975)
)

print(posterior_summary)
write.csv(posterior_summary, "routput/posterior_summary.csv", row.names = FALSE)