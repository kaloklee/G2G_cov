# test-g2g-varying.R
# Test script for G2G time-varying covariates model

# Clean environment
rm(list = ls())

# Load required packages
library(survival)
library(tidyverse)

# Load package functions
# Option 1: If developing the package
devtools::load_all()

# Option 2: If package is installed
# library(G2Gcov)

#######################################################
################# Import and Prepare Data #############
#######################################################

# Import data
kb <- read.csv("../data-raw/kb_covars.csv")

# Find the index of the last instance of each person
last_instance_index <- tapply(seq_len(nrow(kb)), kb$id, tail, n = 1)
print(last_instance_index)

# Initialize Event column
kb$Event <- 0

# Create the Event column
kb$Event[last_instance_index] <- ifelse(kb$censor[last_instance_index] == 0, 1, 0)

# View the data
head(kb)

################ Data Manipulation #####################

# Aggregate censor by id
id_df <- aggregate(censor ~ id, data = kb, FUN = max)

# Join and create status variable
kb2 <- kb %>%
  select(-censor) %>%
  inner_join(id_df, by = "id")

# Create status: 0 = censored, 1 = event
kb2$status <- (kb2$censor + 1) %% 2

# Check data
head(kb2, 50)
min(kb2$status)
table(kb2$status)

##########################################################
############## Fit G2G Model #############################
##########################################################

# Start timer
ptm <- proc.time()

# Fit the G2G model with time-varying covariates
# Note: Using the new function name
test_sol <- G2G_varying_MLE(
  fo = Surv(week, status) ~ coupon + anyp,
  data = kb2,
  subject = "id"
)

# Calculate elapsed time
elapsed_time <- proc.time() - ptm
cat("Elapsed time:", elapsed_time["elapsed"], "seconds\n")

##########################################################
############## Display Results ###########################
##########################################################

# Check convergence
cat("\nConvergence code:", test_sol$convergence, "\n")
if (test_sol$convergence == 0) {
  cat("✓ Model converged successfully\n\n")
} else {
  cat("✗ Warning: Model may not have converged properly\n\n")
}

# Display parameter estimates
cat("=== Parameter Estimates ===\n")
print(test_sol$par)

# Display standard errors
cat("\n=== Standard Errors ===\n")
print(test_sol$par_stderr)

# Display confidence intervals
cat("\n=== 95% Confidence Intervals ===\n")
results_table <- data.frame(
  Parameter = c("r", "alpha", "coupon", "anyp"),
  Estimate = test_sol$par,
  SE = test_sol$par_stderr,
  Lower_95 = test_sol$par_lower,
  Upper_95 = test_sol$par_upper
)
print(results_table)

# Display negative log-likelihood
cat("\nNegative log-likelihood:", test_sol$value, "\n")

##########################################################
############## Additional Diagnostics ####################
##########################################################

# Check for any issues
if (any(is.na(test_sol$par))) {
  cat("\n✗ Warning: Some parameters are NA\n")
}

if (any(!is.finite(test_sol$par_stderr))) {
  cat("\n✗ Warning: Some standard errors are not finite\n")
}

# Check if parameters are within CI (they should be)
within_ci <- (test_sol$par >= test_sol$par_lower) & 
             (test_sol$par <= test_sol$par_upper)
if (all(within_ci)) {
  cat("\n✓ All parameters are within their 95% confidence intervals\n")
} else {
  cat("\n✗ Warning: Some parameters are outside their confidence intervals\n")
}

##########################################################
############## Compare Models ############################
##########################################################

# Fit model with only one covariate for comparison
cat("\n=== Fitting Reduced Models for Comparison ===\n")

# Model with only coupon
fit_coupon <- G2G_varying_MLE(
  fo = Surv(week, status) ~ coupon,
  data = kb2,
  subject = "id"
)

# Model with only anyp
fit_anyp <- G2G_varying_MLE(
  fo = Surv(week, status) ~ anyp,
  data = kb2,
  subject = "id"
)

# Compare log-likelihoods
cat("\nModel Comparison (Log-Likelihoods):\n")
cat("Full model (coupon + anyp):", -test_sol$value, "\n")
cat("Coupon only:", -fit_coupon$value, "\n")
cat("Anyp only:", -fit_anyp$value, "\n")

# Likelihood ratio tests
lrt_coupon <- 2 * (-test_sol$value - (-fit_anyp$value))
lrt_anyp <- 2 * (-test_sol$value - (-fit_coupon$value))

cat("\nLikelihood Ratio Tests:\n")
cat("LRT for coupon:", lrt_coupon, 
    "(p-value:", pchisq(lrt_coupon, df = 1, lower.tail = FALSE), ")\n")
cat("LRT for anyp:", lrt_anyp, 
    "(p-value:", pchisq(lrt_anyp, df = 1, lower.tail = FALSE), ")\n")

##########################################################
############## Save Results ##############################
##########################################################

# Optionally save results
# saveRDS(test_sol, "g2g_varying_results.rds")
# write.csv(results_table, "g2g_varying_estimates.csv", row.names = FALSE)

cat("\n=== Test Complete ===\n")
