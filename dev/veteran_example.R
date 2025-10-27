## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----load-data----------------------------------------------------------------
devtools::load_all()
library(survival)

# Load the veteran dataset
data(veteran)
head(veteran)


## ----fit-model----------------------------------------------------------------
# Fit G2G model with age and Karnofsky score as covariates
fit <- G2G_static(Surv(time, status) ~ age + karno, data = veteran)

# View parameter estimates
fit$par


## ----results------------------------------------------------------------------
# Create a table of estimates with confidence intervals
results <- data.frame(
  Parameter = c("r", "alpha"),
  Estimate = fit$par,
  SE = fit$par_stderr,
  Lower_95 = fit$par_lower,
  Upper_95 = fit$par_upper
)

knitr::kable(results, digits = 3)


## ----visualization, fig.width=7, fig.height=4---------------------------------
# Simulate from the estimated gamma distribution
set.seed(123)
gamma_samples <- rgamma(10000, shape = fit$par[1], rate = fit$par[2])
p_samples <- 1 - exp(-gamma_samples)

# Plot density
hist(p_samples, breaks = 50, probability = TRUE,
     main = "Distribution of Event Probability (p)",
     xlab = "p", col = "lightblue", border = "white")

