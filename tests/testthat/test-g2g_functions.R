# Tests for utility functions, input validation, G2G MLE, and print method

# --- Utility: log1mexp ---

test_that("log1mexp handles boundary at 0", {
  expect_equal(G2Gcov:::log1mexp(0), -Inf)
})

test_that("log1mexp handles -Inf input", {
  expect_equal(G2Gcov:::log1mexp(-Inf), 0)
})

test_that("log1mexp is mathematically correct at log(0.5)", {
  expect_equal(G2Gcov:::log1mexp(log(0.5)), log(0.5), tolerance = 1e-10)
})

test_that("log1mexp returns NaN for positive input", {
  expect_true(is.nan(G2Gcov:::log1mexp(1)))
})

# --- Utility: logdiffexp ---

test_that("logdiffexp returns -Inf when a == b", {
  expect_equal(G2Gcov:::logdiffexp(1, 1), -Inf)
  expect_equal(G2Gcov:::logdiffexp(-Inf, -Inf), -Inf)
})

test_that("logdiffexp returns NaN when a < b", {
  expect_true(is.nan(G2Gcov:::logdiffexp(0, 1)))
})

test_that("logdiffexp is correct for log(exp(0) - exp(-Inf)) == 0", {
  expect_equal(G2Gcov:::logdiffexp(0, -Inf), 0, tolerance = 1e-10)
})

test_that("logdiffexp is correct for log(exp(log2) - exp(0)) == 0", {
  expect_equal(G2Gcov:::logdiffexp(log(2), 0), 0, tolerance = 1e-10)
})

# --- G2G_static_MLE ---

test_that("G2G_static_MLE returns correct structure (veteran data)", {
  skip_if_not_installed("survival")
  library(survival)

  fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)

  expect_type(fit, "list")
  expect_equal(fit$convergence, 0)
  expect_named(fit$par, c("r (shape)", "alpha (rate)", "age", "karno"))
  expect_true(fit$par[["r (shape)"]] > 0)
  expect_true(fit$par[["alpha (rate)"]] > 0)
  # Standard errors are NA when the Hessian is singular; otherwise positive
  expect_true(all(fit$par_stderr > 0, na.rm = TRUE))
  expect_true(all(fit$par_upper > fit$par_lower, na.rm = TRUE))
})

test_that("G2G_static_MLE CI bounds bracket the point estimates", {
  skip_if_not_installed("survival")
  library(survival)

  fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)

  # Exponentiated params: upper/lower already on original scale for first two
  expect_true(all(fit$par_upper >= fit$par))
  expect_true(all(fit$par >= fit$par_lower))
})

# --- G2G_varying_MLE ---

test_that("G2G_varying_MLE returns correct structure (kb_data)", {
  kb <- kb_data
  kb$status <- 1L - kb$censor

  fit <- G2G_varying_MLE("Surv(week, status) ~ coupon + anyp",
                          data = kb, subject = "id")

  expect_type(fit, "list")
  expect_equal(fit$convergence, 0)
  expect_named(fit$par, c("r (shape)", "alpha (rate)", "coupon", "anyp"))
  expect_true(fit$par[["r (shape)"]] > 0)
  expect_true(fit$par[["alpha (rate)"]] > 0)
  expect_true(all(is.finite(fit$par)))
  # Standard errors are NA when the Hessian is singular; otherwise positive
  expect_true(all(fit$par_stderr > 0, na.rm = TRUE))
})

test_that("G2G_varying_MLE negative log-likelihood is finite", {
  kb <- kb_data
  kb$status <- 1L - kb$censor

  fit <- G2G_varying_MLE("Surv(week, status) ~ coupon + anyp",
                          data = kb, subject = "id")

  expect_true(is.finite(fit$value))
  expect_true(fit$value > 0)
})

# --- Input validation: G2G_varying_MLE ---

test_that("G2G_varying_MLE errors on non-data-frame input", {
  expect_error(
    G2G_varying_MLE("Surv(week, status) ~ coupon", data = list(), subject = "id"),
    "'data' must be a data frame"
  )
})

test_that("G2G_varying_MLE errors when subject column is missing", {
  kb <- kb_data; kb$status <- 1L - kb$censor
  expect_error(
    G2G_varying_MLE("Surv(week, status) ~ coupon", data = kb, subject = "person"),
    "not found in 'data'"
  )
})

test_that("G2G_varying_MLE errors when formula column is missing", {
  kb <- kb_data; kb$status <- 1L - kb$censor
  expect_error(
    G2G_varying_MLE("Surv(week, status) ~ nonexistent", data = kb, subject = "id"),
    "not found in 'data'"
  )
})

test_that("G2G_varying_MLE errors when no events present", {
  kb <- kb_data
  kb$status <- 0L          # force all censored
  expect_error(
    G2G_varying_MLE("Surv(week, status) ~ coupon", data = kb, subject = "id"),
    "No events found"
  )
})

test_that("G2G_varying_MLE errors on NA values in key columns", {
  kb <- kb_data; kb$status <- 1L - kb$censor
  kb$coupon[1] <- NA
  expect_error(
    G2G_varying_MLE("Surv(week, status) ~ coupon", data = kb, subject = "id"),
    "NA values found"
  )
})

# --- Input validation: G2G_static_MLE ---

test_that("G2G_static_MLE errors on non-data-frame input", {
  expect_error(
    G2G_static_MLE(Surv(time, status) ~ x, data = list()),
    "'data' must be a data frame"
  )
})

test_that("G2G_static_MLE errors when formula column is missing", {
  skip_if_not_installed("survival")
  library(survival)
  expect_error(
    G2G_static_MLE(Surv(time, status) ~ nonexistent, data = veteran),
    "not found in 'data'"
  )
})

test_that("G2G_static_MLE errors on NA in formula columns", {
  skip_if_not_installed("survival")
  library(survival)
  vet2 <- veteran; vet2$age[1] <- NA
  expect_error(
    G2G_static_MLE(Surv(time, status) ~ age + karno, data = vet2),
    "NA values found"
  )
})

# --- print.G2Gcov ---

test_that("print.G2Gcov returns the fit object invisibly", {
  kb <- kb_data; kb$status <- 1L - kb$censor
  fit <- G2G_varying_MLE("Surv(week, status) ~ coupon + anyp",
                          data = kb, subject = "id")
  expect_s3_class(fit, "G2Gcov")
  out <- capture.output(result <- print(fit))
  expect_identical(result, fit)
  expect_true(any(grepl("G2G Model Fit", out)))
  expect_true(any(grepl("Log-likelihood", out)))
  expect_true(any(grepl("r \\(shape\\)", out)))
})
