# Tests for kb_data and time-varying MLE behaviour

test_that("kb_data has expected columns and no NAs", {
  expect_s3_class(kb_data, "data.frame")
  expect_true(all(c("id", "week", "num_weeks", "censor", "coupon", "anyp") %in%
                    names(kb_data)))
  expect_false(any(is.na(kb_data)))
})

test_that("kb_data censor column is binary", {
  expect_true(all(kb_data$censor %in% c(0L, 1L)))
})

test_that("kb_data is in person-period format (one row per id-week)", {
  dups <- duplicated(kb_data[, c("id", "week")])
  expect_false(any(dups))
})

test_that("kb_data week values are positive integers", {
  expect_true(all(kb_data$week >= 1L))
  expect_true(all(kb_data$week == as.integer(kb_data$week)))
})

test_that("G2G_varying_MLE survival values are in [0, 1] on kb_data", {
  kb <- kb_data
  kb$status <- 1L - kb$censor

  fit <- G2G_varying_MLE("Surv(week, status) ~ coupon + anyp",
                          data = kb, subject = "id")

  par  <- list(r = fit$par[["r (shape)"]], alpha = fit$par[["alpha (rate)"]],
               beta = fit$par[-(1:2)])
  X    <- as.matrix(kb[, c("coupon", "anyp")])
  Ct   <- exp(as.numeric(X %*% par$beta))
  csum <- ave(Ct, kb$id, FUN = cumsum)
  S_t  <- (1 + csum / par$alpha)^(-par$r)

  expect_true(all(S_t >= 0 - 1e-10))
  expect_true(all(S_t <= 1 + 1e-10))
})

test_that("G2G_varying_MLE survival is monotone non-increasing per subject", {
  kb <- kb_data
  kb$status <- 1L - kb$censor

  fit <- G2G_varying_MLE("Surv(week, status) ~ coupon + anyp",
                          data = kb, subject = "id")

  par  <- list(r = fit$par[["r (shape)"]], alpha = fit$par[["alpha (rate)"]],
               beta = fit$par[-(1:2)])
  X    <- as.matrix(kb[, c("coupon", "anyp")])
  Ct   <- exp(as.numeric(X %*% par$beta))
  csum <- ave(Ct, kb$id, FUN = cumsum)
  S_t  <- (1 + csum / par$alpha)^(-par$r)

  monotone <- tapply(S_t, kb$id, function(v) all(diff(v) <= 1e-10))
  expect_true(all(unlist(monotone)))
})
