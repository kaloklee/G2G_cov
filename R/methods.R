#' Print a Fitted G2G Model
#'
#' @description
#' Displays a formatted summary of a G2G model fit returned by
#' \code{\link{G2G_static_MLE}} or \code{\link{G2G_varying_MLE}}, including
#' convergence status, log-likelihood, AIC, and a parameter table with
#' standard errors and 95\% Wald confidence intervals.
#'
#' @param x A fitted model object of class \code{"G2Gcov"}.
#' @param digits Integer. Number of significant digits to print. Default 4.
#' @param ... Further arguments passed to or from other methods (unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{G2G_static_MLE}}, \code{\link{G2G_varying_MLE}}
#' @examples
#' \donttest{
#'   data(kb_data)
#'   kb_data$status <- 1L - kb_data$censor
#'   fit <- G2G_varying_MLE("Surv(week, status) ~ coupon + anyp",
#'                           data = kb_data, subject = "id")
#'   print(fit)
#' }
#' @export
print.G2Gcov <- function(x, digits = 4, ...) {
  conv_msg <- if (x$convergence == 0L) {
    "Converged (code 0)"
  } else {
    paste0("Did NOT converge (code ", x$convergence,
           ") -- interpret estimates with caution")
  }

  k   <- length(x$par)
  aic <- round(2 * k + 2 * x$value, digits)
  ll  <- round(-x$value, digits)

  cat("G2G Model Fit\n")
  cat(strrep("-", 52), "\n")
  cat(sprintf("%-18s %s\n", "Convergence:",    conv_msg))
  cat(sprintf("%-18s %s\n", "Log-likelihood:", ll))
  cat(sprintf("%-18s %s\n", "AIC:",            aic))
  cat(sprintf("%-18s %d\n", "Parameters:",     k))
  cat("\nCoefficients:\n")

  tab <- data.frame(
    Estimate  = round(x$par,        digits),
    Std.Error = round(x$par_stderr, digits),
    Lower.95  = round(x$par_lower,  digits),
    Upper.95  = round(x$par_upper,  digits),
    row.names = names(x$par),
    check.names = FALSE
  )
  print(tab)

  if (anyNA(x$par_stderr))
    cat("\nNote: some standard errors are NA (Hessian is singular at optimum).\n")

  invisible(x)
}
