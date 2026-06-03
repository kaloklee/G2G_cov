#' Fit G2G Model with Static Covariates
#'
#' @param par Parameter vector
#' @param df Data frame with time, status, and covariates
#' @return Negative log-likelihood value
#' @keywords internal
G2G_static_LL <- function(par, df) {
  r <- exp(par[1])
  alpha <- exp(par[2])
  beta <- par[-c(1:2)]
  
  y <- df[, 1]
  status <- df[, 2]
  X <- as.matrix(df[, -c(1:2)])
  
  # Uncensored likelihood
  uncen <- y[which(status == 1)]
  C_u <- exp(X[which(status == 1), ] %*% beta)
  LL_uncen <- logdiffexp(-r * log(1 + C_u * (uncen - 1) / alpha), 
                         -r * log(1 + C_u * uncen / alpha))
  
  # Censored likelihood
  cen <- y[which(status == 0)]
  C_c <- exp(X[which(status == 0), ] %*% beta)
  LL_cen <- -r * log(1 + C_c * cen / alpha)
  
  return(-sum(LL_uncen) - sum(LL_cen))
}

#' Fit a G2G Model with Time-Invariant Covariates
#'
#' @description
#' Fits a Grassia(II)-Geometric (G2G) survival model with time-invariant
#' covariates using maximum likelihood estimation via BFGS. Each subject is
#' represented by a single row. The survivor function at discrete period
#' \eqn{t} is
#' \eqn{S(t) = (1 + t \cdot \exp(\mathbf{x}'\boldsymbol{\beta}) /
#' \alpha)^{-r}},
#' where \eqn{r} and \eqn{\alpha} are the shape and rate of the Gamma
#' mixing distribution and \eqn{\boldsymbol{\beta}} are covariate
#' log-hazard coefficients.
#'
#' @param formula A formula of the form
#'   \code{Surv(time, status) ~ x1 + x2 + ...} where \code{time} is the
#'   observed period, \code{status} is the event indicator (1 = event,
#'   0 = censored), and \code{x1}, \code{x2}, etc. are time-invariant
#'   covariates.
#' @param data A data frame containing all variables referenced in
#'   \code{formula}. One row per subject.
#'
#' @return A list containing the optimization results from
#'   \code{\link[stats]{optim}} with the following additional elements:
#'   \item{par}{Named numeric vector of parameter estimates. The first two
#'     elements are \code{r (shape)} and \code{alpha (rate)}, returned on
#'     the natural (exponentiated) scale. Remaining elements are covariate
#'     log-hazard coefficients.}
#'   \item{par_stderr}{Standard errors derived from the observed Hessian.
#'     \code{NA} when the Hessian is singular or near-singular.}
#'   \item{par_upper}{Upper bounds of 95\% Wald confidence intervals.}
#'   \item{par_lower}{Lower bounds of 95\% Wald confidence intervals.}
#'   \item{value}{Negative log-likelihood at the optimum.}
#'   \item{convergence}{Convergence code from \code{optim} (0 = success).}
#'   \item{hessian}{Observed Hessian matrix at the optimum.}
#'
#' @seealso \code{\link{G2G_varying_MLE}} for the time-varying extension;
#'   \code{\link{plot_p0}} to visualize the fitted baseline propensity
#'   distribution.
#'
#' @importFrom stats optim model.matrix model.frame model.response
#' @examples
#' \donttest{
#'   if (requireNamespace("survival", quietly = TRUE)) {
#'     library(survival)
#'     fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
#'     print(fit$par)
#'     print(fit$par_stderr)
#'   }
#' }
#' @export
G2G_static_MLE <- function(formula, data) {
  
  # Extract variables from formula
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  X <- model.matrix(formula, data = data)[, -1, drop = FALSE]
  
  # Create data frame for optimization
  df <- data.frame(time = y[, 1], status = y[, 2], X)
  
  # Run optimization
  solution <- optim(
    par = c(log(0.5), log(0.5), rep(0, ncol(X))),
    fn = G2G_static_LL,
    df = df,
    method = "BFGS",
    hessian = TRUE
  )
  # Construct parameter names
  par_names <- c("r (shape)", "alpha (rate)", colnames(X))

  names(solution$par) <- par_names
  

  # Calculate standard errors and confidence intervals
  hess_inv <- tryCatch(solve(solution$hessian),
                       error = function(e) NULL)
  solution$par_stderr <- if (is.null(hess_inv)) {
    rep(NA_real_, length(solution$par))
  } else {
    sqrt(diag(hess_inv))
  }
  solution$par_upper <- solution$par + 1.96 * solution$par_stderr
  solution$par_lower <- solution$par - 1.96 * solution$par_stderr
  
  # Exponentiate first two parameters
  exp_idx <- 1:2
  solution$par[exp_idx] <- exp(solution$par[exp_idx])
  solution$par_upper[exp_idx] <- exp(solution$par_upper[exp_idx])
  solution$par_lower[exp_idx] <- exp(solution$par_lower[exp_idx])
  
  return(solution)
}
