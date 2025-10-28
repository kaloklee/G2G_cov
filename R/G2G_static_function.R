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

#' Fit G2G Model with Static Covariates
#'
#' @param formula A formula object (e.g., Surv(time, status) ~ x1 + x2)
#' @param data A data frame containing the variables
#' @return Optimization results with parameter estimates and standard errors
#' @export
#' @importFrom stats aggregate optim as.formula ave model.matrix
#' @examples
#' # Example with veteran dataset
#' # library(survival)
#' # fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
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
  
  # Calculate standard errors and confidence intervals
  solution$par_stderr <- sqrt(diag(solve(solution$hessian)))
  solution$par_upper <- solution$par + 1.96 * solution$par_stderr
  solution$par_lower <- solution$par - 1.96 * solution$par_stderr
  
  # Exponentiate first two parameters
  exp_idx <- 1:2
  solution$par[exp_idx] <- exp(solution$par[exp_idx])
  solution$par_upper[exp_idx] <- exp(solution$par_upper[exp_idx])
  solution$par_lower[exp_idx] <- exp(solution$par_lower[exp_idx])
  
  return(solution)
}
