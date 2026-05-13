# G2G visualizations

.detect_model_type <- function(data, subject) {
  if (any(table(data[[subject]]) > 1L)) "varying" else "static"
}

.g2g_params <- function(solution) {
  list(
    r     = unname(solution$par[1]),
    alpha = unname(solution$par[2]),
    beta  = unname(solution$par[-(1:2)])
  )
}

.build_panel_varying <- function(fo, data, subject) {
  time_name   <- all.vars(as.formula(fo))[1]
  status_name <- all.vars(as.formula(fo))[2]
  X <- stats::model.matrix(~ ., data = data[, c(all.vars(as.formula(fo))[-c(1:2)]), drop = FALSE])
  data.frame(
    id     = data[[subject]],
    time   = data[[time_name]],
    status = data[[status_name]],
    X      = X[, -1, drop = FALSE]
  )
}

.build_panel_static <- function(fo, data, subject) {
  time_name   <- all.vars(as.formula(fo))[1]   
  status_name <- all.vars(as.formula(fo))[2]   
  X <- stats::model.matrix(fo, data = data)[, -1, drop = FALSE]
  
  split_idx <- split(seq_len(nrow(data)), data[[subject]])
  rows <- lapply(split_idx, function(ix) {
    i <- ix[1]
    Ti <- as.integer(data[[time_name]][i]); if (is.na(Ti) || Ti < 1L) return(NULL)
    data.frame(
      id     = data[[subject]][i],
      time   = seq_len(Ti),
      status = c(rep(0L, Ti - 1L), as.integer(data[[status_name]][i])),
      X      = X[rep(i, Ti), , drop = FALSE],
      row.names = NULL
    )
  })
  out <- do.call(rbind, rows)
  names(out)[-(1:3)] <- colnames(X)
  out
}

.augment_survival <- function(solution, panel_df) {
  par <- .g2g_params(solution)
  X   <- as.matrix(panel_df[, -(1:3), drop = FALSE])
  eta <- as.numeric(drop(X %*% par$beta))
  Ct  <- exp(eta)
  ord <- order(panel_df$id, panel_df$time)
  S_t <- (1 + stats::ave(Ct[ord], panel_df$id[ord], FUN = cumsum) / par$alpha)^(-par$r)
  panel_df$survival <- S_t[order(order(ord))]
  panel_df
}

.fit_g2g <- function(fo, data, subject, model_type = c("auto","varying","static")) {
  model_type <- match.arg(model_type)
  if (model_type == "auto") model_type <- .detect_model_type(data, subject)
  sol <- if (model_type == "varying") {
    G2G_varying_MLE(fo, data = data, subject = subject)
  } else {
    G2G_static_MLE(fo, data = data)
  }
  list(solution = sol, model_type = model_type)
}


.km_survival <- function(df_long) {
  last <- stats::aggregate(time ~ id, df_long, max)
  names(last)[2] <- "t_last"
  times <- sort(unique(df_long$time))
  at_risk <- sapply(times, function(t) sum(last$t_last >= t))
  d_t <- stats::aggregate(status ~ time, df_long, sum)
  tab <- merge(data.frame(time = times, n = at_risk), d_t, by = "time", all.x = TRUE)
  tab$status[is.na(tab$status)] <- 0
  tab$S_KM <- cumprod(1 - tab$status / pmax(tab$n, 1))
  tab[, c("time", "S_KM", "n", "d" = "status")]
}

.actual_retention <- function(km_tab) {
  data.frame(time = km_tab$time, R_actual = 1 - km_tab$d / pmax(km_tab$n, 1))
}

.model_retention <- function(aug) {
  last <- stats::aggregate(time ~ id, aug, max)
  names(last)[2] <- "t_last"
  merge_df <- merge(aug, last, by = "id")
  at_risk  <- subset(merge_df, t_last >= time)
  
  ag <- stats::aggregate(retention ~ time, at_risk, mean)
  names(ag) <- c("time", "R_model")
  ag
}

.augment_surv_pred <- function(solution, panel_df) {
  par <- .g2g_params(solution)
  X   <- as.matrix(panel_df[, -(1:3), drop = FALSE])
  eta <- as.numeric(drop(X %*% par$beta))
  Ct  <- exp(eta)
  
  ord <- order(panel_df$id, panel_df$time)
  md  <- panel_df[ord, , drop = FALSE]
  
  csum  <- stats::ave(Ct[ord], md$id, FUN = cumsum)
  S_t   <- (1 + csum / par$alpha)^(-par$r)
  S_tm1 <- stats::ave(S_t, md$id, FUN = function(v) c(1, head(v, -1)))
  
  hazard_t    <- pmax((S_tm1 - S_t) / pmax(S_tm1, 1e-12), 0)
  retention_t <- pmax(S_t / pmax(S_tm1, 1e-12), 0)
  
  md$survival   <- S_t
  md$hazard     <- hazard_t
  md$retention  <- retention_t
  md
}


# Suppress R CMD check notes for ggplot2 non-standard evaluation and
# data-frame columns referenced inside subset() / aes().
utils::globalVariables(c("time", "S_KM", "S_model", "R_actual", "R_model", "t_last"))

#' @importFrom utils head
NULL

# Baseline p0 (gamma)

#' Plot the Baseline Propensity Distribution of a Fitted G2G Model
#'
#' @description
#' Visualises the gamma distribution over baseline churn propensity implied
#' by the fitted \code{r} (shape) and \code{alpha} (rate) parameters.
#'
#' @param solution Object returned by \code{\link{G2G_varying_MLE}} or
#'   \code{\link{G2G_static_MLE}}.
#' @param probs Numeric vector of length 2 giving the lower and upper
#'   probability bounds for the shaded central band.
#'   Default is \code{c(0.05, 0.95)}.
#' @param n Integer grid size for the density curve. Default is \code{400}.
#'
#' @return A \code{ggplot} object.
#' @export

plot_p0 <- function(solution, probs = c(0.05, 0.95), n = 400) {
  par <- .g2g_params(solution)
  lo  <- stats::qgamma(probs[1], shape = par$r, rate = par$alpha)
  hi  <- stats::qgamma(probs[2], shape = par$r, rate = par$alpha)
  span <- hi - lo
  x <- seq(max(0, lo - 0.25 * span), hi + 0.25 * span, length.out = n)
  y <- stats::dgamma(x, shape = par$r, rate = par$alpha)
  mean_ <- par$r / par$alpha
  mode_ <- if (par$r > 1) (par$r - 1) / par$alpha else NA_real_
  
  df <- data.frame(x = x, y = y)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_line() +
    ggplot2::geom_area(data = subset(df, x >= lo & x <= hi), alpha = 0.25) +
    ggplot2::geom_vline(xintercept = mean_) +
    { if (!is.na(mode_)) ggplot2::geom_vline(xintercept = mode_, linetype = 2) } +
    ggplot2::labs(
      title = expression(paste("Baseline ", p[0], " ~ Gamma(r, ", alpha, ")")),
      subtitle = sprintf("mean = %.3g; %s", mean_,
                         if (is.na(mode_)) "mode = NA" else sprintf("mode = %.3g", mode_)),
      x = "propensity", y = "density"
    )
  p
}

# In-sample survivor (s(t))

#' In-Sample Survivor Curve for a Fitted G2G Model
#'
#' @description
#' Fits the G2G model on the supplied training data and plots the mean
#' predicted survivor function S(t) over the training period.
#'
#' @param fo Formula of the form \code{Surv(time, status) ~ covariates}.
#' @param train Data frame used for fitting: person-period format for
#'   time-varying models, one row per subject for static models.
#' @param subject Character string naming the subject ID column.
#' @param model_type One of \code{"auto"} (default), \code{"varying"}, or
#'   \code{"static"}. When \code{"auto"}, the type is detected from the
#'   data structure.
#'
#' @return A \code{ggplot} object.
#' @export
plot_survival_insample_g2g <- function(fo, train, subject,
                                       model_type = c("auto","varying","static")) {
  fit <- .fit_g2g(fo, train, subject, model_type)
  panel <- if (fit$model_type == "varying")
    .build_panel_varying(fo, train, subject) else .build_panel_static(fo, train, subject)
  aug <- .augment_survival(fit$solution, panel)
  ms  <- stats::aggregate(survival ~ time, data = aug, FUN = mean)
  names(ms) <- c("time","S_model")
  
  ggplot2::ggplot(ms, ggplot2::aes(time, S_model)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = sprintf("G2G survivor (in-sample, %s)", fit$model_type),
                  x = "time", y = "S(t)")
}

# Hold-out survivor (s(t))

#' Hold-Out Survivor Curve for a Fitted G2G Model
#'
#' @description
#' Fits the G2G model on the training data and plots the mean predicted
#' survivor function S(t) evaluated on the hold-out test period.
#'
#' @param fo Formula of the form \code{Surv(time, status) ~ covariates}.
#' @param train Data frame used for fitting.
#' @param test Data frame covering the hold-out period.
#' @param subject Character string naming the subject ID column.
#' @param model_type One of \code{"auto"} (default), \code{"varying"}, or
#'   \code{"static"}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_survival_holdout_g2g <- function(fo, train, test, subject,
                                      model_type = c("auto","varying","static")) {
  fit <- .fit_g2g(fo, train, subject, model_type)
  
  full <- rbind(train, test)
  
  panel_full <- if (fit$model_type == "varying")
    .build_panel_varying(fo, full, subject) else .build_panel_static(fo, full, subject)
  
  aug_full <- .augment_survival(fit$solution, panel_full)
  
  time_var <- all.vars(as.formula(fo))[1]
  split_at <- max(train[[time_var]], na.rm = TRUE)
  
  aug_te <- subset(aug_full, time > split_at)
  
  ms <- stats::aggregate(survival ~ time, data = aug_te, FUN = mean)
  names(ms) <- c("time","S_model")
  
  ggplot2::ggplot(ms, ggplot2::aes(time, S_model)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = sprintf("G2G survivor (hold-out, %s)", fit$model_type),
                  x = "time", y = "S(t)")
}


#' Paper-Style Actual-vs-Model Survival and Retention Plot
#'
#' @description
#' Fits the G2G model on the training portion of the data (rows where
#' \code{time <= split_at}) and produces two plots comparing the
#' Kaplan-Meier actual curve against G2G model predictions: one for
#' percent surviving and one for period-over-period retention rate.
#' An optional grouping variable produces faceted panels.
#'
#' @param fo Formula of the form \code{Surv(time, status) ~ covariates}.
#' @param data Full panel data frame containing both training and test rows.
#' @param subject Character string naming the subject ID column.
#' @param split_at Numeric scalar. Observations with time at or below this
#'   value are used for training; the remainder form the hold-out set.
#' @param group Optional character string naming a grouping column for
#'   faceted plots. Use \code{NULL} (default) for an overall-only plot.
#' @param model_type One of \code{"auto"} (default), \code{"varying"}, or
#'   \code{"static"}.
#'
#' @return A named list with elements \code{p_survival} (a \code{ggplot}
#'   of percent surviving), \code{p_retention} (a \code{ggplot} of
#'   retention rate), and \code{solution} (the fitted model object).
#' @export
plot_insample_holdout_paper <- function(fo, data, subject, split_at, group = NULL,
                                        model_type = c("auto","varying","static")) {
  model_type <- match.arg(model_type)
  
  time_var <- all.vars(as.formula(fo))[1]
  train <- data[data[[time_var]] <= split_at, , drop = FALSE]
  fit <- .fit_g2g(fo, train, subject, model_type)
  
  panel_full <- if (fit$model_type == "varying")
    .build_panel_varying(fo, data, subject) else .build_panel_static(fo, data, subject)
  aug_full <- .augment_surv_pred(fit$solution, panel_full)
  
  if (!is.null(group)) {
    grp <- if (is.factor(data[[group]])) data[[group]] else factor(data[[group]])
    gf  <- unique(data.frame(id = data[[subject]], group = grp))
    aug_full <- merge(aug_full, gf, by = "id", all.x = TRUE)
    
    ms_model <- stats::aggregate(survival ~ group + time, aug_full, mean)
    names(ms_model)[3] <- "S_model"
    
    km_list <- lapply(split(aug_full[, c("id","time","status","group")], aug_full$group), .km_survival)
    for (g in names(km_list)) km_list[[g]]$group <- g
    km_df <- do.call(rbind, km_list)
    
    ar_list <- lapply(split(km_df, km_df$group), .actual_retention)
    for (g in names(ar_list)) ar_list[[g]]$group <- g
    ar_df <- do.call(rbind, ar_list)
    
    mr_list <- lapply(split(aug_full, aug_full$group), .model_retention)
    for (g in names(mr_list)) mr_list[[g]]$group <- g
    mr_df <- do.call(rbind, mr_list)
    
    surv_df <- merge(km_df[, c("time","S_KM","group")], ms_model, by = c("time","group"), all = TRUE)
    ret_df  <- merge(ar_df, mr_df, by = c("time","group"), all = TRUE)
  } else {
    ms_model <- stats::aggregate(survival ~ time, aug_full, mean)
    names(ms_model) <- c("time","S_model")
    km_df  <- .km_survival(aug_full[, c("id","time","status")])
    surv_df <- merge(km_df[, c("time","S_KM")], ms_model, by = "time", all = TRUE)
    
    ar_df <- .actual_retention(km_df)
    mr_df <- .model_retention(aug_full)
    ret_df <- merge(ar_df, mr_df, by = "time", all = TRUE)
  }
  
  
  # Plots
  p_surv <- ggplot2::ggplot(surv_df, ggplot2::aes(x = time)) +
    { if (!is.null(group)) ggplot2::facet_wrap(~ group, nrow = 1) } +
    ggplot2::geom_line(ggplot2::aes(y = 100 * S_KM, linetype = "Actual"), linewidth = 0.9) +
    ggplot2::geom_line(ggplot2::aes(y = 100 * S_model, linetype = "G2G"), linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = split_at, linetype = 2) +
    ggplot2::scale_linetype_manual(values = c("Actual" = "solid", "G2G" = "dashed"), name = NULL) +
    ggplot2::labs(x = "Tenure (time)", y = "% Surviving") +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right")
  
  p_ret <- ggplot2::ggplot(ret_df, ggplot2::aes(x = time)) +
    { if (!is.null(group)) ggplot2::facet_wrap(~ group, nrow = 1) } +
    ggplot2::geom_line(ggplot2::aes(y = R_actual, linetype = "Actual"), linewidth = 0.9) +
    ggplot2::geom_line(ggplot2::aes(y = R_model,  linetype = "G2G"),    linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = split_at, linetype = 2) +
    ggplot2::scale_linetype_manual(values = c("Actual" = "solid", "G2G" = "dashed"), name = NULL) +
    ggplot2::labs(x = "Time", y = "Retention Rate") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right")
  
  list(p_survival = p_surv, p_retention = p_ret, solution = fit$solution)
}


