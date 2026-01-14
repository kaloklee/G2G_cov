# ==============================================================================
# G2G visualizations (baseline p0, in-sample S(t), hold-out S(t))
# Works with both time-varying and static fitters via auto-detection.
# ==============================================================================

# ---- utilities ---------------------------------------------------------------

# Detect data shape: multiple rows per id => varying (long), else static (wide)
.detect_model_type <- function(data, subject) {
  if (any(table(data[[subject]]) > 1L)) "varying" else "static"
}

# Extract fitted parameters (your MLEs already exp() r and alpha)
.g2g_params <- function(solution) {
  list(
    r     = unname(solution$par[1]),
    alpha = unname(solution$par[2]),
    beta  = unname(solution$par[-(1:2)])
  )
}

# Build per-(id,time) panel for the VARYING case (already long)
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

# Build per-(id,time) panel for the STATIC case (expand each row to 1..T_i)
.build_panel_static <- function(fo, data, subject) {
  time_name   <- all.vars(as.formula(fo))[1]   # total time T_i
  status_name <- all.vars(as.formula(fo))[2]   # event at T_i (0/1)
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

# Compute model survivor on a panel using fitted (r, alpha, beta)
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

# Fit either model (auto unless you force it)
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

# Discrete-time Kaplan–Meier from long person-period data (one event max/id)
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

# Retention (actual) per time = 1 - d_t / n_t
.actual_retention <- function(km_tab) {
  data.frame(time = km_tab$time, R_actual = 1 - km_tab$d / pmax(km_tab$n, 1))
}

# Model retention per time = 1 - mean(p_t | at risk)
.model_retention <- function(aug) {
  last <- stats::aggregate(time ~ id, aug, max)
  names(last)[2] <- "t_last"
  merge_df <- merge(aug, last, by = "id")
  at_risk  <- subset(merge_df, t_last >= time)
  ag <- stats::aggregate(pred ~ time, at_risk, mean)
  names(ag) <- c("time", "R_model_compl")
  ag$R_model <- 1 - ag$R_model_compl
  ag[, c("time", "R_model")]
}

# When we also need per-period event prob p_t (for retention)
.augment_surv_pred <- function(solution, panel_df) {
  par <- .g2g_params(solution)
  X   <- as.matrix(panel_df[, -(1:3), drop = FALSE])
  eta <- as.numeric(drop(X %*% par$beta))
  Ct  <- exp(eta)
  
  ord      <- order(panel_df$id, panel_df$time)
  md       <- panel_df[ord, , drop = FALSE]
  csum     <- stats::ave(Ct[ord], md$id, FUN = cumsum)
  S_t      <- (1 + csum / par$alpha)^(-par$r)
  S_tm1    <- stats::ave(S_t, md$id, FUN = function(v) c(1, head(v, -1)))
  p_t      <- pmax(S_tm1 - S_t, 0)
  
  md$survival <- S_t
  md$pred     <- p_t
  md
}

# ---- PLOT 1: baseline p0 (Gamma) --------------------------------------------

#' Plot baseline p0 distribution implied by a fitted G2G model
#' @param solution object returned by G2G_[varying|static]_MLE()
#' @param probs    inner probability band to shade (default 5%–95%)
#' @param n        grid size
#' @return ggplot object
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

# ---- PLOT 2: in-sample survivor S(t) -----------------------------------------

#' In-sample survivor: fit on TRAIN and plot model S(t) on TRAIN
#' @param fo      formula Surv(time,status) ~ covariates
#' @param train   person-period (varying) or per-id (static) training data
#' @param subject id column name
#' @param model_type "auto" (default), "varying", or "static"
#' @return ggplot object
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

# ---- PLOT 3: hold-out survivor S(t) ------------------------------------------

#' Hold-out survivor: fit on TRAIN and plot model S(t) on TEST
#' @param fo,train,test,subject,model_type see above
#' @return ggplot object
#' @export
plot_survival_holdout_g2g <- function(fo, train, test, subject,
                                      model_type = c("auto","varying","static")) {
  fit <- .fit_g2g(fo, train, subject, model_type)
  panel_te <- if (fit$model_type == "varying")
    .build_panel_varying(fo, test, subject) else .build_panel_static(fo, test, subject)
  aug <- .augment_survival(fit$solution, panel_te)
  ms  <- stats::aggregate(survival ~ time, data = aug, FUN = mean)
  names(ms) <- c("time","S_model")
  
  ggplot2::ggplot(ms, ggplot2::aes(time, S_model)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = sprintf("G2G survivor (hold-out, %s)", fit$model_type),
                  x = "time", y = "S(t)")
}

#' Paper-style plots: % Surviving and Retention with in-sample/hold-out split
#' @param fo       Surv(time,status) ~ covariates
#' @param data     full panel (train + test rows)
#' @param subject  id column name
#' @param split_at time boundary between train (<=) and test (>)
#' @param group    optional grouping column name; NULL => overall only
#' @return list(p_survival, p_retention, solution)
#' @export
plot_insample_holdout_paper <- function(fo, data, subject, split_at, group = NULL,
                                        model_type = c("auto","varying","static")) {
  model_type <- match.arg(model_type)
  
  # Fit on in-sample subset (auto uses varying for your kb_data)
  time_var <- all.vars(as.formula(fo))[1]
  train <- data[data[[time_var]] <= split_at, , drop = FALSE]
  fit <- .fit_g2g(fo, train, subject, model_type)
  
  # Build full panel and compute model survival + per-period p_t
  panel_full <- if (fit$model_type == "varying")
    .build_panel_varying(fo, data, subject) else .build_panel_static(fo, data, subject)
  aug_full <- .augment_surv_pred(fit$solution, panel_full)
  
  # Actual survival (KM) overall or by group
  if (!is.null(group)) {
    grp <- if (is.factor(data[[group]])) data[[group]] else factor(data[[group]])
    gf  <- unique(data.frame(id = data[[subject]], group = grp))
    aug_full <- merge(aug_full, gf, by = "id", all.x = TRUE)
    
    # Survival: Actual (KM) & Model (mean S) per group
    ms_model <- stats::aggregate(survival ~ group + time, aug_full, mean)
    names(ms_model)[3] <- "S_model"
    
    km_list <- lapply(split(aug_full[, c("id","time","status","group")], aug_full$group), .km_survival)
    for (g in names(km_list)) km_list[[g]]$group <- g
    km_df <- do.call(rbind, km_list)
    
    # Retention: Actual (1 - d/n) and Model (1 - mean p_t | at risk) per group
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








# testing

load("data/kb_data.rda")        
kb_data$status <- 1L - kb_data$censor   

fo <- Surv(week, status) ~ coupon + anyp
cut_week <- 12
train <- subset(kb_data, week <= cut_week)
test  <- subset(kb_data, week >  cut_week)

fit <- .fit_g2g(fo, train, subject = "id")     
plot_p0(fit$solution)

plot_survival_insample_g2g(fo, train, subject = "id")

plot_survival_holdout_g2g(fo, train, test, subject = "id")

kb_data$segment <- ifelse(kb_data$anyp >= median(kb_data$anyp, na.rm = TRUE),
                          "High End", "Regular")
res <- plot_insample_holdout_paper(fo, kb_data, subject = "id",
                                   split_at = cut_week, group = "segment")
print(res$p_survival)
print(res$p_retention)
