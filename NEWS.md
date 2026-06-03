# G2Gcov 0.1.0

* Initial CRAN release.
* `G2G_static_MLE()`: maximum likelihood estimation for G2G models with
  time-invariant covariates. Includes input validation with informative error
  messages for missing columns, NA values, and degenerate data.
* `G2G_varying_MLE()`: maximum likelihood estimation for G2G models with
  time-varying covariates in person-period format. Includes the same input
  validation as the static version plus checks for binary status values,
  positive time periods, and the presence of at least one event.
* `print.G2Gcov()`: S3 print method for fitted model objects. Displays
  convergence status, log-likelihood, AIC, and a formatted parameter table
  with standard errors and 95% Wald confidence intervals.
* `plot_p0()`: visualize the implied Gamma distribution over baseline churn
  propensity.
* `plot_survival_insample_g2g()`: in-sample G2G survivor curve.
* `plot_survival_holdout_g2g()`: hold-out G2G survivor curve fitted on
  training data.
* `plot_insample_holdout_paper()`: paper-style actual-vs-model survival and
  retention plots with optional group faceting.
* Includes `kb_data`, a sample person-period dataset used in the package
  vignettes.
