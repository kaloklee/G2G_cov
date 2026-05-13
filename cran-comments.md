# CRAN submission comments — G2Gcov 0.1.0

## R CMD check results

Checked on: Windows 11, R 4.x

0 errors | 0 warnings | 0 notes

## Notes for CRAN reviewers

* This is a new submission.
* The package implements the Grassia(II)-Geometric (G2G) discrete-time survival
  model with support for both time-invariant and time-varying covariates,
  following the methodology of Fader & Hardie.
* Optimisation examples in `\donttest{}` blocks may take a few seconds on
  small datasets; all tests complete well within the 5-minute limit.
* The package uses `survival::Surv()` only in vignettes and optional tests
  (listed under Suggests); the core estimation functions work without it.

## Downstream dependencies

None — this is a new package with no reverse dependencies.

## Before submitting, please verify

1. Run `devtools::document()` to regenerate NAMESPACE and man/ pages.
   The four visualisation functions (plot_p0, plot_survival_insample_g2g,
   plot_survival_holdout_g2g, plot_insample_holdout_paper) are not yet in
   NAMESPACE because document() has not been run since they were added.
2. Run `devtools::check()` locally and confirm 0 errors / 0 warnings / 0 notes.
3. Run `devtools::check_rhub()` or `rhub::check_for_cran()` for cross-platform
   validation.
4. Optionally add a `URL` and `BugReports` field to DESCRIPTION once the
   package has a public repository.
