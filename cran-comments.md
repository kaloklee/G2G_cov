# CRAN submission comments — G2Gcov 0.1.0

## R CMD check results

Checked on: Windows 11 x64, R 4.5.2

0 errors | 0 warnings | 0 notes

## Platform testing

- Windows 11 x64 (R 4.5.2): 0 errors | 0 warnings | 0 notes
- Cross-platform checks via `devtools::check_win_devel()` and
  `rhub::check_for_cran()` should be run before final submission.

## Notes for CRAN reviewers

* This is a new submission.
* The package implements the Grassia(II)-Geometric (G2G) discrete-time survival
  model with support for both time-invariant and time-varying covariates,
  following the methodology of Fader & Hardie (1997).
* All \donttest{} examples run successfully; they are wrapped in \donttest{}
  because numerical optimisation takes several seconds per call.
* During BFGS optimisation, R may internally produce "NaNs produced" warnings
  as the optimizer explores regions of the parameter space where the
  log-likelihood is undefined. These are expected, benign, and are now
  suppressed within the package using withCallingHandlers(). They will not
  appear to end users.
* The package uses survival::Surv() only in vignettes and optional tests
  (listed under Suggests). Core estimation functions do not require it.
* URL and BugReports fields point to https://github.com/kaloklee/G2G_cov.

## Downstream dependencies

None — this is a new package with no reverse dependencies.
