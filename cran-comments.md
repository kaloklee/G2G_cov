# CRAN submission comments — G2Gcov 0.1.0

## R CMD check results

0 errors | 0 warnings | 0 notes

## Platform testing

All platforms returned `Status: OK` with 0 errors, 0 warnings, 0 notes.

| Platform                        | R version | Result |
|---------------------------------|-----------|--------|
| Windows 11 x64 (local)          | R 4.5.2   | OK     |
| Linux x86_64 Ubuntu 24.04 LTS   | R-devel   | OK     |
| macOS arm64 (macos-latest)      | R-devel   | OK     |
| Windows (windows-latest)        | R-devel   | OK     |

Linux, macOS, and Windows R-devel checks were run via R-hub v2 GitHub Actions
(`rhub::rhub_check()`). Local Windows check run via `devtools::check()`.

## Notes for CRAN reviewers

* This is a new submission.
* The package implements the Grassia(II)-Geometric (G2G) discrete-time survival
  model with support for both time-invariant and time-varying covariates,
  following the methodology of Fader & Hardie (1997).
* All \donttest{} examples run successfully; they are wrapped in \donttest{}
  because numerical optimisation takes several seconds per call.
* During BFGS optimisation, R may internally produce "NaNs produced" warnings
  as the optimizer explores regions of the parameter space where the
  log-likelihood is undefined. These are expected, benign, and are suppressed
  within the package using withCallingHandlers(). They do not appear to end
  users and did not produce any warnings in R CMD check across any platform.
* The package uses survival::Surv() only in vignettes and optional tests
  (listed under Suggests). Core estimation functions do not require it.
* URL and BugReports fields point to https://github.com/kaloklee/G2G_cov.

## Downstream dependencies

None — this is a new package with no reverse dependencies.
