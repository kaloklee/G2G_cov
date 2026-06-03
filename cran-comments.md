# CRAN submission comments — G2Gcov 0.1.0

## R CMD check results

Checked on: Windows 11 x64, R 4.5.2

0 errors | 0 warnings | 1 note

The single note is:

    checking for future file timestamps ... NOTE
    unable to verify current time

This is a known Windows-specific network issue where R cannot reach a time
server to validate file timestamps. It is not reproducible on Linux/macOS and
does not reflect any problem with the package itself. CRAN's own checks do not
flag this as an issue.

## Platform testing

- Windows 11 x64 (R 4.5.2): 0 errors | 0 warnings | 1 note (timestamp, see above)
- Recommend also running rhub::check_for_cran() before final submission for
  Linux and macOS coverage.

## Notes for CRAN reviewers

* This is a new submission.
* The package implements the Grassia(II)-Geometric (G2G) discrete-time survival
  model with support for both time-invariant and time-varying covariates,
  following the methodology of Fader & Hardie (1997).
* Optimisation examples are wrapped in \donttest{} as they take several seconds;
  all examples complete successfully when run.
* The package uses survival::Surv() only in vignettes and optional tests
  (listed under Suggests). Core estimation functions do not require it.
* The single Suggests package 'dplyr' is used only within vignette code for
  data manipulation demonstrations.

## Downstream dependencies

None — this is a new package with no reverse dependencies.
