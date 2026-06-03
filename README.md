# G2Gcov

**G2Gcov** implements the Grassia(II)-Geometric (G2G) distribution for
discrete-time survival analysis with both time-invariant and time-varying
covariates, following the methodology of Fader and Hardie (1997).

The G2G model accounts for unobserved heterogeneity across subjects through a
Gamma-distributed frailty. Covariates enter through a log-linear link on the
cumulative hazard, making it straightforward to assess whether dynamic factors
— such as marketing interventions or clinical treatment changes — accelerate or
delay events over time.

---

## Installation

Install the released version from CRAN:

```r
install.packages("G2Gcov")
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("your-org/G2G_cov")
```

---

## Quick Start

### Time-invariant covariates

```r
library(G2Gcov)
library(survival)

data(veteran, package = "survival")

fit <- G2G_static_MLE(
  Surv(time, status) ~ age + karno,
  data = veteran
)

print(fit$par)          # r, alpha, covariate coefficients
print(fit$par_stderr)   # standard errors
plot_p0(fit)            # baseline propensity distribution
```

### Time-varying covariates

```r
library(G2Gcov)
library(survival)

data(kb_data, package = "G2Gcov")
kb_data$status <- 1L - kb_data$censor   # 1 = event, 0 = still active

fit <- G2G_varying_MLE(
  fo      = Surv(week, status) ~ coupon + anyp,
  data    = kb_data,
  subject = "id"
)

print(fit$par)
plot_p0(fit)

# Actual vs. model survival and retention with train/test split
out <- plot_insample_holdout_paper(
  fo       = Surv(week, status) ~ coupon + anyp,
  data     = kb_data,
  subject  = "id",
  split_at = 12L
)
print(out$p_survival)
print(out$p_retention)
```

---

## Data Format

Time-varying models require **person-period (long) format**: one row per
subject per time period.

| Column    | Description |
|-----------|-------------|
| `id`      | Subject identifier (any type) |
| `time`    | Period number (1, 2, 3, …) |
| `status`  | 1 in the period the event occurred; 0 otherwise |
| `x1, x2`  | Time-varying (or static) covariates |

For time-invariant models, one row per subject is sufficient.

A useful tutorial on constructing person-period datasets can be found at
<https://www.rensvandeschoot.com/tutorials/discrete-time-survival/>.

---

## Main Functions

| Function | Description |
|----------|-------------|
| `G2G_static_MLE()` | MLE for G2G model with time-invariant covariates |
| `G2G_varying_MLE()` | MLE for G2G model with time-varying covariates |
| `plot_p0()` | Plot the Gamma distribution over baseline propensity |
| `plot_survival_insample_g2g()` | In-sample mean survivor curve |
| `plot_survival_holdout_g2g()` | Hold-out mean survivor curve |
| `plot_insample_holdout_paper()` | Actual vs. model survival and retention |

---

## Vignettes

```r
vignette("time-invariant-introduction", package = "G2Gcov")
vignette("time-varying-introduction",   package = "G2Gcov")
```

---

## References

Fader, P.S. and Hardie, B.G.S. (1997). How to project customer retention.
*Journal of Interactive Marketing*, 11(1), 61–73.

Fader, P.S. and Hardie, B.G.S. (2007). How to project customer retention
revisited: The role of duration dependence.
*Journal of Interactive Marketing*, 21(1), 76–90.

Lee, K.L., Fader, P.S. and Hardie, B.G.S. (2007). How to project patient
persistency. *Foresight: The International Journal of Applied Forecasting*,
Issue 8, 31–35.

---

## Contact

Questions or comments: kaloklee@gmail.com
