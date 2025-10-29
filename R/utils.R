#' Numerically Stable Log(1 - exp(x))
#'
#' Computes log(1 - exp(x)) in a numerically stable way following the 
#' algorithm of Mächler (2012).
#'
#' @param x Numeric vector of values (should be <= 0 for valid results)
#' @return Numeric vector. Returns -Inf when x == 0 and NaN when x > 0
#' @keywords internal
#' @references 
#' Mächler, M. (2012). Accurately Computing log(1 - exp(-|a|)). 
#' \url{https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf}
log1mexp <- function(x) {
  ifelse(
    x > -0.6931472,
    log(-expm1(x)),
    log1p(-exp(x))
  )
}


#' Numerically Stable Log(exp(a) - exp(b))
#'
#' Computes log(exp(a) - exp(b)) in a numerically stable way, avoiding 
#' overflow and underflow issues. Follows the algorithm of Mächler (2012).
#'
#' @param a Numeric vector representing the first value
#' @param b Numeric vector representing the second value
#' @return Numeric vector. Returns -Inf when a == b (including a == b == -Inf),
#'   and NaN when a < b or a == Inf
#' @keywords internal
#' @references 
#' Mächler, M. (2012). Accurately Computing log(1 - exp(-|a|)). 
#' \url{https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf}
logdiffexp <- function(a, b) {
  ifelse(
    (a < Inf) & (a > b),
    a + log1mexp(b - a),
    ifelse((a < Inf) & (a == b), -Inf, NaN)
  )
}
