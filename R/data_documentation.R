# R/data_documentation.R

#' KB dataset
#'
#' @description
#' Cleaned dataset created by the script `data-raw/prepare_kb_data.R`.
#' The raw input is read from `data-raw/kb_covars.csv` and column names
#' are converted to lower case for consistency.
#'
#' @format
#' A data frame (tibble) with at least the following columns:
#' \describe{
#'   \item{id}{Integer subject identifier.}
#'   \item{week}{Integer week index for the observation.}
#'   \item{num_weeks}{Total weeks observed for the subject (integer).}
#'   \item{censor}{Binary indicator for event observation.}
#'   \item{coupon}{Coupon stock covariate that measures the aggregate effects
#'                  of distributing coupons to potential customers.}
#'   \item{anyp}{Promotions covariate that measures the percent of stores in 
#'                the market with trade promotions within a given week.}
#' }
#'
#' @source
#' Internal raw CSV: \file{data-raw/kb_covars.csv}, processed by
#' \file{data-raw/prepare_kb_data.R}. 
#'
#' @examples
#' head(kb_data)
"kb_data"
