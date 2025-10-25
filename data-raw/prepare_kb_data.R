# data-raw/prepare_kb_data.R
# Script to generate and save package datasets
# Run manually via devtools::load_all() or devtools::document()

# ---- Setup ----
library(readr)      # for read_csv
library(dplyr)      # for cleaning
library(usethis)    # for use_data()

# ---- Read raw data ----
kb_raw <- read_csv("data-raw/kb_covars.csv")

# ---- Clean / Transform ----
kb_data <- kb_raw %>%
  rename_with(tolower)

# ---- Verify ----
stopifnot(!any(is.na(kb_data)))

# ---- Save ----
usethis::use_data(kb_data, overwrite = TRUE)
