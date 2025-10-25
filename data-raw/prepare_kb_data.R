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
kb_covars <- kb_raw %>%
  rename_with(tolower)

# ---- Print Columns and Data Types ----
message("--- Column Names (kb_data) ---")
print(names(kb_covars))
message("\n--- Data Types (kb_data) ---")
print(str(kb_covars))
message("---------------------------------\n")

# ---- Verify ----
stopifnot(!any(is.na(kb_covars)))

# ---- Save ----
usethis::use_data(kb_covars, overwrite = TRUE)
