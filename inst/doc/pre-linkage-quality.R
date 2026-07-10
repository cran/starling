## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  message  = FALSE,
  warning  = FALSE
)

## ----preflight----------------------------------------------------------------
library(starling)
data(cases_notifiable)
data(vax_air)

audit <- preflight(
  data1        = cases_notifiable,
  data2        = vax_air,
  linkage_vars = c("lettername1", "lettername2", "dob", "medicare10"),
  id_col1      = "id_var",
  id_col2      = "id_var",
  date_cols    = c("dob", "onset_date"),
  medicare_col = "medicare10",
  verbose      = TRUE
)

## ----audit-structure----------------------------------------------------------
# Completeness table
head(audit$completeness)

# Duplicate ID counts
audit$duplicates

# Flags raised (empty = all clear)
audit$flags

## ----medicare-----------------------------------------------------------------
cases_checked <- check_medicare(
  cases_notifiable,
  medicare_col = "medicare10",
  output_col   = "medicare_valid",
  verbose      = TRUE
)

# Distribution of flags
table(cases_checked$medicare_valid, useNA = "always")

## ----fix-medicare-------------------------------------------------------------
cases_checked$medicare10 <- ifelse(
  cases_checked$medicare_valid == 1L,
  cases_checked$medicare10,
  NA_character_
)

## ----flock--------------------------------------------------------------------
# Single-field block (gender only — broadest, lowest specificity)
cases_blocked <- flock(
  cases_checked,
  block1_vars    = "gender",         # block1: gender only
  block2_vars    = "gender",         # block2: same here; use composite in production
  block3_vars    = "postcode",       # block3: postcode (finer)
  birth_year_col = "dob"             # derives birth_year column for composite use
)

vax_blocked <- flock(
  vax_air,
  block1_vars    = "gender",
  block3_vars    = "postcode",
  birth_year_col = "dob"
)

# Inspect blocking key distributions
table(cases_blocked$block1)
head(sort(table(cases_blocked$block3), decreasing = TRUE))

## ----multi-pass, eval = FALSE-------------------------------------------------
# linked_broad <- murmuration(cases_blocked, vax_blocked,
#   blocking_var = "block1", ...)
# linked_fine  <- murmuration(cases_blocked, vax_blocked,
#   blocking_var = "block3", ...)
# 
# # Union, keeping the highest-scoring link per case
# linked_all <- dplyr::bind_rows(linked_broad, linked_fine) |>
#   dplyr::group_by(id_var.x) |>
#   dplyr::slice_max(weights, n = 1, with_ties = FALSE) |>
#   dplyr::ungroup()

## ----summary, eval = FALSE----------------------------------------------------
# library(starling)
# data(cases_notifiable); data(vax_air)
# 
# # 1. Audit
# preflight(cases_notifiable, vax_air,
#           linkage_vars = c("lettername1", "lettername2", "dob", "medicare10"),
#           medicare_col = "medicare10")
# 
# # 2. Fix Medicare
# cases <- check_medicare(cases_notifiable)
# cases$medicare10 <- ifelse(cases$medicare_valid == 1L,
#                             cases$medicare10, NA_character_)
# 
# # 3. Block
# cases <- flock(cases, block1_vars = "gender", birth_year_col = "dob")
# vax   <- flock(vax_air, block1_vars = "gender", birth_year_col = "dob")
# 
# # 4. Link
# linked <- murmuration(cases, vax,
#   linkage_type    = "v2c",
#   event_date      = "onset_date",
#   id_var          = "id_var",
#   blocking_var    = "block1",
#   compare_vars    = c("lettername1", "lettername2", "dob", "medicare10"),
#   threshold_value = 17)

## ----session------------------------------------------------------------------
sessionInfo()

