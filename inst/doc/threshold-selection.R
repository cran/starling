## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  message  = FALSE,
  warning  = FALSE,
  fig.width  = 8,
  fig.height = 5
)

## ----perch-demo---------------------------------------------------------------
library(starling)

# Simulate a scored pairs object (bimodal: non-matches low, matches high)
set.seed(20260624L)
n_nonmatch <- 800
n_match    <- 200
pairs_pred <- data.frame(
  weights = c(
    rnorm(n_nonmatch, mean = 5,  sd = 3),   # non-match cluster
    rnorm(n_match,    mean = 20, sd = 3)    # match cluster
  )
)

# Run the sensitivity sweep
results <- perch(
  pairs_pred    = pairs_pred,
  n_records_df1 = 250L,         # size of the primary dataset
  thresholds    = seq(5, 28, by = 1),
  report        = TRUE,
  plot          = FALSE          # set TRUE in interactive session
)

## ----plot, fig.cap = "Linkage weight distribution. The threshold line should sit in the valley between the two clusters."----
murmuration_plot(pairs_pred, threshold = 17, show_density = FALSE,
                 palette = "sch")

## ----perch-inside, eval = FALSE-----------------------------------------------
# linked <- murmuration(
#   df1                  = cases_blocked,
#   df2                  = vax_blocked,
#   linkage_type         = "v2c",
#   event_date           = "onset_date",
#   id_var               = "id_var",
#   blocking_var         = "block1",
#   compare_vars         = c("lettername1", "lettername2", "dob", "medicare10"),
#   threshold_value      = 17,
#   perch_before_linking = TRUE     # <-- prints the perch() table mid-linkage
# )

## ----session------------------------------------------------------------------
sessionInfo()

