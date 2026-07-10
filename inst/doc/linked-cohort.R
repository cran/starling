## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  message  = FALSE,
  warning  = FALSE,
  eval     = FALSE,      # chunks are illustrative: reclin2 not available at build time
  fig.width  = 8,
  fig.height = 5
)

## ----load---------------------------------------------------------------------
# library(starling)
# 
# data(cases_notifiable)
# data(vax_air)
# 
# cat("Cases linelist:      ", nrow(cases_notifiable), "records\n")
# cat("Vaccination register:", nrow(vax_air),           "records\n")
# cat("True matches (known):", sum(!is.na(cases_notifiable$true_link_id)), "\n")

## ----preflight----------------------------------------------------------------
# audit <- preflight(
#   data1        = cases_notifiable,
#   data2        = vax_air,
#   linkage_vars = c("lettername1", "lettername2", "dob", "medicare10"),
#   id_col1      = "id_var",
#   id_col2      = "id_var",
#   date_cols    = c("dob", "onset_date"),
#   medicare_col = "medicare10"
# )

## ----medicare-----------------------------------------------------------------
# # Validate cases
# cases_checked <- check_medicare(cases_notifiable,
#                                 medicare_col = "medicare10",
#                                 output_col   = "medicare_valid",
#                                 verbose      = TRUE)
# 
# # Confirm AIR Medicare numbers are all valid
# vax_checked <- check_medicare(vax_air,
#                               medicare_col = "medicare10",
#                               output_col   = "medicare_valid",
#                               verbose      = TRUE)
# 
# # Replace corrupted Medicare numbers with NA before linkage
# # so they don't negatively score a true match
# cases_checked$medicare10 <- ifelse(
#   cases_checked$medicare_valid == 1L,
#   cases_checked$medicare10,
#   NA_character_
# )

## ----flock--------------------------------------------------------------------
# # Extract birth year for composite blocking
# cases_blocked <- flock(cases_checked,
#                        block1_vars    = "gender",
#                        block2_vars    = "gender",
#                        block3_vars    = "postcode",
#                        birth_year_col = "dob")
# 
# vax_blocked   <- flock(vax_checked,
#                        block1_vars    = "gender",
#                        block2_vars    = "gender",
#                        block3_vars    = "postcode",
#                        birth_year_col = "dob")
# 
# # Summary of blocking key distributions
# cat("block1 (gender) — unique values in cases:",
#     dplyr::n_distinct(cases_blocked$block1), "\n")
# cat("block3 (postcode) — unique values in cases:",
#     dplyr::n_distinct(cases_blocked$block3), "\n")

## ----perch-standalone, eval = FALSE-------------------------------------------
# # This would require running the EM model first —
# # see the murmuration() call below which does this in one step.
# # For standalone use on a pre-scored pairs object:
# 
# # pairs <- reclin2::pair_blocking(cases_blocked, vax_blocked, "block1")
# # reclin2::compare_pairs(pairs,
# #   on = c("lettername1", "lettername2", "dob", "medicare10"),
# #   default_comparator = reclin2::jaro_winkler(0.9), inplace = TRUE)
# # m          <- reclin2::problink_em(
# #   ~ lettername1 + lettername2 + dob + medicare10, data = pairs)
# # pairs_pred <- predict(m, pairs = pairs, add = TRUE)
# #
# # perch(pairs_pred, n_records_df1 = nrow(cases_blocked),
# #       thresholds = seq(8, 25, by = 1))

## ----murmuration--------------------------------------------------------------
# linked <- murmuration(
#   df1                  = cases_blocked,
#   df2                  = vax_blocked,
#   linkage_type         = "v2c",
#   event_date           = "onset_date",
#   id_var               = "id_var",
#   blocking_var         = "block1",
#   compare_vars         = c("lettername1", "lettername2", "dob", "medicare10"),
#   threshold_value      = 17,
#   perch_before_linking = FALSE,   # set TRUE in interactive sessions to inspect
#   days_allowed_before_event = 14,
#   clean_eggs           = TRUE
# )
# 
# cat("Linked records:      ", nrow(linked), "\n")
# cat("With vaccination:    ", sum(!is.na(linked$vax_date_1)), "\n")
# cat("Without vaccination: ", sum( is.na(linked$vax_date_1)), "\n")

## ----plot, fig.cap = "Linkage weight distribution. The threshold (dashed line) should sit in the valley between the two score clusters."----
# # The linked output retains the weights column when clean_eggs = TRUE
# # For the visualisation, we use the weights from the linked data frame
# if ("weights" %in% names(linked)) {
#   murmuration_plot(linked, threshold = 17, show_density = FALSE,
#                    palette = "sch")
# }

## ----validate-----------------------------------------------------------------
# # Recall: proportion of true matches recovered
# true_positives <- sum(
#   !is.na(linked$true_link_id) &
#   !is.na(linked$id_var_df2) &
#   linked$true_link_id == linked$id_var_df2,
#   na.rm = TRUE
# )
# total_true_matches <- sum(!is.na(cases_notifiable$true_link_id))
# recall <- true_positives / total_true_matches
# 
# # Precision: proportion of accepted links that are true matches
# total_links <- sum(!is.na(linked$id_var_df2))
# precision   <- true_positives / total_links
# 
# cat(sprintf("Recall:    %.1f%%  (%d / %d true matches recovered)\n",
#             recall * 100, true_positives, total_true_matches))
# cat(sprintf("Precision: %.1f%%  (%d / %d links are true matches)\n",
#             precision * 100, true_positives, total_links))
# cat(sprintf("F1 score:  %.3f\n",
#             2 * precision * recall / (precision + recall)))

## ----workflow-summary, eval = FALSE-------------------------------------------
# library(starling)
# data(cases_notifiable); data(vax_air)
# 
# # 1. Pre-linkage audit
# preflight(cases_notifiable, vax_air,
#           linkage_vars = c("lettername1", "lettername2", "dob", "medicare10"),
#           medicare_col = "medicare10")
# 
# # 2. Medicare validation — replace invalid numbers with NA
# cases <- check_medicare(cases_notifiable)
# cases$medicare10 <- ifelse(cases$medicare_valid == 1L, cases$medicare10, NA_character_)
# 
# # 3. Blocking variables
# cases <- flock(cases, block1_vars = "gender", birth_year_col = "dob")
# vax   <- flock(vax_air, block1_vars = "gender", birth_year_col = "dob")
# 
# # 4. Link (perch_before_linking = TRUE in interactive sessions)
# linked <- murmuration(cases, vax,
#   linkage_type    = "v2c",
#   event_date      = "onset_date",
#   id_var          = "id_var",
#   blocking_var    = "block1",
#   compare_vars    = c("lettername1", "lettername2", "dob", "medicare10"),
#   threshold_value = 17)
# 
# # 5. Pass to mudnester or bowerbird for downstream analysis

## ----session, eval = TRUE-----------------------------------------------------
sessionInfo()

