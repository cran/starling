## ----setup, include = FALSE---------------------------------------------------
# The worked data (the clutch datasets) live in mudnester; if it is not
# installed the vignette stops cleanly after this point.
has_mudnester <- requireNamespace("mudnester", quietly = TRUE) &&
  utils::packageVersion("mudnester") >= "0.9.0"   # the clutch datasets arrived in 0.9.0
knitr::opts_chunk$set(
  fig.retina = 1,   # 1x figures: keeps the built HTML small for CRAN
  collapse   = TRUE,
  comment    = "#>",
  message    = FALSE,
  warning    = FALSE,
  fig.width  = 8,
  fig.height = 4.5,
  eval       = has_mudnester
)

## ----needs-mudnester, echo = FALSE, eval = !has_mudnester, results = "asis"----
# cat("This vignette uses the clutch datasets from the mudnester package. Install mudnester (>= 0.9.0) to run it.")
# knitr::knit_exit()

## ----load---------------------------------------------------------------------
library(starling)
library(mudnester)   # clean_the_nest() and the clutch_* datasets
library(dplyr)
library(ggplot2)

aviary_cols <- c(satin = "#325D88", chestnut = "#86251C", wattle = "#D8A31A",
                 evergreen = "#08403B", mud = "#A67D4C", slate = "#3E3F3A")

## ----prepare------------------------------------------------------------------
cases_clean <- clutch_cases |>
  mutate(onset_imputed = is.na(onset_date),
         onset_date    = coalesce(onset_date, specimen_date)) |>
  clean_the_nest(
    data_type = "cases", drop_eggs = TRUE,
    id_var = "notification_id", diagnosis = "condition",
    lettername1 = "first_name", lettername2 = "surname", dob = "date_of_birth",
    gender = "sex", postcode = "postcode", medicare = "medicare_number",
    onset_date = "onset_date", keep_vars = c("truth_id", "onset_imputed")
  )

## ----prepare-air--------------------------------------------------------------
air_clean <- clutch_air |>   # all antigens: every AIR person, for grading identity links
  clean_the_nest(
    data_type = "vaccination", lie_nest_flat = TRUE, drop_eggs = TRUE,
    id_var = "air_id", lettername1 = "given_name", lettername2 = "family_name",
    dob = "date_of_birth", gender = "sex", postcode = "postcode",
    medicare = "medicare_card", medicare_irn = "medicare_irn",
    vax_type = "antigen", vax_date = "date_of_service", keep_vars = "truth_id"
  )

# Invalid Medicare numbers become NA in both comparison fields, then
# block1: gender + birth year (used throughout); block2: birth year alone,
# for the comparison that scores gender instead of blocking on it
cases_blocked <- cases_clean |>
  check_medicare(verbose = FALSE) |>
  mutate(medicare10 = if_else(medicare_valid %in% 0L, NA_character_, medicare10),
         medicare09 = if_else(medicare_valid %in% 0L, NA_character_, medicare09)) |>
  flock(block1_vars = c("gender", "birth_year"), block2_vars = "birth_year",
        birth_year_col = "dob")

air_blocked <- air_clean |>
  check_medicare(verbose = FALSE) |>
  mutate(medicare10 = if_else(medicare_valid %in% 0L, NA_character_, medicare10),
         medicare09 = if_else(medicare_valid %in% 0L, NA_character_, medicare09)) |>
  flock(block1_vars = c("gender", "birth_year"), block2_vars = "birth_year",
        birth_year_col = "dob")

n_true <- cases_blocked |>
  filter(truth_id %in% air_blocked$truth_id) |>
  nrow()
c(cases = nrow(cases_blocked), air_people = nrow(air_blocked),
  cases_with_a_true_air_record = n_true)

## ----inspect-4----------------------------------------------------------------
review_4 <- murmuration(
  data_x             = cases_blocked,   # every notification
  data_y             = air_blocked,     # every AIR person
  linkage_type       = "cohort2event",
  event_type         = "vaccination",
  id_var             = "id_var",
  blocking_var       = "block1",        # gender + birth year
  compare_vars       = c("lettername1", "lettername2", "dob", "medicare10"),
  threshold_value    = 8,               # only marks the pairs; nothing is decided yet
  inspect_matches    = TRUE,            # stop after scoring; event_date not needed
  inspect_extra_vars = c("truth_id", "medicare09", "medicare_irn")  # context, not scored
)

medicare_levels <- c("Same Medicare number",
                     "Same account, card issue digit differs",
                     "Different Medicare number",
                     "Medicare missing on one or both records")

review_4 <- review_4 |>
  mutate(
    pair_type = if_else(truth_id.x == truth_id.y, "True match", "Non-match"),
    medicare_status = case_when(
      is.na(medicare10.x) | is.na(medicare10.y) ~ medicare_levels[4],
      medicare10.x == medicare10.y              ~ medicare_levels[1],
      medicare09.x == medicare09.y              ~ medicare_levels[2],
      TRUE                                      ~ medicare_levels[3]
    ),
    medicare_status = factor(medicare_status, levels = medicare_levels)
  )

c(candidate_pairs = nrow(review_4),
  true_pairs      = sum(review_4$pair_type == "True match"),
  max_weight      = round(max(review_4$weights), 2))

## ----perch-table--------------------------------------------------------------
sweep_4 <- perch(review_4,
                 n_records_primary = nrow(cases_blocked),
                 thresholds        = seq(4, 20, by = 1),
                 report            = FALSE,
                 plot              = FALSE)

sweep_4 |>
  select(threshold, n_above, pct_above, n_clerical, link_rate, reference) |>
  knitr::kable(caption = "Accepted pairs and clerical burden by threshold (4 comparison variables)")

clerical_8  <- sweep_4 |> filter(threshold == 8)  |> pull(n_clerical)
clerical_16 <- sweep_4 |> filter(threshold == 16) |> pull(n_clerical)

## ----perch-plot, fig.height = 5-----------------------------------------------
perch(review_4,
      n_records_primary = nrow(cases_blocked),
      thresholds        = seq(4, 20, by = 0.5),
      report            = FALSE)

## ----plot, fig.cap = "Candidate-pair weights by blocking stratum (gender + birth year; the busiest strata shown, the rest pooled), with the threshold at 8."----
murmuration_plot(review_4, max_levels = 10, palette = "sch")

## ----medicare-split, fig.height = 7, fig.cap = "Weights of true matches and non-matches, by Medicare agreement. Solid line: threshold 8; dashed line: the default 16. Counts on a log scale."----
review_4 |>
  as.data.frame() |>
  filter(weights > 0) |>
  ggplot(aes(x = weights, fill = pair_type)) +
  geom_histogram(binwidth = 0.5, position = "identity", alpha = 0.75) +
  geom_vline(xintercept = 8, colour = aviary_cols[["chestnut"]]) +
  geom_vline(xintercept = 16, linetype = "dashed", colour = aviary_cols[["chestnut"]]) +
  scale_y_log10() +
  scale_fill_manual(values = c("True match" = aviary_cols[["evergreen"]],
                               "Non-match" = aviary_cols[["wattle"]])) +
  facet_wrap(~ medicare_status, ncol = 1) +
  labs(x = "Linkage weight", y = "Pairs (log scale)", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

## ----medicare-clusters--------------------------------------------------------
review_4 |>
  as.data.frame() |>
  mutate(is_true = pair_type == "True match") |>
  group_by(medicare_status, .drop = FALSE) |>
  summarise(
    true_pairs        = sum(is_true),
    true_median_wt    = round(median(weights[is_true]), 1),
    true_5th_pct_wt   = round(unname(quantile(weights[is_true], 0.05)), 1),
    true_ge_16_pct    = round(100 * mean(weights[is_true] >= 16), 1),
    non_matches_ge_8  = sum(!is_true & weights >= 8),
    non_matches_ge_12 = sum(!is_true & weights >= 12)
  ) |>
  knitr::kable(caption = "Three clusters of true matches (NA: no pairs in that group)")

## ----grade-4------------------------------------------------------------------
# The highest-weighted candidate per case is the link at every threshold it
# clears, so grading a threshold is just counting top candidates above it.
top_4 <- review_4 |>
  as.data.frame() |>
  arrange(desc(weights)) |>
  distinct(id_var.x, .keep_all = TRUE) |>
  transmute(id_var.x, weights, correct = truth_id.x == truth_id.y)

grade_4 <- tibble(threshold = seq(5, 18, by = 0.5)) |>
  cross_join(top_4) |>
  filter(weights >= threshold) |>
  group_by(threshold) |>
  summarise(links = n(), correct = sum(correct)) |>
  mutate(precision = correct / links,
         recall    = correct / n_true,
         f1        = 2 * precision * recall / (precision + recall))

best_4    <- grade_4 |> slice_max(f1, n = 1, with_ties = FALSE) |> pull(threshold)
at_best_4 <- grade_4 |> filter(threshold == best_4)
at_16_4   <- grade_4 |> filter(threshold == 16)

grade_4 |>
  filter(threshold %in% c(5, 6, 7, 7.5, 8, 9, 10, 11, 12, 14, 16, 17)) |>
  knitr::kable(digits = 3,
               caption = "Precision, recall and F1 by threshold (4 comparison variables)")

## ----grade-plot, fig.height = 4-----------------------------------------------
grade_long <- bind_rows(
  grade_4 |> transmute(threshold, measure = "Precision", value = precision),
  grade_4 |> transmute(threshold, measure = "Recall",    value = recall),
  grade_4 |> transmute(threshold, measure = "F1",        value = f1)
)

ggplot(grade_long, aes(x = threshold, y = value, colour = measure)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = best_4, colour = aviary_cols[["evergreen"]]) +
  geom_vline(xintercept = 16, colour = aviary_cols[["chestnut"]], linetype = "dashed") +
  annotate("text", x = best_4, y = 0.55, label = paste("best F1 at", best_4),
           hjust = -0.1, colour = aviary_cols[["evergreen"]]) +
  annotate("text", x = 16, y = 0.55, label = "default 16", hjust = 1.1,
           colour = aviary_cols[["chestnut"]]) +
  scale_colour_manual(values = c(Precision = aviary_cols[["satin"]],
                                 Recall = aviary_cols[["wattle"]],
                                 F1 = aviary_cols[["slate"]])) +
  labs(x = "Threshold", y = NULL, colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

## ----nvars-score--------------------------------------------------------------
# threshold_value is left at its default (16) in every call
review_3 <- murmuration(
  data_x             = cases_blocked,
  data_y             = air_blocked,
  linkage_type       = "cohort2event",
  event_type         = "vaccination",
  id_var             = "id_var",
  blocking_var       = "block1",       # gender + birth year
  compare_vars       = c("lettername1", "lettername2", "dob"),
  inspect_matches    = TRUE,           # score only
  inspect_extra_vars = "truth_id"      # for grading
)
review_4_m09 <- murmuration(
  data_x             = cases_blocked,
  data_y             = air_blocked,
  linkage_type       = "cohort2event",
  event_type         = "vaccination",
  id_var             = "id_var",
  blocking_var       = "block1",
  compare_vars       = c("lettername1", "lettername2", "dob", "medicare09"),
  inspect_matches    = TRUE,
  inspect_extra_vars = "truth_id"
)
review_5_pc <- murmuration(
  data_x             = cases_blocked,
  data_y             = air_blocked,
  linkage_type       = "cohort2event",
  event_type         = "vaccination",
  id_var             = "id_var",
  blocking_var       = "block1",
  compare_vars       = c("lettername1", "lettername2", "dob", "medicare10", "postcode"),
  inspect_matches    = TRUE,
  inspect_extra_vars = c("truth_id", "medicare09")
)
review_5_m09 <- murmuration(
  data_x             = cases_blocked,
  data_y             = air_blocked,
  linkage_type       = "cohort2event",
  event_type         = "vaccination",
  id_var             = "id_var",
  blocking_var       = "block1",
  compare_vars       = c("lettername1", "lettername2", "dob", "medicare09", "postcode"),
  inspect_matches    = TRUE,
  inspect_extra_vars = c("truth_id", "medicare10")
)
review_5_sex <- murmuration(
  data_x             = cases_blocked,
  data_y             = air_blocked,
  linkage_type       = "cohort2event",
  event_type         = "vaccination",
  id_var             = "id_var",
  blocking_var       = "block2",       # birth year only
  compare_vars       = c("lettername1", "lettername2", "dob", "medicare10", "gender"),
  inspect_matches    = TRUE,
  inspect_extra_vars = "truth_id"
)

## ----nvars-grade--------------------------------------------------------------
set_names <- c("3: names + DOB",
               "4: + medicare10",
               "4: + medicare09",
               "5: + medicare10 + postcode",
               "5: + medicare09 + postcode",
               "5: + medicare10 + gender (year blocks)")

# every scored pair from the six runs, labelled by variable set
all_pairs <- bind_rows(
  review_3     |> as.data.frame() |> mutate(variables = set_names[1]),
  review_4     |> as.data.frame() |> mutate(variables = set_names[2]),
  review_4_m09 |> as.data.frame() |> mutate(variables = set_names[3]),
  review_5_pc  |> as.data.frame() |> mutate(variables = set_names[4]),
  review_5_m09 |> as.data.frame() |> mutate(variables = set_names[5]),
  review_5_sex |> as.data.frame() |> mutate(variables = set_names[6])
) |>
  select(variables, id_var.x, truth_id.x, truth_id.y, weights) |>
  mutate(variables = factor(variables, levels = set_names))

set_sizes <- all_pairs |>
  group_by(variables) |>
  summarise(candidate_pairs = n(), max_weight = round(max(weights), 1))

# the top candidate per case, in each set
top_pairs <- all_pairs |>
  arrange(desc(weights)) |>
  distinct(variables, id_var.x, .keep_all = TRUE) |>
  mutate(correct = truth_id.x == truth_id.y)

thresholds <- tibble(threshold = seq(2, 24, by = 0.5))

curves <- top_pairs |>
  cross_join(thresholds) |>
  filter(weights >= threshold) |>
  group_by(variables, threshold) |>
  summarise(links = n(), correct = sum(correct), .groups = "drop") |>
  right_join(cross_join(tibble(variables = factor(set_names, levels = set_names)), thresholds),
             by = c("variables", "threshold")) |>
  mutate(links     = coalesce(links, 0L),       # no pair above this threshold
         correct   = coalesce(correct, 0L),
         precision = if_else(links > 0, correct / links, NA_real_),
         recall    = correct / n_true,
         f1        = 2 * precision * recall / (precision + recall)) |>
  arrange(variables, threshold)

best_by_set <- curves |>
  filter(!is.na(f1)) |>
  group_by(variables) |>
  slice_max(f1, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(variables,
            best_f1_threshold = threshold,
            recall_at_best    = round(100 * recall, 1),
            precision_at_best = round(100 * precision, 1))

nvars_table <- set_sizes |>
  left_join(curves |>
              filter(threshold == 16) |>
              transmute(variables,
                        recall_at_16    = round(100 * recall, 1),
                        precision_at_16 = round(100 * precision, 1)),
            by = "variables") |>
  left_join(curves |>
              filter(threshold == 17) |>
              transmute(variables, recall_at_17 = round(100 * recall, 1)),
            by = "variables") |>
  left_join(best_by_set, by = "variables") |>
  arrange(variables)

knitr::kable(nvars_table, format.args = list(big.mark = ","),
             caption = "Recall and precision (%) at 16 and 17, and at each set's best-F1 threshold")

## ----nvars-plot, fig.height = 6.5, fig.cap = "Recall and precision by threshold for each comparison-variable set, keeping the best-weighted candidate per case. Dashed line: the default 16."----
curves_long <- bind_rows(
  curves |> transmute(variables, threshold, measure = "Recall (%)",    value = 100 * recall),
  curves |> transmute(variables, threshold, measure = "Precision (%)", value = 100 * precision)
) |>
  filter(!is.na(value))

ggplot(curves_long, aes(x = threshold, y = value, colour = variables)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 16, linetype = "dashed", colour = aviary_cols[["slate"]]) +
  facet_wrap(~ measure, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = unname(aviary_cols[c("slate", "satin", "evergreen",
                                                    "wattle", "mud", "chestnut")])) +
  labs(x = "Threshold", y = NULL, colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 2))

## ----nvars-medicare-----------------------------------------------------------
# the true pairs of the two five-variable (postcode) runs, by Medicare agreement
true_5 <- bind_rows(
  review_5_pc  |> as.data.frame() |> mutate(medicare_field = "medicare10"),
  review_5_m09 |> as.data.frame() |> mutate(medicare_field = "medicare09")
) |>
  filter(truth_id.x == truth_id.y) |>
  mutate(
    medicare_status = case_when(
      is.na(medicare10.x) | is.na(medicare10.y) ~ medicare_levels[4],
      medicare10.x == medicare10.y              ~ medicare_levels[1],
      medicare09.x == medicare09.y              ~ medicare_levels[2],
      TRUE                                      ~ medicare_levels[3]
    ),
    medicare_status = factor(medicare_status, levels = medicare_levels)
  )

by_status_10 <- true_5 |>
  filter(medicare_field == "medicare10") |>
  group_by(medicare_status, .drop = FALSE) |>
  summarise(true_pairs           = n(),
            median_wt_medicare10 = round(median(weights), 1),
            pct_ge_16_medicare10 = round(100 * mean(weights >= 16), 1))

by_status_09 <- true_5 |>
  filter(medicare_field == "medicare09") |>
  group_by(medicare_status, .drop = FALSE) |>
  summarise(median_wt_medicare09 = round(median(weights), 1),
            pct_ge_16_medicare09 = round(100 * mean(weights >= 16), 1))

medicare_by_field <- by_status_10 |>
  left_join(by_status_09, by = "medicare_status") |>
  mutate(pct_of_true_pairs = round(100 * true_pairs / sum(true_pairs)))

knitr::kable(medicare_by_field,
             caption = "True pairs by Medicare agreement: five variables (names, DOB, Medicare, postcode)")

median_reissued <- medicare_by_field |> filter(medicare_status == medicare_levels[2]) |> pull(median_wt_medicare10)
median_missing  <- medicare_by_field |> filter(medicare_status == medicare_levels[4]) |> pull(median_wt_medicare10)
pct_agrees      <- medicare_by_field |> filter(medicare_status == medicare_levels[1]) |> pull(pct_of_true_pairs)
pct_reissued    <- medicare_by_field |> filter(medicare_status == medicare_levels[2]) |> pull(pct_of_true_pairs)
pct_missing     <- medicare_by_field |> filter(medicare_status == medicare_levels[4]) |> pull(pct_of_true_pairs)

## ----perch-inside, eval = FALSE-----------------------------------------------
# linked <- murmuration(
#   data_x               = cases_blocked,
#   data_y               = air_covid_blocked,  # for a vaccine question: ONE vaccine type
#   linkage_type         = "cohort2event",
#   event_type           = "vaccination",
#   event_date           = "onset_date",
#   id_var               = "id_var",
#   blocking_var         = "block1",
#   compare_vars         = c("lettername1", "lettername2", "dob", "medicare10"),
#   threshold_value      = 8,
#   perch_before_linking = TRUE,     # <-- prints the perch() table mid-linkage
#   keep_pairs_pred      = TRUE      # <-- keeps the pairs for later re-exploration
# )

## ----reexplore, eval = FALSE--------------------------------------------------
# murmuration_plot(linked, threshold = 9)                     # a different cutoff
# perch(attr(linked, "pairs_pred"), n_records_primary = nrow(cases_blocked),
#       thresholds = seq(4, 20, by = 1))                       # re-sweep, no re-fit

## ----family-flag--------------------------------------------------------------
review_4 <- review_4 |>
  mutate(
    # above the threshold, same card, but a different first name AND date of birth
    family_flag = threshold &
      (medicare10.x == medicare10.y) %in% TRUE &
      (dob.x != dob.y) %in% TRUE &
      (lettername1.x != lettername1.y) %in% TRUE,
    # ...and both records carry the IRN
    irn_known = family_flag & !is.na(medicare_irn.x) & !is.na(medicare_irn.y)
  )

n_family_flag <- sum(review_4$family_flag)
n_irn_known   <- sum(review_4$irn_known)

review_4 |>
  as.data.frame() |>
  filter(family_flag) |>
  select(lettername1.x, lettername1.y, lettername2.x, lettername2.y,
         dob.x, dob.y, medicare_irn.x, medicare_irn.y, weights) |>
  knitr::kable(digits = 1,
               caption = "Above-threshold pairs sharing a Medicare card but not a first name or DOB")

## ----family-override----------------------------------------------------------
review_plain <- review_4    # kept unchanged for comparison below

review_4 <- review_4 |>
  mutate(
    match_override = case_when(
      irn_known   ~ if_else(medicare_irn.x == medicare_irn.y, "match", "non-match"),
      family_flag ~ if_else(truth_id.x == truth_id.y, "match", "non-match"),  # the reviewer
      TRUE        ~ match_override
    ),
    override_reason = case_when(
      irn_known   ~ "Shared Medicare card: decided on the IRN",
      family_flag ~ "Shared Medicare card: reviewed for family member",
      TRUE        ~ NA_character_
    )
  )

review_4 |>
  as.data.frame() |>
  filter(!is.na(override_reason)) |>
  count(override_reason, match_override, name = "pairs") |>
  knitr::kable()

## ----best-candidate-----------------------------------------------------------
# accepted = a reviewed match, or above the threshold and not reviewed;
# rank within each case (reviewed match first, then weight) and keep the first
runner_up_rows <- review_4 |>
  as.data.frame() |>
  mutate(row = row_number()) |>
  filter(match_override %in% "match" | (is.na(match_override) & threshold)) |>
  arrange(id_var.x, desc(match_override %in% "match"), desc(weights)) |>
  filter(duplicated(id_var.x)) |>
  pull(row)

review_4 <- review_4 |>
  mutate(
    runner_up       = row_number() %in% runner_up_rows,
    match_override  = if_else(runner_up, "non-match", match_override),
    override_reason = if_else(runner_up, "Lower-ranked of several accepted candidates",
                              override_reason)
  )

c(cases_with_several_candidates = n_distinct(review_4$id_var.x[review_4$runner_up]),
  pairs_overridden = sum(!is.na(review_4$match_override)))

## ----resume-------------------------------------------------------------------
linked_plain <- murmuration(
  data_x            = cases_blocked,   # the SAME frames that built review_4
  data_y            = air_blocked,
  linkage_type      = "cohort2event",
  event_type        = "vaccination",
  event_date        = "onset_date",
  id_var            = "id_var",
  threshold_value   = 8,
  resume_from_pairs = review_plain     # no overrides: the plain threshold-8 linkage
)
linked_reviewed <- murmuration(
  data_x            = cases_blocked,
  data_y            = air_blocked,
  linkage_type      = "cohort2event",
  event_type        = "vaccination",
  event_date        = "onset_date",
  id_var            = "id_var",
  threshold_value   = 8,
  resume_from_pairs = review_4         # with the reviewer's decisions
)

# grade both against the truth
air_truth <- air_blocked |>
  select(id_var_y = id_var, truth_id_air = truth_id)

bind_rows(
  linked_plain    |> mutate(linkage = "Threshold 8, no review"),
  linked_reviewed |> mutate(linkage = "Threshold 8, reviewed")
) |>
  left_join(air_truth, by = "id_var_y") |>
  group_by(linkage) |>
  summarise(links   = sum(!is.na(truth_id_air)),
            correct = sum(truth_id_air == truth_id, na.rm = TRUE)) |>
  mutate(precision = round(100 * correct / links, 1)) |>
  knitr::kable(caption = "Precision (%) before and after review")

## ----session------------------------------------------------------------------
sessionInfo()

