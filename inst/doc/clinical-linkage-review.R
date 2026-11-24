## ----setup, include = FALSE---------------------------------------------------
# The worked data (the clutch datasets) live in mudnester; if it is not
# installed the vignette stops cleanly after this point. Phonetic blocking
# needs the phonics package (Suggests); those chunks are skipped without it.
has_mudnester <- requireNamespace("mudnester", quietly = TRUE) &&
  utils::packageVersion("mudnester") >= "0.9.0"   # the clutch datasets arrived in 0.9.0
has_phonics   <- requireNamespace("phonics", quietly = TRUE)
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

## ----contacts-----------------------------------------------------------------
clutch_contacts |>
  count(exposure_setting, name = "contacts") |>
  knitr::kable(caption = "Where the contacts were exposed")

## ----prep---------------------------------------------------------------------
contacts_clean <- clean_the_nest(
  clutch_contacts,
  data_type         = "cohort",
  drop_eggs         = TRUE,
  id_var            = "contact_id",
  lettername1       = "given_name",
  lettername2       = "family_name",
  dob               = "date_of_birth",
  gender            = "sex",
  postcode          = "postcode",
  medicare          = "medicare_number",
  cohort_entry_date = "exposure_date",
  keep_vars         = c("truth_id", "exposure_setting", "pep_type")
)

air_clean <- clean_the_nest(
  clutch_air,                  # every antigen: a person frame for the identity link
  data_type       = "vaccination",
  lie_nest_flat   = TRUE,
  drop_the_na_vax = FALSE,     # keep registered people with no vaccine recorded
  drop_eggs       = TRUE,
  id_var          = "air_id",
  lettername1     = "given_name",
  lettername2     = "family_name",
  dob             = "date_of_birth",
  gender          = "sex",
  postcode        = "postcode",
  medicare        = "medicare_card",
  medicare_irn    = "medicare_irn",
  vax_type        = "antigen",
  vax_date        = "date_of_service",
  keep_vars       = "truth_id"
)

# Step 1: audit (see vignette("pre-linkage-quality") for what each flag means)
audit <- preflight(contacts_clean, air_clean,
                   linkage_vars = c("lettername1", "lettername2", "dob", "medicare10"),
                   id_col1      = "id_var",
                   id_col2      = "id_var",
                   date_cols    = "dob",
                   medicare_col = "medicare10",
                   verbose      = FALSE)
knitr::kable(audit$completeness, caption = "Completeness of the linkage variables")

# Step 2: Medicare validation -- invalid numbers become NA
contacts_checked <- contacts_clean |>
  check_medicare(verbose = FALSE) |>
  mutate(medicare10 = if_else(medicare_valid %in% 0L, NA_character_, medicare10))
air_checked <- air_clean |>
  check_medicare(verbose = FALSE) |>
  mutate(medicare10 = if_else(medicare_valid %in% 0L, NA_character_, medicare10))

# Step 3: blocking on gender + birth year
contacts_blocked <- contacts_checked |>
  flock(block1_vars = c("gender", "birth_year"), birth_year_col = "dob")
air_blocked <- air_checked |>
  flock(block1_vars = c("gender", "birth_year"), birth_year_col = "dob")

## ----given-names--------------------------------------------------------------
# each contact next to the same person's AIR record (truth_id: measurement only)
contacts_clean |>
  inner_join(air_clean, by = "truth_id", suffix = c("_contact", "_air")) |>
  summarise(
    first_given_name_agrees_pct = round(100 * mean(lettername1_contact == lettername1_air,
                                                   na.rm = TRUE), 1),
    all_given_names_agree_pct   = round(100 * mean(lettername1_full_contact == lettername1_full_air,
                                                   na.rm = TRUE), 1)
  ) |>
  knitr::kable()

## ----inspect------------------------------------------------------------------
review <- murmuration(
  data_x             = contacts_blocked,  # one row per contact
  data_y             = air_blocked,       # one row per registered AIR person
  linkage_type       = "cohort2event",
  event_type         = "vaccination",
  id_var             = "id_var",
  blocking_var       = "block1",          # gender + birth year
  compare_vars       = c("lettername1", "lettername2", "dob"),  # Medicare mostly absent
  threshold_value    = 7,                 # provisional: only marks the pairs
  inspect_matches    = TRUE,              # stop after scoring
  inspect_extra_vars = c("postcode", "medicare10", "truth_id")  # context, not scored
)

n_accepted_pairs    <- sum(review$threshold)
n_accepted_contacts <- n_distinct(review$id_var.x[review$threshold])

## ----plot-review, fig.cap = "Contact-to-AIR candidate-pair weights, run directly on the inspect_matches output."----
murmuration_plot(review, max_levels = 10)

## ----perch-review-------------------------------------------------------------
sweep <- perch(review,
               n_records_primary = nrow(contacts_blocked),
               thresholds        = seq(4, 11, by = 0.5),
               report            = FALSE,
               plot              = FALSE)

sweep |>
  select(threshold, n_above, n_clerical, link_rate) |>
  knitr::kable(caption = "Accepted pairs and clerical burden by threshold")

## ----suspects-----------------------------------------------------------------
# accepted candidates per contact
accepted_counts <- review |>
  as.data.frame() |>
  filter(threshold) |>
  count(id_var.x, name = "n_accepted")

multi_ids  <- accepted_counts |> filter(n_accepted > 1) |> pull(id_var.x)
linked_ids <- accepted_counts |> pull(id_var.x)

review <- review |>
  mutate(competing = threshold & id_var.x %in% multi_ids)   # one of several accepted

review |>
  as.data.frame() |>
  filter(competing) |>
  select(id_var.x, lettername1.x, lettername1.y, lettername2.x, lettername2.y,
         dob.x, dob.y, postcode.x, postcode.y, weights) |>
  arrange(id_var.x, desc(weights)) |>
  knitr::kable(digits = 1, caption = "Contacts with more than one accepted AIR candidate")

n_jan1 <- contacts_blocked |>
  filter(format(dob, "%m-%d") == "01-01") |>
  nrow()

## ----near-miss----------------------------------------------------------------
review <- review |>
  mutate(near_miss = weights >= 6 & !threshold & !(id_var.x %in% linked_ids))

review |>
  as.data.frame() |>
  filter(near_miss) |>
  select(id_var.x, lettername1.x, lettername1.y, lettername2.x, lettername2.y,
         dob.x, dob.y, weights) |>
  knitr::kable(digits = 1, caption = "Near misses: weight 6-7, contact otherwise unlinked")

## ----override-----------------------------------------------------------------
review <- review |>
  mutate(in_review = competing | near_miss)
review_untouched <- review    # kept unchanged for comparison in Step 10

review <- review |>
  mutate(match_override = if_else(
    in_review,
    if_else(truth_id.x == truth_id.y, "match", "non-match"),  # the follow-up call
    match_override
  ))

review |>
  as.data.frame() |>
  count(match_override, name = "pairs") |>
  knitr::kable()

## ----reason-column------------------------------------------------------------
review <- review |>
  mutate(override_reason = case_when(
    in_review & threshold  ~ "Several accepted AIR candidates: identity confirmed at follow-up call",
    in_review & !threshold ~ "Near miss (weight 6-7): identity confirmed at follow-up call",
    TRUE                   ~ NA_character_
  ))

review |>
  as.data.frame() |>
  filter(!is.na(override_reason)) |>
  count(override_reason, match_override, name = "pairs") |>
  knitr::kable()

## ----finish-------------------------------------------------------------------
linked_reviewed <- murmuration(
  data_x            = contacts_blocked,     # the SAME frames that built review
  data_y            = air_blocked,
  linkage_type      = "cohort2event",
  event_type        = "vaccination",
  event_date        = "cohort_entry_date",  # the exposure date
  id_var            = "id_var",
  threshold_value   = 7,                    # used wherever match_override is NA
  vax_window        = list(days_before = 0, days_after = 0, lookback_days = Inf),
  resume_from_pairs = review                # the reviewed table
)

## ----validate-----------------------------------------------------------------
linked_naive <- murmuration(
  data_x            = contacts_blocked,
  data_y            = air_blocked,
  linkage_type      = "cohort2event",
  event_type        = "vaccination",
  event_date        = "cohort_entry_date",
  id_var            = "id_var",
  threshold_value   = 7,
  vax_window        = list(days_before = 0, days_after = 0, lookback_days = Inf),
  resume_from_pairs = review_untouched      # no overrides: the plain threshold-7 linkage
)

n_true <- contacts_blocked |>
  filter(truth_id %in% air_blocked$truth_id) |>
  nrow()

# the truth_id behind each matched AIR record -- used only to grade the links
air_truth <- air_blocked |>
  select(id_var_y = id_var, truth_id_air = truth_id)

review_summary <- bind_rows(
  linked_naive    |> mutate(linkage = "Threshold 7 only"),
  linked_reviewed |> mutate(linkage = "Threshold 7 + review")
) |>
  left_join(air_truth, by = "id_var_y") |>
  group_by(linkage) |>
  summarise(links   = sum(!is.na(truth_id_air)),
            correct = sum(truth_id_air == truth_id, na.rm = TRUE)) |>
  mutate(precision_pct = round(100 * correct / links, 1),
         recall_pct    = round(100 * correct / n_true, 1))

knitr::kable(review_summary, caption = paste("Contacts with a true AIR record:", n_true))

n_pairs_reviewed    <- sum(review$in_review)
n_contacts_reviewed <- n_distinct(review$id_var.x[review$in_review])

## ----pass2, eval = has_mudnester && has_phonics-------------------------------
contacts_phon <- contacts_checked |>
  flock(block1_vars     = c("lettername2_soundex", "gender"),
        phonetic_vars   = c("lettername1", "lettername2"),
        phonetic_method = "soundex")
air_phon <- air_checked |>
  flock(block1_vars     = c("lettername2_soundex", "gender"),
        phonetic_vars   = c("lettername1", "lettername2"),
        phonetic_method = "soundex")

linked_pass2 <- murmuration(
  data_x          = contacts_phon,
  data_y          = air_phon,
  linkage_type    = "cohort2event",
  event_type      = "vaccination",
  event_date      = "cohort_entry_date",
  id_var          = "id_var",
  blocking_var    = "block1",            # Soundex(surname) + gender
  compare_vars    = c("lettername1", "lettername2", "dob"),
  threshold_value = 6,                   # chosen from this pass's own sweep (below)
  vax_window      = list(days_before = 0, days_after = 0, lookback_days = Inf),
  keep_pairs_pred = TRUE
)

sweep_sdx <- perch(attr(linked_pass2, "pairs_pred"),
                   n_records_primary = nrow(contacts_phon),
                   thresholds        = seq(4, 11, by = 1),
                   report            = FALSE,
                   plot              = FALSE)

sweep_sdx |>
  select(threshold, n_above, n_clerical) |>
  knitr::kable(caption = "Second pass (Soundex block): accepted pairs by threshold")

## ----pass2-gain, eval = has_mudnester && has_phonics--------------------------
first_pass_ids <- linked_reviewed |>
  filter(!is.na(id_var_y)) |>
  pull(id_var)

pass2_new <- linked_pass2 |>
  filter(!(id_var %in% first_pass_ids), !is.na(id_var_y)) |>
  select(id_var, id_var_y, truth_id) |>
  left_join(air_truth, by = "id_var_y") |>
  mutate(correct = truth_id_air == truth_id)

tibble(
  measure = c("Contacts unlinked after pass one", "Linked by pass two",
              "of which correct"),
  value   = c(nrow(contacts_blocked) - length(first_pass_ids),
              nrow(pass2_new), sum(pass2_new$correct))
) |>
  knitr::kable()

## ----deterministic, eval = has_mudnester && has_phonics-----------------------
det_exact <- murmuration(
  data_x       = contacts_phon,
  data_y       = air_phon,
  linkage_type = "cohort2event",
  event_type   = "vaccination",
  method       = "deterministic",    # every compare_vars field must match exactly
  event_date   = "cohort_entry_date",
  id_var       = "id_var",
  blocking_var = "block1",           # required, though exact matching ignores it
  compare_vars = c("lettername1", "lettername2", "dob"),
  vax_window   = list(days_before = 0, days_after = 0, lookback_days = Inf)
)
det_phonetic <- murmuration(
  data_x       = contacts_phon,
  data_y       = air_phon,
  linkage_type = "cohort2event",
  event_type   = "vaccination",
  method       = "deterministic",
  event_date   = "cohort_entry_date",
  id_var       = "id_var",
  blocking_var = "block1",
  compare_vars = c("lettername1_soundex", "lettername2_soundex", "dob"),  # the codes match exactly
  vax_window   = list(days_before = 0, days_after = 0, lookback_days = Inf)
)

bind_rows(
  det_exact       |> mutate(key = "Exact: first name + last name + DOB [1]"),
  det_phonetic    |> mutate(key = "Soundex(first name) + Soundex(last name) + DOB"),
  linked_reviewed |> mutate(key = "Probabilistic, reviewed (Steps 4-10)")
) |>
  left_join(air_truth, by = "id_var_y") |>
  group_by(key) |>
  summarise(links   = sum(!is.na(truth_id_air)),
            correct = sum(truth_id_air == truth_id, na.rm = TRUE)) |>
  knitr::kable(caption = paste("Contacts with a true AIR record:", n_true))

## ----mmr-status---------------------------------------------------------------
final <- linked_reviewed
if (has_phonics) {
  final <- bind_rows(
    linked_reviewed |> filter(!is.na(id_var_y)),               # linked in pass one
    linked_pass2    |> filter(!(id_var %in% first_pass_ids))   # the rest, from pass two
  )
}

# ONE vaccine type: MMR and MMRV doses only
mmr_doses <- clutch_air |>
  filter(grepl("^Measles", antigen)) |>
  distinct(id_var_y = air_id, date_of_service)

mmr_before_exposure <- final |>
  select(id_var, id_var_y, cohort_entry_date) |>
  inner_join(mmr_doses, by = "id_var_y", relationship = "many-to-many") |>
  filter(date_of_service <= cohort_entry_date) |>
  count(id_var, name = "mmr_doses")

final <- final |>
  left_join(mmr_before_exposure, by = "id_var") |>
  mutate(mmr_doses  = coalesce(mmr_doses, 0L),
         air_status = case_when(is.na(id_var_y) ~ "Not found on the AIR",
                                mmr_doses == 0  ~ "On the AIR, no MMR before exposure",
                                mmr_doses == 1  ~ "1 MMR dose",
                                TRUE            ~ "2+ MMR doses"))

final |>
  count(air_status, pep_type, name = "contacts") |>
  knitr::kable(caption = "AIR status at exposure by post-exposure management recorded")
final |>
  count(air_status, vaccination_status, name = "contacts") |>
  knitr::kable(caption = "MMR status against murmuration()'s any-vaccine vaccination_status")

## ----mmr-truth----------------------------------------------------------------
mmr_truth <- clutch_contacts |>
  select(id_var = contact_id, truth_id, exposure_date) |>
  inner_join(clutch_air |>
               filter(grepl("^Measles", antigen)) |>
               distinct(truth_id, date_of_service),
             by = "truth_id", relationship = "many-to-many") |>
  filter(date_of_service <= exposure_date) |>
  count(id_var, name = "mmr_truth")

final |>
  left_join(mmr_truth, by = "id_var") |>
  mutate(mmr_truth        = coalesce(mmr_truth, 0L),
         air_status_truth = case_when(!(truth_id %in% clutch_air$truth_id) ~ "Not found on the AIR",
                                      mmr_truth == 0 ~ "On the AIR, no MMR before exposure",
                                      mmr_truth == 1 ~ "1 MMR dose",
                                      TRUE           ~ "2+ MMR doses")) |>
  count(linked = air_status, truth = air_status_truth, name = "contacts") |>
  knitr::kable(caption = "AIR status from the linkage against the truth")

## ----workflow-summary, eval = FALSE-------------------------------------------
# library(starling)
# library(mudnester)
# library(dplyr)
# 
# # 1-3: clean, audit, Medicare validation, blocking
# contacts <- clutch_contacts |>
#   clean_the_nest(data_type = "cohort", cohort_entry_date = "exposure_date", ...) |>
#   check_medicare() |>
#   flock(block1_vars = c("gender", "birth_year"), birth_year_col = "dob")
# 
# # every registered person for the identity link (the exception to one-vaccine-at-a-time)
# air <- clutch_air |>
#   clean_the_nest(data_type = "vaccination", lie_nest_flat = TRUE,
#                  drop_the_na_vax = FALSE, ...) |>
#   check_medicare() |>
#   flock(block1_vars = c("gender", "birth_year"), birth_year_col = "dob")
# 
# # 4. Stop after scoring, before any threshold decision is applied
# review <- murmuration(data_x = contacts, data_y = air,
#   linkage_type = "cohort2event", event_type = "vaccination", id_var = "id_var",
#   blocking_var = "block1", compare_vars = c("lettername1", "lettername2", "dob"),
#   threshold_value = 7, inspect_matches = TRUE)
# 
# # 5-6. Look at the shape of the problem, and quantify it
# murmuration_plot(review)
# perch(review, n_records_primary = nrow(contacts))
# 
# # 7-8. Inspect competing candidates and near misses; record decisions with reasons
# review <- review |>
#   mutate(match_override  = case_when(... ~ "match", ... ~ "non-match", TRUE ~ match_override),
#          override_reason = case_when(... ~ "..."))
# 
# # 9. Finish -- same event-window logic as a normal murmuration() call
# linked <- murmuration(data_x = contacts, data_y = air,
#   linkage_type = "cohort2event", event_type = "vaccination",
#   event_date = "cohort_entry_date", id_var = "id_var", threshold_value = 7,
#   vax_window = list(days_before = 0, days_after = 0, lookback_days = Inf),
#   resume_from_pairs = review)
# 
# # 11. Optional second pass (Soundex block) for contacts still unlinked
# # 12. Count MMR doses before exposure from the MEASLES-ONLY extract, joined on id_var_y
# mmr <- clutch_air |>
#   filter(grepl("^Measles", antigen)) |>
#   distinct(id_var_y = air_id, date_of_service)

## ----session------------------------------------------------------------------
sessionInfo()

