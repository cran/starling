## ----setup, include = FALSE---------------------------------------------------
# The worked data (the clutch datasets) live in mudnester, a suggested
# package: without it the vignette stops after the notice below.
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
# cat("This vignette uses the clutch datasets from the mudnester package.",
#     "Install mudnester (>= 0.9.0) to run it.")
# knitr::knit_exit()

## ----load---------------------------------------------------------------------
library(starling)
library(mudnester)   # clean_the_nest() and the clutch_* datasets
library(dplyr)
library(ggplot2)

aviary_cols <- c(satin = "#325D88", chestnut = "#86251C", wattle = "#D8A31A",
                 evergreen = "#08403B", mud = "#A67D4C", slate = "#3E3F3A")

## ----raw----------------------------------------------------------------------
clutch_cases |>
  select(notification_id, condition, first_name, surname, date_of_birth, sex,
         medicare_number) |>
  head(6) |>
  knitr::kable()

clutch_air |>
  select(air_id, given_name, family_name, date_of_birth, sex, medicare_card,
         medicare_irn, date_of_service, antigen) |>
  head(6) |>
  knitr::kable()

## ----clean--------------------------------------------------------------------
cases_raw <- clutch_cases |>
  mutate(
    onset_imputed = is.na(onset_date),                  # keep a flag
    onset_date    = coalesce(onset_date, specimen_date) # fall back to the specimen date
  )

cases_clean <- clean_the_nest(
  cases_raw,
  data_type   = "cases",
  drop_eggs   = TRUE,
  id_var      = "notification_id",
  diagnosis   = "condition",
  lettername1 = "first_name",
  lettername2 = "surname",
  dob         = "date_of_birth",
  gender      = "sex",
  postcode    = "postcode",
  medicare    = "medicare_number",
  onset_date  = "onset_date",
  keep_vars   = c("truth_id", "onset_imputed")
)

## ----clean-air----------------------------------------------------------------
air_clean <- clean_the_nest(
  clutch_air,                    # all antigens: a person frame for identifier checks
  data_type     = "vaccination",
  lie_nest_flat = TRUE,          # one row per person: vax_date_1, vax_type_1, ...
  drop_eggs     = TRUE,
  id_var        = "air_id",
  lettername1   = "given_name",
  lettername2   = "family_name",
  dob           = "date_of_birth",
  gender        = "sex",
  postcode      = "postcode",
  medicare      = "medicare_card",
  medicare_irn  = "medicare_irn",  # the IRN is stored apart from the card number
  vax_type      = "antigen",
  vax_date      = "date_of_service",
  keep_vars     = "truth_id"
)

c(cases = nrow(cases_clean), air_people = nrow(air_clean))

## ----names--------------------------------------------------------------------
air_names <- clutch_air |>
  distinct(air_id, given_name, family_name) |>
  inner_join(air_clean |> select(air_id = id_var, lettername1, lettername1_full, lettername2),
             by = "air_id")

bind_rows(
  air_names |> filter(grepl(" ", trimws(given_name)))  |> head(3),  # compound given names
  air_names |> filter(grepl(" ", trimws(family_name))) |> head(2),  # spaced family names
  air_names |> filter(grepl("[-']", family_name))      |> head(3)   # hyphens, apostrophes
) |>
  select(-air_id) |>
  knitr::kable(caption = "Raw and cleaned AIR names")

## ----names-agreement----------------------------------------------------------
# each notification next to the same person's AIR record (truth_id: measurement only)
same_person <- cases_clean |>
  inner_join(air_clean, by = "truth_id", suffix = c("_case", "_air"))

name_agreement <- same_person |>
  summarise(
    l1      = mean(lettername1_case == lettername1_air, na.rm = TRUE),
    l1_full = mean(lettername1_full_case == lettername1_full_air, na.rm = TRUE),
    l2      = mean(lettername2_case == lettername2_air, na.rm = TRUE)
  )
name_agreement <- tibble(
  field     = c("lettername1 (first given name)",
                "lettername1_full (all given names)",
                "lettername2 (family name)"),
  agree_pct = round(100 * c(name_agreement$l1, name_agreement$l1_full,
                            name_agreement$l2), 1)
)
knitr::kable(name_agreement,
             caption = "Same person, notification vs AIR: % of cleaned names agreeing exactly")

## ----preflight----------------------------------------------------------------
audit <- preflight(
  data1        = cases_clean,
  data2        = air_clean,
  linkage_vars = c("lettername1", "lettername2", "dob", "medicare10"),
  id_col1      = "id_var",
  id_col2      = "id_var",
  date_cols    = c("dob", "onset_date"),
  medicare_col = "medicare10",
  verbose      = FALSE
)

knitr::kable(audit$completeness, digits = 1,
             caption = "Completeness of the linkage variables")
knitr::kable(audit$duplicates, digits = 1,
             caption = "Identifier uniqueness")
audit$flags

## ----medicare-----------------------------------------------------------------
cases_checked <- cases_clean |>
  check_medicare(medicare_col = "medicare10", output_col = "medicare_valid",
                 verbose = FALSE)
air_checked <- air_clean |>
  check_medicare(medicare_col = "medicare10", output_col = "medicare_valid",
                 verbose = FALSE)

medicare_summary <- bind_rows(
  cases_checked |> transmute(dataset = "clutch_cases", medicare_valid),
  air_checked   |> transmute(dataset = "clutch_air",   medicare_valid)
) |>
  mutate(status = case_when(medicare_valid %in% 1L ~ "Valid",
                            medicare_valid %in% 0L ~ "Fails checksum",
                            TRUE                   ~ "Missing")) |>
  count(dataset, status) |>
  group_by(dataset) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  ungroup()

knitr::kable(medicare_summary)

## ----medicare-plot, fig.height = 3--------------------------------------------
ggplot(medicare_summary, aes(x = pct, y = dataset, fill = status)) +
  geom_col() +
  scale_fill_manual(values = c("Valid" = aviary_cols[["evergreen"]],
                               "Fails checksum" = aviary_cols[["chestnut"]],
                               "Missing" = "grey75")) +
  labs(x = "Percent of records", y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

## ----fix-medicare-------------------------------------------------------------
cases_checked <- cases_checked |>
  mutate(medicare10 = if_else(medicare_valid %in% 0L, NA_character_, medicare10))
air_checked <- air_checked |>
  mutate(medicare10 = if_else(medicare_valid %in% 0L, NA_character_, medicare10))

## ----household-cards----------------------------------------------------------
air_checked |>
  filter(!is.na(medicare10)) |>
  count(medicare10, name = "people_on_card") |>   # people per card...
  count(people_on_card, name = "n_cards") |>       # ...and cards of each size
  knitr::kable(caption = "AIR people sharing a 10-digit Medicare number")

## ----irn-availability---------------------------------------------------------
detail_levels <- c("Card number + IRN", "Card number only", "No valid number")

medicare_detail <- bind_rows(
  cases_checked |> transmute(dataset = "clutch_cases", medicare10, medicare_irn),
  air_checked   |> transmute(dataset = "clutch_air",   medicare10, medicare_irn)
) |>
  mutate(detail = case_when(!is.na(medicare10) & !is.na(medicare_irn) ~ detail_levels[1],
                            !is.na(medicare10)                        ~ detail_levels[2],
                            TRUE                                      ~ detail_levels[3]),
         detail = factor(detail, levels = rev(detail_levels))) |>
  count(dataset, detail, .drop = FALSE) |>
  group_by(dataset) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  ungroup() |>
  arrange(desc(dataset), desc(detail))

pct_cases_with_irn <- medicare_detail |>
  filter(dataset == "clutch_cases", detail == "Card number + IRN") |>
  pull(pct)

knitr::kable(medicare_detail, caption = "How much of the Medicare number each source records")

## ----irn-plot, fig.height = 2.8, fig.cap = "Medicare detail recorded, by source. Invalid numbers have already been set to NA."----
ggplot(medicare_detail, aes(x = pct, y = dataset, fill = detail)) +
  geom_col() +
  scale_fill_manual(values = c("Card number + IRN" = aviary_cols[["satin"]],
                               "Card number only"  = aviary_cols[["wattle"]],
                               "No valid number"   = "grey75"),
                    breaks = detail_levels) +
  labs(x = "Percent of records", y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

## ----irn-households-----------------------------------------------------------
cases_with_irn <- cases_checked |>
  filter(!is.na(medicare10), !is.na(medicare_irn)) |>
  select(medicare10, medicare_irn, truth_id)
air_with_irn <- air_checked |>
  filter(!is.na(medicare10), !is.na(medicare_irn)) |>
  select(medicare10, medicare_irn, truth_id)

cases_with_irn |>
  inner_join(air_with_irn, by = "medicare10", suffix = c("_case", "_air"),
             relationship = "many-to-many") |>
  mutate(irn    = if_else(medicare_irn_case == medicare_irn_air, "IRN agrees", "IRN differs"),
         person = if_else(truth_id_case == truth_id_air, "Same person", "Different person")) |>
  count(irn, person, name = "pairs") |>
  knitr::kable(caption = "Notification-AIR pairs on the same card, IRN recorded on both")

## ----medicare-digits----------------------------------------------------------
medicare_agreement <- same_person |>
  filter(!is.na(medicare10_case), !is.na(medicare10_air)) |>
  summarise(m08 = mean(medicare08_case == medicare08_air),
            m09 = mean(medicare09_case == medicare09_air),
            m10 = mean(medicare10_case == medicare10_air))
medicare_agreement <- tibble(
  field = c("medicare08 (account)",
            "medicare09 (account + check digit)",
            "medicare10 (+ card issue digit)"),
  agree_pct = round(100 * c(medicare_agreement$m08, medicare_agreement$m09,
                            medicare_agreement$m10), 1)
)
knitr::kable(medicare_agreement,
             caption = "Same person, both sources recorded a number: % agreeing")

## ----sex-coding---------------------------------------------------------------
clutch_hospital |> count(sex)
clutch_air      |> count(sex)

air_blocks <- clutch_air |>
  distinct(block = paste(sex, format(date_of_birth, "%Y")))

clutch_hospital |>
  distinct(block = paste(sex, format(birth_date, "%Y"))) |>
  semi_join(air_blocks, by = "block") |>
  nrow()   # zero: no hospital record would ever meet an AIR record

## ----sex-fix------------------------------------------------------------------
hosp_raw <- clutch_hospital |>
  mutate(sex = case_when(sex == "Male"   ~ "M",
                         sex == "Female" ~ "F"))

hosp_raw |>
  distinct(block = paste(sex, format(birth_date, "%Y"))) |>
  semi_join(air_blocks, by = "block") |>
  nrow()

## ----flock--------------------------------------------------------------------
cases_blocked <- cases_checked |>
  flock(block1_vars    = "gender",
        block2_vars    = c("gender", "birth_year"),
        block3_vars    = "postcode",
        birth_year_col = "dob")
air_blocked <- air_checked |>
  flock(block1_vars    = "gender",
        block2_vars    = c("gender", "birth_year"),
        block3_vars    = "postcode",
        birth_year_col = "dob")

## ----block-evaluation---------------------------------------------------------
# candidate pairs = sum over shared keys of (records with the key) x (records with the key)
pairs_block1 <- cases_blocked |>
  filter(!is.na(block1)) |>
  count(block1, name = "n_x") |>
  inner_join(air_blocked |> count(block1, name = "n_y"), by = "block1") |>
  summarise(pairs = sum(as.numeric(n_x) * n_y)) |>
  pull(pairs)

pairs_block2 <- cases_blocked |>
  filter(!is.na(block2)) |>
  count(block2, name = "n_x") |>
  inner_join(air_blocked |> count(block2, name = "n_y"), by = "block2") |>
  summarise(pairs = sum(as.numeric(n_x) * n_y)) |>
  pull(pairs)

pairs_block3 <- cases_blocked |>
  filter(!is.na(block3)) |>
  count(block3, name = "n_x") |>
  inner_join(air_blocked |> count(block3, name = "n_y"), by = "block3") |>
  summarise(pairs = sum(as.numeric(n_x) * n_y)) |>
  pull(pairs)

# blocking recall: the share of true pairs that land in the same block
recall <- cases_blocked |>
  select(truth_id, x1 = block1, x2 = block2, x3 = block3) |>
  inner_join(air_blocked |> select(truth_id, y1 = block1, y2 = block2, y3 = block3),
             by = "truth_id") |>
  summarise(block1 = mean(x1 == y1, na.rm = TRUE),
            block2 = mean(x2 == y2, na.rm = TRUE),
            block3 = mean((x3 == y3) %in% TRUE))   # a missing postcode never matches

blocking_table <- tibble(
  block           = c("block1: gender", "block2: gender + birth year", "block3: postcode"),
  candidate_pairs = c(pairs_block1, pairs_block2, pairs_block3),
  blocking_recall = round(100 * c(recall$block1, recall$block2, recall$block3), 1)
)
knitr::kable(blocking_table, format.args = list(big.mark = ","))

## ----phonetic, eval = has_mudnester && has_phonics----------------------------
# block1: Soundex(surname) + gender; block2: gender + birth year, for comparison
cases_phon <- cases_checked |>
  flock(block1_vars     = c("lettername2_soundex", "gender"),
        block2_vars     = c("gender", "birth_year"),
        phonetic_vars   = "lettername2",
        phonetic_method = "soundex",
        birth_year_col  = "dob")
air_phon <- air_checked |>
  flock(block1_vars     = c("lettername2_soundex", "gender"),
        block2_vars     = c("gender", "birth_year"),
        phonetic_vars   = "lettername2",
        phonetic_method = "soundex",
        birth_year_col  = "dob")

cases_phon |>
  distinct(lettername2, lettername2_soundex) |>
  head(8) |>
  knitr::kable()

pairs_soundex <- cases_phon |>
  filter(!is.na(block1)) |>
  count(block1, name = "n_x") |>
  inner_join(air_phon |> count(block1, name = "n_y"), by = "block1") |>
  summarise(pairs = sum(as.numeric(n_x) * n_y)) |>
  pull(pairs)

phonetic_recall <- cases_phon |>
  select(truth_id, sdx_x = block1, year_x = block2) |>
  inner_join(air_phon |> select(truth_id, sdx_y = block1, year_y = block2),
             by = "truth_id") |>
  mutate(same_soundex   = sdx_x == sdx_y,
         same_birthyear = year_x == year_y) |>
  summarise(year   = mean(same_birthyear, na.rm = TRUE),
            sdx    = mean(same_soundex, na.rm = TRUE),
            either = mean(same_birthyear | same_soundex, na.rm = TRUE))

tibble(
  strategy        = c("gender + birth year", "Soundex(surname) + gender",
                      "either (two passes, union)"),
  candidate_pairs = c(pairs_block2, pairs_soundex, pairs_block2 + pairs_soundex),
  blocking_recall = round(100 * c(phonetic_recall$year, phonetic_recall$sdx,
                                  phonetic_recall$either), 1)
) |>
  knitr::kable(format.args = list(big.mark = ","))

## ----multi-pass, eval = FALSE-------------------------------------------------
# # For a vaccine question, the AIR side holds ONE vaccine type (see the callout above)
# air_covid_phon <- clutch_air |>
#   filter(antigen == "COVID-19") |>
#   clean_the_nest(data_type = "vaccination", lie_nest_flat = TRUE, ...) |>
#   check_medicare() |>
#   flock(block1_vars     = c("lettername2_soundex", "gender"),
#         block2_vars     = c("gender", "birth_year"),
#         phonetic_vars   = "lettername2",
#         phonetic_method = "soundex",
#         birth_year_col  = "dob")
# 
# # pass 1: block on gender + birth year
# linked_year <- murmuration(
#   data_x          = cases_phon,
#   data_y          = air_covid_phon,
#   linkage_type    = "cohort2event",
#   event_type      = "vaccination",
#   event_date      = "onset_date",
#   id_var          = "id_var",
#   blocking_var    = "block2",             # gender + birth year
#   compare_vars    = c("lettername1", "lettername2", "dob", "medicare10"),
#   threshold_value = 8
# )
# 
# # pass 2: block on Soundex(surname) + gender
# linked_sdx <- murmuration(
#   data_x          = cases_phon,
#   data_y          = air_covid_phon,
#   linkage_type    = "cohort2event",
#   event_type      = "vaccination",
#   event_date      = "onset_date",
#   id_var          = "id_var",
#   blocking_var    = "block1",             # Soundex(surname) + gender
#   compare_vars    = c("lettername1", "lettername2", "dob", "medicare10"),
#   threshold_value = 8
# )
# 
# # keep the best-scoring link per case across both passes
# linked_all <- bind_rows(linked_year, linked_sdx) |>
#   arrange(desc(threshold), desc(weights)) |>
#   distinct(id_var, .keep_all = TRUE)

## ----summary, eval = FALSE----------------------------------------------------
# library(starling)
# library(mudnester)
# library(dplyr)
# 
# # 1. Clean both sources the same way (impute missing onset dates first)
# cases <- clutch_cases |>
#   mutate(onset_date = coalesce(onset_date, specimen_date)) |>
#   clean_the_nest(data_type = "cases", id_var = "notification_id", ...)
# 
# # ... and the AIR for ONE vaccine type
# air <- clutch_air |>
#   filter(antigen == "COVID-19") |>
#   clean_the_nest(data_type = "vaccination", lie_nest_flat = TRUE,
#                  medicare = "medicare_card", medicare_irn = "medicare_irn", ...)
# 
# # 2. Audit
# preflight(cases, air,
#           linkage_vars = c("lettername1", "lettername2", "dob", "medicare10"),
#           medicare_col = "medicare10")
# 
# # 3. Medicare: invalid numbers become NA; 4. Blocks (harmonise coding first)
# cases <- cases |>
#   check_medicare() |>
#   mutate(medicare10 = if_else(medicare_valid %in% 0L, NA_character_, medicare10)) |>
#   flock(block1_vars = c("gender", "birth_year"), birth_year_col = "dob")
# air <- air |>
#   check_medicare() |>
#   mutate(medicare10 = if_else(medicare_valid %in% 0L, NA_character_, medicare10)) |>
#   flock(block1_vars = c("gender", "birth_year"), birth_year_col = "dob")
# 
# # 5. Link -- see vignette("linked-cohort")

## ----session------------------------------------------------------------------
sessionInfo()

