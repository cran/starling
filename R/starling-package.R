#' @keywords internal
"_PACKAGE"

#' @import dplyr
#' @import tidyr
#' @importFrom magrittr %>%
#' @importFrom rlang := .data sym
#' @importFrom lubridate ymd
#' @importFrom stringr str_extract str_replace_all
#' @importFrom janitor clean_names
#' @importFrom datawizard data_rename
#' @importFrom reclin2 pair_blocking compare_pairs jaro_winkler problink_em select_threshold link
#' @importFrom stats median na.omit predict quantile reformulate weights
NULL

# Suppress R CMD check notes for variables used in dplyr operations
utils::globalVariables(c(
  ".", "dob_yr", "last_vax_type", "threshold", "admission_date",
  "id_var.x", "id_var.y", "total_admissions", "onset_adm_diff",
  "discharge_date", "icu_date", "icu_hours", "los", "hospital",
  "icd_code", "drg", "diagnosis_description", "dialysis",
  "first_admission_date", "dates", "total_valid_vax", "types",
  "vaccination_status", "last_vax_date", "onset_date", "event_date"
))
