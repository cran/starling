#' @name murmuration
#' @title Links case, hospital or vaccination datasets
#' @description Machine learning data linkage. The murmuration command will link diagnostic registry data (cases or linelist) to hospitalization and immunization records (e.g. Australian Immunization Register).
#'
#' @details A murmuration is a shape-shifting flock of thousands of starlings all flying in synch with each other.
#' Murmuration means that each bird must be linked (through observation of their movements) to approximately seen other birds to achieve the beautiful sky art that moves through the sky.
#' Make sure that you do not have the same variables (other than linkage variables e.g. letternames, DOB, gender) in both datasets.
#' Always make sure your date columns are properly formatted using as_date, or as.Date.
#' For example, if both datasets have date of death, choose the dataset with the highest confidence, and drop out the date of death from the other dataset.
#' If the dataset is being linked to a hospitalization dataset, the difference in time between onset_date and admission_date will be used to identify related hospitalizations.
#' The user can filter out unrelated hospitalizations using diagnostic-related groups or ICD-10 codes separately, prior to linkage.
#' Classic workflow would be:
#' \enumerate{
#'  \item \code{\link{clean_the_nest}} to clean and prep data for linkage. Pay close attention to your linkage variables (letternames, date of birth, medicare number, gender and/or postcode), and ensure all dates are formatted as dates.
#'  \item \code{\link{murmuration}} with linkage_type="v2c" to link cases to vaccination data.
#'  \item \code{\link{murmuration}} with linkage_type="v2h" to link a v2c dataset to hospitalization data. Or skip linking to case data, and just build a v2h dataset for test-negative case-control studies.
#'  \item \code{\link{murmuration}} with linkage_type="v2e" to link event linelists (flight manifests, outbreak investigations) to vaccination data.
#'  \item \code{\link{preening}} to prettify the dataframe prepping it for exploration, analysis and presentation. Great to use with \code{gtsummary::tbl_summary()}.
#' }
#' @note Ensure there are no missing vaccination dates in vaccination dataset prior to murmuration. Murmuration requires complete vaccination data (equal date and type columns per observation) to achieve correct matching of vaccination columns.
#' If there are too few variables to match on, then matching will not work well. For example, if you have first name, last name and date of birth, and a very large dataset (Immunization Register), then the scoring will not differentiate true from false matches. Consider deterministic linkage when there is a paucity of information to use to derive linkage scores.
#'
#' @param df1 This is a dataframe object, cleaned using clean_build_nest, and would often represent the base, or "x" dataset (when doing left joins).
#' Typically this would be a dataset of cases, have enough data to create linkages, and have onset dates.
#' @param df2 This is a dataframe object, cleaned using clean_build_nest, and would often represent the admissions or vaccination dataset ("y" dataset when doing left joins).
#' Typically this would have enough data to create linkages, and include either admission data or vaccination event data (e.g. Australian Immunization Register).
#' @param linkage_type Parameter name. Either "c2h", for linkage of cases to hospital admissions data (default).
#' "v2c" for linkage of cases to vaccination datasets.
#' "v2h" for linkage of hospitalizations to vaccination history (e.g. building a dataset for test-negative case-control studies). Use "v2h" if you want to link a "v2c" dataset to a hospitalization dataset.
#' "v2e" for linkage of event participants (flight manifest, outbreak linelist) to vaccination history to determine vaccination status at time of event.
#' If using linking to a vaccination dataset, must use single row per person dataset. If you have multiple vaccines per person, run it through the clean_the_nest command with "lie_nest_flat" option set to TRUE.
#' @param onset_date Variable name for onset date (used in c2h and v2c linkage types). Should be present in df1.
#' @param event_date A date object (e.g., ymd('2024-12-15')) representing when the event occurred. Required for v2e linkage type. All valid vaccinations must occur before this date.
#' @param id_var Variable name (e.g. "id")This is critical for data-linkage and the base dataset is the dataset you would left join onto (e.g. the "x" dataset). Cannot have missing data, or the observation will be lost in the linking process.
#' @param blocking_var Variable name (e.g. "block1"). Choice of blocking variable. You can create your own. Up to three blocking vars are created in the past
#' @param compare_vars Vector of variables. Used to compare variables between each dataset and calculate the string score differences. Typically names, dates of births and medicare/social security numbers.
#' @param threshold_value Numeric (e.g. "12"), default is 12. This represents the threshold above which you decide that the linkage is true or false.
#' The higher the number, the higher the specificity of your linkages (compare_vars match more exactly).
#' The lower the threshold, the more sensitive you are to selecting matches, at the expense of specificity. Default is 12, and arbitrarily chosen.
#' @param days_allowed_before_event Numeric (e.g. "7"). How much time you choose to allow prior to the onset_date of a disease-related admission for a c2h dataset (see linkage type).
#' For c2h linkages, this represents the lower limit of the window for disease related admissions.
#' For v2h datasets this represents the minimum time between latest vaccination date and admission date to be considered a valid vaccination dose.
#' For a v2c datasets this represents the minimum time between latest vaccination date and onset_date to be considered a valid vaccination dose.
#' For v2e datasets this represents the minimum time between latest vaccination date and event_date to be considered a valid vaccination dose.
#' For example, if you choose seven days, then you are allowing for an admission to occur up to seven days prior to the diagnosis, which means the disease was diagnosed while an inpatient.
#' @param days_allowed_after_event Numeric (e.g. "30"). How much time you choose to allow after the onset of a disease related admission. Upper limit of window for disease related admissions.
#' For example, if you choose 30 days, then you are allowing for a disease-related admission to occur up to 30 days after the diagnosis, which means the disease was diagnosed very close to or prior to the admission.
#' @param one_row_per_person Logical (TRUE or FALSE) with the default being TRUE. It will take multiple admissions per person, and create a series of variables prefixed with "first_", such as "first_admission_date",
#' and put into a single row all admission events, and create a series of variables suffixed with "s", such as "admission_dates". Will work with single admissions per person.
#' @param clean_eggs Logical (TRUE or FALSE) with the default being TRUE. Drops all the .y variables that are duplicates of the second dataset (df2), and keeps the variables and removes the .x from df1. If you leave this on, many,
#' if not most variables will have ".x" or ".y" attached to them (e.g. gender) and thus keep this as TRUE for default, and FALSE if you want to check the linkages are true and working.
#' @param days_between_onset_death Numeric (e.g. "30"). If you have put a date of death into the clean_build_nest command (which will rename it to dod), then the command will find disease related dates of death.
#' This is chosen number of days between onset and death for a disease-related death. Often this may be 30 days for SARS-CoV-2 or can be much longer for HIV. If you don't want an upper limit, use "9999".
#' @param last_follow_up represents a date (input as ymd(2024-11-22)) that represents last follow-up. This could be the latest admission date of a dataset. Used for calculating survival time.
#' @return A linked dataset with some new variables.
#' @export
#' @examples
#' \donttest{
#' # Example 1: Link cases to vaccination history
#' # First, clean the datasets to standardize column names
#' dx_clean <- clean_the_nest(dx_data,
#'   data_type = "cases",
#'   id_var = "identity",
#'   lettername1 = "first_name",
#'   lettername2 = "surname",
#'   dob = "date_of_birth",
#'   gender = "gender",
#'   postcode = "postcode",
#'   medicare = "medicare_no",
#'   diagnosis = "disease_name")
#'
#' vax_clean <- clean_the_nest(vax_data,
#'   data_type = "vaccination",
#'   id_var = "patient_id",
#'   lettername1 = "firstname",
#'   lettername2 = "last_name",
#'   dob = "birth_date",
#'   gender = "gender",
#'   postcode = "postcode",
#'   medicare = "medicare_number",
#'   vax_type = "vaccine_delivered",
#'   vax_date = "service_date")
#'
#' # Now link cases to vaccination history
#' df1 <- murmuration(dx_clean, vax_clean,
#'   linkage_type = "v2c",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob"),
#'   clean_eggs = FALSE)
#'
#' # Example 2: Link hospitalization data to vaccination history
#' hosp_clean <- clean_the_nest(hosp_data,
#'   data_type = "hospital",
#'   id_var = "patient_id",
#'   lettername1 = "firstname",
#'   lettername2 = "last_name",
#'   dob = "birth_date",
#'   gender = "sex",
#'   postcode = "zip_codes",
#'   medicare = "medicare_number",
#'   admission_date = "date_of_admission",
#'   discharge_date = "date_of_discharge")
#'
#' df2 <- murmuration(hosp_clean, vax_clean,
#'   linkage_type = "v2c",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "medicare10", "dob"),
#'   clean_eggs = FALSE,
#'   one_row_per_person = TRUE)
#'
#' # Example 3: Link flight manifest to vaccination history
#' manifest_clean <- clean_the_nest(manifest_data,
#'   data_type = "cases",
#'   id_var = "passenger_id",
#'   lettername1 = "first_name",
#'   lettername2 = "surname",
#'   dob = "date_of_birth",
#'   gender = "gender")
#'
#' df_flight <- murmuration(manifest_clean, vax_clean,
#'   linkage_type = "v2e",
#'   event_date = as.Date("2024-03-15"),
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob"),
#'   days_allowed_before_event = 14,
#'   clean_eggs = FALSE)
#'
#' # Example 4: Link outbreak linelist to vaccination history
#' linelist_clean <- clean_the_nest(linelist_data,
#'   data_type = "cases",
#'   id_var = "case_id",
#'   lettername1 = "first_name",
#'   lettername2 = "surname",
#'   dob = "date_of_birth",
#'   gender = "gender",
#'   postcode = "postcode",
#'   medicare = "medicare_no",
#'   onset_date = "onset_date")
#'
#' df_outbreak <- murmuration(linelist_clean, vax_clean,
#'   linkage_type = "v2e",
#'   event_date = as.Date("2024-06-01"),
#'   blocking_var = "postcode",
#'   compare_vars = c("lettername1", "lettername2", "dob", "medicare10"),
#'   days_allowed_before_event = 7,
#'   clean_eggs = FALSE)
#' }
murmuration <- function(df1,
                        df2,
                        linkage_type = "c2h",
                        onset_date=NULL,
                        event_date=NULL,
                        id_var = id_var,
                        blocking_var,
                        compare_vars,
                        threshold_value=12,
                        days_allowed_before_event=7,
                        days_allowed_after_event=14,
                        one_row_per_person=TRUE,
                        clean_eggs=TRUE,
                        days_between_onset_death=30,
                        last_follow_up=NULL) {


  if ("weights" %in% colnames(df1)) {
    df1 <- df1 %>% select(-weights)
  }
  if ("weights" %in% colnames(df2)) {
    df2 <- df2 %>% select(-weights)
  }
  if ("threshold" %in% colnames(df1)) {
    df1 <- df1 %>% select(-threshold)
  }
  if ("threshold" %in% colnames(df2)) {
    df2 <- df2 %>% select(-threshold)
  }

  pairs <- pair_blocking(df1, df2, blocking_var)

  if (identical(pairs$.x, pairs$.y)==TRUE) {
    stop("Blocking variable too specific, matching 1:1. Choose broader/looser blocking variable like postcode, or gender, rather than postcode + year of birth. <(*)")
  }

  compare_pairs(pairs, on = compare_vars, default_comparator = jaro_winkler(0.9), inplace = TRUE)

  m <- problink_em(reformulate(compare_vars), data = pairs)

  pairs_pred <- predict(m, pairs = pairs, add = TRUE)

  pairs_thresh <- select_threshold(pairs_pred, "threshold", score = "weights", threshold = threshold_value)

  df <- link(pairs_thresh, selection = TRUE, keep_from_pairs = c("weights", "threshold"), all = TRUE)

  df <- df %>% filter(!is.na(id_var.x)) # This removes persons that are created that do not have an admission ("x" data entry) but have a Y data entry

  #___ Case to Hospitalization (Refactored to match v2h) ___

  if (linkage_type == "c2h") {

    if (!"onset_date" %in% colnames(df1)) {
      stop("You need an onset_date var, and your df1 does not have it. Run your data through clean_the_nest first <(*)" )
    }
    if (!"admission_date" %in% colnames(df2)) {
      stop("You need an admission_date var, and your df2 does not have it. Run your data through clean_the_nest first <(*)")
    }

    # Step 1: Invalidate admission data based on linkage quality and timing
    df <- df %>%
      mutate(
        onset_adm_diff = admission_date - onset_date,
        admission_date = case_when(
          threshold == FALSE ~ NA_Date_,
          admission_date < onset_date - days_allowed_before_event ~ NA_Date_,
          admission_date > onset_date + days_allowed_after_event ~ NA_Date_,
          TRUE ~ admission_date
        ),
        # Invalidate all admission-related variables when admission_date is NA
        across(c(ends_with(".y"), starts_with("admission"),
                 any_of(c("discharge_date", "onset_adm_diff", "los", "icu",
                          "icu_date", "icu_hours", "hospital", "icd_code",
                          "diagnosis_description", "drg", "dod"))),
               ~ case_when(is.na(admission_date) ~ NA, TRUE ~ .))
      )

    # Step 2: Select the single best record for each person
    df <- df %>%
      group_by(id_var.x) %>%
      arrange(desc(threshold), admission_date, .by_group = TRUE) %>%
      slice(1) %>%
      ungroup() %>%
      rename(id_var_df2 = id_var.y)

    # Step 3: Calculate admission counts and outcomes
    df <- df %>%
      mutate(
        total_admissions = if_else(!is.na(admission_date), 1L, 0L),
        adm_no = total_admissions
      )

    # Step 4: Handle multiple admissions per person (if one_row_per_person == FALSE)
    # NOTE: Since we're doing slice(1) above, we only ever have 1 admission per person now
    # If you want to preserve multiple admissions, you'll need to adjust Step 2

    if (one_row_per_person == TRUE) {

      if ("onset_adm_diff" %in% colnames(df)) {
        df <- df %>% mutate(
          first_onset_adm_diff = onset_adm_diff,
          all_onset_adm_diffs = as.character(onset_adm_diff)
        )
      }

      if ("admission_date" %in% colnames(df)) {
        df <- df %>% mutate(
          first_admission_date = admission_date,
          last_admission_date = admission_date,
          all_admission_dates = as.character(admission_date)
        )
      }

      if ("discharge_date" %in% colnames(df)) {
        df <- df %>% mutate(
          first_discharge_date = discharge_date,
          last_discharge_date = discharge_date,
          all_discharge_dates = as.character(discharge_date)
        )
      }

      if ("icu_date" %in% colnames(df)) {
        df <- df %>% mutate(first_icu_date = icu_date)
      }

      if ("icu_hours" %in% colnames(df)) {
        df <- df %>% mutate(total_icu_hours = icu_hours)
      }

      if ("los" %in% colnames(df)) {
        df <- df %>% mutate(
          total_los = as.numeric(los),
          all_los = as.character(los)
        )
      }

      if ("hospital" %in% colnames(df)) {
        df <- df %>% mutate(
          first_hospital = hospital,
          all_hospitals = as.character(hospital)
        )
      }

      if ("icd_code" %in% colnames(df)) {
        df <- df %>% mutate(
          first_icd_code = icd_code,
          all_icd_codes = as.character(icd_code)
        )
      }

      if ("drg" %in% colnames(df)) {
        df <- df %>% mutate(
          first_drg = drg,
          all_drgs = as.character(drg)
        )
      }

      if ("diagnosis_description" %in% colnames(df)) {
        df <- df %>% mutate(
          first_diagnosis_description = diagnosis_description,
          all_diag_desc = as.character(diagnosis_description)
        )
      }

      if ("dialysis_outcome" %in% colnames(df)) {
        df <- df %>% mutate(dialysis_outcome = dialysis)
      }

      if (!is.null(last_follow_up)) {
        df <- df %>% mutate(
          survtime = case_when(
            !is.na(dod) ~ as.numeric(dod - onset_date),
            TRUE ~ as.numeric(last_follow_up - onset_date)
          )
        )
      }

      # Create outcome variables
      df <- df %>%
        mutate(
          admission_outcome = factor(
            if_else(!is.na(first_admission_date), 1L, 0L),
            levels = c(0, 1),
            labels = c("No Admission", "Admission"),
            ordered = TRUE
          )
        )

      if ("icu_date" %in% colnames(df)) {
        df <- df %>% mutate(
          icu_outcome = factor(
            if_else(!is.na(icu_date), 1L, 0L),
            levels = c(0, 1),
            labels = c("No ICU Admission", "ICU Admission"),
            ordered = TRUE
          )
        )
      } else if ("total_icu_hours" %in% colnames(df)) {
        df <- df %>% mutate(
          icu_outcome = factor(
            case_when(
              total_icu_hours > 0 ~ 1L,
              TRUE ~ 0L
            ),
            levels = c(0, 1),
            labels = c("No ICU Admission", "ICU Admission"),
            ordered = TRUE
          )
        )
      }

      if ("dod" %in% colnames(df)) {
        df <- df %>% mutate(
          death_outcome = factor(
            case_when(
              !is.na(dod) & dod <= onset_date + days_between_onset_death ~ 1L,
              TRUE ~ 0L
            ),
            levels = c(0, 1),
            labels = c("Alive", "Died"),
            ordered = TRUE
          )
        )
      }

      # Remove intermediate columns
      df <- df %>% select(-any_of(c("adm_no", "admission_date", "onset_adm_diff",
                                    "los", "hospital", "icd_code",
                                    "diagnosis_description", "drg")))
    }  # End one_row_per_person == TRUE
  }  # End c2h


  #___ Vaccine to Case (Refactored to match v2h) ___

  if (linkage_type == "v2c") {

    # Step 1: Invalidate vaccination data based on linkage quality and timing
    df <- df %>%
      mutate(
        across(starts_with('vax_date_'), ~ case_when(
          threshold == FALSE ~ NA_Date_,
          . + days_allowed_before_event > onset_date ~ NA_Date_,
          TRUE ~ .
        )),
        across(c(starts_with('vax_type_'), ends_with(".y")), ~ case_when(
          threshold == FALSE ~ NA,
          TRUE ~ .
        ))
      )

    # Step 2: Select the single best record for each person
    df <- df %>%
      group_by(id_var.x) %>%
      arrange(desc(threshold), .by_group = TRUE) %>%
      slice(1) %>%
      ungroup()

    # Step 3: Calculate first and last vaccination dates
    date_cols <- grep("^vax_date_\\d+$", names(df), value = TRUE)
    type_cols <- grep("^vax_type_\\d+$", names(df), value = TRUE)

    # ____ Validate vaccination doses by numeric pairing
    for (dc in date_cols) {
      # Extract the numeric suffix from the date column name
      x <- str_extract(dc, "\\d+$")
      tc <- paste0("vax_type_", x)

      # Only update if the matching type column exists
      if (tc %in% type_cols) {
        df <- df %>%
          mutate(
            !!tc := if_else(is.na(.data[[dc]]), NA_character_, .data[[tc]])
          )
      }}

    # ____ Compute first / last vax & type based on valid dates
    df <- df %>%
      rowwise() %>% mutate(
        dates = list(c_across(all_of(date_cols))),
        types = list(c_across(all_of(type_cols))),
        total_valid_vax = sum(!is.na(dates)),

        first_vax_date = if (total_valid_vax > 0) min(dates, na.rm = TRUE) else as.Date(NA),
        first_vax_type = if (total_valid_vax > 0) types[which.min(dates)] else NA_character_,

        last_vax_date = if (total_valid_vax > 0) max(dates, na.rm = TRUE) else as.Date(NA),
        last_vax_type = if (total_valid_vax > 0) types[which.max(dates)] else NA_character_)

  }


  #___ Vaccine to Hospitalization (Original - already good) ___

  if (linkage_type == "v2h") {

    # Step 1: Invalidate vaccination data based on linkage quality and timing
    df <- df %>%
      mutate(
        across(starts_with('vax_date_'), ~ case_when(
          threshold == FALSE ~ NA_Date_,
          . + days_allowed_before_event > admission_date ~ NA_Date_,
          TRUE ~ .
        )),
        across(c(starts_with('vax_type_'), ends_with(".y")), ~ case_when(
          threshold == FALSE ~ NA,
          TRUE ~ .
        ))
      )

    # Step 2: Select the single best record for each person
    df <- df %>%
      group_by(id_var.x) %>%
      arrange(desc(threshold), .by_group = TRUE) %>%
      slice(1) %>%
      ungroup()

    # Step 3: Calculate vaccination status and outcomes                               _______________________________________________________________________________________
    # Section updated 06/11/2025 - vax dates that were filtered out by days_allowed_before_event
    date_cols <- grep("^vax_date_\\d+$", names(df), value = TRUE)                     # were not updated in the corresponding vax type columns, leading to a mismatch in the
    type_cols <- grep("^vax_type_\\d+$", names(df), value = TRUE)                     # number of vax dates and types recorded, and a misrepresentation of last_vax_type. This
    # code was applied to each of the other dataset matching methods.
    # ____ Validate vaccination doses by numeric pairing - necessary to ensure
    # the appropriate corresponding vax_type_ column is filtered out.
    for (dc in date_cols) {
      # Extract the numeric suffix from the date column name
      x <- str_extract(dc, "\\d+$")
      # Match it to the vax_type_ column with the same suffix
      tc <- paste0("vax_type_", x)

      # Only update if the matching type column exists
      if (tc %in% type_cols) {
        df <- df %>%
          mutate(
            # Replace vax_type_ column with NA if the corresponding vax_date column is NA
            !!tc := if_else(is.na(.data[[dc]]), NA_character_, .data[[tc]])
          )
      }}

    # ____ Compute first / last vax & type based on valid dates
    df <- df %>%
      rowwise() %>% mutate(
        dates = list(c_across(all_of(date_cols))),
        types = list(c_across(all_of(type_cols))),
        total_valid_vax = sum(!is.na(dates)),

        first_vax_date = if (total_valid_vax > 0) min(dates, na.rm = TRUE) else as.Date(NA),
        first_vax_type = if (total_valid_vax > 0) types[which.min(dates)] else NA_character_,

        last_vax_date = if (total_valid_vax > 0) max(dates, na.rm = TRUE) else as.Date(NA),
        last_vax_type = if (total_valid_vax > 0) types[which.max(dates)] else NA_character_,

        # Vaccination status
        vaccination_status = if (total_valid_vax > 0) "Vaccinated" else "Unvaccinated",
        vaccination_status_num = if_else(vaccination_status == "Unvaccinated", 0, 1),

        tsv = if (total_valid_vax > 0) as.numeric(admission_date - last_vax_date) else 0
      ) %>%
      select(-dates, -types) %>% # clean up temporary columns
      ungroup()
  }
  #                                                                                     ______________________________________________________________________________________

  #___ Vaccine to Event (New mode for flight manifests, outbreak linelists, etc.) ___

  if (linkage_type == "v2e") {

    if (is.null(event_date)) {
      stop("You need an event_date parameter for v2e linkage type. This should be a date (e.g., ymd('2024-12-15')) representing when the event occurred. <(*)")
    }

    # Step 1: Invalidate vaccination data based on linkage quality and timing
    # Only vaccines delivered before the event are valid
    df <- df %>%
      mutate(
        across(starts_with('vax_date_'), ~ case_when(
          threshold == FALSE ~ NA_Date_,
          . + days_allowed_before_event > event_date ~ NA_Date_,
          TRUE ~ .
        )),
        across(c(starts_with('vax_type_'), ends_with(".y")), ~ case_when(
          threshold == FALSE ~ NA,
          TRUE ~ .
        ))
      )

    # Step 2: Select the single best record for each person
    df <- df %>%
      group_by(id_var.x) %>%
      arrange(desc(threshold), .by_group = TRUE) %>%
      slice(1) %>%
      ungroup()

    # Step 3: Calculate vaccination status and outcomes at time of event
    date_cols <- grep("^vax_date_\\d+$", names(df), value = TRUE)
    type_cols <- grep("^vax_type_\\d+$", names(df), value = TRUE)

    # ____ Validate vaccination doses by numeric pairing
    for (dc in date_cols) {
      # Extract the numeric suffix from the date column name
      x <- str_extract(dc, "\\d+$")
      tc <- paste0("vax_type_", x)

      # Only update if the matching type column exists
      if (tc %in% type_cols) {
        df <- df %>%
          mutate(
            !!tc := if_else(is.na(.data[[dc]]), NA_character_, .data[[tc]])
          )
      }}

    # ____ Compute first / last vax & type based on valid dates
    df <- df %>%
      rowwise() %>% mutate(
        dates = list(c_across(all_of(date_cols))),
        types = list(c_across(all_of(type_cols))),
        total_valid_vax = sum(!is.na(dates)),

        first_vax_date = if (total_valid_vax > 0) min(dates, na.rm = TRUE) else as.Date(NA),
        first_vax_type = if (total_valid_vax > 0) types[which.min(dates)] else NA_character_,

        last_vax_date = if (total_valid_vax > 0) max(dates, na.rm = TRUE) else as.Date(NA),
        last_vax_type = if (total_valid_vax > 0) types[which.max(dates)] else NA_character_,

        # Vaccination status
        vaccination_status = if (total_valid_vax > 0) "Vaccinated" else "Unvaccinated",
        vaccination_status_num = if_else(vaccination_status == "Unvaccinated", 0, 1),

        tsv = if (total_valid_vax > 0) as.numeric(event_date - last_vax_date) else 0,

        event_date = event_date     # Add event_date to output for reference
      ) %>%
      select(-dates, -types) %>% # clean up temporary columns
      ungroup()
  }


  # Clean up column names
  if (clean_eggs == TRUE) {
    df <- df %>% select(-ends_with(".y"))
    df <- df %>% rename_with(~str_remove(., '[.]x$'))
  }

  # Reorder columns
  df <- df %>%
    select(
      any_of(c("lettername1_lettername2_dob", "lettername1_lettername2_dob.x",
               "onset_date", "admission_outcome", "first_admission_date",
               "first_onset_adm_diff", "icu_outcome", "death_outcome", "dod",
               "diagnosis", "drg", "icd_code", "id_var", "id_var_df2",
               "weights", "threshold")),
      ends_with("_vax"),
      starts_with('vax_date_'),
      starts_with('vax_type'),
      starts_with('all_'),
      everything()
    )

  return(df)
}
