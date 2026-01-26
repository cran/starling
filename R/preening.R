#' @name preening
#' @title Prettification of infectious diseases datasets
#' @description Prettifies your dataset in preparation for data exploration and presenting tables. Adds variable labels and creates a series of age and time categories for analysis.
#' Just list the dataframe and it let it clean your variables and create exploratory variables. Use it as late in the workflow as possible, but, can be used at anytime.
#'
#' Classic workflow would be:
#' \enumerate{
#'  \item \code{\link{clean_the_nest}} to clean and prep data for linkage. Pay close attention to your linkage variables (letternames, date of birth, medicare number, gender and/or postcode), and ensure all dates are formatted as dates.
#'  \item \code{\link{murmuration}} to link cases to vaccination data (named here "c2v").
#'  \item \code{\link{murmuration}} to link c2v to hospitalization data (named here c2v2h). Of note, you can skip linking the vaccination dataset.
#'  \item \code{\link{preening}} to prettify the dataframe prepping it for exploration, analysis and presentation. Great to use with \code{gtsummary::tbl_summary()}.
#' }
#' @param df The dataset as a dataframe, which can be a case notifications dataset (infections), hospital admissions or vaccination dataset.
#' @param create_age_categories Logical. If TRUE (default), creates 21 standardized age category variables. Requires an 'age' variable in the dataset.
#' @param create_temporal_vars Logical. If TRUE (default), creates temporal variables (ISO weeks, quarters, months) for date columns.
#' @param calculate_age Logical. If TRUE (default), attempts to calculate age from dob if age variable is missing.
#' @param age_reference_date Character. Column name to use as reference date for age calculation if age is missing.
#'   If NULL (default), uses first available from: event_date, onset_date, admission_date, first_vax_date, last_vax_date, vax_date_*.
#'
#' @return The output is a dataframe with variable labels (useful for making pretty tables and graphics), and creates several age categories and time categories (month-year, quarter-year etc.)
#'
#' @details
#' This function enhances infectious disease datasets by:
#' \itemize{
#'   \item Adding descriptive variable labels for cleaner tables and graphics
#'   \item Creating comprehensive temporal variables (ISO weeks, quarters, months) from date fields
#'   \item Generating 21 standardized age category variables for flexible analysis
#'   \item Calculating age from date of birth if not already present
#'   \item Adding useful derived variables for epidemiological analysis
#' }
#'
#' IMPORTANT - Date Format Requirements:
#'
#' All date columns MUST be in R's Date format before using this function.
#' The function expects dates to already be properly formatted and will error
#' with a clear message if they are not.
#'
#' Common date conversions:
#' \itemize{
#'   \item From character: data$dob <- as.Date(data$dob, format = "\%Y-\%m-\%d")
#'   \item From character (alternative): data$dob <- lubridate::ymd(data$dob)
#'   \item From Excel dates: data$dob <- as.Date(data$dob, origin = "1899-12-30")
#'   \item Always check: class(data$dob) should return "Date"
#' }
#'
#' If you receive an error like "column must be in Date format", convert your
#' date columns first, then run preening().
#'
#' Age Categorization: If create_age_categories = TRUE and an 'age' variable exists (or can be calculated),
#' the function creates 21 standardized age category variables with nomenclature age[x]cat where x indicates
#' the number of categories:
#'
#' \describe{
#'   \item{age2cat}{2 categories: Pediatric vs Adult (<18, 18+)}
#'   \item{age3cat}{3 categories: Child, Adult, Older Adult (<18, 18-64, 65+)}
#'   \item{age4cat}{4 categories: Infant/Child, Young Adult, Adult, Older Adult (<5, 5-17, 18-64, 65+)}
#'   \item{age5cat}{5 categories: Standard public health categories (0-4, 5-17, 18-64, 65-74, 75+)}
#'   \item{age6cat}{6 categories: Granular infant categories (<1, 1-4, 5-17, 18-64, 65-74, 75+)}
#'   \item{age7cat}{7 categories: Fine pediatric cuts (<1, 1, 2-4, 5-11, 12-17, 18-64, 65+)}
#'   \item{age8cat}{8 categories: Infant subcategories (<3mo, 3-5mo, 6-11mo, 1-4, 5-17, 18-64, 65-74, 75+)}
#'   \item{age9cat}{9 categories: Monthly infant categories (<1mo, 1mo, 2-5mo, 6-11mo, 1-4, 5-17, 18-64, 65-74, 75+)}
#'   \item{age10cat}{10 categories: Decade bands (0-4, 5-9, 10-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+)}
#'   \item{age11cat}{11 categories: Fine older adult categories (0-4, 5-17, 18-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80-89, 90-99, 100+)}
#'   \item{age12cat}{12 categories: Detailed pediatric + adult decades (<1, 1-4, 5-9, 10-14, 15-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+)}
#'   \item{age13cat}{13 categories: Very fine infant + standard adult (<1mo, 1mo, 2mo, 3-5mo, 6-11mo, 1, 2-4, 5-11, 12-17, 18-39, 40-64, 65-79, 80+)}
#'   \item{age14cat}{14 categories: ABS-like with fine elderly (0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80-84, 85-89, 90+)}
#'   \item{age15cat}{15 categories: Vaccine schedule aligned (<2mo, 2-3mo, 4-5mo, 6-11mo, 1, 2-3, 4, 5-11, 12-17, 18-49, 50-64, 65-74, 75-84, 85-94, 95+)}
#'   \item{age16cat}{16 categories: Granular pediatric + 10-year adult bands (<1, 1, 2, 3, 4, 5-9, 10-14, 15-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80-89, 90+)}
#'   \item{age17cat}{17 categories: WHO/UNICEF standard with extensions (<1mo, 1-5mo, 6-11mo, 1, 2-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80-89, 90+)}
#'   \item{age18cat}{18 categories: Standard 5-year bands (census/ABS style) (0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, 80-84, 85+)}
#'   \item{age19cat}{19 categories: Extended 5-year bands with fine elderly (0-4, 5-9, ..., 80-84, 85-89, 90+)}
#'   \item{age20cat}{20 categories: Monthly up to 12 months + standard thereafter (<1mo, 1mo, 2mo, 3mo, 4mo, 5mo, 6mo, 7mo, 8mo, 9mo, 10mo, 11mo, 1-4, 5-17, 18-39, 40-64, 65-74, 75-84, 85-94, 95+)}
#'   \item{age21cat}{21 categories: Comprehensive life course categories (<1mo, 1-2mo, 3-5mo, 6-11mo, 1, 2-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-74, 75-84, 85+)}
#' }
#'
#' @export
preening <- function(df,
                     create_age_categories = TRUE,
                     create_temporal_vars = TRUE,
                     calculate_age = TRUE,
                     age_reference_date = NULL) {

  # Check for required packages
  required_pkgs <- c("dplyr", "lubridate")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required but not installed. ",
                  "Please install it with: install.packages('", pkg, "')"),
           call. = FALSE)
    }
  }

  # Create a function that sets labels without using tiny_labelled
  set_column_label <- function(data, column_name, label_text) {
    if (column_name %in% colnames(data)) {
      attr(data[[column_name]], "label") <- label_text
    }
    return(data)
  }

  # Calculate age if requested and missing
  if (calculate_age && !"age" %in% colnames(df)) {
    df <- calculate_age_if_missing(df, age_reference_date)
  }

  # Identity metadata
  df <- set_column_label(df, "id_var", "ID")
  df <- set_column_label(df, "event_id_var", "Event ID")
  df <- set_column_label(df, "lettername1", "First Name")
  df <- set_column_label(df, "lettername2", "Last Name")
  df <- set_column_label(df, "dob", "Date of Birth")
  df <- set_column_label(df, "postcode", "Post Code")

  # Medicare information
  df <- set_column_label(df, "medicare", "Medicare Number")
  df <- set_column_label(df, "medicare11", "Medicare (11 digits)")
  df <- set_column_label(df, "medicare10", "Medicare (10 digits)")

  # Demographic information
  df <- set_column_label(df, "gender", "Gender")
  df <- set_column_label(df, "fn", "First Nations Status")
  df <- set_column_label(df, "diagnosis", "Pathogen")

  # Create temporal variables if requested
  if (create_temporal_vars) {
    # Onset date and derived variables
    if ("onset_date" %in% colnames(df)) {
      df <- set_column_label(df, "onset_date", "Onset Date")
      df <- df %>%
        dplyr::mutate(
          onset_iso_week = lubridate::epiweek(onset_date),
          onset_year = lubridate::year(onset_date),
          onset_year_iso_week = format(onset_date, format = "%Y-W%V"),
          onset_q = lubridate::quarter(onset_date, type = "year.quarter"),
          onset_q_number = lubridate::quarter(onset_date),
          onset_m = format(onset_date, format = "%b-%Y"),
          onset_m_number = format(onset_date, format = "%m"),
          onset_m_human = format(onset_date, format = "%Y-%b")
        )

      # Set labels for onset derived variables
      df <- set_column_label(df, "onset_iso_week", "Onset CDC ISO Week")
      df <- set_column_label(df, "onset_year", "Onset Year")
      df <- set_column_label(df, "onset_year_iso_week", "Onset Year CDC ISO Week")
      df <- set_column_label(df, "onset_q", "Onset Quarter")
      df <- set_column_label(df, "onset_q_number", "Onset Quarter Number")
      df <- set_column_label(df, "onset_m", "Onset Month")
      df <- set_column_label(df, "onset_m_number", "Onset Month Number")
      df <- set_column_label(df, "onset_m_human", "Onset Month (Human Readable)")
    }

    # Admission date and derived variables
    if ("admission_date" %in% colnames(df)) {
      df <- set_column_label(df, "admission_date", "Admission Date")
      df <- df %>%
        dplyr::mutate(
          admission_iso_week = lubridate::epiweek(admission_date),
          admission_year_iso_week = format(admission_date, format = "%Y-W%V"),
          admission_y = lubridate::year(admission_date),
          admission_q = lubridate::quarter(admission_date, type = "year.quarter"),
          admission_q_number = lubridate::quarter(admission_date),
          admission_m = format(admission_date, format = "%Y-%m"),
          admission_m_number = format(admission_date, format = "%m"),
          admission_m_human = format(admission_date, format = "%Y-%b")
        )

      # Set labels for admission derived variables
      df <- set_column_label(df, "admission_iso_week", "Admission CDC ISO Week")
      df <- set_column_label(df, "admission_year_iso_week", "Admission Year CDC ISO Week")
      df <- set_column_label(df, "admission_y", "Admission Year")
      df <- set_column_label(df, "admission_q", "Admission Quarter")
      df <- set_column_label(df, "admission_q_number", "Admission Quarter Number")
      df <- set_column_label(df, "admission_m", "Admission Month")
      df <- set_column_label(df, "admission_m_number", "Admission Month Number")
      df <- set_column_label(df, "admission_m_human", "Admission Month (Human Readable)")
    }

    # Event date temporal variables
    if ("event_date" %in% colnames(df)) {
      df <- set_column_label(df, "event_date", "Event Date")
      df <- df %>%
        dplyr::mutate(
          event_iso_week = lubridate::epiweek(event_date),
          event_year = lubridate::year(event_date),
          event_year_iso_week = format(event_date, format = "%Y-W%V"),
          event_q = lubridate::quarter(event_date, type = "year.quarter"),
          event_q_number = lubridate::quarter(event_date),
          event_m = format(event_date, format = "%b-%Y"),
          event_m_number = format(event_date, format = "%m"),
          event_m_human = format(event_date, format = "%Y-%b")
        )

      # Set labels for event derived variables
      df <- set_column_label(df, "event_iso_week", "Event CDC ISO Week")
      df <- set_column_label(df, "event_year", "Event Year")
      df <- set_column_label(df, "event_year_iso_week", "Event Year CDC ISO Week")
      df <- set_column_label(df, "event_q", "Event Quarter")
      df <- set_column_label(df, "event_q_number", "Event Quarter Number")
      df <- set_column_label(df, "event_m", "Event Month")
      df <- set_column_label(df, "event_m_number", "Event Month Number")
      df <- set_column_label(df, "event_m_human", "Event Month (Human Readable)")
    }
  }

  # Admission information
  df <- set_column_label(df, "admission_outcome", "Admission Outcome")

  # Discharge information
  df <- set_column_label(df, "discharge_date", "Discharge Date")
  df <- set_column_label(df, "all_discharge_dates", "All Discharge Dates")

  # Length of stay
  df <- set_column_label(df, "los", "Length of Stay")
  df <- set_column_label(df, "all_los", "All Length of Stays")

  # Hospital information
  df <- set_column_label(df, "hospital", "Admission Location")
  df <- set_column_label(df, "all_hospitals", "All Hospital Locations")

  # Diagnosis codes
  df <- set_column_label(df, "icd_code", "ICD Code")
  df <- set_column_label(df, "all_icd_codes", "All ICD Codes")
  df <- set_column_label(df, "diagnosis_description", "Diagnosis Description")
  df <- set_column_label(df, "all_diag_desc", "All Diagnosis Descriptions")
  df <- set_column_label(df, "drg", "DRG Code")
  df <- set_column_label(df, "all_drgs", "All DRGs")

  # ICU and dialysis information
  df <- set_column_label(df, "icu_date", "ICU Admission Date")
  df <- set_column_label(df, "icu_outcome", "ICU Admission")
  df <- set_column_label(df, "icu_hours", "Hours in ICU")
  df <- set_column_label(df, "dialysis", "Dialysis")
  df <- set_column_label(df, "dialysis_outcome", "Dialysis Event")
  df <- set_column_label(df, "onset_adm_diff", "Days between Onset and Admission")
  df <- set_column_label(df, "first_onset_adm_diff", "Days between First Onset and Admission")
  df <- set_column_label(df, "all_onset_adm_diffs", "Days between All Onsets and Admissions")

  # First admission information
  df <- set_column_label(df, "first_admission_date", "First Admission Date")
  df <- set_column_label(df, "last_admission_date", "Last Admission Date")
  df <- set_column_label(df, "all_admission_dates", "All Admission Dates")
  df <- set_column_label(df, "first_admission_iso_week", "First Admission CDC ISO Week")
  df <- set_column_label(df, "first_admission_year_iso_week", "First Admission Year CDC ISO Week")
  df <- set_column_label(df, "first_admission_q", "First Admission Quarter")
  df <- set_column_label(df, "first_admission_m", "First Admission Month")
  df <- set_column_label(df, "first_discharge_date", "First Discharge Date")
  df <- set_column_label(df, "last_discharge_date", "Last Discharge Date")
  df <- set_column_label(df, "first_los", "First Length of Stay")
  df <- set_column_label(df, "first_hospital", "First Admission Location")
  df <- set_column_label(df, "first_icd_code", "First ICD Code")
  df <- set_column_label(df, "first_diagnosis_description", "First Diagnosis Description")
  df <- set_column_label(df, "first_drg", "First DRG Code")
  df <- set_column_label(df, "first_icu_date", "First ICU Admission Date")
  df <- set_column_label(df, "first_icu_outcome", "First ICU Admission")
  df <- set_column_label(df, "first_icu_hours", "First Hours in ICU")
  df <- set_column_label(df, "first_dialysis", "First Dialysis")
  df <- set_column_label(df, "first_dialysis_outcome", "First Dialysis Event")

  # Other information
  df <- set_column_label(df, "genomics", "Genomic Subtype")
  df <- set_column_label(df, "dod", "Date of Death")
  df <- set_column_label(df, "death_outcome", "Life Status")
  df <- set_column_label(df, "vax_type", "Vaccine Type")
  df <- set_column_label(df, "vax_date", "Vaccination Date")
  df <- set_column_label(df, "last_follow_up", "Last Follow-Up")
  df <- set_column_label(df, "last_vax", "Last Vaccination Date")
  df <- set_column_label(df, "tsv", "Time Since Last Vaccination")
  df <- set_column_label(df, "vaccination_status", "Vaccination Status at Time of Event")
  df <- set_column_label(df, "vaccination_status_num", "Vaccination Status at Time of Event (Numeric)")
  df <- set_column_label(df, "total_admissions", "Total Admissions")
  df <- set_column_label(df, "total_los", "Total Length of Stay")
  df <- set_column_label(df, "weights", "Linkage Score")
  df <- set_column_label(df, "threshold", "Threshold Met")
  df <- set_column_label(df, "id_var_df2", "ID (Linked Data)")

  # Create comprehensive age categories if requested
  if (create_age_categories && "age" %in% colnames(df)) {
    df <- set_column_label(df, "age", "Age")
    df <- create_comprehensive_age_categories(df)

    # Set labels for all age categories
    df <- set_column_label(df, "age2cat", "Age: 2 Categories")
    df <- set_column_label(df, "age3cat", "Age: 3 Categories")
    df <- set_column_label(df, "age4cat", "Age: 4 Categories")
    df <- set_column_label(df, "age5cat", "Age: 5 Categories")
    df <- set_column_label(df, "age6cat", "Age: 6 Categories")
    df <- set_column_label(df, "age7cat", "Age: 7 Categories")
    df <- set_column_label(df, "age8cat", "Age: 8 Categories")
    df <- set_column_label(df, "age9cat", "Age: 9 Categories")
    df <- set_column_label(df, "age10cat", "Age: 10 Categories (Decades)")
    df <- set_column_label(df, "age11cat", "Age: 11 Categories")
    df <- set_column_label(df, "age12cat", "Age: 12 Categories")
    df <- set_column_label(df, "age13cat", "Age: 13 Categories (Fine Infant)")
    df <- set_column_label(df, "age14cat", "Age: 14 Categories (ABS-like)")
    df <- set_column_label(df, "age15cat", "Age: 15 Categories (Vaccine Schedule)")
    df <- set_column_label(df, "age16cat", "Age: 16 Categories")
    df <- set_column_label(df, "age17cat", "Age: 17 Categories (WHO/UNICEF)")
    df <- set_column_label(df, "age18cat", "Age: 18 Categories (Census/ABS)")
    df <- set_column_label(df, "age19cat", "Age: 19 Categories (Extended 5-year)")
    df <- set_column_label(df, "age20cat", "Age: 20 Categories (Monthly Infant)")
    df <- set_column_label(df, "age21cat", "Age: 21 Categories (Life Course)")
  } else if (create_age_categories && !"age" %in% colnames(df)) {
    message("Age categories not created: 'age' variable not found in dataset.")
  }

  # Set an attribute on the data frame that it has been preprocessed
  attr(df, "preened") <- TRUE

  # Summary message
  message(paste0(
    "\n--- Preening Complete ---\n",
    "Rows processed: ", nrow(df), "\n",
    "Variables labeled: ", sum(sapply(df, function(x) !is.null(attr(x, "label")))), "\n",
    if (create_age_categories && "age" %in% colnames(df)) {
      "Age categories: 21 variables created\n"
    } else {
      ""
    },
    if (create_temporal_vars) {
      paste0("Temporal variables created for: ",
             paste(intersect(c("onset_date", "admission_date", "event_date"), names(df)),
                   collapse = ", "), "\n")
    } else {
      ""
    }
  ))

  return(df)
}


#' Calculate Age if Missing
#'
#' Internal function to calculate age from date of birth if not present
#'
#' @param data A data frame
#' @param reference_date Character. Column name to use as reference date.
#'   If NULL, uses first available from standard date columns.
#'
#' @return Data frame with age variable added if possible
#' @keywords internal
calculate_age_if_missing <- function(data, reference_date = NULL) {

  # Look for dob variable
  dob_col <- find_column_preening(data, c("dob", "date_of_birth", "birth_date", "birthdate"))

  if (is.null(dob_col)) {
    message("Cannot calculate age: no date of birth variable found.")
    return(data)
  }

  # Check if dob column is in Date format
  if (!inherits(data[[dob_col]], "Date")) {
    stop(paste0("The '", dob_col, "' column must be in Date format. ",
                "Please convert it to Date format before using preening(). ",
                "Example: data$", dob_col, " <- as.Date(data$", dob_col, ")"),
         call. = FALSE)
  }

  # Determine reference date
  if (!is.null(reference_date)) {
    if (!reference_date %in% names(data)) {
      warning(paste0("Specified reference_date '", reference_date,
                     "' not found in data. Attempting to find alternative."))
      reference_date <- NULL
    } else {
      date_col <- reference_date
    }
  }

  if (is.null(reference_date)) {
    # Look for reference date in priority order
    date_col <- find_column_preening(data, c(
      "event_date",
      "onset_date",
      "admission_date",
      "first_vax_date",
      "last_vax_date"
    ))

    # If not found, look for any vax_date_* column
    if (is.null(date_col)) {
      vax_date_cols <- grep("^vax_date_\\d+$", names(data), value = TRUE, ignore.case = TRUE)
      if (length(vax_date_cols) > 0) {
        date_col <- vax_date_cols[1]
      }
    }
  }

  if (is.null(date_col)) {
    message("Cannot calculate age: no reference date variable found.")
    return(data)
  }

  # Check if reference date column is in Date format
  if (!inherits(data[[date_col]], "Date")) {
    stop(paste0("The '", date_col, "' column must be in Date format. ",
                "Please convert it to Date format before using preening(). ",
                "Example: data$", date_col, " <- as.Date(data$", date_col, ")"),
         call. = FALSE)
  }

  message(paste0("Calculating age from: ", date_col, " - ", dob_col))

  # Calculate age directly from Date columns
  data <- data %>%
    dplyr::mutate(
      age = as.numeric(lubridate::time_length(
        lubridate::interval(.data[[dob_col]], .data[[date_col]]),
        "years"
      ))
    )

  return(data)
}


#' Create Comprehensive Age Categories
#'
#' Internal function to create all 21 age category variables
#'
#' @param data A data frame with an 'age' variable
#'
#' @return Data frame with age category variables added
#' @keywords internal
create_comprehensive_age_categories <- function(data) {

  if (!"age" %in% names(data)) {
    warning("Cannot create age categories: 'age' variable not found.")
    return(data)
  }

  age_years <- data$age

  data <- data %>%
    dplyr::mutate(
      # 2 categories: Pediatric vs Adult
      age2cat = dplyr::case_when(
        age_years < 18 ~ "<18",
        age_years >= 18 ~ "18+",
        TRUE ~ NA_character_
      ),

      # 3 categories: Child, Adult, Older Adult
      age3cat = dplyr::case_when(
        age_years < 18 ~ "<18",
        age_years >= 18 & age_years < 65 ~ "18-64",
        age_years >= 65 ~ "65+",
        TRUE ~ NA_character_
      ),

      # 4 categories: Infant/Child, Young Adult, Adult, Older Adult
      age4cat = dplyr::case_when(
        age_years < 5 ~ "<5",
        age_years >= 5 & age_years < 18 ~ "5-17",
        age_years >= 18 & age_years < 65 ~ "18-64",
        age_years >= 65 ~ "65+",
        TRUE ~ NA_character_
      ),

      # 5 categories: Standard public health categories
      age5cat = dplyr::case_when(
        age_years < 5 ~ "0-4",
        age_years >= 5 & age_years < 18 ~ "5-17",
        age_years >= 18 & age_years < 65 ~ "18-64",
        age_years >= 65 & age_years < 75 ~ "65-74",
        age_years >= 75 ~ "75+",
        TRUE ~ NA_character_
      ),

      # 6 categories: Granular infant categories
      age6cat = dplyr::case_when(
        age_years < 1 ~ "<1",
        age_years >= 1 & age_years < 5 ~ "1-4",
        age_years >= 5 & age_years < 18 ~ "5-17",
        age_years >= 18 & age_years < 65 ~ "18-64",
        age_years >= 65 & age_years < 75 ~ "65-74",
        age_years >= 75 ~ "75+",
        TRUE ~ NA_character_
      ),

      # 7 categories: Fine pediatric cuts
      age7cat = dplyr::case_when(
        age_years < 1 ~ "<1",
        age_years >= 1 & age_years < 2 ~ "1",
        age_years >= 2 & age_years < 5 ~ "2-4",
        age_years >= 5 & age_years < 12 ~ "5-11",
        age_years >= 12 & age_years < 18 ~ "12-17",
        age_years >= 18 & age_years < 65 ~ "18-64",
        age_years >= 65 ~ "65+",
        TRUE ~ NA_character_
      ),

      # 8 categories: Infant subcategories
      age8cat = dplyr::case_when(
        age_years < 0.25 ~ "<3mo",
        age_years >= 0.25 & age_years < 0.5 ~ "3-5mo",
        age_years >= 0.5 & age_years < 1 ~ "6-11mo",
        age_years >= 1 & age_years < 5 ~ "1-4",
        age_years >= 5 & age_years < 18 ~ "5-17",
        age_years >= 18 & age_years < 65 ~ "18-64",
        age_years >= 65 & age_years < 75 ~ "65-74",
        age_years >= 75 ~ "75+",
        TRUE ~ NA_character_
      ),

      # 9 categories: Monthly infant categories
      age9cat = dplyr::case_when(
        age_years < 1/12 ~ "<1mo",
        age_years >= 1/12 & age_years < 2/12 ~ "1mo",
        age_years >= 2/12 & age_years < 6/12 ~ "2-5mo",
        age_years >= 6/12 & age_years < 12/12 ~ "6-11mo",
        age_years >= 1 & age_years < 5 ~ "1-4",
        age_years >= 5 & age_years < 18 ~ "5-17",
        age_years >= 18 & age_years < 65 ~ "18-64",
        age_years >= 65 & age_years < 75 ~ "65-74",
        age_years >= 75 ~ "75+",
        TRUE ~ NA_character_
      ),

      # 10 categories: Decade bands
      age10cat = dplyr::case_when(
        age_years < 5 ~ "0-4",
        age_years >= 5 & age_years < 10 ~ "5-9",
        age_years >= 10 & age_years < 20 ~ "10-19",
        age_years >= 20 & age_years < 30 ~ "20-29",
        age_years >= 30 & age_years < 40 ~ "30-39",
        age_years >= 40 & age_years < 50 ~ "40-49",
        age_years >= 50 & age_years < 60 ~ "50-59",
        age_years >= 60 & age_years < 70 ~ "60-69",
        age_years >= 70 & age_years < 80 ~ "70-79",
        age_years >= 80 ~ "80+",
        TRUE ~ NA_character_
      ),

      # 11 categories: Fine older adult categories
      age11cat = dplyr::case_when(
        age_years < 5 ~ "0-4",
        age_years >= 5 & age_years < 18 ~ "5-17",
        age_years >= 18 & age_years < 30 ~ "18-29",
        age_years >= 30 & age_years < 40 ~ "30-39",
        age_years >= 40 & age_years < 50 ~ "40-49",
        age_years >= 50 & age_years < 60 ~ "50-59",
        age_years >= 60 & age_years < 70 ~ "60-69",
        age_years >= 70 & age_years < 80 ~ "70-79",
        age_years >= 80 & age_years < 90 ~ "80-89",
        age_years >= 90 & age_years < 100 ~ "90-99",
        age_years >= 100 ~ "100+",
        TRUE ~ NA_character_
      ),

      # 12 categories: Detailed pediatric + adult decades
      age12cat = dplyr::case_when(
        age_years < 1 ~ "<1",
        age_years >= 1 & age_years < 5 ~ "1-4",
        age_years >= 5 & age_years < 10 ~ "5-9",
        age_years >= 10 & age_years < 15 ~ "10-14",
        age_years >= 15 & age_years < 20 ~ "15-19",
        age_years >= 20 & age_years < 30 ~ "20-29",
        age_years >= 30 & age_years < 40 ~ "30-39",
        age_years >= 40 & age_years < 50 ~ "40-49",
        age_years >= 50 & age_years < 60 ~ "50-59",
        age_years >= 60 & age_years < 70 ~ "60-69",
        age_years >= 70 & age_years < 80 ~ "70-79",
        age_years >= 80 ~ "80+",
        TRUE ~ NA_character_
      ),

      # 13 categories: Very fine infant + standard adult
      age13cat = dplyr::case_when(
        age_years < 1/12 ~ "<1mo",
        age_years >= 1/12 & age_years < 2/12 ~ "1mo",
        age_years >= 2/12 & age_years < 3/12 ~ "2mo",
        age_years >= 3/12 & age_years < 6/12 ~ "3-5mo",
        age_years >= 6/12 & age_years < 12/12 ~ "6-11mo",
        age_years >= 1 & age_years < 2 ~ "1",
        age_years >= 2 & age_years < 5 ~ "2-4",
        age_years >= 5 & age_years < 12 ~ "5-11",
        age_years >= 12 & age_years < 18 ~ "12-17",
        age_years >= 18 & age_years < 40 ~ "18-39",
        age_years >= 40 & age_years < 65 ~ "40-64",
        age_years >= 65 & age_years < 80 ~ "65-79",
        age_years >= 80 ~ "80+",
        TRUE ~ NA_character_
      ),

      # 14 categories: ABS-like with fine elderly
      age14cat = dplyr::case_when(
        age_years < 5 ~ "0-4",
        age_years >= 5 & age_years < 10 ~ "5-9",
        age_years >= 10 & age_years < 15 ~ "10-14",
        age_years >= 15 & age_years < 20 ~ "15-19",
        age_years >= 20 & age_years < 25 ~ "20-24",
        age_years >= 25 & age_years < 30 ~ "25-29",
        age_years >= 30 & age_years < 40 ~ "30-39",
        age_years >= 40 & age_years < 50 ~ "40-49",
        age_years >= 50 & age_years < 60 ~ "50-59",
        age_years >= 60 & age_years < 70 ~ "60-69",
        age_years >= 70 & age_years < 80 ~ "70-79",
        age_years >= 80 & age_years < 85 ~ "80-84",
        age_years >= 85 & age_years < 90 ~ "85-89",
        age_years >= 90 ~ "90+",
        TRUE ~ NA_character_
      ),

      # 15 categories: Vaccine schedule aligned
      age15cat = dplyr::case_when(
        age_years < 2/12 ~ "<2mo",
        age_years >= 2/12 & age_years < 4/12 ~ "2-3mo",
        age_years >= 4/12 & age_years < 6/12 ~ "4-5mo",
        age_years >= 6/12 & age_years < 12/12 ~ "6-11mo",
        age_years >= 1 & age_years < 2 ~ "1",
        age_years >= 2 & age_years < 4 ~ "2-3",
        age_years >= 4 & age_years < 5 ~ "4",
        age_years >= 5 & age_years < 12 ~ "5-11",
        age_years >= 12 & age_years < 18 ~ "12-17",
        age_years >= 18 & age_years < 50 ~ "18-49",
        age_years >= 50 & age_years < 65 ~ "50-64",
        age_years >= 65 & age_years < 75 ~ "65-74",
        age_years >= 75 & age_years < 85 ~ "75-84",
        age_years >= 85 & age_years < 95 ~ "85-94",
        age_years >= 95 ~ "95+",
        TRUE ~ NA_character_
      ),

      # 16 categories: Granular pediatric + 10-year adult bands
      age16cat = dplyr::case_when(
        age_years < 1 ~ "<1",
        age_years >= 1 & age_years < 2 ~ "1",
        age_years >= 2 & age_years < 3 ~ "2",
        age_years >= 3 & age_years < 4 ~ "3",
        age_years >= 4 & age_years < 5 ~ "4",
        age_years >= 5 & age_years < 10 ~ "5-9",
        age_years >= 10 & age_years < 15 ~ "10-14",
        age_years >= 15 & age_years < 20 ~ "15-19",
        age_years >= 20 & age_years < 30 ~ "20-29",
        age_years >= 30 & age_years < 40 ~ "30-39",
        age_years >= 40 & age_years < 50 ~ "40-49",
        age_years >= 50 & age_years < 60 ~ "50-59",
        age_years >= 60 & age_years < 70 ~ "60-69",
        age_years >= 70 & age_years < 80 ~ "70-79",
        age_years >= 80 & age_years < 90 ~ "80-89",
        age_years >= 90 ~ "90+",
        TRUE ~ NA_character_
      ),

      # 17 categories: WHO/UNICEF standard with extensions
      age17cat = dplyr::case_when(
        age_years < 1/12 ~ "<1mo",
        age_years >= 1/12 & age_years < 6/12 ~ "1-5mo",
        age_years >= 6/12 & age_years < 12/12 ~ "6-11mo",
        age_years >= 1 & age_years < 2 ~ "1",
        age_years >= 2 & age_years < 5 ~ "2-4",
        age_years >= 5 & age_years < 10 ~ "5-9",
        age_years >= 10 & age_years < 15 ~ "10-14",
        age_years >= 15 & age_years < 20 ~ "15-19",
        age_years >= 20 & age_years < 25 ~ "20-24",
        age_years >= 25 & age_years < 30 ~ "25-29",
        age_years >= 30 & age_years < 40 ~ "30-39",
        age_years >= 40 & age_years < 50 ~ "40-49",
        age_years >= 50 & age_years < 60 ~ "50-59",
        age_years >= 60 & age_years < 70 ~ "60-69",
        age_years >= 70 & age_years < 80 ~ "70-79",
        age_years >= 80 & age_years < 90 ~ "80-89",
        age_years >= 90 ~ "90+",
        TRUE ~ NA_character_
      ),

      # 18 categories: Standard 5-year bands (census/ABS style)
      age18cat = dplyr::case_when(
        age_years < 5 ~ "0-4",
        age_years >= 5 & age_years < 10 ~ "5-9",
        age_years >= 10 & age_years < 15 ~ "10-14",
        age_years >= 15 & age_years < 20 ~ "15-19",
        age_years >= 20 & age_years < 25 ~ "20-24",
        age_years >= 25 & age_years < 30 ~ "25-29",
        age_years >= 30 & age_years < 35 ~ "30-34",
        age_years >= 35 & age_years < 40 ~ "35-39",
        age_years >= 40 & age_years < 45 ~ "40-44",
        age_years >= 45 & age_years < 50 ~ "45-49",
        age_years >= 50 & age_years < 55 ~ "50-54",
        age_years >= 55 & age_years < 60 ~ "55-59",
        age_years >= 60 & age_years < 65 ~ "60-64",
        age_years >= 65 & age_years < 70 ~ "65-69",
        age_years >= 70 & age_years < 75 ~ "70-74",
        age_years >= 75 & age_years < 80 ~ "75-79",
        age_years >= 80 & age_years < 85 ~ "80-84",
        age_years >= 85 ~ "85+",
        TRUE ~ NA_character_
      ),

      # 19 categories: Extended 5-year bands with fine elderly
      age19cat = dplyr::case_when(
        age_years < 5 ~ "0-4",
        age_years >= 5 & age_years < 10 ~ "5-9",
        age_years >= 10 & age_years < 15 ~ "10-14",
        age_years >= 15 & age_years < 20 ~ "15-19",
        age_years >= 20 & age_years < 25 ~ "20-24",
        age_years >= 25 & age_years < 30 ~ "25-29",
        age_years >= 30 & age_years < 35 ~ "30-34",
        age_years >= 35 & age_years < 40 ~ "35-39",
        age_years >= 40 & age_years < 45 ~ "40-44",
        age_years >= 45 & age_years < 50 ~ "45-49",
        age_years >= 50 & age_years < 55 ~ "50-54",
        age_years >= 55 & age_years < 60 ~ "55-59",
        age_years >= 60 & age_years < 65 ~ "60-64",
        age_years >= 65 & age_years < 70 ~ "65-69",
        age_years >= 70 & age_years < 75 ~ "70-74",
        age_years >= 75 & age_years < 80 ~ "75-79",
        age_years >= 80 & age_years < 85 ~ "80-84",
        age_years >= 85 & age_years < 90 ~ "85-89",
        age_years >= 90 ~ "90+",
        TRUE ~ NA_character_
      ),

      # 20 categories: Monthly up to 12 months + standard thereafter
      age20cat = dplyr::case_when(
        age_years < 1/12 ~ "<1mo",
        age_years >= 1/12 & age_years < 2/12 ~ "1mo",
        age_years >= 2/12 & age_years < 3/12 ~ "2mo",
        age_years >= 3/12 & age_years < 4/12 ~ "3mo",
        age_years >= 4/12 & age_years < 5/12 ~ "4mo",
        age_years >= 5/12 & age_years < 6/12 ~ "5mo",
        age_years >= 6/12 & age_years < 7/12 ~ "6mo",
        age_years >= 7/12 & age_years < 8/12 ~ "7mo",
        age_years >= 8/12 & age_years < 9/12 ~ "8mo",
        age_years >= 9/12 & age_years < 10/12 ~ "9mo",
        age_years >= 10/12 & age_years < 11/12 ~ "10mo",
        age_years >= 11/12 & age_years < 12/12 ~ "11mo",
        age_years >= 1 & age_years < 5 ~ "1-4",
        age_years >= 5 & age_years < 18 ~ "5-17",
        age_years >= 18 & age_years < 40 ~ "18-39",
        age_years >= 40 & age_years < 65 ~ "40-64",
        age_years >= 65 & age_years < 75 ~ "65-74",
        age_years >= 75 & age_years < 85 ~ "75-84",
        age_years >= 85 & age_years < 95 ~ "85-94",
        age_years >= 95 ~ "95+",
        TRUE ~ NA_character_
      ),

      # 21 categories: Comprehensive life course categories
      age21cat = dplyr::case_when(
        age_years < 1/12 ~ "<1mo",
        age_years >= 1/12 & age_years < 3/12 ~ "1-2mo",
        age_years >= 3/12 & age_years < 6/12 ~ "3-5mo",
        age_years >= 6/12 & age_years < 12/12 ~ "6-11mo",
        age_years >= 1 & age_years < 2 ~ "1",
        age_years >= 2 & age_years < 5 ~ "2-4",
        age_years >= 5 & age_years < 10 ~ "5-9",
        age_years >= 10 & age_years < 15 ~ "10-14",
        age_years >= 15 & age_years < 20 ~ "15-19",
        age_years >= 20 & age_years < 25 ~ "20-24",
        age_years >= 25 & age_years < 30 ~ "25-29",
        age_years >= 30 & age_years < 35 ~ "30-34",
        age_years >= 35 & age_years < 40 ~ "35-39",
        age_years >= 40 & age_years < 45 ~ "40-44",
        age_years >= 45 & age_years < 50 ~ "45-49",
        age_years >= 50 & age_years < 55 ~ "50-54",
        age_years >= 55 & age_years < 60 ~ "55-59",
        age_years >= 60 & age_years < 65 ~ "60-64",
        age_years >= 65 & age_years < 75 ~ "65-74",
        age_years >= 75 & age_years < 85 ~ "75-84",
        age_years >= 85 ~ "85+",
        TRUE ~ NA_character_
      )
    )

  # Convert all age categories to factors with proper ordering
  age_cat_vars <- paste0("age", 2:21, "cat")
  for (var in age_cat_vars) {
    if (var %in% names(data)) {
      unique_levels <- sort(unique(data[[var]][!is.na(data[[var]])]))
      data[[var]] <- factor(data[[var]], levels = unique_levels, ordered = TRUE)
    }
  }

  message("Created 21 comprehensive age category variables.")

  return(data)
}


#' Find Column by Pattern Matching
#'
#' Internal helper function to find a column matching a list of patterns
#'
#' @param data A data frame
#' @param patterns Character vector of column name patterns to search for
#'
#' @return Name of the first matching column, or NULL if none found
#' @keywords internal
find_column_preening <- function(data, patterns) {
  for (pattern in patterns) {
    # Try exact match first
    if (pattern %in% names(data)) {
      return(pattern)
    }

    # Try case-insensitive match
    matches <- grep(paste0("^", pattern, "$"), names(data),
                    ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      return(matches[1])
    }
  }

  return(NULL)
}
