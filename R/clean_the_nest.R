#' @name clean_the_nest
#' @title Clean datasets and establishes common variable name nomenclature
#' @description Cleans three dataset types and prepares them for data-linkage. This command is the first step in creating the datasets for analysis. Building a solid "nest" is akin to building a solid foundation for future work. Of note, Starlings are cavity nesters,
#' meaning that they prefer to build their homes inside holes and crevices. This command is meant to work with diagnosis datasets (linelists like Notifiable Conditions registers) and, hospitalization datasets (administrative datasets), and vaccination datasets.
#' This command is used to prepare datasets for linkage with \code{\link{murmuration}}. There are no mandatory variables to include. However, a dataset of infections would include at minimum an onset date (date of diagnosis),
#' a dataset of admissions would include admission dates, and a dataset of vaccinations would include dates of vaccination and type of vaccines.
#' All of the datasets should include information that would allow for data-linkage, such as first name, last name, date of birth, address etc etc.
#'
#' Classic workflow would be:
#' \enumerate{
#'  \item \code{\link{clean_the_nest}} to clean and prep data for linkage. Pay close attention to your linkage variables (letternames, date of birth, medicare number, gender and/or postcode), and ensure all dates are formatted as dates.
#'  \item \code{\link{murmuration}} to link cases to vaccination data (named here "c2v").
#'  \item \code{\link{murmuration}} to link c2v to hospitalization data (named here c2v2h). Of note, you can skip linking the vaccination dataset.
#'  \item \code{\link{preening}} to prettify the dataframe prepping it for exploration, analysis and presentation. Great to use with \code{gtsummary::tbl_summary()}.
#' }
#' @param data The dataset, which can be a case notifications dataset (infections), hospital admissions or vaccination dataset (must pre-specify if it is a vaccinations dataset). Make sure dates are in date format.
#' @param drop_eggs This effectively drops the variables that are not being used. May turn this off if you need lots of extra information, but certainly good for the early stages of an analysis. Enables a lean dataset.
#' @param data_type Three options: "vaccination", "hospital", or "cases". The key information required is that for linkage, and the vaccination events. No age or age categories will be calculated if it is a vaccination dataset.
#' @param lie_nest_flat Takes a long vaccination dataset (like Australian Immunization Register; 1 or more rows per person) and turns it into a wide dataset - one row per person
#' @param drop_the_na_vax Drops (removes) vaccines that are listed as having no names.
#' @param keep_vars Vector list of variables. Variables in a vector list with quotation marks, as it will be used in a select statement.
#' @param id_var Any format as long as unique to individual. This is important This ID variable is critical. Must ensure for case data that it only has one row per person, or first infection only. Identifies the multiple rows associated with a person who has multiple vaccines, admissions or infections. Cannot have missing data, or the observation will be lost in the linking process.
#' @param event_id_var Any format as long as unique for the whole dataset. This represents the ID of the vaccination event, or the hospitalization event, which MUST be distinct. A person (id_var) can have multiple events (event_id). Some datasets will surprise you with multiple entries for the same admission.
#' @param diagnosis Character format. The column with the infectious disease diagnosis listed. e.g. COVID-19, SARS-CoV-2, RSV, Influenza.
#' @param lettername1 Character format. First Name variable. If there is a second first name (some cases this might be a middle name), it will be removed during cleaning. All non-alphanumeric characters will be removed and everything becomes lower case.
#' @param lettername2 Character format. Last name variable. All non-alphanumeric characters will be removed and everything becomes lower case. Two part last names will be kept.
#' @param dob Date format. The date of birth (make sure dates are in date format).
#' @param age Numeric format. Include age only if it has been pre-specified in the dataset, and you don't want it re-calculated.
#' @param gender Character format. Pay close attention that your genders are in a similar format for data-linkage - "F", vs "0" vs "Female". This is left up to the user to clean.
#' @param medicare Numeric format. Medicare number. A medicare number with 9, 10 and 11 numbers will have been created. In Australia, the 10th number represents the card ID, and the 11th number represents the person ID. A family or individual will get a new card id (10th digit) every time their card expires.
#' @param postcode Numeric format. Post code of person with no restriction on the number of digits.
#' @param fn Character format. First Nations Status.
#' @param latitude Numeric format. Latitude of address. Not explicitly required for linkage.
#' @param longitude Numeric format. Longitude of address. Not explicitly required for linkage.
#' @param onset_date Date format. Onset date of the illness. Commonly the date of diagnosis (date of the lab test or date of the first symptom). Must be in date format.
#' @param admission_date Date format. Admission date variable. Typically, this should be later than the date of onset, but there are times when the disease is diagnosed in hospital.
#' @param discharge_date Date format. Discharge date variable. This date should be later than the date of admission.
#' @param icu_date Date format. ICU admission date preferably. Typically, this should be later than the date of onset and admission, but there are times when the disease is diagnosed in ICU.
#' @param genomics Character format. Genomics variable. Can be variant of SARS-CoV-2, or similarly the Hepatitis A.
#' @param icd_code Character format. ICD code variable for the admission. No pre-specified format required.
#' @param diagnosis_description Character format. Written description of the ICD code. For ease of understanding what the ICD codes mean, not a critical variable.
#' @param drg Character format. Diagnostic related group variable for the admission. No pre-specified format required.
#' @param dod Date format. Variable representing date of death. Must only have one date of death chosen (in diagnosis dataset or hospitalization dataset, not both). If dod selected is from the hospitalization dataset, it will be deleted for persons without an admission.
#' @param died Variable representing death, best use 0 and 1.
#' @param hospital Hospital identifier. Typically name of the hospital.
#' @param icu_hours ICU hours. Hours spent in ICU. Should be numeric.
#' @param dialysis Dialysis indicator (0/1).
#' @param vax_type Character format. Variable that indicates the vaccine type, brand, or antigen
#' @param vax_date Date format. Variable that indicates the vaccination event date. Make sure is in date format, and arranged in order of dates you would like it to appear when it goes to wide format. For example, if it is not in order, \code{vax_date_1} (an output variable) may be the latest vaccination date, instead of the first.
#' @param lag Numeric format. Number of days to add to the vaccination event date. Useful to define when a person reaches peak immunity post-vaccination. For COVID-19 this is often thought to be 14 days. Default lag is zero days.
#' @return The output is a dataframe that is cleaned and could be ready for machine learning data-linkage.
#' @export
#' @examples
#' # Basic usage of clean_the_nest.
#' # Use this to set up for datalinkage using the murmuration command and then cleaning with preening
#' data(dx_data)
#' df_diag <- clean_the_nest(dx_data, drop_eggs=TRUE, data_type = "cases",
#'   id_var ="identity",
#'   diagnosis = "disease_name",
#'   lettername1 = "first_name",
#'   lettername2 = "surname",
#'   dob = "date_of_birth",
#'   medicare = "medicare_no",
#'   gender = "gender",
#'   postcode="postcode",
#'   fn="indigenous_status",
#'   onset_date = "diagnosis_date")
#'
#' data(hosp_data)
#' df_hosp <- clean_the_nest(hosp_data, drop_eggs=TRUE,
#'   data_type = "hospital",
#'   id_var ="patient_id",
#'   lettername1 = "firstname",
#'   lettername2 = "last_name",
#'   dob = "birth_date",
#'   medicare = "medicare_number",
#'   gender = "sex",
#'   postcode="zip_codes",
#'   fn="cultural_heritage",
#'   icd_code = "icd_codes",
#'   admission_date = "date_of_admission",
#'   discharge_date = "date_of_discharge")
#'
#' data(vax_data)
#' df_vax <- clean_the_nest(data = vax_data,
#'   data_type = "vaccination",
#'   lie_nest_flat=TRUE,
#'   id_var = "patient_id",
#'   lettername1="firstname",
#'   lettername2="last_name",
#'   dob="birth_date",
#'   medicare="medicare_number",
#'   gender = "gender",
#'   postcode = "postcode",
#'   vax_type = "vaccine_delivered",
#'   vax_date = "service_date")


clean_the_nest <- function(data,
                           id_var=NULL,
                           event_id_var=NULL,
                           drop_eggs=FALSE,
                           data_type = NULL,
                           lie_nest_flat=FALSE,
                           drop_the_na_vax=TRUE,
                           keep_vars=NULL,
                           diagnosis=NULL,
                           lettername1=NULL,
                           lettername2=NULL,
                           dob=NULL,
                           age=NULL,
                           medicare=NULL,
                           postcode=NULL,
                           gender=NULL,
                           fn=NULL,
                           latitude=NULL,
                           longitude=NULL,
                           onset_date=NULL,
                           vax_type=NULL,
                           vax_date=NULL,
                           lag=0,
                           admission_date=NULL,
                           discharge_date=NULL,
                           hospital=NULL,
                           icd_code=NULL,
                           diagnosis_description=NULL,
                           drg=NULL,
                           icu_date=NULL,
                           icu_hours=NULL,
                           dialysis=NULL,
                           genomics=NULL,
                           dod=NULL,
                           died=NULL) {

  col_names <- c("id_var", "event_id_var", "diagnosis", "lettername1", "lettername2", "dob", "age", "medicare", "postcode", "gender", "fn",
                 "latitude", "longitude", "onset_date", "vax_code", "vax_type", "vax_disease", "vax_date", "admission_date", "discharge_date", "hospital",
                 "icd_code", "diagnosis_description", "drg", "icu_date", "icu_hours", "icu_days", "dialysis", "genomics", "dod", "died")


  if (is.null(data_type)){
    stop("No data type selected. Please specify data_type as 'cases', 'hospital', or 'vaccination'.")
  }

  if (!is.null(id_var)){
    data <- data_rename(data, {{ id_var }}, "id_var")

    if (data_type=="cases"){
      if (!all(duplicated(data$id_var)==FALSE)){
        stop("Your id_var is not unique, and there are some duplicates. Use groupby and sequence to select out your first case only (This function is currently working for first cases only at this stage.).")
      }
    }
  }







  if (!is.null(event_id_var)){
    data <- data_rename(data, {{ event_id_var }}, "event_id_var")

    if (data_type=="hospital" || data_type=="vaccination"){
      if (!all(duplicated(data$event_id_var)==FALSE)){
        stop("Your event_id_var is not unique, and there are some duplicates. All events must be unique.")
      }
    }
  }


  if (!is.null(diagnosis)){
    data <- data_rename(data, {{ diagnosis }}, "diagnosis")
  }

  if (is.null(lettername1)){
    warning("No first name variable selected. First name, last name, date of birth, and gender are recommended for data linkage.")
  } else {
    data <- data_rename(data, {{ lettername1 }}, "lettername1")
    data <- data %>% dplyr::mutate(lettername1 = tolower(lettername1), #to only take the first name, and not middle name
                                   lettername1 = str_replace_all(lettername1, "[^[:alnum:]]", ""))
  }

  if (is.null(lettername2)){
    warning("No last name variable selected. First name, last name, date of birth, and gender are recommended for data linkage.")
  } else {
    data <- data_rename(data, {{ lettername2 }}, "lettername2")
    data <- data %>% dplyr::mutate(lettername2 = tolower(lettername2),
                                   lettername2 = str_replace_all(lettername2, "[^[:alnum:]]", ""))
  }

  if (is.null(dob)){
    warning("No date of birth variable selected. First name, last name, date of birth, and gender are recommended for data linkage.")
  } else {
    data <- data_rename(data, {{ dob }}, "dob")
    if (!inherits(data$dob, "Date")) {
      stop("Date of birth variable is not in Date format.")
    }
    data <- data %>% dplyr::mutate(dob_yr = lubridate::year(dob),
                                   dob_mo = lubridate::month(dob),
                                   dob_day = lubridate::day(dob))
  }

  if (!is.null(lettername1) & !is.null(lettername2)){
    data <- data %>% dplyr::mutate(lettername1_lettername2 = paste(lettername1, lettername2))
  }

  if (!is.null(lettername1) & !is.null(lettername2) & !is.null(dob)){
    data <- data %>% dplyr::mutate(lettername1_lettername2_dob = paste(lettername1, lettername2, dob))
  }


  if (!is.null(postcode)){
    data <- data_rename(data, {{ postcode }}, "postcode")
  }


  if (!is.null(medicare)){
    data <- data_rename(data, {{ medicare }}, "medicare")
    data <- data %>% dplyr::mutate(medicare09 = substr(medicare, 1, 09),
                                   medicare10 = substr(medicare, 1, 10),
                                   medicare11 = str_replace_all(medicare, " ", ""))
  }

  if (!is.null(gender)){
    data <- data_rename(data, {{ gender }}, "gender")
  }

  if (!is.null(fn)){
    data <- data_rename(data, {{ fn }}, "fn")
  }

  if (!is.null(onset_date)) {
    date_col_name <- rlang::as_name(rlang::enquo(onset_date))

    if (!inherits(data[[date_col_name]], "Date")) {
      stop("Date of onset variable is not in proper Date format (not a Date class).")
    }

    if (inherits(data[[date_col_name]], "tiny_labelled")) {
      warning("onset_date has tiny_labelled class in addition to Date class")
    }

    data <- data_rename(data, {{ onset_date }}, "onset_date")
  }

  if (is.null(age) & !is.null(dob) & !is.null(onset_date)){
    data <- data %>% dplyr::mutate(age = as.numeric(round((onset_date - dob)/365.25, digits = 1), na.rm = TRUE))
  }

  if (!is.null({{ age }})){
    data <- data_rename(data, {{ age }} , "age")
  }


  if (!is.null(admission_date)){
    date_col_name <- rlang::as_name(rlang::enquo(admission_date))

    if (!inherits(data[[date_col_name]], "Date")) {
      stop("Admission date variable is not in proper Date format (not a Date class).")
    }

    if (inherits(data[[date_col_name]], "tiny_labelled")) {
      warning("Admission date has tiny_labelled class in addition to Date class")
    }

    data <- data_rename(data, {{ admission_date }}, "admission_date")

    if (any(is.na(data[[date_col_name]]))) {
      data <- data %>%
        dplyr::mutate(admission_outcome = factor(
          case_when(
            is.na(!!sym(date_col_name)) ~ 0,
            TRUE ~ 1
          ),
          levels = c(0, 1),
          labels = c("No Admission", "Admission")
        ))
      message("Admission outcome created: missing admission_date = 'No Admission', non-missing = 'Admission'")
    } else {
      data <- data %>%
        dplyr::mutate(admission_outcome = factor(1, levels = c(0, 1), labels = c("No Admission", "Admission")))
    }
  }

  if (!is.null(discharge_date)){
    data <- data_rename(data, {{ discharge_date }}, "discharge_date")
  }

  if (!is.null(admission_date) & !is.null(discharge_date)){
    data <- data %>% dplyr::mutate(los = as.numeric(discharge_date - admission_date))
  }

  if (!is.null(hospital)){
    data <- data_rename(data, {{ hospital }}, "hospital")
  }

  if (!is.null(icd_code)){
    data <- data_rename(data, {{ icd_code }}, "icd_code")
  }

  if (!is.null(diagnosis_description)){
    data <- data_rename(data, {{ diagnosis_description }}, "diagnosis_description")
  }

  if (!is.null(drg)){
    data <- data_rename(data, {{ drg }}, "drg")
  }

  if (!is.null(icu_hours)){
    if (!is.numeric(data[[icu_hours]])) {
      warning("ICU Hours is not numeric. Please check prior to linkage.")
    }
    data <- data %>%
      dplyr::rename(icu_hours = {{ icu_hours }}) %>%
      dplyr::mutate(icu_days = round(icu_hours / 24, 1))
  }



  if (!is.null(icu_date)){
    data <- data_rename(data, {{ icu_date }}, "icu_date")
    data <- data %>% dplyr::mutate(icu_outcome = factor(case_when(is.na(icu_date) ~ 0,
                                                                  TRUE == "Yes" ~ 1), levels = c("0", "1"),
                                                        labels = c("No ICU Admission", "ICU Admission"), ordered=TRUE))
  }

  if (is.null(icu_date) & !is.null(icu_hours)){
    data <- data_rename(data, {{ icu_hours }}, "icu_hours")
    data <- data %>% dplyr::mutate(icu_outcome = factor(case_when(icu_hours>0 ~ 1,
                                                                  TRUE ~ 0), levels = c("0", "1"),
                                                        labels = c("No ICU Admission", "ICU Admission"), ordered=TRUE))
  }


  if (!is.null(dialysis)){
    data <- data_rename(data, {{ dialysis }}, "dialysis")
  }

  if (!is.null(genomics)){
    data <- data_rename(data, {{ genomics }}, "genomics")
  }

  if (!is.null(dod)){
    data <- data_rename(data, {{ dod }}, "dod")
    if (!inherits(data$dod, "Date")) {
      stop("Date of death variable is not in Date format.")
    }

    if (!is.null(dob)){
      data <- data %>% mutate(age_at_death = as.numeric(round((dod - dob)/365.25, 1)))
    }
  }

  if (is.null(died) & !is.null(dod)){
    data <- data %>% dplyr::mutate(death_outcome = factor(case_when(is.na(dod) ~ 0,
                                                                    !is.na(dod) ~ 1), levels = c("0", "1"), labels = c("Alive", "Died")))
  }

  if (!is.null(data$gender) & !is.null(data$postcode) & !is.null(data$dob_yr)){
    data <- data %>% dplyr::mutate(block1 = paste(gender, postcode, dob_yr),
                                   block2 = paste(postcode, dob_yr),
                                   block3 = paste(gender, dob_yr))
  } else if (!is.null(data$dob_yr) & !is.null(data$postcode) ) {
    data <- data %>% dplyr::mutate(block2 = paste(postcode, dob_yr))
  } else if (!is.null(data$dob_yr) & !is.null(data$gender) ) {
    data <- data %>% dplyr::mutate(block3 = paste(gender, dob_yr))
  }




  if (!is.null(vax_type)){
    data <- data_rename(data, {{ vax_type }}, "vax_type")
  }


  if (!is.null(vax_date)){
    data <- data_rename(data, {{ vax_date }}, "vax_date")
    data <- data %>% mutate(vax_date = vax_date + lag)
  }

  if (lie_nest_flat == TRUE) {

    if (data_type == "hospital") {
      stop("You need to specify that this is vaccination data, not hospital data.")
    }


    if (is.null(id_var) | is.null(vax_type) | is.null(vax_date)) {
      stop("You need id_var, vax_type and vax_date to be specified. Do not need to have more than one vaccine name, and names can be generic.")
    }

    if (is.null(event_id_var)){
      data <- data %>% dplyr::mutate(event_id_var = 1:n())
    }

    data <- data %>%
      group_by(id_var) %>%
      dplyr::arrange(vax_date, .by_group = TRUE) %>%
      dplyr::mutate(seq = 1:n(),
                    total_vax = max(seq),
                    first_vax_date = min(vax_date),
                    last_vax_date = max(vax_date),
                    first_vax_type = case_when(seq == 1 ~ vax_type),
                    last_vax_type = case_when(seq == total_vax ~ vax_type)) %>%
      fill(last_vax_type, .direction = "up")

    vax_df <- data %>%
      pivot_wider(
        id_cols = c("id_var"),
        names_from = c("seq"),
        values_from = c("vax_date", "vax_type"),
        values_fill = NA) %>% clean_names()

    data <- data %>% filter(row_number() == 1) %>% select(-starts_with("vax_"))
    data <- left_join(data, vax_df, by="id_var") %>% ungroup() %>% select(-seq, -event_id_var)
  }


  if (drop_eggs==TRUE) {
    data <- data %>%
      dplyr::select(starts_with(c("letter","age", "medicare", "onset", "block")), any_of(col_names), starts_with('vax_'), any_of("los"), contains("admiss"), any_of(keep_vars))
  }

  return(data)
}
