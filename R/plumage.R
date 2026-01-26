#' @name plumage
#' @title Identify Chronic Conditions Using ICD-10-AM U-Codes
#'
#' @description
#' Analyzes a hospitalization dataset to identify chronic conditions based on ICD-10-AM U-codes.
#' Like identifying a bird by its distinctive plumage (feathers), this function identifies patients
#' by their chronic condition patterns. Creates binary indicators for each condition and calculates
#' total condition counts by category.
#'
#' @param df A data frame containing hospitalization records
#' @param icd_column Character string specifying the name of the column containing ICD-10-AM codes
#' @param prefix Optional character string to prefix all output column names (default: NULL)
#' @param decimal Logical indicating whether to match U-codes with decimal points (TRUE, default)
#'   or without decimal points (FALSE). When TRUE, matches "U78.1" format; when FALSE, matches "U781" format.
#' @param drop_eggs Logical indicating whether to drop individual condition columns and retain only
#'   summary columns. Default is FALSE.
#'
#' @details
#' This function identifies chronic conditions from ICD-10-AM U-codes (Australian modification codes
#' for chronic conditions). The function recognizes the following conditions:
#'
#' **Metabolic/Endocrine:**
#' * U78.1: Obesity
#' * U78.2: Cystic fibrosis
#'
#' **Mental Health:**
#' * U79.1: Dementia
#' * U79.2: Schizophrenia
#' * U79.3: Depression
#' * U79.4: Intellectual/developmental disability
#'
#' **Neurological:**
#' * U80.1: Parkinson's disease
#' * U80.2: Multiple sclerosis
#' * U80.3: Epilepsy
#' * U80.4: Cerebral palsy
#' * U80.5: Paralysis
#'
#' **Cardiovascular:**
#' * U82.1: Ischaemic heart disease
#' * U82.2: Heart failure
#' * U82.3: Hypertension
#'
#' **Respiratory:**
#' * U83.1: Emphysema
#' * U83.2: COPD
#' * U83.3: Asthma
#' * U83.4: Bronchiectasis
#' * U83.5: Respiratory failure
#'
#' **Gastrointestinal:**
#' * U84.1: Crohn's disease
#' * U84.2: Ulcerative colitis
#' * U84.3: Liver failure
#'
#' **Musculoskeletal:**
#' * U86.1: Rheumatoid arthritis
#' * U86.2: Osteoarthritis
#' * U86.3: Systemic lupus erythematosus
#' * U86.4: Osteoporosis
#'
#' **Renal:**
#' * U87.1: Chronic kidney disease
#'
#' **Congenital:**
#' * U88.1: Spina bifida
#' * U88.2: Down syndrome
#'
#' The function searches for these codes within the specified ICD column and creates binary
#' indicators for each condition. It also calculates summary measures including total conditions
#' overall and by disease category.
#'
#' Note: Cystic fibrosis (U78.2) is counted in both metabolic and respiratory categories.
#'
#' @return
#' Returns the input data frame with additional columns:
#' * Binary indicators (0/1) for each chronic condition, optionally prefixed (unless drop_eggs = TRUE)
#' * total_conditions: Sum of all identified conditions
#' * total_metabolic_conditions: Sum of metabolic/endocrine conditions
#' * total_mental_health_conditions: Sum of mental health conditions
#' * total_neurological_conditions: Sum of neurological conditions
#' * total_cardiovascular_conditions: Sum of cardiovascular conditions
#' * total_respiratory_conditions: Sum of respiratory conditions
#' * total_gastrointestinal_conditions: Sum of gastrointestinal conditions
#' * total_musculoskeletal_conditions: Sum of musculoskeletal conditions
#' * total_renal_conditions: Sum of renal conditions
#' * total_congenital_conditions: Sum of congenital conditions
#' * conditions_category: Factor with levels "0", "1", "2", "3+" based on total_conditions
#'
#' @examples
#' # Create sample hospitalization data
#' hospital_data <- data.frame(
#'   patient_id = 1:4,
#'   icd_codes = c(
#'     "K29.70",
#'     "U78.1, U83.2, U82.3",
#'     "U79.3, U83.3",
#'     "U80.1, U86.2"
#'   )
#' )
#'
#' # Identify chronic conditions with decimal format (default)
#' results1 <- plumage(hospital_data, "icd_codes")
#'
#' # View category summaries
#' results1[, c("patient_id", "total_conditions",
#'              "total_cardiovascular_conditions",
#'              "total_respiratory_conditions")]
#'
#' # Identify chronic conditions without decimal format
#' results2 <- plumage(hospital_data, "icd_codes", decimal = FALSE)
#'
#' # Identify chronic conditions with prefix
#' results3 <- plumage(hospital_data, "icd_codes", prefix = "chronic_")
#'
#' # Keep only summary columns, drop individual conditions
#' results4 <- plumage(hospital_data, "icd_codes", drop_eggs = TRUE)
#'
#' @importFrom dplyr mutate across select all_of sym case_when
#' @importFrom stringr str_detect
#' @export
plumage <- function(df, icd_column, prefix = NULL, decimal = TRUE, drop_eggs = FALSE) {

  # Input validation
  if (!inherits(df, "data.frame")) {
    stop("Input 'df' must be a data frame")
  }

  if (!icd_column %in% names(df)) {
    stop("Column '", icd_column, "' not found in data frame")
  }

  if (!is.logical(decimal)) {
    stop("Parameter 'decimal' must be logical (TRUE or FALSE)")
  }

  if (!is.logical(drop_eggs)) {
    stop("Parameter 'drop_eggs' must be logical (TRUE or FALSE)")
  }

  if (!is.null(prefix) && !is.character(prefix)) {
    stop("Parameter 'prefix' must be NULL or a character string")
  }

  # Create a list of all U-codes with their descriptions
  # Base codes without formatting
  base_codes <- list(
    obesity = c("U78", "1"),
    cystic_fibrosis = c("U78", "2"),
    dementia = c("U79", "1"),
    schizophrenia = c("U79", "2"),
    depression = c("U79", "3"),
    intellectual_dev = c("U79", "4"),
    parkinsons = c("U80", "1"),
    multiple_sclerosis = c("U80", "2"),
    epilepsy = c("U80", "3"),
    cerebral_palsy = c("U80", "4"),
    paralysis = c("U80", "5"),
    ihd = c("U82", "1"),
    heart_failure = c("U82", "2"),
    hypertension = c("U82", "3"),
    emphysema = c("U83", "1"),
    copd = c("U83", "2"),
    asthma = c("U83", "3"),
    bronchiectasis = c("U83", "4"),
    respiratory_failure = c("U83", "5"),
    crohns = c("U84", "1"),
    ulcerative_colitis = c("U84", "2"),
    liver_failure = c("U84", "3"),
    rheumatoid_arthritis = c("U86", "1"),
    osteoarthritis = c("U86", "2"),
    lupus = c("U86", "3"),
    osteoporosis = c("U86", "4"),
    kidney_disease = c("U87", "1"),
    spina_bifida = c("U88", "1"),
    downs = c("U88", "2")
  )

  # Format patterns based on decimal parameter
  ucode_patterns <- lapply(base_codes, function(x) {
    if (decimal) {
      paste0(x[1], "\\.", x[2])  # Format: U78\.1
    } else {
      paste0(x[1], x[2])  # Format: U781
    }
  })

  # Define disease category groupings (without prefix)
  category_groups <- list(
    metabolic = c("obesity", "cystic_fibrosis"),
    mental_health = c("dementia", "schizophrenia", "depression", "intellectual_dev"),
    neurological = c("parkinsons", "multiple_sclerosis", "epilepsy", "cerebral_palsy", "paralysis"),
    cardiovascular = c("ihd", "heart_failure", "hypertension"),
    respiratory = c("emphysema", "copd", "asthma", "bronchiectasis", "respiratory_failure", "cystic_fibrosis"),
    gastrointestinal = c("crohns", "ulcerative_colitis", "liver_failure"),
    musculoskeletal = c("rheumatoid_arthritis", "osteoarthritis", "lupus", "osteoporosis"),
    renal = c("kidney_disease"),
    congenital = c("spina_bifida", "downs")
  )

  # Add prefix to condition names if specified
  if (!is.null(prefix)) {
    names(ucode_patterns) <- paste0(prefix, names(ucode_patterns))

    # Update category groups with prefix
    category_groups <- lapply(category_groups, function(conditions) {
      paste0(prefix, conditions)
    })
  }

  # Create indicator variables for each U-code
  result <- df %>%
    mutate(
      across(
        .cols = all_of(icd_column),
        .names = "{.col}",
        .fns = identity
      )
    )

  # Add a column for each U-code
  for (condition in names(ucode_patterns)) {
    result <- result %>%
      mutate(
        !!condition := as.integer(str_detect(!!sym(icd_column), ucode_patterns[[condition]]))
      )
  }

  # Create names for summary columns (with prefix if specified)
  total_name <- if (!is.null(prefix)) paste0(prefix, "total_conditions") else "total_conditions"
  category_name <- if (!is.null(prefix)) paste0(prefix, "conditions_category") else "conditions_category"

  # Add total conditions count
  result <- result %>%
    mutate(
      !!total_name := rowSums(select(., all_of(names(ucode_patterns))))
    )

  # Add category-specific totals
  for (category in names(category_groups)) {
    category_total_name <- if (!is.null(prefix)) {
      paste0(prefix, "total_", category, "_conditions")
    } else {
      paste0("total_", category, "_conditions")
    }

    result <- result %>%
      mutate(
        !!category_total_name := rowSums(select(., all_of(category_groups[[category]])))
      )
  }

  # Add categorical variable based on total_conditions
  result <- result %>%
    mutate(
      !!category_name := factor(
        case_when(
          !!sym(total_name) == 0 ~ "0",
          !!sym(total_name) == 1 ~ "1",
          !!sym(total_name) == 2 ~ "2",
          !!sym(total_name) >= 3 ~ "3+"
        ),
        levels = c("0", "1", "2", "3+")
      )
    )

  # Drop individual condition columns if requested
  if (drop_eggs) {
    result <- result %>%
      select(-all_of(names(ucode_patterns)))
  }

  return(result)
}
