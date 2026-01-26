#' @name tweet
#' @title Examine and summarize variables in a dataset
#' @description Provides a comprehensive summary of variables in a dataset after cleaning with clean_the_nest. This function
#' examines variables by type, providing appropriate statistics for numeric, date, factor, and character variables.
#' For numeric variables, it shows min/max values, quartiles and missing data counts. For date variables, it displays
#' the date range and percentage of non-missing values. For factor and character variables, it shows the number of
#' unique levels, frequency of top levels, and missing data counts.
#'
#' @param data The dataset, typically output from clean_the_nest function
#' @param select_vars Optional vector of variable names to examine. If NULL, all variables will be summarized.
#' @param top_n Number of top categories to display for factor and character variables. Default is 3.
#' @param sort_by How to sort variables in the output. Options are "name" (alphabetical) or "type" (grouped by data type). Default is "type".
#' @return A data frame with one row per variable, containing variable name, type, missingness, and type-specific statistics.
#' @export
#' @examples
#' # basic usage of tweet after clean_the_nest
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
#' # Examine all variables in the cleaned dataset
#' summary_df <- tweet(df_diag)
#'
#' # Examine only specific variables
#' summary_df_subset <- tweet(df_diag, select_vars = c("age", "gender", "onset_date"))
#'
#' # Show more categories for factor variables
#' summary_df_detailed <- tweet(df_diag, top_n = 5)

tweet <- function(data,
                  select_vars = NULL,
                  top_n = 3,
                  sort_by = "type") {

  # Validate inputs
  if (!inherits(data, "data.frame")) {
    stop("Input must be a data frame or tibble")
  }

  if (!is.null(select_vars)) {
    missing_vars <- select_vars[!select_vars %in% names(data)]
    if (length(missing_vars) > 0) {
      warning("The following variables were not found in the dataset: ",
              paste(missing_vars, collapse = ", "))
    }
    var_names <- intersect(select_vars, names(data))
  } else {
    var_names <- names(data)
  }

  if (length(var_names) == 0) {
    stop("No valid variables to summarize")
  }

  # Check sort_by parameter
  if (!sort_by %in% c("name", "type")) {
    warning("sort_by must be 'name' or 'type'. Using default 'type'")
    sort_by <- "type"
  }

  # Function to classify variable types
  get_var_type <- function(var) {
    if (inherits(var, "Date")) {
      return("date")
    } else if (inherits(var, "factor") || inherits(var, "ordered")) {
      return("factor")
    } else if (is.numeric(var)) {
      return("numeric")
    } else if (is.logical(var)) {
      return("logical")
    } else {
      return("character")
    }
  }

  # Initialize result list
  result_list <- list()

  # Process each variable
  for (var_name in var_names) {
    var_data <- data[[var_name]]
    var_type <- get_var_type(var_data)

    # Calculate missingness
    n_missing <- sum(is.na(var_data))
    pct_missing <- round(100 * n_missing / nrow(data), 1)

    # Base row with common info - initialize ALL possible columns with NA
    row_data <- list(
      variable = var_name,
      type = var_type,
      n_total = nrow(data),
      n_missing = n_missing,
      pct_missing = pct_missing,
      # Numeric columns
      min = NA_real_,
      max = NA_real_,
      mean = NA_real_,
      median = NA_real_,
      q1 = NA_real_,
      q3 = NA_real_,
      # Date columns
      date_min = NA_character_,
      date_max = NA_character_,
      date_range_days = NA_integer_,
      n_valid = NA_integer_,
      pct_valid = NA_real_,
      # Categorical columns
      n_unique = NA_integer_,
      top_values = NA_character_
    )

    # Type-specific summaries
    if (var_type == "numeric") {
      if (!all(is.na(var_data))) {
        row_data$min <- min(var_data, na.rm = TRUE)
        row_data$max <- max(var_data, na.rm = TRUE)
        row_data$mean <- mean(var_data, na.rm = TRUE)
        row_data$median <- median(var_data, na.rm = TRUE)
        row_data$q1 <- quantile(var_data, 0.25, na.rm = TRUE, names = FALSE)
        row_data$q3 <- quantile(var_data, 0.75, na.rm = TRUE, names = FALSE)
      }
    }
    else if (var_type == "date") {
      if (!all(is.na(var_data))) {
        valid_dates <- sum(!is.na(var_data))
        pct_valid <- round(100 * valid_dates / nrow(data), 1)

        date_min <- min(var_data, na.rm = TRUE)
        date_max <- max(var_data, na.rm = TRUE)
        date_range_days <- as.integer(date_max - date_min)

        row_data$date_min <- as.character(date_min)
        row_data$date_max <- as.character(date_max)
        row_data$date_range_days <- date_range_days
        row_data$n_valid <- valid_dates
        row_data$pct_valid <- pct_valid
      } else {
        row_data$n_valid <- 0L
        row_data$pct_valid <- 0
      }
    }
    else if (var_type %in% c("factor", "character", "logical")) {
      if (!all(is.na(var_data))) {
        # Convert to character for consistent handling
        if (var_type == "factor") {
          char_data <- as.character(var_data)
        } else {
          char_data <- var_data
        }

        # Count unique values
        unique_values <- na.omit(unique(char_data))
        n_unique <- length(unique_values)

        # Get value counts
        value_counts <- table(char_data, useNA = "no")
        sorted_counts <- sort(value_counts, decreasing = TRUE)

        # Create top values string
        display_n <- min(top_n, n_unique)

        if (n_unique == 0) {
          top_values_str <- NA_character_
        } else if (n_unique <= 10) {
          # For few categories, show all with counts
          top_values_str <- paste(
            names(sorted_counts), " (", sorted_counts, ")",
            sep = "", collapse = "; "
          )
        } else {
          # For many categories, show top N
          top_values <- names(sorted_counts)[1:display_n]
          top_counts <- sorted_counts[1:display_n]
          top_values_str <- paste(
            top_values, " (", top_counts, ")",
            sep = "", collapse = "; "
          )
        }

        row_data$n_unique <- n_unique
        row_data$top_values <- top_values_str
      } else {
        row_data$n_unique <- 0L
      }
    }

    result_list[[var_name]] <- row_data
  }

  # Convert list to data frame using bind_rows (handles different columns)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for tweet(). Please install it.",
         call. = FALSE)
  }

  result_df <- dplyr::bind_rows(result_list)

  # Determine variable types for sorting
  var_types <- sapply(data[var_names], get_var_type)
  result_df$type <- var_types[result_df$variable]

  # Sort if requested
  if (sort_by == "type") {
    result_df <- result_df[order(result_df$type, result_df$variable), ]
  } else {
    result_df <- result_df[order(result_df$variable), ]
  }

  # Reset row names after sorting
  rownames(result_df) <- NULL

  return(result_df)
}
