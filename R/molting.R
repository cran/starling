#' @name molting
#' @title Molt: De-identify a Dataset with Hash-based Relinking
#'
#' @description
#' Like a bird molting its feathers for new plumage, this function removes
#' identifiable information and replaces it with a unique hash for each row.
#' It returns both the de-identified dataset and a lookup table for relinking.
#' Age category variables (age2cat, age3cat, etc.) are automatically retained.
#'
#' @param data A data frame to be de-identified.
#' @param id_cols An optional character vector of column names to use for
#'   creating the hash. If NULL (the default), the function will use the
#'   PII columns it automatically detects.
#' @param pii_patterns An optional character vector of regular expression
#'   patterns used to detect PII columns for removal. The default list
#'   includes common identifiers.
#' @param additional_pii_cols An optional character vector of specific column
#'   names to remove as PII, in addition to those detected by pattern matching.
#'   Useful for adding dataset-specific identifiers without modifying patterns.
#' @param hash_method The hashing algorithm to use. Options include "sha256"
#'   (default), "md5", "sha1", "sha512", "crc32", "xxhash32", "xxhash64",
#'   "murmur32", "spookyhash", or "blake3". See ?digest::digest for details.
#' @param hash_col_name A string for the name of the new hash column.
#'   Defaults to "row_hash".
#' @param return_lookup Logical. If TRUE (default), returns a list containing
#'   both the de-identified data and a lookup table. If FALSE, returns only
#'   the de-identified data frame.
#' @param seed An optional integer seed for reproducible hashing with certain
#'   algorithms. Defaults to NULL.
#'
#' @return If return_lookup = TRUE (default), a list with two elements:
#'   \itemize{
#'     \item \code{deidentified}: The de-identified data frame with hash column
#'     \item \code{lookup}: A data frame containing only the identifier columns
#'           and the hash for relinking
#'   }
#'   If return_lookup = FALSE, returns only the de-identified data frame.
#'
#' @importFrom digest digest
#' @importFrom dplyr mutate across select all_of rowwise ungroup c_across relocate n_distinct
#' @importFrom stringr str_subset regex
#' @importFrom rlang :=
#'
#' @details
#' The function identifies PII columns based on pattern matching, creates a
#' unique hash for each row based on the concatenated identifier values, and
#' returns both a de-identified dataset and a secure lookup table.
#'
#' Age category variables (variables matching the pattern "age\\d+cat" such as
#' age2cat, age5cat, age10cat, etc.) are automatically retained in the
#' de-identified dataset as they are not considered directly identifying.
#'
#' Security Note: The lookup table contains sensitive information and should be
#' stored securely with appropriate access controls. Consider encrypting this
#' file if storing to disk.
#'
#' @examples
#' # Create sample data
#' patient_data <- data.frame(
#'   patient_name = c("John Doe", "Jane Smith"),
#'   dob = as.Date(c("1980-01-01", "1975-05-15")),
#'   mrn = c("12345", "67890"),
#'   age5cat = factor(c("18-64", "18-64")),
#'   diagnosis = c("Condition A", "Condition B"),
#'   lab_value = c(120, 95)
#' )
#'
#' # Basic de-identification (age categories automatically retained)
#' result <- suppressMessages(molting(patient_data))
#' names(result$deidentified)  # Check column names
#' head(result$deidentified, 2)  # View de-identified data
#'
#' # Use different hash method
#' result_md5 <- suppressMessages(
#'   molting(patient_data, hash_method = "md5")
#' )
#'
#' # Return only de-identified data (no lookup table)
#' deidentified_only <- suppressMessages(
#'   molting(patient_data, return_lookup = FALSE)
#' )
#'
#' # Add specific columns to PII removal
#' result_custom <- suppressMessages(
#'   molting(patient_data, additional_pii_cols = c("study_id"))
#' )
#'
#' # Specify custom identifier columns for hashing
#' result_ids <- suppressMessages(
#'   molting(patient_data, id_cols = c("mrn", "dob"))
#' )
#'
#' @export
molting <- function(data,
                    id_cols = NULL,
                    pii_patterns = NULL,
                    additional_pii_cols = NULL,
                    hash_method = "sha256",
                    hash_col_name = "row_hash",
                    return_lookup = TRUE,
                    seed = NULL) {

  # Check for required packages
  required_pkgs <- c("dplyr", "stringr", "digest")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed. ",
           "Please install it with: install.packages('", pkg, "')",
           call. = FALSE)
    }
  }

  # Validate input
  if (!inherits(data, "data.frame")) {
    stop("Input 'data' must be a data frame", call. = FALSE)
  }

  # Validate hash_method
  valid_methods <- c("md5", "sha1", "crc32", "sha256", "sha512",
                     "xxhash32", "xxhash64", "murmur32", "spookyhash", "blake3")
  if (!hash_method %in% valid_methods) {
    stop("Invalid hash_method. Must be one of: ",
         paste(valid_methods, collapse = ", "),
         call. = FALSE)
  }

  # Define default PII patterns if not provided
  if (is.null(pii_patterns)) {
    pii_patterns <- c(
      "name", "address", "street",
      "phone", "email", "contact", "dob", "birth", "ssn",
      "mrn", "urn", "medicare", "patient_id", "subject_id",
      "participant", "surname", "firstname", "lastname",
      "_id$"  # Matches any column ending in "_id"
    )
  }

  # Find PII columns to remove based on pattern matching
  pii_search_pattern <- paste(pii_patterns, collapse = "|")
  cols_to_remove <- stringr::str_subset(
    names(data),
    stringr::regex(pii_search_pattern, ignore_case = TRUE)
  )

  # Identify age category variables to preserve (age2cat, age5cat, etc.)
  age_cat_pattern <- "^age\\d+cat$"
  age_cat_vars <- grep(age_cat_pattern, cols_to_remove, value = TRUE, ignore.case = TRUE)

  # Remove age category variables from the removal list
  if (length(age_cat_vars) > 0) {
    cols_to_remove <- setdiff(cols_to_remove, age_cat_vars)
    message("Preserving age category variables: ",
            paste(age_cat_vars, collapse = ", "))
  }

  # Add any additional specific column names provided by user
  if (!is.null(additional_pii_cols)) {
    # Validate that additional columns exist in the data
    missing_additional_cols <- setdiff(additional_pii_cols, names(data))
    if (length(missing_additional_cols) > 0) {
      warning("The following columns in `additional_pii_cols` are not in the dataset and will be ignored: ",
              paste(missing_additional_cols, collapse = ", "))
    }

    # Add existing additional columns to removal list
    existing_additional_cols <- intersect(additional_pii_cols, names(data))
    cols_to_remove <- unique(c(cols_to_remove, existing_additional_cols))

    if (length(existing_additional_cols) > 0) {
      message("Adding user-specified PII columns: ",
              paste(existing_additional_cols, collapse = ", "))
    }
  }

  if (length(cols_to_remove) == 0) {
    warning("No potential PII columns were detected based on the patterns. ",
            "Returning original data unchanged.")
    if (return_lookup) {
      return(list(deidentified = data, lookup = NULL))
    } else {
      return(data)
    }
  }

  # Determine which columns to use for hashing
  if (is.null(id_cols)) {
    cols_for_hash <- cols_to_remove
    message("Using auto-detected PII columns for hashing: ",
            paste(cols_for_hash, collapse = ", "))
  } else {
    cols_for_hash <- id_cols
    missing_cols <- setdiff(cols_for_hash, names(data))
    if (length(missing_cols) > 0) {
      stop("The following specified `id_cols` are not in the dataset: ",
           paste(missing_cols, collapse = ", "),
           call. = FALSE)
    }
    message("Using specified columns for hashing: ",
            paste(cols_for_hash, collapse = ", "))
  }

  # Create hash for each row
  data_with_hash <- data %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(cols_for_hash), as.character)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      !!hash_col_name := digest::digest(
        object = paste(dplyr::c_across(dplyr::all_of(cols_for_hash)),
                       collapse = "_"),
        algo = hash_method,
        seed = seed
      )
    ) %>%
    dplyr::ungroup()

  # Create de-identified dataset
  deidentified_data <- data_with_hash %>%
    dplyr::select(-dplyr::all_of(cols_to_remove)) %>%
    dplyr::relocate(!!hash_col_name)

  # Create lookup table if requested
  if (return_lookup) {
    lookup_data <- data_with_hash %>%
      dplyr::select(!!hash_col_name, dplyr::all_of(cols_to_remove)) %>%
      dplyr::relocate(!!hash_col_name)

    # Check for hash collisions
    n_rows <- nrow(lookup_data)
    n_unique_hashes <- dplyr::n_distinct(lookup_data[[hash_col_name]])

    if (n_rows != n_unique_hashes) {
      warning("Hash collision detected! ", n_rows, " rows resulted in ",
              n_unique_hashes, " unique hashes. ",
              "Consider using a different hash_method or adding more ",
              "identifier columns.")
    }
  }

  # Summary message
  message(
    "\n--- Molting Complete ---\n",
    "Hash algorithm: ", hash_method, "\n",
    "Rows processed: ", nrow(data), "\n",
    "Columns removed: ", length(cols_to_remove), " (",
    paste(cols_to_remove, collapse = ", "), ")\n",
    "Columns retained: ", ncol(deidentified_data) - 1, "\n",
    if (return_lookup) {
      paste0("Lookup table created: ", nrow(lookup_data), " rows\n",
             "IMPORTANT: Store the lookup table securely!")
    } else {
      "Lookup table not created (return_lookup = FALSE)"
    }
  )

  # Return results
  if (return_lookup) {
    return(list(
      deidentified = deidentified_data,
      lookup = lookup_data
    ))
  } else {
    return(deidentified_data)
  }
}
