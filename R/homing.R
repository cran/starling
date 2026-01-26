#' Homing: Relink De-identified Data Using Lookup Table
#'
#' Like homing pigeons finding their way back, this function relinks
#' de-identified data with original identifiers using the lookup table
#' created by molting().
#'
#' @param deidentified_data A de-identified data frame containing a hash column
#'   (typically the output from molting()$deidentified).
#' @param lookup_table The lookup table data frame that maps anonymous hash values
#'   back to original identifiers. Created by molting(), it contains the hash column
#'   plus all removed identifier columns (names, dates of birth, medical record
#'   numbers, etc.). This serves as the secure "key" for relinking de-identified
#'   data back to real identities. Each row maps one hash to one set of identifiers.
#'   Typically obtained as molting()$lookup. Example structure for a dataset that
#'   had patient_name, dob, and mrn removed: row_hash | patient_name | dob | mrn.
#' @param hash_col_name The name of the hash column used for linking. Must exist
#'   in both deidentified_data and lookup_table. Defaults to "row_hash".
#' @param keep_hash Logical. If TRUE (default), keeps the hash column in the
#'   relinked data. If FALSE, removes it after relinking.
#'
#' @return A data frame with the original identifiers merged back in.
#'
#' @examples
#' \dontrun{
#' # First, de-identify the data
#' result <- molting(patient_data)
#'
#' # Later, relink when needed
#' original_data <- homing(
#'   result$deidentified,
#'   result$lookup
#' )
#' }
#'
#' @export
homing <- function(deidentified_data,
                   lookup_table,
                   hash_col_name = "row_hash",
                   keep_hash = TRUE) {

  # Check for required package
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.", call. = FALSE)
  }

  # Validate inputs
  if (!hash_col_name %in% names(deidentified_data)) {
    stop(paste0("Hash column '", hash_col_name,
                "' not found in deidentified_data."), call. = FALSE)
  }

  if (!hash_col_name %in% names(lookup_table)) {
    stop(paste0("Hash column '", hash_col_name,
                "' not found in lookup_table."), call. = FALSE)
  }

  # Perform the join
  relinked <- deidentified_data %>%
    dplyr::left_join(lookup_table, by = hash_col_name)

  # Check for unmatched records
  unmatched <- sum(is.na(relinked[[names(lookup_table)[2]]]))
  if (unmatched > 0) {
    warning(paste0(unmatched, " row(s) in deidentified_data could not be ",
                   "matched to the lookup table."))
  }

  # Remove hash column if requested
  if (!keep_hash) {
    relinked <- relinked %>%
      dplyr::select(-!!hash_col_name)
  }

  message(paste0("Homing complete. ", nrow(relinked), " rows processed, ",
                 nrow(relinked) - unmatched, " successfully relinked."))

  return(relinked)
}
