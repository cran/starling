#' starling: Probabilistic Record Linkage for Public Health Surveillance
#'
#' @description
#' **Bird note**: Starlings form murmurations — thousands of individual birds moving
#' as one coordinated, fluid shape with no central conductor, each responding only
#' to its nearest neighbours until a single coherent pattern emerges from the flock.
#' That is the visual the package is named for: \code{murmuration()} takes two
#' separate, uncoordinated sets of records and lets pairwise local comparisons
#' resolve into one coherent linked dataset. The supporting functions —
#' \code{flock()}, \code{preflight()}, \code{check_medicare()}, and
#' \code{murmuration_plot()} — are the pre-flight preparation and post-flight
#' review that make the murmuration trustworthy.
#'
#' @details
#' ## Core workflow
#'
#' \enumerate{
#'   \item \strong{\code{\link{preflight}()}} — Run a battery of pre-linkage data
#'     quality checks: completeness, duplicates, date plausibility, Medicare
#'     validity, name quality, and factor-level consistency across both datasets.
#'   \item \strong{\code{\link{check_medicare}()}} — Validate Medicare numbers using
#'     the Services Australia Modulus 10 weighted checksum (can also be called
#'     standalone or from \code{preflight()}).
#'   \item \strong{\code{\link{flock}()}} — Generate blocking variables (single-field,
#'     phonetic, or composite) to partition records before linkage.
#'   \item \strong{\code{\link{murmuration}()}} — Link a cohort to a dated event
#'     stream (\code{linkage_type = "cohort2event"} with an \code{event_type} of
#'     \code{"vaccination"}, \code{"hospitalisation"}, or \code{"case"}), using
#'     either \code{method = "probabilistic"} (Fellegi-Sunter EM via \pkg{reclin2})
#'     or \code{method = "deterministic"} (exact-key matching, no scoring). Of the
#'     historical four-code taxonomy, case-to-hospitalisation (\code{"c2h"}) and
#'     vaccination-to-case (\code{"v2c"}) still work but are deprecated in favour of
#'     \code{cohort2event}; vaccination-to-hospitalisation (\code{"v2h"}) and
#'     vaccination-to-event (\code{"v2e"}) are defunct as of 2.0.0 (use
#'     \code{cohort2event} with \code{event_type = "vaccination"}, and a fixed
#'     \code{event_date} for an event linelist).
#'   \item \strong{\code{\link{murmuration_plot}()}} — Plot the linkage weight
#'     distribution to confirm the threshold is placed at the natural valley between
#'     match and non-match clusters before finalising. Only applicable to
#'     \code{method = "probabilistic"} output — deterministic linkage has no score
#'     to plot.
#' }
#'
#' ## Threshold guidance
#'
#' See \code{\link{murmuration_plot}} for a discussion of threshold selection.
#' \code{murmuration()}'s default \code{threshold_value} is \strong{16}, calibrated
#' for five comparison variables (first name, surname, date of birth, Medicare number,
#' and postcode or gender) with an agreeing Medicare number. With fewer variables, or
#' Medicare numbers that are missing or from an older card, true matches score lower:
#' always plot the distribution first, and see
#' \code{vignette("threshold-selection", package = "starling")}.
#'
#' @keywords internal
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "-- starling ", utils::packageVersion("starling"), " ",
    paste(rep("-", 40), collapse = ""), "\n",
    "Probabilistic record linkage for public health surveillance.\n",
    "  murmuration()      link cohort to event stream (cohort2event; legacy c2h/v2c)\n",
    "  flock()            generate blocking variables\n",
    "  check_medicare()   validate Medicare checksums\n",
    "  preflight()        pre-linkage data quality audit\n",
    "  murmuration_plot() inspect weight distribution + threshold\n",
    "  perch()            sweep candidate thresholds (AIHW/PHRN benchmarks)\n",
    paste(rep("-", 55), collapse = "")
  )
}
