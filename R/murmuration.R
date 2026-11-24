#' @name murmuration
#' @title Links case, hospital, outbreak, event linelists, vaccination datasets

# Suppress R CMD check notes for unquoted column names used in dplyr pipelines
utils::globalVariables(c(
  "weights", "threshold", ".data",
  "admission_date", "total_admissions", "first_admission_date",
  "dates", "types", "total_valid_vax", "vaccination_status", "last_vax_date",
  "x_pos", "match_grp", "n_above", "n_clerical", "link_rate",
  "pct_above", "pct_clerical", ".starling_det_matched"
))
#' @description
#' **Bird note**: A murmuration is one of nature's most arresting phenomena: tens
#' of thousands of starlings moving as a single fluid shape across the sky, with
#' no conductor and no central instruction — each bird responding only to its
#' nearest neighbours until a coherent, shifting pattern emerges from the whole
#' flock. \code{murmuration()} works the same way: it takes two separate,
#' uncoordinated datasets and lets pairwise local comparisons between records
#' resolve, through the Fellegi-Sunter EM algorithm, into one coherent linked
#' result. No single record "knows" about the others — the shape emerges from
#' the scoring.
#'
#' Probabilistic record linkage of two surveillance datasets using the
#' Fellegi-Sunter framework. Links diagnostic case linelists, hospital admission
#' records, outbreak linelists, and event manifests to vaccination history
#' (Australian Immunisation Register) or to each other. Wraps the full
#' \pkg{reclin2} pipeline — blocking, comparison, EM scoring, threshold
#' selection, and post-linkage window filtering — in a single practitioner-
#' facing call.
#'
#' @details
#' Make sure that you do not have the same variables (other than linkage
#' variables e.g. letternames, DOB, gender) in both datasets. Always make sure
#' your date columns are properly formatted using \code{as.Date()}. For example,
#' if both datasets have a date of death, choose the dataset with the highest
#' confidence and drop the date of death from the other before linking. If
#' linking to a hospitalisation dataset, the difference between
#' \code{event_date} and \code{admission_date} is used to identify
#' disease-related admissions; unrelated hospitalisations can be filtered
#' separately using ICD-10 codes or AR-DRG codes prior to linkage.
#'
#' The recommended pre-linkage workflow is:
#' \enumerate{
#'  \item \code{mudnester::clean_the_nest()} to clean and prep data for linkage. Pay close attention to your linkage variables (letternames, date of birth, medicare number, gender and/or postcode), and ensure all dates are formatted as dates.
#'  \item \code{linkage_quality()} to run a structured battery of pre-linkage quality checks across both datasets (completeness, duplicates, Medicare validity, date plausibility, gender and DOB year distribution).
#'  \item \code{\link{check_medicare}} to validate Medicare check digits and flag invalid numbers before they enter comparison scoring.
#'  \item \code{flock()} to create up to three blocking variables at different specificity levels (coarse, medium, fine) — pass one to the \code{blocking_var} argument below.
#'  \item \code{\link{murmuration}} with \code{linkage_type = "cohort2event", event_type = "vaccination"} to link cases (or any cohort -- hospitalisations for a test-negative design, a birth cohort) to vaccination data, anchored on a date column of data_x.
#'  \item \code{\link{murmuration}} with \code{linkage_type = "cohort2event", event_type = "hospitalisation"} to link cases to hospital admissions.
#'  \item \code{\link{murmuration}} with \code{linkage_type = "cohort2event", event_type = "vaccination"} and a single \code{Date} as \code{event_date} to link event linelists (flight manifests, outbreak investigations) to vaccination history as at that date.
#'  \item \code{\link{murmuration_plot}}, called directly on this function's own return value (pass \code{keep_pairs_pred = TRUE} to the \code{\link{murmuration}} call above first — default \code{FALSE} skips building this table), to visualise the full linkage score distribution (jittered by \code{blocking_var}) and confirm the threshold before accepting the linked dataset. No need to re-derive the pairs table via \pkg{reclin2} yourself.
#'  \item \code{mudnester::preening()} to prettify the dataframe prepping it for exploration, analysis and presentation. Great to use with \code{gtsummary::tbl_summary()}.
#' }
#'
#' ## On threshold selection
#' The \code{threshold_value} is a Fellegi-Sunter log-likelihood ratio score and is
#' \strong{dataset-specific} — there is no universal correct value. The default of 16 is
#' calibrated against SCPHU datasets with 5 comparison variables (first name, last name,
#' date of birth, Medicare number and gender or postcode) where the Medicare number is
#' present and current on both sides. Use \code{\link{murmuration_plot}} (after calling
#' this function with \code{keep_pairs_pred = TRUE}) to inspect the bimodal score
#' distribution from your specific dataset and identify the natural valley between the
#' match and non-match peaks. Typical working ranges: conservative (high specificity)
#' about 20-25; balanced (default) about 14-18; sensitive (high recall)
#' about 8-13.
#'
#' \strong{Reminder}: 16 is calibrated for 5 \code{compare_vars} with an agreeing Medicare
#' number -- it is not a safe default otherwise. Fewer variables (e.g. AIR-style
#' vaccination extracts with no address field, often down to 3) structurally cap the
#' achievable weight below what 16 assumes, which reads as "zero links" even when the
#' underlying data is fine. A Medicare number that is missing on one side, or recorded
#' from an older card (the 10th digit is the card issue number and changes every time a
#' card is replaced), removes most of that field's weight, so those true matches form a
#' second, lower cluster that 16 can cut off entirely -- see
#' \code{vignette("threshold-selection", package = "starling")}, which measures recall
#' for 3, 4 and 5 comparison variables. Comparing on \code{medicare09} (account plus
#' check digit) instead of \code{medicare10} avoids the card-issue penalty. Treat a change
#' in \code{compare_vars} -- or in how complete your Medicare field is -- as a reason to
#' re-plot and re-choose the threshold, not just a reason to re-run the same call.
#'
#' ## Phonetic-key deterministic matching
#' \code{method = "deterministic"} requires every \code{compare_vars} field to match
#' \emph{exactly} — by design, it has no tolerance for typos or transcription variants
#' the way \code{method = "probabilistic"} does. A common way to add tolerance back in,
#' without moving to probabilistic scoring at all, is to match on a \strong{phonetic
#' code} instead of the literal string: "Smith"/"Smyth" or "Catherine"/"Kathryn"
#' produce the same Soundex/Metaphone code even though the literal strings differ, so
#' an exact match on the \emph{code} tolerates exactly this kind of spelling variation
#' while remaining fully deterministic — no scoring, no threshold, same as any other
#' deterministic key.
#'
#' This needs \strong{no change to \code{murmuration()} itself} — \code{\link{flock}()}
#' already generates the phonetic codes required, via its \code{phonetic_vars}/
#' \code{phonetic_method = "soundex"} or \code{"metaphone"} arguments, appending a
#' \code{<col>_soundex}/\code{<col>_metaphone} column. Prepare both datasets before
#' calling this function, then compare on the phonetic column instead of the literal
#' name field:
#' \preformatted{
#' data_x_phon <- flock(data_x, phonetic_vars = c("lettername1", "lettername2"),
#'                    phonetic_method = "metaphone")
#' data_y_phon <- flock(data_y, phonetic_vars = c("lettername1", "lettername2"),
#'                    phonetic_method = "metaphone")
#'
#' murmuration(data_x_phon, data_y_phon,
#'   method       = "deterministic",
#'   compare_vars = c("lettername1_metaphone", "lettername2_metaphone", "dob"),
#'   ...)
#' }
#' Double Metaphone is generally preferred over Soundex when initial-letter/sound
#' variation is common, since Soundex keeps the first letter of the string literal —
#' see \code{\link{flock}}'s own documentation for the full comparison.
#'
#' This trades precision for recall in a specific, controllable way: it will now count
#' "Smith"/"Smyth" as identical, but it will \emph{also} count some genuinely different
#' names that happen to share a phonetic code as identical (both Soundex and Metaphone
#' have known false-positive pairs). Reserve it for fields where transcription
#' variation is the dominant error mode you're trying to catch, not as a general
#' replacement for exact matching across every field — most linkages benefit from
#' phonetic tolerance on \emph{name} fields specifically, while keeping
#' \code{dob}/\code{postcode}/\code{medicare} as ordinary exact matches. See Example 11
#' below for a complete worked version.
#'
#' ## Linking one cohort to several event streams
#' A common design is linking a single origin dataset (a cohort, a dialysis register, a
#' birth cohort) against several \emph{separate} event-stream datasets in independent
#' calls — one call per vaccine type, for instance, rather than one call with several
#' events mixed together. Every call produces the same generic column names
#' (\code{weights}, \code{vaccination_status}, \code{total_vax}, ...) regardless of which
#' event stream it linked, which collide the moment more than one call's result needs to
#' sit in the same data frame. \code{prefix}/\code{suffix} solve this at the source — see
#' their documentation below and Example 9 — rather than requiring a manual
#' \code{select()}/\code{rename()} pass on each call's output before joining. Keep each
#' call independent (against the origin dataset directly, never chaining one call's
#' result into the next as \code{data_x}) and combine the results afterwards with
#' \code{dplyr::left_join()} on \code{id_var}, which \code{prefix}/\code{suffix} never
#' touch.
#'
#' ## Troubleshooting: zero (or unexpectedly few) links
#' Five real, previously-encountered causes, not a theoretical list. The first three
#' affect the \emph{score} a genuine pair receives; the fourth stops candidate pairs
#' from being formed at all, upstream of scoring; the fifth is specific to
#' vaccination-type linkage. Run \code{\link{preflight}} before linking — it catches
#' most of causes 1 and 3 up front. After linking, \code{keep_pairs_pred = TRUE} plus
#' \code{\link{murmuration_plot}} (or a bare \code{summary(attr(x, "pairs_pred")$weights)})
#' tells you which kind of zero you have: a \strong{scoring} problem (plenty of candidate
#' pairs, none score high enough) or a \strong{blocking} problem (almost no candidates
#' exist to begin with) — the fix is different for each.
#' \enumerate{
#'   \item \strong{Letter-name fields swapped between \code{data_x} and \code{data_y}.} If one
#'     dataset's \code{mudnester::clean_the_nest()} call maps \code{lettername1} to first name
#'     and \code{lettername2} to last name, and the other dataset has them the other way
#'     around, every comparison is first-name-to-last-name and vice versa — structurally
#'     never a match, however clean the underlying data actually is. Each side looks
#'     correct in isolation (it's internally consistent), so this doesn't surface until
#'     the two \code{clean_the_nest()} calls are checked side by side, field mapping by
#'     field mapping.
#'   \item \strong{Too few \code{compare_vars} for the \code{threshold_value} in use.} The
#'     maximum achievable Fellegi-Sunter weight rises with the number of comparison
#'     fields — a threshold calibrated against 5 \code{compare_vars} (this function's
#'     default of \code{16} assumes that many) is often unreachable with 3 or 4,
#'     since AIR-style vaccination extracts frequently have no address field to compare.
#'     The console output naming this exact combination is "Pairs selected
#'     (threshold=TRUE): 0" against a large candidate-pair count — that pairing means the
#'     threshold, not the data, is the problem. Lower \code{threshold_value} and
#'     re-inspect with \code{\link{murmuration_plot}} rather than concluding the linkage
#'     itself failed.
#'   \item \strong{Inconsistent cleaning/encoding between \code{data_x} and \code{data_y}.} Two
#'     records that are obviously the same person to a human don't match if one side has
#'     \code{gender = "M"/"F"} and the other \code{"Male"/"Female"}, or one side's
#'     postcode is numeric and the other's is a zero-padded string. This degrades both
#'     \code{compare_vars} scoring (cause 2 above) and \code{blocking_var} alignment
#'     (cause 4 below) at once, since both are built from the same cleaned fields.
#'   \item \strong{\code{blocking_var} doesn't mean the same thing on both sides.}
#'     Blocking keys are literal string concatenations — a block built as
#'     \code{gender + postcode} on \code{data_x} and \code{gender + dob_yr} on \code{data_y}
#'     (or the "same" block built by two \code{\link{flock}()} calls whose inputs quietly
#'     drifted apart over time) never places a true pair in the same block, so they are
#'     never even compared, let alone scored. This is a \emph{candidate-formation}
#'     problem, not a scoring one — \code{keep_pairs_pred = TRUE} shows a near-zero
#'     candidate-pair count, not just a low-scoring one. A \code{blocking_var} that's too
#'     fine-grained on its own (e.g. exact \code{dob} with no coarser fallback) produces
#'     the identical symptom for a different reason: every block ends up with exactly one
#'     row per side, so there's nothing for a true pair to be compared against.
#'   \item \strong{Vaccination data left in long (one-row-per-dose) format.}
#'     \code{clean_the_nest(data_type = "vaccination")} does not pivot to
#'     \code{vax_date_1}/\code{vax_type_1}, ... on its own — \code{lie_nest_flat = TRUE}
#'     is required. Without it, the data passes through unpivoted with no error from
#'     \code{clean_the_nest()} itself; the failure only surfaces here, silently, as a
#'     vaccination-type (\code{event_type = "vaccination"}) linkage that never finds the
#'     columns it looks for.
#' }
#' @note Ensure there are no missing vaccination dates in vaccination dataset prior to murmuration. Murmuration requires complete vaccination data (equal date and type columns per observation) to achieve correct matching of vaccination columns.
#' If there are too few variables to match on, then matching will not work well. For example, if you have first name, last name and date of birth, and a very large dataset (Immunization Register), then the scoring will not differentiate true from false matches. Consider deterministic linkage when there is a paucity of information to use to derive linkage scores.
#'
#' @param data_x This is a dataframe object, cleaned using clean_the_nest, and would often represent the base, or "x" dataset (when doing left joins).
#' Typically this would be a dataset of cases, have enough data to create linkages, and have event dates (e.g., onset_date).
#' @param data_y This is a dataframe object, cleaned using clean_the_nest, and would often represent the admissions or vaccination dataset ("y" dataset when doing left joins).
#' Typically this would have enough data to create linkages, and include either admission data or vaccination event data (e.g. Australian Immunization Register).
#' \code{data_x}/\code{data_y} deliberately match the \code{.x}/\code{.y} suffixes
#' \pkg{reclin2}'s own \code{link()} applies to any colliding column name (including
#' \code{id_var} itself — see \code{id_var}\code{_y} in \code{@return} below): the naming
#' here isn't arbitrary, it's the same x/y distinction the join engine itself already
#' uses internally.
#' @param linkage_type The linkage to perform. As of starling 1.1.0 the preferred
#'   value is \code{"cohort2event"} — a single, explicit linkage concept: link a
#'   base linelist of people (a "cohort" — cases, hospitalisations, a birth cohort,
#'   an outbreak linelist, a flight manifest) to a stream of dated events (a "case"
#'   diagnosis, a "hospitalisation", or a "vaccination"), keeping events that fall in
#'   a time window relative to an anchor. Pair it with \code{event_type} to name the
#'   event stream. The window may sit \emph{before}, \emph{during}, or \emph{after}
#'   the cohort's follow-up, governed by \code{vax_window} (event-relative) or
#'   \code{cohort_window} (calendar).
#'
#'   \strong{Legacy codes}: \code{"c2h"} and \code{"v2c"} still work and route
#'   identically, but emit a deprecation warning nudging to the collapsed form:
#'   \code{"c2h"} (cases->hospitalisations) becomes \code{event_type =
#'   "hospitalisation"}; \code{"v2c"} (->vaccination history) becomes
#'   \code{event_type = "vaccination"}. \code{"v2h"} and \code{"v2e"} are
#'   \strong{defunct} (an error naming the replacement): \code{"v2h"} was
#'   \code{event_type = "vaccination"} with hospitalisations as data_x, and
#'   \code{"v2e"} is \code{event_type = "vaccination"} with a scalar \code{Date}
#'   passed to \code{event_date}. If \code{NULL} (default), \code{"c2h"} is assumed
#'   for backwards compatibility. When linking to a vaccination dataset, use a single
#'   row per person (run the vaccination data through \code{clean_the_nest()} with
#'   \code{lie_nest_flat = TRUE}).
#' @param event_type Required when \code{linkage_type = "cohort2event"}: names the
#'   event stream in data_y. One of \code{"vaccination"}, \code{"hospitalisation"}
#'   (\code{"hospitalization"} also accepted), or \code{"case"}. Ignored (with the
#'   routing inferred from the code) when a legacy \code{linkage_type} is supplied.
#' @param method Linkage engine. \code{"probabilistic"} (default) uses the
#'   Fellegi-Sunter EM pipeline (\pkg{reclin2}) with fuzzy field comparison and a
#'   score threshold — appropriate when identifiers are imperfect and you need
#'   weighted, partial-agreement matching. \code{"deterministic"} skips EM scoring
#'   entirely and links only records whose \code{compare_vars} match \emph{exactly}
#'   (all fields, no partial credit). Deterministic linkage is the right choice when
#'   you have a clean, reliable composite key (e.g. full name + date of birth, as in
#'   Smoll et al. 2023) and want reproducible, explainable, threshold-free links.
#'   Under \code{"deterministic"}: \code{weights} is \code{NA} (no composite score
#'   exists), \code{threshold_value}/\code{perch_before_linking} are ignored, and a
#'   non-unique key in data_y emits a message rather than silently resolving. All
#'   linkage types and both window mechanisms work under either method.
#' @param event_date The anchor date used to define the linkage window and determine valid vaccinations or related admissions.
#'   \itemize{
#'     \item Usually: a character string naming any date column present in data_x.
#'       Any date column produced by \code{mudnester::clean_the_nest()} is valid here — not just \code{"onset_date"} or
#'       \code{"admission_date"}. In birth-cohort or coverage studies, pass \code{"cohort_entry_date"} to anchor
#'       the vaccination window from the start of follow-up, or \code{"cohort_exit_date"} to bound from the end.
#'       Examples: \code{"onset_date"}, \code{"admission_date"}, \code{"cohort_entry_date"}, \code{"cohort_exit_date"}.
#'     \item For a fixed-date event (\code{event_type = "vaccination"} only): a single \code{Date} object (e.g., \code{as.Date("2024-12-15")}) representing a fixed
#'       event date shared across all records (e.g. a flight, an outbreak event).
#'   }
#'   This parameter is critical for determining valid vaccinations (must occur before \code{event_date - days_allowed_before_event})
#'   and for \code{"c2h"} linkages, determining disease-related admissions.
#' @param id_var Variable name (e.g. "id"). This is critical for data-linkage and the base dataset is the dataset you would left join onto (e.g. the "x" dataset). Cannot have missing data, or the observation will be lost in the linking process.
#' @param blocking_var Variable name (e.g. "block2"). Choice of blocking variable.
#' Use \code{flock()} to automatically construct up to three blocking variables
#' (coarse, medium, fine) from standardised linkage fields — it produces columns named
#' \code{block1}, \code{block2}, \code{block3} by default. You can also create your own.
#' Choose a broader/looser variable (e.g. gender) for small datasets or low data quality,
#' and a finer variable (e.g. postcode + DOB year) for large, high-quality datasets.
#' Required unless \code{resume_from_pairs} is used, in which case blocking has already
#' happened and this is ignored.
#' @param compare_vars Vector of variables. Used to compare variables between each dataset and calculate the string score differences. Typically names, dates of births and medicare/social security numbers.
#' Required unless \code{resume_from_pairs} is used, in which case comparing has already
#' happened and this is ignored.
#' @param threshold_value Numeric, default 16. The Fellegi-Sunter log-likelihood
#' ratio score at or above which a pair is classified as a match. The score is
#' \strong{dataset-specific} — there is no universal correct value.
#' Use \code{\link{murmuration_plot}} to inspect your score distribution and select a
#' threshold that sits in the natural valley between the match and non-match peaks — pass
#' this function's own return value straight into it (call this function with
#' \code{keep_pairs_pred = TRUE} first; see that argument), no need to re-derive the pairs table.
#' Higher values increase specificity (fewer false links); lower values increase sensitivity
#' (fewer missed true matches). \strong{The default of 16 assumes 5 \code{compare_vars}
#' with an agreeing Medicare number} (see "On threshold selection" above): with 3 or 4
#' fields, or an incomplete Medicare field, the valley usually sits much lower. Re-derive
#' the threshold via \code{\link{murmuration_plot}} whenever \code{compare_vars}' length
#' changes — don't carry a threshold forward from a differently-sized comparison set (see
#' Troubleshooting, cause 2, below).
#' @param perch_before_linking Logical. If \code{TRUE}, calls \code{\link{perch}} on the
#' scored pairs immediately after the EM algorithm fits and before the threshold is applied,
#' printing a sensitivity table showing match counts and clerical burden across the score
#' range. In interactive sessions the plot is also displayed. In non-interactive sessions
#' (e.g. Quarto render, batch jobs) only the table is printed and execution continues
#' automatically. Default \code{FALSE}.
#' @param inspect_matches Logical. If \code{TRUE}, stops immediately after blocking and
#'   EM scoring and returns the full candidate-pairs table for manual review, instead of
#'   proceeding to threshold selection, \code{link()}, or any \code{linkage_type}-specific
#'   event-window filtering. Requires \code{method = "probabilistic"}. Unlike the plain
#'   \code{pairs_pred} attached elsewhere, this table has the actual \code{id_var} and
#'   \code{compare_vars} field values joined in from both sides (suffixed \code{.x}/\code{.y}),
#'   so a specific pair can be read and judged, not just indexed. A mechanical
#'   \code{threshold} column (\code{weights >= threshold_value}) and an empty
#'   \code{match_override} column (\code{NA}/\code{"match"}/\code{"non-match"}, for you to
#'   fill in on individual pairs after inspection) are included. \code{event_date} is not
#'   required when this is \code{TRUE}, since the event-window stage that needs it is never
#'   reached. Pass the result straight to \code{\link{murmuration_plot}} to see the score
#'   distribution. Default \code{FALSE}.
#' @param inspect_extra_vars Character vector of additional column names (present in both
#'   \code{data_x} and \code{data_y}) to join into the \code{inspect_matches} review table beyond
#'   \code{id_var} and \code{compare_vars} -- e.g. a field you want visible for context while
#'   reviewing but did not include in the linkage key itself. Ignored unless
#'   \code{inspect_matches = TRUE}. Default \code{NULL}.
#' @param resume_from_pairs A reviewed pairs table -- normally the object returned by
#'   \code{murmuration(inspect_matches = TRUE)}, optionally with its \code{match_override}
#'   column edited (\code{NA}/\code{"match"}/\code{"non-match"}) -- to resume the linkage
#'   from, skipping blocking, comparing, and EM-fitting entirely. \code{blocking_var} and
#'   \code{compare_vars} are not needed in this mode. \code{data_x}/\code{data_y} must be the
#'   exact same data frames used to produce it (a row-count check guards against an
#'   obvious mismatch, but cannot catch every case) -- \code{.x}/\code{.y} are row indices
#'   into them. Mutually exclusive with \code{inspect_matches}. \code{threshold_value}
#'   still applies here, to any pair whose \code{match_override} is \code{NA} -- it does
#'   not have to match what was used when \code{inspect_matches} built this table.
#'   Default \code{NULL}.
#' @param vax_window A named list defining the vaccination validity window relative
#'   to \code{event_date}. Used for every vaccination linkage
#'   (\code{event_type = "vaccination"}). Contains up to three elements:
#'   \describe{
#'     \item{\code{days_before}}{Minimum days \emph{before} \code{event_date} a vaccine
#'       must have been administered to count as prior valid exposure. Vaccines administered
#'       within this buffer (i.e. fewer than \code{days_before} days before the event)
#'       are excluded because they may not yet have conferred protection. Default \code{14}.
#'       Set to \code{0} if you want all prior vaccinations regardless of recency
#'       (e.g. birth-cohort studies where protection is assumed immediate).}
#'     \item{\code{days_after}}{Maximum days \emph{after} \code{event_date} to retain
#'       a vaccination record. Default \code{0} (no post-event vaccinations retained).
#'       Set to a positive value to capture post-event vaccination -- useful for studying
#'       vaccination uptake following a diagnosis or intervention (e.g. 365 to capture
#'       all vaccinations in the year after diagnosis). Post-event doses are kept
#'       in the \code{vax_date_<n>} columns like any other valid dose; compare them
#'       with \code{event_date} downstream if pre- and post-event doses need
#'       separating.}
#'     \item{\code{lookback_days}}{Maximum days \emph{before} \code{event_date} to look
#'       back for valid vaccinations. Default \code{Inf} (no upper lookback limit --
#'       any prior vaccination regardless of age counts). Set to \code{365} to restrict
#'       to vaccinations in the past year only (appropriate for annual influenza VE
#'       studies where a vaccine from five years ago is not the relevant exposure).}
#'   }
#'   Example: \code{vax_window = list(days_before = 14, days_after = 0, lookback_days = 365)}
#'   retains only vaccinations that occurred 14 or more days before the event and within
#'   the past year, excluding both too-recent and too-old doses.
#'
#'   \strong{Overridden by \code{cohort_window} when both are supplied.} For
#'   \code{"c2h"} linkages, use \code{days_allowed_before_event} and
#'   \code{days_allowed_after_event} instead (those govern the admission window, not
#'   vaccination).
#'
#' @param cohort_window An optional named list defining a \emph{calendar-based}
#'   observation window that \strong{replaces} \code{vax_window} when supplied.
#'   Use this for cohort studies, coverage analyses, or any design where the
#'   observation window is defined at the study level rather than relative to
#'   individual event dates:
#'   \describe{
#'     \item{\code{entry}}{A \code{Date} object (or character coercible to Date via
#'       \code{as.Date()}) giving the start of the observation window. Vaccinations
#'       before this date are excluded. Example: \code{as.Date("2024-01-01")}.}
#'     \item{\code{exit}}{A \code{Date} object giving the end of the observation
#'       window. Vaccinations after this date are excluded. Example:
#'       \code{as.Date("2024-12-31")}.}
#'   }
#'   When \code{cohort_window} is supplied, \code{vax_window} is ignored entirely
#'   and a message is emitted confirming which window is active. Useful when the
#'   event date is unknown, unreliable, or irrelevant to the vaccination exposure
#'   definition (e.g. a prevalence survey, a screening programme evaluation, or a
#'   birth-cohort study anchored to a fixed follow-up period rather than individual
#'   birth dates). Can be used with any \code{linkage_type}.
#'
#' @param days_allowed_before_event Numeric. For \code{"c2h"} linkages only: the
#'   lower bound of the disease-related admission window -- how many days before
#'   \code{event_date} an admission may still be considered disease-related. Default
#'   \code{7}. Not used for vaccination linkage types; use \code{vax_window$days_before}
#'   for those.
#' @param days_allowed_after_event Numeric. For \code{"c2h"} linkages only: the upper
#'   bound of the disease-related admission window -- how many days after
#'   \code{event_date} an admission may still be considered disease-related. Default
#'   \code{14}. Not used for vaccination linkage types; use \code{vax_window$days_after}
#'   for those.
#' @param one_row_per_person Logical (TRUE or FALSE) with the default being TRUE. It will take multiple admissions per person, and create a series of variables prefixed with "first_", such as "first_admission_date",
#' and put into a single row all admission events, and create a series of variables suffixed with "s", such as "admission_dates". Will work with single admissions per person.
#' @param clean_eggs Logical (TRUE or FALSE) with the default being TRUE. Drops all the .y variables that are duplicates of the second dataset (data_y), and keeps the variables and removes the .x from data_x. If you leave this on, many,
#' if not most variables will have ".x" or ".y" attached to them (e.g. gender) and thus keep this as TRUE for default, and FALSE if you want to check the linkages are true and working.
#' @param prefix,suffix Character or \code{NULL} (default \code{NULL} for both --
#'   every existing call is unaffected unless you supply one). Text added to every
#'   column \code{murmuration()} contributes from \code{data_y} or computes about the
#'   link itself (\code{weights}, \code{threshold}, and every derived column the
#'   relevant \code{linkage_type} block produces -- e.g. \code{vaccination_status},
#'   \code{total_vax}, \code{first_vax_date}, \code{first_vax_type},
#'   \code{last_vax_date}, \code{last_vax_type}, a colliding identifier from
#'   \code{data_y}'s side). \code{data_x}'s own columns are never touched -- \code{id_var}
#'   in particular stays plain and stable, which is what makes it a reliable join
#'   key afterwards. If both are supplied, prefix is applied first, then suffix.
#'
#'   Built for linking one origin dataset against several separate event-stream
#'   datasets in independent calls -- e.g. a cohort against flu, COVID, RSV, and
#'   shingles vaccination registries -- where every call would otherwise produce
#'   the same generic column names (\code{weights}, \code{vaccination_status}, ...)
#'   and collide the moment you try to combine the results:
#'   \preformatted{
#'   flu_linked   <- murmuration(cohort, flu_vax,   ..., suffix = "_flu")
#'   covid_linked <- murmuration(cohort, covid_vax, ..., suffix = "_covid")
#'   rsv_linked   <- murmuration(cohort, rsv_vax,   ..., suffix = "_rsv")
#'
#'   master <- cohort |>
#'     dplyr::left_join(dplyr::select(flu_linked,   id_var, ends_with("_flu")),   by = "id_var") |>
#'     dplyr::left_join(dplyr::select(covid_linked, id_var, ends_with("_covid")), by = "id_var") |>
#'     dplyr::left_join(dplyr::select(rsv_linked,   id_var, ends_with("_rsv")),   by = "id_var")
#'   }
#'   No manual \code{select()}/\code{rename()} step is needed between each
#'   \code{murmuration()} call and the join, and each call stays fully independent
#'   -- none of them ever needs a previous call's result as input, so there is
#'   nothing accumulated for column names to collide with in the first place.
#' @param days_between_onset_death Numeric (e.g. "30"). If you have put a date of death into the clean_the_nest command (which will rename it to dod), then the command will find disease related dates of death.
#' This is chosen number of days between event_date and death for a disease-related death. Often this may be 30 days for SARS-CoV-2 or can be much longer for HIV. If you don't want an upper limit, use "9999".
#' @param last_follow_up Represents a date (input as ymd(2024-11-22)) that represents last follow-up. This could be the latest admission date of a dataset. Used for calculating survival time.
#' @param keep_pairs_pred Logical. If \code{TRUE}, attach the full candidate-pairs
#'   scoring table as \code{attr(x, "pairs_pred")} on the returned data frame (see
#'   \code{@return} below) -- needed for \code{\link{murmuration_plot}}'s score
#'   distribution, or for pulling the table out directly. Default \code{FALSE}:
#'   OPT-IN, not opt-out. This table is every candidate pair blocking produced,
#'   before thresholding -- confirmed in practice to reach tens of millions of
#'   rows for a single call against a large \code{data_y}, and it is serialised in
#'   full by \code{saveRDS()} along with everything else on the returned object,
#'   not just held in memory for the duration of the call. Leave at \code{FALSE}
#'   unless you specifically need \code{murmuration_plot()} or
#'   \code{attr(x, "pairs_pred")} for THIS call -- re-run just that one call with
#'   \code{keep_pairs_pred = TRUE} when you do, rather than defaulting it on
#'   everywhere. Does not affect \code{inspect_matches = TRUE}, which already
#'   returns the full pairs table directly as its return value -- an
#'   already-explicit request, not governed by this argument.
#' @return A linked dataset with some new variables. When \code{method = "probabilistic"},
#'   linkage metadata travels with the returned data frame as
#'   attributes, so \code{\link{murmuration_plot}} can be called directly on it with no
#'   need to re-derive anything via \pkg{reclin2} -- except \code{pairs_pred} itself,
#'   which is only attached when \code{keep_pairs_pred = TRUE} (default \code{FALSE},
#'   see that argument above):
#'   \describe{
#'     \item{\code{attr(x, "pairs_pred")}}{Only present when \code{keep_pairs_pred = TRUE}
#'       was passed to this call; \code{NULL} otherwise (the new default). The full
#'       candidate-pairs scoring table (every
#'       compared pair, not just the ones that cleared \code{threshold_value}), with the
#'       realised \code{blocking_var} value joined on for each pair.}
#'     \item{\code{attr(x, "blocking_var")}}{The \code{blocking_var} value actually used.}
#'     \item{\code{attr(x, "threshold_value")}}{The \code{threshold_value} actually used.}
#'     \item{\code{attr(x, "murmuration_method")}}{\code{"probabilistic"} or
#'       \code{"deterministic"}.}
#'   }
#'   Under \code{method = "deterministic"}, \code{pairs_pred} is \code{NULL} and
#'   \code{blocking_var}/\code{threshold_value} are \code{NA}: there is no Fellegi-Sunter
#'   score to inspect (deterministic linkage is threshold-free by design -- see
#'   \code{method} above).
#'
#'   \strong{When \code{inspect_matches = TRUE}}, the return shape is different: instead of
#'   a linked dataset, you get the candidate-pairs table itself (see \code{inspect_matches}
#'   above) -- there is no linked \code{df} yet, because threshold selection and any
#'   event-window filtering have not run. The same four attributes are still attached
#'   directly to this table (not nested under a separate \code{pairs_pred} attribute, since
#'   the returned object already \emph{is} that table).
#'
#'   \strong{\code{<id_var>_y}} (corrected this revision -- previously the fixed literal
#'   \code{"id_var_df2"}, and \code{"c2h"}-only; now dynamically named the same way
#'   \code{id_x_col}/\code{id_y_col} are internally, and consistent across all four
#'   \code{linkage_type} routes): the matched \code{data_y} record's own \code{id_var}
#'   value, one row per person after the \code{one_row_per_person} collapse -- e.g.
#'   \code{id_var = "identity"} produces a column literally named \code{identity_y}.
#'   \code{NA} for a \code{data_x} row with no above-threshold match. This is the only
#'   trace of \emph{which specific data_y record} a person matched to that survives
#'   \code{clean_eggs = TRUE} (the default) -- every other \code{.y}-suffixed column
#'   from \code{data_y} is dropped by that step, by design. Keep it if you need an
#'   audit trail back to a specific vaccination dose or hospital admission record;
#'   otherwise it's safe to ignore alongside the domain-specific derived columns
#'   (\code{vaccination_status}, \code{first_admission_date}, ...).
#'
#' @references
#' \strong{Deterministic linkage (\code{method = "deterministic"}) -- Applied.} Smoll NR,
#' Al Imam MH, Shulz C, Booy R, Khandaker G. The effectiveness of vaccination for
#' preventing hospitalisation with COVID-19 in regional Queensland: a data linkage study.
#' \emph{Med J Aust} 2023;219(4):162--165. \doi{10.5694/mja2.52019}. This study's linkage key --
#' full first name (middle names omitted), full last name, and date of birth -- across
#' PCR results (NOCS), hospital admissions, and the Australian Immunisation Register is
#' exactly what \code{method = "deterministic"} reproduces directly: pass
#' \code{compare_vars = c("lettername1", "lettername2", "dob")} to recreate this study's
#' linkage recipe (see Example 3).
#'
#' \strong{Probabilistic linkage (\code{method = "probabilistic"}) -- Precedent.} The same
#' Smoll et al. 2023 study is also the methodological precedent for
#' \code{method = "probabilistic"}: the Fellegi-Sunter approach generalises that study's
#' exact-key design to settings where identifiers are too incomplete or inconsistent for
#' an all-fields-exact match to find every true pair.
#'
#' \strong{Applied.} Roughan WH, Cavanagh-Ochiai B, Barnekow D, Lottkowitz K, Schafer K,
#' Camadoo L, Smoll NR. Infant respiratory syncytial virus (RSV) immunisation coverage in
#' an Australian regional area, 1 May 2024 - 31 October 2025. \emph{Commun Dis Intell}
#' 2026;50. \doi{10.33321/cdi.2026.50.031}. About 10\% of this study's birth-register
#' records carried identifier or date-of-birth errors that ruled out an exact match to
#' hospital and immunisation records; these were resolved with weighted-score matching
#' across multiple data points -- the situation \code{method = "probabilistic"} together
#' with \code{\link{murmuration_plot}} is built for, where a clean composite key (as in
#' the Smoll et al. 2023 design above) cannot be assumed and the score distribution should
#' be inspected before a threshold is chosen. The birth-cohort linkage in Example 4 below
#' mirrors this study's maternal-neonatal dyad plus vaccination-register design; Example 10
#' reproduces its two-pathway (nirsevimab/Abrysvo) linkage directly.
#'
#' @importFrom dplyr mutate filter select group_by ungroup arrange slice
#' @importFrom dplyr bind_rows n_distinct case_when if_else starts_with ends_with
#' @importFrom dplyr any_of all_of across rowwise c_across rename rename_with
#' @importFrom dplyr everything desc left_join
#' @importFrom stringr str_extract str_remove
#' @importFrom lubridate ymd
#' @importFrom reclin2 pair_blocking compare_pairs jaro_winkler problink_em
#' @importFrom reclin2 select_threshold link
#' @importFrom rlang sym :=
#' @importFrom stats predict reformulate
#' @export
#' @examples
#' # All eleven examples below are kept as \dontrun{} rather than \donttest{}
#' # (aviary ecosystem plan \S8.9.2): each illustrates a distinct calling
#' # convention (linkage_type/event_type routing, threshold selection,
#' # manual pair review, multi-registry linking) against realistically
#' # multi-hundred-row blocking -- the score-distribution and threshold
#' # behaviour these are actually meant to show does not survive shrinking
#' # to a handful of synthetic rows small enough to embed here, and several
#' # depend on \pkg{reclin2}'s full pair_blocking()/problink_em() pipeline
#' # rather than the simpler deterministic exact-join path. See
#' # \code{\link{murmuration_plot}} for a genuinely runnable \donttest{}
#' # example against this same package's plotting function.
#' \dontrun{
#' # Example 1: Link cases to vaccination history (onset_date as anchor),
#' # using the current preferred interface (linkage_type = "cohort2event" +
#' # event_type) introduced in starling 1.1.0.
#' dx_clean <- clean_the_nest(dx_data,
#'   data_type = "cases",
#'   id_var = "identity",
#'   lettername1 = "first_name",
#'   lettername2 = "surname",
#'   dob = "date_of_birth",
#'   gender = "gender",
#'   postcode = "postcode",
#'   medicare = "medicare_no",
#'   diagnosis = "disease_name",
#'   onset_date = "date_of_onset")
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
#' data_x <- murmuration(dx_clean, vax_clean,
#'   linkage_type = "cohort2event",
#'   event_type   = "vaccination",
#'   event_date = "onset_date",
#'   id_var = "identity",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob"),
#'   clean_eggs = FALSE)
#' # Legacy equivalent (still works, but emits a deprecation warning):
#' #   murmuration(dx_clean, vax_clean, linkage_type = "v2c", ...)
#'
#' # Example 2: Link hospitalization data to vaccination history (admission_date as anchor)
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
#' data_y <- murmuration(hosp_clean, vax_clean,
#'   linkage_type = "cohort2event",
#'   event_type   = "vaccination",
#'   event_date = "admission_date",
#'   id_var = "patient_id",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "medicare10", "dob"),
#'   clean_eggs = FALSE,
#'   one_row_per_person = TRUE)
#' # The hospitalisations are simply the cohort here (data_x); this replaces the
#' # defunct linkage_type = "v2h".
#'
#' # Example 3: Link cases to hospitalisations (onset_date as anchor) --
#' # method = "deterministic" reproducing Smoll et al. 2023's exact linkage
#' # recipe (see @references): full first name, full last name, and date of
#' # birth, linking 1 Jan-31 Mar 2022 NOCS PCR results to Central Queensland
#' # hospital admissions. The paper defines a COVID-19-related admission as
#' # ICD-10-AM B97.2 within 30 days of diagnosis (no pre-diagnosis allowance
#' # is described in the paper, hence days_allowed_before_event = 0 here --
#' # widen this if your own admissions data can genuinely precede diagnosis,
#' # e.g. a PCR test taken during an already-in-progress admission).
#' # Reported result: 47/9,682 adults hospitalised (0.48%); vaccine
#' # effectiveness 69.9% for a primary course alone, 81.8% with a booster.
#' det_linked <- murmuration(dx_clean, hosp_clean,
#'   linkage_type = "cohort2event",
#'   event_type   = "hospitalisation",
#'   method       = "deterministic",
#'   event_date = "onset_date",
#'   id_var = "identity",
#'   blocking_var = "postcode",
#'   compare_vars = c("lettername1", "lettername2", "dob"),
#'   days_allowed_before_event = 0,
#'   days_allowed_after_event = 30,
#'   clean_eggs = FALSE)
#' # Legacy equivalent (still works, but emits a deprecation warning):
#' #   murmuration(dx_clean, hosp_clean, linkage_type = "c2h", method = "deterministic", ...)
#'
#' # Example 4: Birth-cohort study — cohort_window approach, reproducing the
#' # observation window of Roughan et al. 2026 (see @references): all 5,116
#' # infants born at Sunshine Coast University Hospital, 1 May 2024-31 Oct
#' # 2025. Captures any vaccination administered during the defined follow-up
#' # period, regardless of individual event dates. That study resolved its
#' # messier ~10% of birth-register records (identifier/DOB errors) with
#' # weighted-score matching rather than an exact key -- exactly the case
#' # method = "probabilistic" (and murmuration_plot(), Example 7 below) is
#' # built to support. Reported result: 71.7% (3,669/5,116) received valid
#' # RSV immunisation.
#' cohort_clean <- clean_the_nest(birth_cohort_data,
#'   data_type = "cases",
#'   id_var = "baby_id",
#'   lettername1 = "first_name",
#'   lettername2 = "last_name",
#'   dob = "babys_date_of_birth",
#'   gender = "sex")
#'
#' df_cohort <- murmuration(cohort_clean, vax_clean,
#'   linkage_type  = "cohort2event",
#'   event_type    = "vaccination",
#'   event_date    = "dob",
#'   id_var        = "baby_id",
#'   blocking_var  = "gender",
#'   compare_vars  = c("lettername1", "lettername2", "dob"),
#'   cohort_window = list(
#'     entry = as.Date("2024-05-01"),
#'     exit  = as.Date("2025-10-31")
#'   ),
#'   clean_eggs = FALSE)
#'
#' # Example 4b: Annual influenza VE study — 1-year lookback, 14-day buffer.
#' # Only vaccinations in the 12 months before onset, and >= 14 days before onset,
#' # are retained as valid prior exposure.
#' df_flu_ve <- murmuration(cases_clean, vax_clean,
#'   linkage_type = "v2c",
#'   event_date   = "onset_date",
#'   id_var       = "identity",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob", "medicare10"),
#'   vax_window   = list(
#'     days_before   = 14,
#'     days_after    = 0,
#'     lookback_days = 365
#'   ),
#'   clean_eggs = FALSE)
#'
#' # Example 4c: Post-diagnosis vaccination uptake study.
#' # Retains vaccinations up to 1 year AFTER diagnosis as well as valid prior doses;
#' # compare each vax_date_<n> with onset_date to separate pre- from post-event doses.
#' df_post_dx <- murmuration(cases_clean, vax_clean,
#'   linkage_type = "v2c",
#'   event_date   = "onset_date",
#'   id_var       = "identity",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob", "medicare10"),
#'   vax_window   = list(
#'     days_before   = 14,
#'     days_after    = 365,
#'     lookback_days = Inf
#'   ),
#'   clean_eggs = FALSE)
#'
#' # Example 5: Link flight manifest to vaccination history (fixed event Date object)
#' manifest_clean <- clean_the_nest(manifest_data,
#'   data_type = "cases",
#'   id_var = "passenger_id",
#'   lettername1 = "first_name",
#'   lettername2 = "surname",
#'   dob = "date_of_birth",
#'   gender = "gender")
#'
#' df_flight <- murmuration(manifest_clean, vax_clean,
#'   linkage_type = "cohort2event",
#'   event_type   = "vaccination",
#'   event_date = as.Date("2024-03-15"),   # a single Date: the fixed-date event path
#'   id_var = "passenger_id",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob"),
#'   vax_window = list(days_before = 14, days_after = 0, lookback_days = Inf),
#'   clean_eggs = FALSE)
#' # Replaces the defunct linkage_type = "v2e".
#'
#' # Example 6: Link outbreak linelist to vaccination history (fixed event Date object)
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
#'   linkage_type = "cohort2event",
#'   event_type   = "vaccination",
#'   event_date = as.Date("2024-06-01"),
#'   id_var = "case_id",
#'   blocking_var = "postcode",
#'   compare_vars = c("lettername1", "lettername2", "dob", "medicare10"),
#'   vax_window = list(days_before = 7, days_after = 0, lookback_days = Inf),
#'   clean_eggs = FALSE)
#'
#' # Examples 7-9 below continue to use the legacy "v2c" code for brevity, since
#' # the point of each is the surrounding workflow (threshold inspection, manual
#' # review, multi-registry linking) rather than the linkage_type argument itself.
#' # Substitute linkage_type = "cohort2event", event_type = "vaccination" in any
#' # of them if you are starting a fresh script.
#'
#' # Example 7: Inspect the linkage score distribution before committing to a
#' # threshold, in the spirit of Roughan et al. 2026 (@references) — that study
#' # did not assume its composite key would always match exactly, and resolved
#' # its messier ~10% of birth-register records with weighted-score matching
#' # instead.
#' linked <- murmuration(cohort_clean, vax_clean,
#'   linkage_type  = "v2c",
#'   event_date    = "dob",
#'   id_var        = "baby_id",
#'   blocking_var  = "gender",
#'   compare_vars  = c("lettername1", "lettername2", "dob"),
#'   threshold_value = 16,
#'   clean_eggs = FALSE,
#'   keep_pairs_pred = TRUE)   # needed for murmuration_plot()/attr() below --
#'                             # default FALSE would leave pairs_pred as NULL
#'
#' # murmuration_plot() runs directly off this return value — no need to hand-roll
#' # the reclin2 pipeline (pair_blocking() -> compare_pairs() -> problink_em() ->
#' # predict()) just to see the score distribution. threshold_value and
#' # blocking_var are picked up automatically from attr(linked, ...).
#' murmuration_plot(linked)
#'
#' # Try a different threshold without re-running murmuration():
#' murmuration_plot(linked, threshold = 14)
#'
#' # The pairs-level scoring table is also available directly if you need it:
#' pp <- attr(linked, "pairs_pred")
#'
#' # Example 8: Manual review of borderline pairs before finalising. Useful
#' # when murmuration_plot() shows an above-threshold count that doesn't
#' # match independent expectations (e.g. 79 candidate matches where only 76
#' # are plausible) and you need to see WHICH pairs are in question, not just
#' # how many.
#' review <- murmuration(cohort_clean, vax_clean,
#'   linkage_type  = "v2c",
#'   id_var        = "baby_id",
#'   blocking_var  = "gender",
#'   compare_vars  = c("lettername1", "lettername2", "dob"),
#'   threshold_value  = 16,
#'   inspect_matches  = TRUE)   # event_date not needed for this mode
#'
#' murmuration_plot(review)   # same plot, run on the full candidate-pairs table
#'
#' # Pull out the above-threshold pairs with their identifying fields
#' # (lettername1.x/.y, lettername2.x/.y, dob.x/.y, baby_id.x/.y) for a look:
#' suspects <- review[review$threshold, ]
#'
#' # After eyeballing names/DOBs, record a manual decision on specific pairs:
#' review$match_override[review$baby_id.x == "B0042" &
#'                        review$baby_id.y == "V0917"] <- "non-match"
#'
#' # Finish the linkage, respecting that override. data_x/data_y must be the exact
#' # same data frames used above -- .x/.y are row indices into them.
#' linked <- murmuration(cohort_clean, vax_clean,
#'   linkage_type      = "v2c",
#'   id_var            = "baby_id",
#'   event_date        = "dob",
#'   threshold_value   = 16,
#'   resume_from_pairs = review)
#'
#' # Example 9: Link one origin cohort against several separate vaccine
#' # registries (flu, COVID, RSV), each call fully independent -- never
#' # chaining one call's result into the next -- with suffix distinguishing
#' # each call's contributed columns so the results can be joined back onto
#' # the cohort without any manual select()/rename() step in between.
#' flu_linked <- murmuration(dialysis_cohort, flu_vax,
#'   linkage_type = "v2c", event_date = "onset_date", id_var = "id_var",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob"),
#'   suffix = "_flu")
#'
#' covid_linked <- murmuration(dialysis_cohort, covid_vax,
#'   linkage_type = "v2c", event_date = "onset_date", id_var = "id_var",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob"),
#'   suffix = "_covid")
#'
#' rsv_linked <- murmuration(dialysis_cohort, rsv_vax,
#'   linkage_type = "v2c", event_date = "onset_date", id_var = "id_var",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob"),
#'   suffix = "_rsv")
#'
#' # id_var is untouched by suffix on every call, so it's a stable join key
#' master_linkage <- dialysis_cohort |>
#'   dplyr::left_join(dplyr::select(flu_linked,   id_var, ends_with("_flu")),
#'                     by = "id_var") |>
#'   dplyr::left_join(dplyr::select(covid_linked, id_var, ends_with("_covid")),
#'                     by = "id_var") |>
#'   dplyr::left_join(dplyr::select(rsv_linked,   id_var, ends_with("_rsv")),
#'                     by = "id_var")
#'
#' # Example 10: Link a birth cohort to its two separate RSV immunisation
#' # pathways -- neonatal nirsevimab and antenatal Abrysvo -- reproducing the
#' # two-pathway design of Roughan et al. 2026 (@references) as two independent
#' # calls joined via id_var, following the multi-registry pattern in Example 9.
#' # Reported result: of the 71.7% of infants with valid RSV immunisation,
#' # 45.2 percentage points came through nirsevimab only and 26.4 through
#' # Abrysvo only.
#' #
#' # Abrysvo is only valid if given >= 14 days before delivery -- vax_window
#' # captures that floor directly. The paper's additional 28-36-week-gestation
#' # criterion is a property of the vaccination date itself, not of its distance
#' # to event_date, so it is not something vax_window can express: check it as
#' # a separate filter on vax_date/gestational_age after linkage.
#' abrysvo_linked <- murmuration(cohort_clean, abrysvo_vax,
#'   linkage_type = "cohort2event", event_type = "vaccination",
#'   event_date   = "dob",
#'   id_var       = "baby_id",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob"),
#'   vax_window   = list(days_before = 14, days_after = 0, lookback_days = 280),
#'   suffix       = "_abrysvo")
#'
#' # Nirsevimab is only valid if given before hospital discharge; discharge
#' # timing varies per infant, so this is approximated with a short fixed
#' # post-birth window wide enough to cover a typical birth admission --
#' # confirm against your own discharge_date field rather than relying on
#' # this alone.
#' nirsevimab_linked <- murmuration(cohort_clean, nirsevimab_vax,
#'   linkage_type = "cohort2event", event_type = "vaccination",
#'   event_date   = "dob",
#'   id_var       = "baby_id",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob"),
#'   vax_window   = list(days_before = 0, days_after = 10, lookback_days = 0),
#'   suffix       = "_nirsevimab")
#'
#' rsv_coverage <- cohort_clean |>
#'   dplyr::left_join(dplyr::select(abrysvo_linked, baby_id, ends_with("_abrysvo")),
#'                     by = "baby_id") |>
#'   dplyr::left_join(dplyr::select(nirsevimab_linked, baby_id, ends_with("_nirsevimab")),
#'                     by = "baby_id")
#'
#' # Example 11: Deterministic linkage with phonetic tolerance -- flock() first to
#' # generate metaphone codes, then match on those instead of the literal strings
#' # (see "Phonetic-key deterministic matching" above). No murmuration() argument
#' # changes at all -- this is purely a data-preparation technique.
#' dx_phon <- flock(dx_clean,
#'   phonetic_vars   = c("lettername1", "lettername2"),
#'   phonetic_method = "metaphone")
#' hosp_phon <- flock(hosp_clean,
#'   phonetic_vars   = c("lettername1", "lettername2"),
#'   phonetic_method = "metaphone")
#'
#' phon_linked <- murmuration(dx_phon, hosp_phon,
#'   linkage_type = "cohort2event",
#'   event_type   = "hospitalisation",
#'   method       = "deterministic",
#'   event_date   = "onset_date",
#'   id_var       = "identity",
#'   blocking_var = "postcode",
#'   compare_vars = c("lettername1_metaphone", "lettername2_metaphone", "dob"),
#'   clean_eggs   = FALSE)
#' # Catches "Smith"/"Smyth"-style transcription variants that exact-string
#' # deterministic matching (Example 3) would miss -- at the cost of also
#' # conflating rarer, genuinely different names sharing a phonetic code.
#' }
murmuration <- function(data_x,
                        data_y,
                        linkage_type = NULL,
                        event_type   = NULL,
                        method       = "probabilistic",
                        event_date = NULL,
                        id_var,
                        blocking_var       = NULL,
                        compare_vars       = NULL,
                        threshold_value      = 16,
                        perch_before_linking = FALSE,
                        inspect_matches      = FALSE,
                        inspect_extra_vars   = NULL,
                        resume_from_pairs    = NULL,
                        vax_window = list(
                          days_before   = 14,
                          days_after    = 0,
                          lookback_days = Inf
                        ),
                        cohort_window             = NULL,
                        days_allowed_before_event = 7,
                        days_allowed_after_event  = 14,
                        one_row_per_person        = TRUE,
                        clean_eggs                = TRUE,
                        prefix                    = NULL,
                        suffix                    = NULL,
                        days_between_onset_death  = 30,
                        last_follow_up            = NULL,
                        keep_pairs_pred           = FALSE) {

  # --- Resolve method -------------------------------------------------------
  method <- match.arg(method, c("probabilistic", "deterministic"))

  if (isTRUE(inspect_matches) && method == "deterministic") {
    stop(
      "(*)> starling::murmuration() \u2014 inspect_matches = TRUE requires ",
      "method = \"probabilistic\": deterministic linkage has no Fellegi-Sunter ",
      "score to inspect (weights is NA throughout, by design).",
      call. = FALSE
    )
  }

  if (isTRUE(inspect_matches) && !is.null(resume_from_pairs)) {
    stop(
      "(*)> starling::murmuration() \u2014 'inspect_matches' and 'resume_from_pairs' ",
      "are mutually exclusive: the first stops before scoring is applied, the ",
      "second resumes after it. Use one or the other, not both.",
      call. = FALSE
    )
  }

  # blocking_var/compare_vars drive blocking + EM scoring, which never run
  # when resuming from an already-scored pairs table -- so they are only
  # required on a normal (or inspect_matches) run.
  if (is.null(resume_from_pairs)) {
    if (is.null(blocking_var)) {
      stop(
        "(*)> starling::murmuration() \u2014 'blocking_var' must be supplied ",
        "(unless 'resume_from_pairs' is used)."
      )
    }
    if (is.null(compare_vars)) {
      stop(
        "(*)> starling::murmuration() \u2014 'compare_vars' must be supplied ",
        "(unless 'resume_from_pairs' is used)."
      )
    }
  }

  if (!is.null(prefix) && (!is.character(prefix) || length(prefix) != 1 || is.na(prefix))) {
    stop("(*)> starling::murmuration() \u2014 'prefix' must be a single character string, or NULL.")
  }
  if (!is.null(suffix) && (!is.character(suffix) || length(suffix) != 1 || is.na(suffix))) {
    stop("(*)> starling::murmuration() \u2014 'suffix' must be a single character string, or NULL.")
  }

  # --- Resolve linkage nomenclature (collapsed model + legacy aliases) ------
  # The four historical codes (c2h/v2c/v2h/v2e) are all the same operation:
  # link a base linelist of people ("cohort") to a dated event stream, keeping
  # events in a time window relative to an anchor. As of starling 1.1.0 the
  # preferred API expresses this as ONE linkage_type = "cohort2event" plus an
  # event_type naming the stream. "c2h" and "v2c" still work (deprecation
  # warning); "v2h" and "v2e" are defunct as of 2.0.0.
  #
  # Internal routing keys off a single `.route` token: "c2h" (hospitalisation
  # or case events), "v2c" (vaccination events anchored on a date column of
  # data_x) or "v2e" (vaccination events anchored on one fixed Date).
  .legacy_codes  <- c("c2h", "v2c")   # deprecated: warn, then route as before
  .defunct_codes <- c("v2h", "v2e")   # defunct: error naming the replacement

  if (is.null(linkage_type)) {
    linkage_type <- "c2h"  # preserve historical default
  }

  if (identical(linkage_type, "cohort2event")) {
    if (is.null(event_type)) {
      stop(
        "(*)> starling::murmuration() \u2014 linkage_type = \"cohort2event\" requires ",
        "`event_type` to be supplied: one of \"vaccination\", \"hospitalisation\", ",
        "or \"case\".\n",
        "  cohort-to-vaccination  -> event_type = \"vaccination\"\n",
        "  cohort-to-hospital     -> event_type = \"hospitalisation\"\n",
        "  cohort-to-case         -> event_type = \"case\"",
        call. = FALSE
      )
    }
    event_type <- match.arg(
      event_type,
      c("vaccination", "hospitalisation", "hospitalization", "case", "cases")
    )
    # Map the collapsed (linkage_type, event_type) API onto the internal route
    # token. Any cohort (cases, hospitalisations, a birth cohort) links to
    # vaccination through the "v2c" route; the fixed-date-event case is
    # detected from a scalar Date event_date and uses the "v2e" route.
    .route <- switch(
      event_type,
      vaccination      = "v2c",   # cohort -> vaccination stream (event-relative or cohort window)
      hospitalisation  = "c2h",   # cohort -> hospitalisation stream
      hospitalization  = "c2h",
      case             = "c2h",   # cohort -> case/diagnosis stream (same routing as hosp events)
      cases            = "c2h"
    )
    if (event_type == "vaccination" && inherits(event_date, "Date")) {
      .route <- "v2e"
    }
  } else if (linkage_type %in% .defunct_codes) {
    stop(
      "(*)> starling::murmuration() \u2014 linkage_type = \"", linkage_type,
      "\" is defunct (deprecated since starling 1.1.0, removed in 2.0.0).\n",
      if (identical(linkage_type, "v2h"))
        "  Use linkage_type = \"cohort2event\", event_type = \"vaccination\", with the hospitalisations as data_x."
      else
        "  Use linkage_type = \"cohort2event\", event_type = \"vaccination\", event_date = as.Date(\"YYYY-MM-DD\").",
      call. = FALSE
    )
  } else if (linkage_type %in% .legacy_codes) {
    warning(
      "(*)> starling::murmuration() \u2014 linkage_type = \"", linkage_type,
      "\" is deprecated as of starling 1.1.0 and will be removed in a future ",
      "release.\n",
      "  Use linkage_type = \"cohort2event\" with an explicit event_type instead:\n",
      "    \"v2c\" -> linkage_type = \"cohort2event\", event_type = \"vaccination\"\n",
      "    \"c2h\" -> linkage_type = \"cohort2event\", event_type = \"hospitalisation\"\n",
      "  The legacy code still works for now and routes identically.",
      call. = FALSE
    )
    .route <- linkage_type
  } else {
    stop(
      "(*)> starling::murmuration() \u2014 unrecognised linkage_type = \"", linkage_type,
      "\".\n",
      "  Use linkage_type = \"cohort2event\" (with event_type), or one of the ",
      "deprecated legacy codes \"c2h\"/\"v2c\".",
      call. = FALSE
    )
  }

  # From here down, the historical internal variable name is retained so the
  # large per-type blocks below need no edits: linkage_type == the resolved route.
  linkage_type <- .route

  # --- Dynamic id_var.x / id_var.y column name resolution -------------------
  # reclin2::link() names the kept id columns "<id_var>.x" and "<id_var>.y"
  # after whatever string was passed as id_var — NOT literally "id_var.x".
  # Build these names dynamically so the function works for any id_var value
  # (e.g. id_var = "denovo_id" produces "denovo_id.x" / "denovo_id.y").
  if (missing(id_var) || is.null(id_var) || !is.character(id_var) || length(id_var) != 1) {
    stop(
      "(*)> starling::murmuration() \u2014 'id_var' must be supplied as a single character ",
      "string naming the unique identifier column in data_x (e.g. id_var = \"patient_id\")."
    )
  }
  id_x_col <- paste0(id_var, ".x")
  id_y_col <- paste0(id_var, ".y")

  # Dynamic name for the preserved "which df2/data_y record actually matched"
  # column (was the hardcoded literal "id_var_df2" before the data_x/data_y
  # rename -- now built the same dynamic way as id_x_col/id_y_col above, an
  # underscore rather than a dot specifically so clean_eggs's `ends_with(".y")`
  # drop doesn't catch it; see clean_eggs handling below).
  id_y_preserved_col <- paste0(id_var, "_y")

  if ("weights" %in% colnames(data_x)) {
    data_x <- data_x |> select(-weights)
  }
  if ("weights" %in% colnames(data_y)) {
    data_y <- data_y |> select(-weights)
  }
  if ("threshold" %in% colnames(data_x)) {
    data_x <- data_x |> select(-threshold)
  }
  if ("threshold" %in% colnames(data_y)) {
    data_y <- data_y |> select(-threshold)
  }

  # Baseline for prefix/suffix renaming at the very end of the function: any
  # column present in df, at the end, that ISN'T one of data_x's own original
  # columns (or an .x-suffixed copy of one) is "new" -- contributed by data_y or
  # computed by the linkage itself -- and is what prefix/suffix target. This
  # is captured here (post-cleanup, pre-everything-else) rather than derived
  # from a hardcoded list of known output column names, so it stays correct
  # regardless of linkage_type, method, or future changes to what the
  # event-window blocks below compute.
  .original_data_x_names <- names(data_x)

  # Columns only data_y carries (admission fields, dates of death, anything kept
  # via clean_the_nest(keep_vars = )). The hospitalisation route blanks ALL of
  # them when a linked admission falls outside the event window -- previously a
  # fixed list, so extra columns kept the out-of-window admission's values.
  .data_y_only_names <- setdiff(names(data_y), names(data_x))

  # ------------------------------------------------------------------
  # Resolve vaccination window
  #
  # Two modes:
  #   cohort_window supplied -> calendar-based window (entry/exit dates)
  #                             ignores vax_window entirely
  #   cohort_window NULL     -> event-relative window (vax_window list)
  #
  # In both cases we resolve to four concrete scalars used in the
  # per-linkage-type filter blocks below:
  #   .vax_days_before   : minimum days before event (exclusion buffer)
  #   .vax_days_after    : maximum days after event to retain
  #   .vax_lookback_days : maximum days before event to look back
  #   .vax_cohort_entry  : calendar entry date (Date; NA when unused)
  #   .vax_cohort_exit   : calendar exit  date (Date; NA when unused)
  # ------------------------------------------------------------------
  .use_cohort_window <- !is.null(cohort_window)

  if (.use_cohort_window) {

    if (!is.list(cohort_window) ||
        !all(c("entry", "exit") %in% names(cohort_window))) {
      stop(
        "(*)> starling::murmuration() \u2014 'cohort_window' must be a named list with ",
        "elements 'entry' and 'exit' (both Date objects or coercible strings).\n",
        "Example: cohort_window = list(entry = as.Date('2024-01-01'), ",
        "exit = as.Date('2024-12-31'))"
      )
    }

    .vax_cohort_entry <- as.Date(cohort_window$entry)
    .vax_cohort_exit  <- as.Date(cohort_window$exit)

    if (is.na(.vax_cohort_entry) || is.na(.vax_cohort_exit))
      stop(
        "(*)> starling::murmuration() \u2014 'cohort_window$entry' and 'cohort_window$exit' ",
        "must be valid dates. Check format (use as.Date() if passing a string)."
      )
    if (.vax_cohort_exit < .vax_cohort_entry)
      stop(
        "(*)> starling::murmuration() \u2014 'cohort_window$exit' (", .vax_cohort_exit,
        ") must be on or after 'cohort_window$entry' (", .vax_cohort_entry, ")."
      )

    # Unused in cohort-window mode but set to safe values to avoid
    # undefined-variable errors in the filter blocks below
    .vax_days_before   <- 0
    .vax_days_after    <- Inf
    .vax_lookback_days <- Inf

    message(
      "(*)> starling::murmuration() \u2014 cohort_window active: retaining vaccinations ",
      "between ", .vax_cohort_entry, " and ", .vax_cohort_exit, ".\n",
      "  vax_window is ignored."
    )

  } else {

    # Event-relative window from vax_window list
    if (!is.list(vax_window))
      stop(
        "(*)> starling::murmuration() \u2014 'vax_window' must be a named list.\n",
        "Example: vax_window = list(days_before = 14, days_after = 0, ",
        "lookback_days = 365)"
      )

    .vax_days_before   <- vax_window[["days_before"]]
    .vax_days_after    <- vax_window[["days_after"]]
    .vax_lookback_days <- vax_window[["lookback_days"]]

    # Apply defaults for any element not supplied
    if (is.null(.vax_days_before))   .vax_days_before   <- 14
    if (is.null(.vax_days_after))    .vax_days_after    <- 0
    if (is.null(.vax_lookback_days)) .vax_lookback_days <- Inf

    if (!is.numeric(.vax_days_before) || .vax_days_before < 0)
      stop(
        "(*)> starling::murmuration() \u2014 'vax_window$days_before' must be a ",
        "non-negative number. Got: ", .vax_days_before
      )
    if (!is.numeric(.vax_days_after) || .vax_days_after < 0)
      stop(
        "(*)> starling::murmuration() \u2014 'vax_window$days_after' must be a ",
        "non-negative number. Got: ", .vax_days_after
      )
    if (!is.numeric(.vax_lookback_days) || .vax_lookback_days <= 0)
      stop(
        "(*)> starling::murmuration() \u2014 'vax_window$lookback_days' must be a ",
        "positive number or Inf. Got: ", .vax_lookback_days
      )

    # Unused in event-relative mode. NA, not NULL: the vectorised
    # `.use_cohort_window & (d < entry | d > exit)` in the vaccination routes is
    # then FALSE for every row, whereas NULL gives a length-0 condition that
    # case_when() cannot recycle.
    .vax_cohort_entry <- as.Date(NA)
    .vax_cohort_exit  <- as.Date(NA)

    message(
      "(*)> starling::murmuration() \u2014 vax_window active:",
      "\n  days_before   = ", .vax_days_before,
      " (vaccines within ", .vax_days_before, " days of event excluded as too recent)",
      "\n  days_after    = ", .vax_days_after,
      if (.vax_days_after == 0) " (no post-event vaccinations retained)"
      else paste0(" (post-event vaccinations up to ", .vax_days_after, " days retained)"),
      "\n  lookback_days = ", .vax_lookback_days,
      if (is.infinite(.vax_lookback_days)) " (no upper lookback limit)"
      else paste0(" (vaccines older than ", .vax_lookback_days, " days before event excluded)")
    )
  }

  # Validate event_date parameter based on linkage_type.
  # For column-based linkage types, any date column present in data_x is valid -
  # including onset_date, admission_date, cohort_entry_date, and cohort_exit_date.
  # Skipped entirely when inspect_matches = TRUE: that mode stops right after
  # blocking + scoring and never reaches the event-window stage that actually
  # needs event_date, so there is no reason to force the caller to supply one.
  if (!isTRUE(inspect_matches)) {
    if (linkage_type %in% c("c2h", "v2c")) {
      if (is.null(event_date) || !is.character(event_date)) {
        stop(
          "(*)> starling::murmuration() \u2014 event_date must be a character string ",
          "naming a date column in data_x (e.g. 'onset_date', 'admission_date', ",
          "'cohort_entry_date', 'cohort_exit_date')."
        )
      }
      if (!event_date %in% colnames(data_x)) {
        stop(
          "(*)> starling::murmuration() \u2014 Column '", event_date, "' not found in data_x. ",
          "Run your data through clean_the_nest first, and ensure the column you are ",
          "passing as event_date was specified in that call."
        )
      }
    } else if (linkage_type == "v2e") {
      if (is.null(event_date) || !inherits(event_date, "Date")) {
        stop("(*)> starling::murmuration() \u2014 For a fixed-date event, event_date must be a Date object (e.g., as.Date('2024-12-15')).")
      }
    }
  }

  # ------------------------------------------------------------------
  # Linkage engine: resume-from-review, probabilistic (Fellegi-Sunter EM), OR
  # deterministic (exact key match). All three produce the SAME linked `df`
  # shape — an all_x left join carrying `weights` and a logical `threshold`
  # column — so every per-linkage-type block below consumes the result
  # identically and needs no awareness of which path produced it.
  # ------------------------------------------------------------------
  .resumed <- !is.null(resume_from_pairs)

  if (.resumed) {

    # --------------------------------------------------------------
    # Resume from a reviewed pairs table (see inspect_matches above). Blocking,
    # comparing, and EM-fitting are skipped entirely -- resume_from_pairs
    # already carries a `weights` column from a prior probabilistic run.
    # This is deliberately NOT a copy of the event-window logic further
    # below: it produces `df` via the exact same link() call as the normal
    # probabilistic branch, so every linkage_type block after this point
    # runs completely unchanged regardless of which branch built `df`.
    # --------------------------------------------------------------
    required_cols <- c(".x", ".y", "weights")
    missing_cols  <- setdiff(required_cols, names(resume_from_pairs))
    if (length(missing_cols) > 0) {
      stop(
        "(*)> starling::murmuration() \u2014 'resume_from_pairs' is missing required ",
        "column(s): ", paste(missing_cols, collapse = ", "), ".\n",
        "  Pass the object returned by murmuration(inspect_matches = TRUE), ",
        "optionally with 'match_override' edited.",
        call. = FALSE
      )
    }

    # .x/.y are row indices into the SAME data_x/data_y that produced them -- if
    # either has since changed size, they can no longer be trusted to point
    # at the right rows.
    .n_data_x <- attr(resume_from_pairs, "n_data_x", exact = TRUE)
    .n_data_y <- attr(resume_from_pairs, "n_data_y", exact = TRUE)
    if (!is.null(.n_data_x) && nrow(data_x) != .n_data_x) {
      stop(
        "(*)> starling::murmuration() \u2014 'data_x' has ", nrow(data_x), " rows, but ",
        "'resume_from_pairs' was built from a data_x with ", .n_data_x, " rows. Pass ",
        "the SAME data_x used to create it (via inspect_matches = TRUE) \u2014 .x/.y ",
        "are row indices into it.",
        call. = FALSE
      )
    }
    if (!is.null(.n_data_y) && nrow(data_y) != .n_data_y) {
      stop(
        "(*)> starling::murmuration() \u2014 'data_y' has ", nrow(data_y), " rows, but ",
        "'resume_from_pairs' was built from a data_y with ", .n_data_y, " rows. Pass ",
        "the SAME data_y used to create it (via inspect_matches = TRUE) \u2014 .x/.y ",
        "are row indices into it.",
        call. = FALSE
      )
    }

    pairs_pred  <- resume_from_pairs
    pairs_thresh <- pairs_pred

    # Final classification: a manual match_override wins where set; every
    # other pair falls back to weights >= threshold_value (which may be a
    # NEW value supplied here, not necessarily the one used when
    # inspect_matches originally built this table).
    if ("match_override" %in% names(pairs_thresh)) {
      bad_override <- setdiff(unique(pairs_thresh$match_override), c(NA, "match", "non-match"))
      if (length(bad_override) > 0) {
        stop(
          "(*)> starling::murmuration() \u2014 'match_override' contains invalid ",
          "value(s): ", paste(bad_override, collapse = ", "), ". Valid values ",
          "are NA, \"match\", and \"non-match\".",
          call. = FALSE
        )
      }
      n_overridden <- sum(!is.na(pairs_thresh$match_override))
      pairs_thresh$threshold <- ifelse(
        !is.na(pairs_thresh$match_override),
        pairs_thresh$match_override == "match",
        pairs_thresh$weights >= threshold_value
      )
    } else {
      n_overridden <- 0L
      pairs_thresh$threshold <- pairs_thresh$weights >= threshold_value
    }

    message(
      "(*)> starling::murmuration() \u2014 resuming from reviewed pairs:\n",
      sprintf("  Candidate pairs                : %d\n", nrow(pairs_thresh)),
      sprintf("  Manual overrides applied       : %d\n", n_overridden),
      sprintf("  Pairs selected (threshold=TRUE): %d\n", sum(pairs_thresh$threshold, na.rm = TRUE)),
      sprintf("  data_x has                        : %d records", nrow(data_x))
    )

    # data_x/data_y passed explicitly rather than read from attr(pairs, "x"/"y"):
    # a reviewed table edited with dplyr may have lost those attributes, and
    # these are the same frames its .x/.y indices point into (checked above).
    df <- link(pairs_thresh, selection = "threshold", x = data_x, y = data_y,
               keep_from_pairs = c("weights", "threshold"), all_x = TRUE, all_y = FALSE)

  } else if (method == "probabilistic") {

    pairs <- pair_blocking(data_x, data_y, blocking_var)

    if (identical(pairs$.x, pairs$.y) == TRUE) {
      stop("(*)> starling::murmuration() \u2014 Blocking variable too specific, matching 1:1. Choose a broader/looser blocking variable like postcode, or gender, rather than postcode + year of birth.")
    }

    compare_pairs(pairs, on = compare_vars, default_comparator = jaro_winkler(0.9), inplace = TRUE)

    m <- problink_em(reformulate(compare_vars), data = pairs)

    pairs_pred <- predict(m, pairs = pairs, add = TRUE)

    # Attach the blocking variable's realised value to every candidate pair,
    # pulled from data_x via the row index reclin2 stores in `.x`. pair_blocking()
    # only ever generates a pair when the blocking value already agrees between
    # data_x and data_y, so data_x's value is authoritative for both sides of the pair.
    # This lets murmuration_plot() jitter/position pairs by blocking stratum
    # without needing separate access to data_x -- see the `pairs_pred` attribute
    # on the returned data frame, below.
    pairs_pred[[blocking_var]] <- data_x[[blocking_var]][pairs_pred$.x]

    # ------------------------------------------------------------------
    # inspect_matches: stop here and hand back every candidate pair for
    # manual review, before threshold selection, link(), or any
    # linkage_type-specific event-window logic is applied.
    # ------------------------------------------------------------------
    # This exists because a single threshold_value is a blunt instrument --
    # some pairs just above it may be false matches (similar names, different
    # people) and some just below it may be true matches with a typo or a
    # missing field. Neither murmuration_plot() nor perch() can currently
    # show which, because pairs_pred carries no identifying information --
    # only .x/.y row indices, weights, and blocking_var. This block joins the
    # actual field values back in from both sides so a person can be
    # eyeballed, not just a row index.
    if (isTRUE(inspect_matches)) {

      if (length(inspect_extra_vars) > 0) {
        missing_vars <- setdiff(inspect_extra_vars, intersect(names(data_x), names(data_y)))
        if (length(missing_vars) > 0) {
          stop(
            "(*)> starling::murmuration() \u2014 inspect_extra_vars not found in both ",
            "data_x and data_y: ", paste(missing_vars, collapse = ", "),
            call. = FALSE
          )
        }
      }

      # id_var and compare_vars are already guaranteed to exist in both data
      # frames by this point (compare_pairs() above would already have
      # failed otherwise); inspect_extra_vars was just checked above.
      review_vars <- unique(c(id_var, compare_vars, inspect_extra_vars))

      for (v in review_vars) {
        pairs_pred[[paste0(v, ".x")]] <- data_x[[v]][pairs_pred$.x]
        pairs_pred[[paste0(v, ".y")]] <- data_y[[v]][pairs_pred$.y]
      }

      # A mechanical reference classification against threshold_value -- NOT
      # final. This is what murmuration_plot()'s "Above/Below threshold"
      # colouring will use unless overridden below.
      pairs_pred$threshold <- pairs_pred$weights >= threshold_value

      # The reviewer's column. NA (the default for every row) means "use the
      # threshold column as-is"; set individual rows to "match" or
      # "non-match" to override the mechanical cutoff after eyeballing a
      # pair's identifying fields (the .x/.y columns just added above).
      # Consumed by murmuration(resume_from_pairs = ...) via perch().
      pairs_pred$match_override <- NA_character_

      message(
        "(*)> starling::murmuration() \u2014 inspect_matches = TRUE: stopping after ",
        "blocking and scoring,\n",
        "  before threshold selection or any event-window filtering.\n",
        sprintf("  Candidate pairs                 : %d\n", nrow(pairs_pred)),
        sprintf("  Above threshold (%.1f)           : %d\n", threshold_value,
                sum(pairs_pred$threshold, na.rm = TRUE)),
        sprintf("  Identifier columns added (.x/.y): %s\n",
                paste(review_vars, collapse = ", ")),
        "  Inspect with murmuration_plot(); set 'match_override' to \"match\" or\n",
        "  \"non-match\" for specific pairs (NA = use the threshold column as-is).\n",
        "  Then finalise with murmuration(resume_from_pairs = ...)."
      )

      attr(pairs_pred, "blocking_var")       <- blocking_var
      attr(pairs_pred, "threshold_value")    <- threshold_value
      attr(pairs_pred, "murmuration_method") <- "probabilistic"
      attr(pairs_pred, "inspect_matches")    <- TRUE
      # So a later murmuration(resume_from_pairs = ...) call can confirm it's
      # being given the same data_x/data_y that .x/.y actually index into.
      attr(pairs_pred, "n_data_x") <- nrow(data_x)
      attr(pairs_pred, "n_data_y") <- nrow(data_y)

      # Return a plain data.table, not reclin2's "pairs" subclass: dplyr verbs
      # (mutate(), filter(), ...) restore extra attributes only for plain
      # data.frame/data.table objects, so on a "pairs" object they silently
      # drop blocking_var, threshold_value and the data_x/data_y references
      # that murmuration_plot() and resume_from_pairs rely on.
      class(pairs_pred) <- setdiff(class(pairs_pred), "pairs")

      return(pairs_pred)
    }

    if (isTRUE(perch_before_linking)) {
      perch(
        pairs_pred,
        n_records_primary = nrow(data_x),
        report = TRUE,
        plot   = interactive()
      )
    }

    pairs_thresh <- select_threshold(pairs_pred, "threshold", score = "weights", threshold = threshold_value)

    .threshold_note <- if (isTRUE(keep_pairs_pred)) {
      "use murmuration_plot() to review score distribution"
    } else {
      "re-run with keep_pairs_pred = TRUE to enable murmuration_plot()"
    }
    message(
      "(*)> starling::murmuration() \u2014 threshold selection complete:\n",
      sprintf("  Pairs before link              : %d\n", nrow(pairs_thresh)),
      sprintf("  Pairs selected (threshold=TRUE): %d\n", sum(pairs_thresh$threshold, na.rm = TRUE)),
      sprintf("  Unique data_x IDs in pairs        : %d\n", length(unique(pairs_thresh$.x))),
      sprintf("  data_x has                        : %d records\n", nrow(data_x)),
      sprintf("  Threshold applied              : %.1f  (%s)", threshold_value, .threshold_note)
    )

    df <- link(pairs_thresh, selection = "threshold", keep_from_pairs = c("weights", "threshold"), all_x = TRUE, all_y = FALSE)

  } else {  # method == "deterministic"

    df <- .murmuration_deterministic(
      data_x          = data_x,
      data_y          = data_y,
      compare_vars = compare_vars
    )

    # Deterministic linkage has no EM-scored candidate-pairs table -- weights
    # is NA throughout (see .murmuration_deterministic()) so there is nothing
    # meaningful for murmuration_plot() to draw. Explicitly NULL rather than
    # leaving `pairs_pred` undefined, so the attribute-assignment block below
    # (which runs for both methods) can rely on it always existing.
    pairs_pred <- NULL
  }


  message(
    "(*)> starling::murmuration() \u2014 link complete:\n",
    sprintf("  After link       : %d records\n", nrow(df)),
    sprintf("  Records with %s: %d\n", id_x_col, sum(!is.na(df[[id_x_col]]))),
    sprintf("  Records with %s: %d", id_y_col, sum(!is.na(df[[id_y_col]])))
  )

  if (linkage_type == "c2h") {
    df <- df |> filter(!is.na(.data[[id_x_col]]))
  }

  #___ Case to Hospitalization ___

  if (linkage_type == "c2h") {

    df <- df |>
      mutate(
        event_adm_diff = admission_date - .data[[event_date]],
        admission_date = case_when(
          threshold == FALSE ~ as.Date(NA),
          admission_date < .data[[event_date]] - days_allowed_before_event ~ as.Date(NA),
          admission_date > .data[[event_date]] + days_allowed_after_event ~ as.Date(NA),
          TRUE ~ admission_date
        ),
        across(c(ends_with(".y"), starts_with("admission"),
                 any_of(c("discharge_date", "event_adm_diff", "los", "icu",
                          "icu_date", "icu_hours", "hospital", "icd_code",
                          "diagnosis_description", "drg", "dod")),
                 any_of(setdiff(.data_y_only_names, id_var))),
               ~ case_when(is.na(admission_date) ~ NA, TRUE ~ .))
      )

    df <- df |>
      group_by(.data[[id_x_col]]) |>
      arrange(desc(threshold), admission_date, .by_group = TRUE) |>
      slice(1) |>
      ungroup() |>
      rename(!!id_y_preserved_col := !!id_y_col)

    df <- df |>
      mutate(
        total_admissions = if_else(!is.na(admission_date), 1L, 0L),
        adm_no = total_admissions
      )

    if (one_row_per_person == TRUE) {

      if ("event_adm_diff" %in% colnames(df)) {
        df <- df |> mutate(
          first_event_adm_diff = .data[["event_adm_diff"]],
          all_event_adm_diffs = as.character(.data[["event_adm_diff"]])
        )
      }

      if ("admission_date" %in% colnames(df)) {
        df <- df |> mutate(
          first_admission_date = .data[["admission_date"]],
          last_admission_date = .data[["admission_date"]],
          all_admission_dates = as.character(.data[["admission_date"]])
        )
      }

      if ("discharge_date" %in% colnames(df)) {
        df <- df |> mutate(
          first_discharge_date = .data[["discharge_date"]],
          last_discharge_date = .data[["discharge_date"]],
          all_discharge_dates = as.character(.data[["discharge_date"]])
        )
      }

      if ("icu_date" %in% colnames(df)) {
        df <- df |> mutate(first_icu_date = .data[["icu_date"]])
      }

      if ("icu_hours" %in% colnames(df)) {
        df <- df |> mutate(total_icu_hours = .data[["icu_hours"]])
      }

      if ("los" %in% colnames(df)) {
        df <- df |> mutate(
          total_los = as.numeric(.data[["los"]]),
          all_los = as.character(.data[["los"]])
        )
      }

      if ("hospital" %in% colnames(df)) {
        df <- df |> mutate(
          first_hospital = .data[["hospital"]],
          all_hospitals = as.character(.data[["hospital"]])
        )
      }

      if ("icd_code" %in% colnames(df)) {
        df <- df |> mutate(
          first_icd_code = .data[["icd_code"]],
          all_icd_codes = as.character(.data[["icd_code"]])
        )
      }

      if ("drg" %in% colnames(df)) {
        df <- df |> mutate(
          first_drg = .data[["drg"]],
          all_drgs = as.character(.data[["drg"]])
        )
      }

      if ("diagnosis_description" %in% colnames(df)) {
        df <- df |> mutate(
          first_diagnosis_description = .data[["diagnosis_description"]],
          all_diag_desc = as.character(.data[["diagnosis_description"]])
        )
      }

      if ("dialysis_outcome" %in% colnames(df)) {
        df <- df |> mutate(dialysis_outcome = .data[["dialysis"]])
      }

      if (!is.null(last_follow_up)) {
        df <- df |> mutate(
          survtime = case_when(
            !is.na(dod) ~ as.numeric(dod - .data[[event_date]]),
            TRUE ~ as.numeric(last_follow_up - .data[[event_date]])
          )
        )
      }

      df <- df |>
        mutate(
          admission_outcome = factor(
            if_else(!is.na(first_admission_date), 1L, 0L),
            levels = c(0, 1),
            labels = c("No Admission", "Admission"),
            ordered = TRUE
          )
        )

      if ("icu_date" %in% colnames(df)) {
        df <- df |> mutate(
          icu_outcome = factor(
            if_else(!is.na(.data[["icu_date"]]), 1L, 0L),
            levels = c(0, 1),
            labels = c("No ICU Admission", "ICU Admission"),
            ordered = TRUE
          )
        )
      } else if ("total_icu_hours" %in% colnames(df)) {
        df <- df |> mutate(
          icu_outcome = factor(
            case_when(
              .data[["total_icu_hours"]] > 0 ~ 1L,
              TRUE ~ 0L
            ),
            levels = c(0, 1),
            labels = c("No ICU Admission", "ICU Admission"),
            ordered = TRUE
          )
        )
      }

      if ("dod" %in% colnames(df)) {
        df <- df |> mutate(
          death_outcome = factor(
            case_when(
              !is.na(.data[["dod"]]) & .data[["dod"]] <= .data[[event_date]] + days_between_onset_death ~ 1L,
              TRUE ~ 0L
            ),
            levels = c(0, 1),
            labels = c("Alive", "Died"),
            ordered = TRUE
          )
        )
      }

      df <- df |> select(-any_of(c("adm_no", "admission_date", "event_adm_diff",
                                   "los", "hospital", "icd_code",
                                   "diagnosis_description", "drg")))
    }
  }  # End c2h


  #___ Vaccine to Case ___

  if (linkage_type == "v2c") {

    df <- df |>
      mutate(
        across(starts_with('vax_date_'), ~ {
          d <- .
          ev <- .data[[event_date]]
          case_when(
            # Not a match
            threshold == FALSE                              ~ as.Date(NA),
            # Cohort-window mode: outside calendar bounds. Vectorised & (not
            # &&: d is a whole column, and R >= 4.3 errors on && with a
            # vector). The bounds are NA, not NULL, in event-relative mode,
            # so this is FALSE there rather than logical(0).
            .use_cohort_window & (d < .vax_cohort_entry |
                                    d > .vax_cohort_exit)   ~ as.Date(NA),
            # Event-relative mode: inside the pre-event buffer (on or before
            # the event, but fewer than days_before days before it). Doses
            # after the event are left to the days_after rule below --
            # previously `d + days_before > ev` also caught every post-event
            # dose, so days_after could never keep one.
            !.use_cohort_window &
              !is.na(ev) &
              d <= ev & d > ev - .vax_days_before         ~ as.Date(NA),
            # Event-relative mode: beyond lookback limit
            !.use_cohort_window &
              !is.na(ev) & is.finite(.vax_lookback_days) &
              d < ev - .vax_lookback_days                 ~ as.Date(NA),
            # Event-relative mode: after the post-event window
            !.use_cohort_window &
              !is.na(ev) &
              d > ev + .vax_days_after                    ~ as.Date(NA),
            TRUE ~ d
          )
        }),
        across(c(starts_with('vax_type_'), ends_with(".y")), ~ case_when(
          threshold == FALSE ~ NA,
          TRUE ~ .
        ))
      )

    df <- df |>
      group_by(.data[[id_x_col]]) |>
      # Highest-scoring candidate first: household members share surname,
      # postcode and Medicare card, so several can clear the threshold;
      # previously slice(1) kept whichever came first in data_y order.
      arrange(desc(threshold), desc(weights), .by_group = TRUE) |>
      slice(1) |>
      ungroup()

    # Preserve data_y's own matched-record id under the same dynamic name
    # c2h already uses (id_y_preserved_col), renamed BEFORE clean_eggs's global
    # `select(-ends_with(".y"))` step runs -- without this, clean_eggs = TRUE
    # (the default) silently drops the matched vaccination record's own id
    # entirely, with nothing else in the output tying a match back to its
    # specific data_y source row. Matches c2h's own equivalent step exactly.
    df <- df |> rename(!!id_y_preserved_col := !!id_y_col)

    date_cols <- grep("^vax_date_\\d+$", names(df), value = TRUE)
    type_cols <- grep("^vax_type_\\d+$", names(df), value = TRUE)

    for (dc in date_cols) {
      x <- str_extract(dc, "\\d+$")
      tc <- paste0("vax_type_", x)
      if (tc %in% type_cols) {
        df <- df |>
          mutate(!!tc := if_else(is.na(.data[[dc]]), NA_character_, .data[[tc]]))
      }
    }

    df <- df |>
      rowwise() |> mutate(
        dates = list(c_across(all_of(date_cols))),
        types = list(c_across(all_of(type_cols))),
        total_valid_vax = sum(!is.na(dates)),
        first_vax_date = if (total_valid_vax > 0) min(dates, na.rm = TRUE) else as.Date(NA),
        first_vax_type = if (total_valid_vax > 0) types[which.min(dates)] else NA_character_,
        last_vax_date  = if (total_valid_vax > 0) max(dates, na.rm = TRUE) else as.Date(NA),
        last_vax_type  = if (total_valid_vax > 0) types[which.max(dates)] else NA_character_,
        vaccination_status     = if (total_valid_vax > 0) "Vaccinated" else "Unvaccinated",
        vaccination_status_num = if_else(vaccination_status == "Unvaccinated", 0, 1),
        tsv = if (total_valid_vax > 0) as.numeric(.data[[event_date]] - last_vax_date) else 0
      ) |>
      select(-dates, -types) |>
      ungroup()
  }


  # (The former "v2h" block was identical to the vaccination route above apart
  # from its closing columns; "v2h" is defunct and that block is removed.)


  #___ Vaccine to Event ___

  if (linkage_type == "v2e") {

    df <- df |>
      mutate(
        across(starts_with('vax_date_'), ~ {
          d <- .
          ev <- event_date   # v2e: event_date is a scalar Date object
          case_when(
            threshold == FALSE                              ~ as.Date(NA),
            # Cohort-window mode: outside calendar bounds (vectorised; see the
            # v2c block above)
            .use_cohort_window & (d < .vax_cohort_entry |
                                    d > .vax_cohort_exit)   ~ as.Date(NA),
            # Inside the pre-event buffer only; post-event doses are left to
            # the days_after rule (see the v2c block above)
            !.use_cohort_window &
              d <= ev & d > ev - .vax_days_before         ~ as.Date(NA),
            !.use_cohort_window &
              is.finite(.vax_lookback_days) &
              d < ev - .vax_lookback_days                 ~ as.Date(NA),
            !.use_cohort_window &
              d > ev + .vax_days_after                    ~ as.Date(NA),
            TRUE ~ d
          )
        }),
        across(c(starts_with('vax_type_'), ends_with(".y")), ~ case_when(
          threshold == FALSE ~ NA,
          TRUE ~ .
        ))
      )

    df <- df |>
      group_by(.data[[id_x_col]]) |>
      # Highest-scoring candidate first: household members share surname,
      # postcode and Medicare card, so several can clear the threshold;
      # previously slice(1) kept whichever came first in data_y order.
      arrange(desc(threshold), desc(weights), .by_group = TRUE) |>
      slice(1) |>
      ungroup()

    # Preserve data_y's own matched-record id under the same dynamic name
    # c2h already uses (id_y_preserved_col), renamed BEFORE clean_eggs's global
    # `select(-ends_with(".y"))` step runs -- without this, clean_eggs = TRUE
    # (the default) silently drops the matched vaccination record's own id
    # entirely, with nothing else in the output tying a match back to its
    # specific data_y source row. Matches c2h's own equivalent step exactly.
    df <- df |> rename(!!id_y_preserved_col := !!id_y_col)

    date_cols <- grep("^vax_date_\\d+$", names(df), value = TRUE)
    type_cols <- grep("^vax_type_\\d+$", names(df), value = TRUE)

    for (dc in date_cols) {
      x <- str_extract(dc, "\\d+$")
      tc <- paste0("vax_type_", x)
      if (tc %in% type_cols) {
        df <- df |>
          mutate(!!tc := if_else(is.na(.data[[dc]]), NA_character_, .data[[tc]]))
      }
    }

    df <- df |>
      rowwise() |> mutate(
        dates = list(c_across(all_of(date_cols))),
        types = list(c_across(all_of(type_cols))),
        total_valid_vax = sum(!is.na(dates)),
        first_vax_date = if (total_valid_vax > 0) min(dates, na.rm = TRUE) else as.Date(NA),
        first_vax_type = if (total_valid_vax > 0) types[which.min(dates)] else NA_character_,
        last_vax_date  = if (total_valid_vax > 0) max(dates, na.rm = TRUE) else as.Date(NA),
        last_vax_type  = if (total_valid_vax > 0) types[which.max(dates)] else NA_character_,
        vaccination_status     = if (total_valid_vax > 0) "Vaccinated" else "Unvaccinated",
        vaccination_status_num = if_else(vaccination_status == "Unvaccinated", 0, 1),
        tsv        = if (total_valid_vax > 0) as.numeric(event_date - last_vax_date) else 0,
        event_date = event_date
      ) |>
      select(-dates, -types) |>
      ungroup()
  }


  # Clean up column names
  if (clean_eggs == TRUE) {
    df <- df |> select(-ends_with(".y"))
    df <- df |> rename_with(~str_remove(., '[.]x$'))
  }

  # Reorder columns
  df <- df |>
    select(
      any_of(c("lettername1_lettername2_dob", "lettername1_lettername2_dob.x",
               event_date, "admission_outcome", "first_admission_date",
               "first_event_adm_diff", "icu_outcome", "death_outcome", "dod",
               "diagnosis", "drg", "icd_code", id_var, id_y_preserved_col,
               "weights", "threshold")),
      ends_with("_vax"),
      starts_with('vax_date_'),
      starts_with('vax_type'),
      starts_with('all_'),
      everything()
    )

  # --------------------------------------------------------------------
  # Apply prefix/suffix to newly-contributed columns (from data_y, or computed
  # by the linkage itself) -- NOT data_x's own passthrough columns, including
  # id_var, which stays plain and stable so it remains a predictable join
  # key across multiple murmuration() calls onto the same origin dataset.
  # Built for exactly this workflow: linking one cohort against several
  # separate event-stream datasets (e.g. flu/covid/rsv/shingles vaccination
  # registries) in independent calls, then left_join()-ing the results back
  # onto the origin dataset by id_var without every call's "weights",
  # "vaccination_status", etc. colliding under the same generic names.
  # --------------------------------------------------------------------
  # A column is "new" here if it isn't one of data_x's original names (captured
  # in .original_data_x_names before any processing began), and isn't an
  # .x-suffixed copy of one either (relevant only when clean_eggs = FALSE,
  # where an overlapping data_x column survives as "name.x" rather than being
  # renamed back to "name"). This is deliberately not a hardcoded list of
  # known output columns (weights, threshold, total_vax, first_vax_date,
  # ...): it stays correct automatically as linkage_type-specific processing
  # changes, rather than needing a matching update every time it does.
  if (!is.null(prefix) || !is.null(suffix)) {
    current_names   <- names(df)
    is_own_x_suffix <- grepl("\\.x$", current_names) &
      sub("\\.x$", "", current_names) %in% .original_data_x_names
    is_new_col      <- !(current_names %in% .original_data_x_names) & !is_own_x_suffix

    new_names <- current_names
    if (!is.null(prefix)) new_names[is_new_col] <- paste0(prefix, new_names[is_new_col])
    if (!is.null(suffix)) new_names[is_new_col] <- paste0(new_names[is_new_col], suffix)
    names(df) <- new_names
  }

  # --------------------------------------------------------------------
  # Attach pairs-level scoring output as metadata (non-breaking)
  # --------------------------------------------------------------------
  # Per the ecosystem's existing attribute-based side-channel pattern (see
  # mudnester::preening()'s attr(x, "age_scheme_source")), the full
  # candidate-pairs table travels with the returned data frame as an
  # attribute rather than changing murmuration()'s return type. Every
  # existing caller that pipes murmuration()'s output into dplyr,
  # mudnester::preening(), or a further murmuration() call keeps working
  # completely unchanged.
  #
  # CHANGED: keep_pairs_pred now defaults to FALSE -- OPT-IN, not opt-out.
  # pairs_pred is the FULL candidate-pairs table (every pair blocking
  # produced, before thresholding) -- confirmed from a real run to reach
  # into the tens of millions of rows for a single call against a large
  # extract, riding along as a hidden attribute on an otherwise small
  # result. It was previously attached unconditionally on every
  # probabilistic/resumed call regardless of whether anything downstream
  # used it, which is what made it invisible: it doesn't show up in
  # str()/glimpse() of the columns, doesn't change nrow()/ncol(), but does
  # get serialised in full by saveRDS() and does sit in memory for as long
  # as the result object does. Callers that want murmuration_plot()'s score
  # distribution (or attr(x, "pairs_pred") directly) need to pass
  # keep_pairs_pred = TRUE explicitly on that call -- there is no other way
  # to get it now, by design. inspect_matches = TRUE is unaffected either
  # way: that path already returns the full pairs table directly as the
  # function's return value (see above), an already-explicit request, not
  # governed by this argument.
  #
  # attr(df, "pairs_pred") is NULL under method = "deterministic" regardless
  # of keep_pairs_pred -- there is no Fellegi-Sunter score to attach (see
  # the deterministic branch above). murmuration_plot() checks for this and
  # errors informatively rather than drawing a meaningless single-point
  # "distribution" -- the same error it now also gives when keep_pairs_pred
  # was left at its default FALSE.
  #
  # NOTE: attributes set this way are not guaranteed to survive every dplyr
  # verb applied to `df` afterwards (dplyr does not promise to preserve
  # arbitrary custom attributes -- see the same caveat already noted for
  # roost_df metadata, section 8.3 of the ecosystem plan). Call murmuration_plot()
  # on murmuration()'s output directly, or extract what you need immediately:
  #   linked <- murmuration(..., keep_pairs_pred = TRUE)
  #   pp     <- attr(linked, "pairs_pred")   # safe to keep even if `linked`
  #                                          # is later filtered/mutated
  attr(df, "pairs_pred")         <- if (isTRUE(keep_pairs_pred)) pairs_pred else NULL
  attr(df, "blocking_var")       <- if (.resumed || method == "probabilistic") blocking_var else NA_character_
  attr(df, "threshold_value")    <- if (.resumed || method == "probabilistic") threshold_value else NA_real_
  attr(df, "murmuration_method") <- if (.resumed) "probabilistic" else method
  attr(df, "resumed_from_review") <- .resumed

  return(df)
}


# ==============================================================================
# Internal: deterministic exact-key linkage engine
# ==============================================================================

#' Deterministic exact-key linkage (internal)
#'
#' Produces the same linked-data-frame shape as the reclin2 probabilistic
#' pipeline (`pair_blocking` -> ... -> `link(all_x = TRUE, all_y = FALSE)`),
#' but by exact matching on `compare_vars` rather than EM-scored comparison.
#' This lets \code{murmuration(method = "deterministic")} route into every
#' existing per-linkage-type block unchanged: those blocks only require a
#' left-joined `df` carrying a logical `threshold` column and `.x`/`.y`
#' suffixed shared columns.
#'
#' Design decisions (mirroring the plan's \\S4.1 deterministic spec):
#' \itemize{
#'   \item A pair links iff \emph{every} column in \code{compare_vars} matches
#'     exactly (all-TRUE). No fuzzy comparison, no partial credit — that is
#'     what \code{method = "probabilistic"} is for.
#'   \item \code{weights} is set to \code{NA_real_}: deterministic linkage has
#'     no Fellegi-Sunter composite score to report. Downstream code that
#'     filters on \code{weights} must be NA-aware.
#'   \item \code{threshold} is \code{TRUE} for every row that found an exact
#'     match and \code{FALSE} for unmatched data_x rows (kept via the left join),
#'     exactly as the probabilistic path sets it after \code{select_threshold()}.
#'   \item Shared column names get \code{.x}/\code{.y} suffixes; data_y-only
#'     columns (e.g. \code{admission_date}, \code{vax_date_1}) stay unsuffixed,
#'     matching reclin2's behaviour so the per-type blocks find them by bare name.
#' }
#'
#' @param data_x,data_y Cleaned input frames (base and event stream).
#' @param compare_vars Character vector of columns that must match exactly.
#' @return A data frame: all data_x rows left-joined to matching data_y rows, with
#'   `.x`/`.y` suffixes on shared columns, plus `weights` (NA) and `threshold`.
#' @keywords internal
#' @noRd
.murmuration_deterministic <- function(data_x, data_y, compare_vars) {

  # --- Validate that the exact-match key exists in both frames ---------------
  missing_a <- setdiff(compare_vars, names(data_x))
  missing_b <- setdiff(compare_vars, names(data_y))
  if (length(missing_a) > 0 || length(missing_b) > 0) {
    stop(
      "(*)> starling::murmuration() \u2014 method = \"deterministic\" requires every ",
      "column in `compare_vars` to be present in BOTH datasets to form the exact ",
      "match key.\n",
      if (length(missing_a) > 0)
        paste0("  Missing from data_x: ", paste(missing_a, collapse = ", "), "\n") else "",
      if (length(missing_b) > 0)
        paste0("  Missing from data_y: ", paste(missing_b, collapse = ", "), "\n") else "",
      "  Deterministic linkage cannot fall back to fuzzy comparison the way ",
      "method = \"probabilistic\" can. Either supply a key present in both, or ",
      "use method = \"probabilistic\".",
      call. = FALSE
    )
  }

  # --- Warn on duplicate keys (data-quality condition, not silently resolved)-
  # An identical key appearing more than once in data_y will many-to-many expand
  # the join. We surface it rather than silently pick one — the per-type blocks
  # below already do their own one_row_per_person collapsing, but the analyst
  # should know the key is not unique.
  key_b <- do.call(paste, c(data_y[compare_vars], sep = "\r"))
  dup_b <- sum(duplicated(key_b))
  if (dup_b > 0) {
    message(
      "(*)> starling::murmuration() \u2014 deterministic key is not unique in data_y: ",
      dup_b, " row(s) share a key with another. These will produce multiple ",
      "candidate links per person; downstream one_row_per_person handling ",
      "collapses them, but review if unexpected."
    )
  }

  # --- Exact left join (all data_x rows kept), reclin2-compatible suffixes ------
  message(
    "(*)> starling::murmuration() \u2014 deterministic linkage:\n",
    sprintf("  Exact-match key: %s\n", paste(compare_vars, collapse = " + ")),
    sprintf("  data_x has        : %d records", nrow(data_x))
  )

  # Add a sentinel presence column to data_y BEFORE the join so match detection
  # does not depend on how dplyr suffixes the id column. reclin2::link() always
  # keeps <id_var>.x/.y, but a dplyr join only suffixes columns present in BOTH
  # frames — so keying "was a data_y partner found?" off id_y_col is fragile if
  # data_y's id column has a different name. A dedicated sentinel is robust.
  data_y <- data_y
  data_y[[".starling_det_matched"]] <- TRUE

  df <- dplyr::left_join(
    data_x, data_y,
    by     = compare_vars,
    suffix = c(".x", ".y"),
    relationship = "many-to-many"
  )

  # threshold is TRUE where a data_y partner was found (sentinel present after the
  # left join), FALSE for unmatched data_x rows — mirroring the probabilistic
  # path's post-select_threshold() logical column. weights is NA: deterministic
  # linkage has no Fellegi-Sunter composite score.
  df <- df |>
    dplyr::mutate(
      weights   = NA_real_,
      threshold = !is.na(.data[[".starling_det_matched"]])
    )
  df[[".starling_det_matched"]] <- NULL

  df
}
