# Suppress R CMD check notes for unquoted column names in ggplot2 aes() mappings
utils::globalVariables(c("weights", "match_grp", "block", "n_pairs", "text"))

#' Visualise Linkage Weight Distribution Across Blocking Strata
#'
#' @description
#' **Bird note**: When a murmuration peaks, a watcher on the ground sees the flock
#' split momentarily into two layers — the dense, high-altitude core of confirmed
#' companions and the looser, lower-flying fringe of uncertain travellers — before
#' the whole shape resolves again. \code{murmuration_plot()} captures exactly that
#' split, and now shows whether it holds up flock by flock: the linkage weight
#' distribution separates into a high-score cluster of likely true matches, a
#' low-score cluster of likely non-matches, and a threshold line between them —
#' broken out by blocking stratum so a threshold that looks fine pooled across the
#' whole dataset can be checked for strata where it quietly performs worse.
#'
#' Generates a jittered scatter plot of Fellegi-Sunter linkage weights, positioned
#' along the x-axis by the \code{blocking_var} stratum each candidate pair belongs to
#' (jittered within each stratum for readability), with a horizontal line marking the
#' chosen threshold and shaded regions distinguishing likely matches from likely
#' non-matches.
#'
#' @details
#' ## Running directly off murmuration()'s own output
#' As of starling 1.2.0, \code{murmuration(method = "probabilistic")} attaches its
#' full candidate-pairs scoring table to its own return value as
#' \code{attr(x, "pairs_pred")}, along with the \code{blocking_var} and
#' \code{threshold_value} actually used. \code{murmuration_plot()} reads all three
#' straight off that object:
#' \preformatted{
#' linked <- murmuration(data_x, data_y, ..., blocking_var = "gender", threshold_value = 16)
#' murmuration_plot(linked)                 # threshold and blocking_var picked up
#'                                           # automatically
#' murmuration_plot(linked, threshold = 14) # explore a different cutoff without
#'                                           # re-running murmuration()
#' }
#' There is no need to hand-roll the underlying pairs-scoring pipeline yourself just
#' to see the score distribution. \code{data} can also be a bare pairs-level data
#' frame with a \code{weights} column directly (a saved \code{attr(x, "pairs_pred")},
#' or output from the legacy raw pairs-scoring workflow) — see \code{data} below.
#'
#' \code{method = "deterministic"} linkage has no Fellegi-Sunter score to plot
#' (\code{weights} is \code{NA} throughout by design — see \code{\link{murmuration}}):
#' passing its output here is a \code{stop()}, not a silent empty plot.
#'
#' ## How points are sampled
#' With large datasets, plotting every candidate pair is slow and not usually any
#' more informative than a representative sample — but a sample taken uniformly
#' across the whole dataset (or even just within each blocking stratum) can be
#' actively misleading here: above-threshold pairs are often a tiny fraction of all
#' candidate pairs — well under 1\% is common — and a uniform sample can easily miss
#' every single one of them, making a real, if small, separation look like no
#' separation at all. Sampling is therefore asymmetric:
#' \itemize{
#'   \item Every above-threshold pair is kept, in full, no matter how rare —
#'     unless there are more of them than \code{n_sample} itself, an unusual case
#'     in which they are themselves subsampled down to \code{n_sample}.
#'   \item The (usually much larger) below-threshold population is subsampled to
#'     fill whatever budget remains, proportionally \emph{within each blocking
#'     stratum} so small strata aren't crowded out by large ones.
#' }
#' When this asymmetric sampling is active, the plot subtitle says so explicitly
#' ("every above-threshold pair shown; below-threshold sampled") when
#' \code{show_subtitle = TRUE}. Above-threshold points are also drawn larger, more
#' opaque, and as a different shape (a triangle rather than a circle) than
#' below-threshold points, so they stay visible even when very few of them exist
#' against a large background cloud. The header counts (n above/below threshold) are
#' always computed from the complete, unsampled data — only the scatter of
#' individual points is affected by \code{n_sample}.
#'
#' ## Reading the per-stratum layout
#' Strata are ordered left to right by how many above-threshold pairs they
#' contain — most first, tapering to zero — rather than by raw size, so where
#' matches actually concentrate is visible at a glance instead of requiring a
#' scan across every column. Ties (typically among strata with no
#' above-threshold pairs at all) are broken by total pair volume, largest
#' first. "Other" is always last, regardless of its own aggregate count, since
#' it is a heterogeneous pool rather than a single meaningful stratum. If
#' \code{max_levels} forces some strata to be pooled, any stratum containing
#' at least one above-threshold pair is kept individually in preference to a
#' larger zero-match stratum — see \code{max_levels}.
#'
#' Each blocking stratum (e.g. each level of \code{gender}, or each \code{postcode}
#' block) gets its own jittered column of points, all sharing the same y-axis and
#' threshold line. A well-behaved threshold shows the same match/non-match separation
#' in every column; a stratum where points cluster near the threshold line, or where
#' the "likely non-match" cloud creeps above it, is a sign that stratum's identifiers
#' are noisier and may need a different \code{blocking_var}, additional
#' \code{compare_vars}, or a stratum-specific threshold applied via a separate
#' \code{murmuration()} call and a manual \code{dplyr::bind_rows()}.
#'
#' ## Threshold guidance for Australian public health linkage
#' The Fellegi-Sunter weights produced by the EM algorithm are log-likelihood
#' ratios, not fixed absolute scores, so the "right" threshold is
#' dataset-specific. That said, common reference points for Australian population
#' health work (Medicare, 2 names, DOB) are:
#'
#' | Threshold | Practical meaning |
#' |---|---|
#' | 8–12 | High sensitivity, lower specificity — suitable when recall matters more than precision (e.g. small datasets, rare disease surveillance) |
#' | 12–17 | Balanced — typical starting point for surveillance linkage with a complete variable set (the \code{murmuration()} default, 16, assumes 5 comparison variables with an agreeing Medicare number) |
#' | 17–22 | High specificity — suitable when false matches are especially costly (e.g. VE studies where false vaccination attribution biases the estimate) |
#' | > 22 | Very conservative — consider adding clerical review for records in the 17–22 range |
#'
#' The most principled approach is to look at the weight distribution (this plot)
#' and choose the threshold at the **lowest density point** between the two modes,
#' checked across every blocking stratum in turn rather than only in the pooled
#' picture. If the two modes overlap substantially, your variable set may lack
#' discriminating power for this dataset — see \code{\link{preflight}} for
#' diagnostic guidance.
#'
#' @param data The object returned by \code{\link{murmuration}} (preferred) — its
#'   attached \code{attr(x, "pairs_pred")} is used automatically. Also accepts a bare
#'   pairs-level data frame with a \code{weight_col} column directly, for a saved
#'   \code{pairs_pred} object or legacy pairs-scoring output produced outside
#'   \code{murmuration()}.
#' @param threshold Numeric or \code{NULL} (default). \code{NULL} resolves to
#'   \code{attr(data, "threshold_value")} — whatever \code{threshold_value} was
#'   actually passed to \code{murmuration()} — falling back to \code{16} (the
#'   \code{murmuration()} default) with a
#'   \code{message()} if that attribute isn't available (e.g. a bare pairs-level data
#'   frame was supplied).
#' @param blocking_var Character or \code{NULL} (default). Name of the column giving
#'   each pair's blocking stratum, used to position points along the x-axis (jittered
#'   within each stratum). \code{NULL} resolves to \code{attr(data, "blocking_var")};
#'   if that's unavailable either, falls back to a single unstratified column with a
#'   \code{message()} explaining how to enable per-stratum jittering.
#' @param max_levels Integer. If more blocking strata are present than this,
#'   \code{max_levels - 1} are shown individually and the rest are pooled into an
#'   "Other" column, with a \code{warning()}. Any stratum containing at least one
#'   above-threshold pair is kept individually first; remaining slots go to the
#'   highest-volume zero-match strata — see Details, "Reading the per-stratum
#'   layout". This is consistent with the ecosystem's facet-count convention (see
#'   \code{bowerbird::roost_plot(facet_by = )}). Default \code{12}.
#' @param weight_col Name of the column containing linkage weights. Default
#'   \code{"weights"}.
#' @param n_sample Integer. Maximum number of points to plot, in total. Above-threshold
#'   pairs are always kept in full; below-threshold pairs are subsampled
#'   proportionally \emph{within each blocking stratum} to fill whatever budget
#'   remains, so small strata aren't swamped out by large ones — see Details, "How
#'   points are sampled". Plotting millions of points is slow and uninformative; the
#'   below-threshold distribution shape is well-represented by 5,000–10,000 points in
#'   total. Default \code{5000}.
#' @param jitter_width Numeric. Horizontal jitter width within each blocking stratum's
#'   column. Default \code{0.3}.
#' @param point_alpha Numeric (0–1). Point transparency. Default \code{0.35}.
#' @param point_size Numeric. Point size. Default \code{1.2}.
#' @param show_subtitle Logical. Show the subtitle line (threshold, above/below
#'   counts, and the sampling note when applicable)? Default \code{TRUE}.
#' @param show_counts Logical. Label each blocking stratum's column with its pair
#'   count (\code{n = ...})? Default \code{TRUE}.
#' @param palette Character. Colour palette for match/non-match zones:
#'   \code{"sch"} (Sunshine Coast HHS greens — default), \code{"default"}
#'   (teal/coral), or \code{"grey"} (greyscale, publication-safe).
#' @param title Optional plot title string. \code{NULL} uses a default title.
#' @param interactive Logical. Return a \pkg{plotly} object (\code{TRUE}) or a
#'   \pkg{ggplot2} object (\code{FALSE}, default).
#' @param pairs_pred \code{Deprecated}. Pre-1.2.0 name for \code{data}, when it took
#'   only a raw pairs-level object. Positional calls (\code{murmuration_plot(pp, ...)})
#'   are unaffected by the rename and never need this — it exists solely to catch old
#'   calls that named the argument explicitly (\code{pairs_pred = pp}), which emit a
#'   \code{warning()} and are routed to \code{data} for now. Will be removed in a
#'   future release; update calls to use \code{data} instead.
#'
#' @return A \pkg{ggplot2} object (\code{interactive = FALSE}) or a \pkg{plotly}
#'   htmlwidget (\code{interactive = TRUE}).
#'
#' @examples
#' \donttest{
#' # Minimal, fully self-contained example: a bare pairs-level data frame
#' # with just a weights column (e.g. a saved attr(x, "pairs_pred")) -- the
#' # "Also accepts a bare pairs-level data frame" path documented in @param
#' # data above -- so this is a genuine, runnable demonstration that needs
#' # neither murmuration() nor reclin2 to have actually run.
#' pairs_toy <- data.frame(weights = c(2, 4, 5, 18, 21, 25, 3, 19, 11, 16))
#' murmuration_plot(pairs_toy, threshold = 16)
#' }
#'
#' \dontrun{
#' # The following illustrate running directly off murmuration()'s own
#' # output -- kept as \dontrun{} rather than \donttest{} because they need
#' # a real reclin2-scored pairs table (via murmuration(method =
#' # "probabilistic")) to show per-stratum jittering, threshold exploration,
#' # and asymmetric sampling meaningfully; a handful of synthetic rows would
#' # not exercise any of that. The \donttest{} example above already
#' # demonstrates the core plotting function runs correctly end to end.
#'
#' # Preferred workflow: pass murmuration()'s own output straight in
#' linked <- murmuration(dx_clean, vax_clean,
#'   linkage_type = "cohort2event", event_type = "vaccination",
#'   event_date = "onset_date", id_var = "identity",
#'   blocking_var = "gender",
#'   compare_vars = c("lettername1", "lettername2", "dob"),
#'   threshold_value = 16)
#'
#' murmuration_plot(linked)                     # threshold + blocking_var automatic
#' murmuration_plot(linked, threshold = 14)     # explore a different cutoff
#' murmuration_plot(linked, max_levels = 6)     # collapse a fine blocking_var
#' murmuration_plot(linked, palette = "grey")   # greyscale for a journal figure
#' murmuration_plot(linked, interactive = TRUE) # plotly, for Shiny/Quarto HTML
#'
#' # Advanced: the pairs-level table on its own, saved for later
#' pp <- attr(linked, "pairs_pred")
#' murmuration_plot(pp, threshold = 16, blocking_var = "gender")
#' }
#'
#' @seealso \code{\link{murmuration}}, \code{\link{preflight}}
#'
#' @references
#' Roughan WH, Cavanagh-Ochiai B, Barnekow D, Lottkowitz K, Schafer K, Camadoo L,
#' Smoll NR. Infant respiratory syncytial virus (RSV) immunisation coverage in an
#' Australian regional area, 1 May 2024 - 31 October 2025. \emph{Commun Dis Intell}
#' 2026;50. \doi{10.33321/cdi.2026.50.031}. \strong{Applied, illustrative.} About 10\%
#' of this study's birth-register records carried identifier or date-of-birth errors
#' and were resolved with weighted-score matching rather than an exact key —
#' \code{murmuration_plot()} is built for inspecting exactly that situation, checking
#' whether a candidate threshold holds up once messier records are in the mix, rather
#' than accepting a single cutoff on faith. See \code{\link{murmuration}} for the
#' full \code{@references} entry and the deterministic-linkage precedent
#' (Smoll et al. 2023), which has no score distribution to plot here.
#'
#' @importFrom ggplot2 ggplot aes geom_jitter geom_hline geom_rect geom_text annotate
#' @importFrom ggplot2 scale_colour_manual scale_fill_manual labs theme_minimal theme
#' @importFrom ggplot2 element_text element_line element_blank
#' @importFrom utils head
#' @export
murmuration_plot <- function(data,
                             threshold    = NULL,
                             blocking_var = NULL,
                             max_levels   = 12L,
                             weight_col   = "weights",
                             n_sample     = 5000L,
                             jitter_width = 0.3,
                             point_alpha  = 0.35,
                             point_size   = 1.2,
                             show_subtitle = TRUE,
                             show_counts  = TRUE,
                             palette      = "sch",
                             title        = NULL,
                             interactive  = FALSE,
                             pairs_pred   = NULL) {

  # ------------------------------------------------------------------
  # Legacy argument name shim: murmuration_plot() < 1.2.0 took the raw pairs
  # object as `pairs_pred`. Positional calls (the documented usage everywhere
  # pre-1.2.0) are completely unaffected by the rename to `data` -- this only
  # catches callers who used pairs_pred = ... by name.
  # ------------------------------------------------------------------
  if (!missing(pairs_pred) && missing(data)) {
    warning(
      "(*)> starling::murmuration_plot() \u2014 the 'pairs_pred' argument is deprecated ",
      "as of starling 1.2.0; use 'data' instead ('data' now accepts murmuration()'s own ",
      "return value directly, not just a raw pairs-level object).\n",
      "  Using the value passed to 'pairs_pred' for now \u2014 please update your call.",
      call. = FALSE
    )
    data <- pairs_pred
  }

  # ------------------------------------------------------------------
  # Input validation
  # ------------------------------------------------------------------
  if (!inherits(data, "data.frame")) {
    stop(
      "(*)> starling::murmuration_plot() \u2014 'data' must be a data frame or pairs ",
      "object: pass the object returned by murmuration(), or a pairs-level data frame ",
      "with a weights column.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------------
  # Resolve the pairs-level scoring table
  # ------------------------------------------------------------------
  pairs_tbl      <- attr(data, "pairs_pred", exact = TRUE)
  murm_method    <- attr(data, "murmuration_method", exact = TRUE)
  has_pairs_attr <- !is.null(pairs_tbl)

  if (!has_pairs_attr) {

    if (!is.null(murm_method) && identical(murm_method, "deterministic")) {
      stop(
        "(*)> starling::murmuration_plot() \u2014 'data' was produced by ",
        "murmuration(method = \"deterministic\"), which has no Fellegi-Sunter score to ",
        "plot (weights is NA throughout, by design \u2014 deterministic linkage is ",
        "exact-match and threshold-free).\n",
        "  Re-run with method = \"probabilistic\" if you want to inspect a score ",
        "distribution before choosing a threshold.",
        call. = FALSE
      )
    }

    # No pairs_pred attribute: treat `data` itself as the pairs-level table
    # (the legacy raw pairs-scoring workflow, or any bare pairs-shaped data frame).
    # Whether it actually has a usable weight_col -- present at all, and not
    # all-NA -- is checked uniformly below, same as for a "real" pairs_pred.
    pairs_tbl <- data
  }

  if (!weight_col %in% names(pairs_tbl)) {
    stop(
      "(*)> starling::murmuration_plot() \u2014 Column '", weight_col, "' not found in pairs_pred.\n",
      "Available columns: ", paste(utils::head(names(pairs_tbl), 20), collapse = ", "),
      if (has_pairs_attr) "" else paste0(
        "\n  If this is meant to be murmuration()'s output, pass the object returned ",
        "directly by murmuration(method = \"probabilistic\") \u2014 no ",
        "attr(data, \"pairs_pred\") was found on it."
      )
    )
  }

  # ------------------------------------------------------------------
  # Resolve threshold: explicit argument > murmuration() attribute > default
  # ------------------------------------------------------------------
  if (is.null(threshold)) {
    threshold <- attr(data, "threshold_value", exact = TRUE)
    if (is.null(threshold) || is.na(threshold)) {
      threshold <- 16
      message(
        "(*)> starling::murmuration_plot() \u2014 No threshold supplied and none found on ",
        "'data' (attr(data, \"threshold_value\")). Using the murmuration() default of 16 ",
        "\u2014 pass threshold explicitly to override."
      )
    }
  }
  if (!is.numeric(threshold) || length(threshold) != 1) {
    stop("(*)> starling::murmuration_plot() \u2014 'threshold' must be a single numeric value.")
  }

  # ------------------------------------------------------------------
  # Resolve blocking_var: explicit argument > murmuration() attribute > none
  # ------------------------------------------------------------------
  if (is.null(blocking_var)) {
    candidate <- attr(data, "blocking_var", exact = TRUE)
    if (!is.null(candidate) && !is.na(candidate) && candidate %in% names(pairs_tbl)) {
      blocking_var <- candidate
    }
  } else if (!blocking_var %in% names(pairs_tbl)) {
    stop(
      "(*)> starling::murmuration_plot() \u2014 Column '", blocking_var, "' not found in the ",
      "pairs-level table.\nAvailable columns: ",
      paste(utils::head(names(pairs_tbl), 20), collapse = ", ")
    )
  }

  use_blocks <- !is.null(blocking_var)

  # Only worth a message when `data` came from murmuration() itself and its
  # blocking_var attribute is unexpectedly missing/invalid -- a bare/legacy
  # pairs object was never going to carry blocking info, so unstratified
  # plotting there is simply the expected default, not a fallback worth
  # flagging every time.
  if (!use_blocks && has_pairs_attr) {
    message(
      "(*)> starling::murmuration_plot() \u2014 No blocking variable found on 'data' (its ",
      "attr(data, \"blocking_var\") is missing or invalid). Plotting all pairs in a ",
      "single column. Supply blocking_var explicitly to jitter by blocking stratum ",
      "instead."
    )
  }

  palette <- match.arg(palette, c("sch", "default", "grey"))

  # ------------------------------------------------------------------
  # Colour palette definitions
  # ------------------------------------------------------------------
  colours <- switch(palette,
    sch = list(
      match    = "#08403B",   # SCH Evergreen
      nonmatch = "#C8A96E",   # SCH Gold/Warm
      zone_match    = "#08403B22",
      zone_nonmatch = "#C8A96E22",
      threshold_line = "#D64045",
      threshold_text = "#D64045"
    ),
    default = list(
      match    = "#1B7B8A",
      nonmatch = "#E07070",
      zone_match    = "#1B7B8A22",
      zone_nonmatch = "#E0707022",
      threshold_line = "#333333",
      threshold_text = "#333333"
    ),
    grey = list(
      match    = "#333333",
      nonmatch = "#999999",
      zone_match    = "#33333315",
      zone_nonmatch = "#99999915",
      threshold_line = "#000000",
      threshold_text = "#000000"
    )
  )

  # ------------------------------------------------------------------
  # Extract weights + block membership as plain vectors (works identically
  # whether pairs_tbl is a data.frame, tibble, or a data.table) and
  # drop missing weights
  # ------------------------------------------------------------------
  weight_vals <- pairs_tbl[[weight_col]]

  if (use_blocks) {
    block_vals_raw <- as.character(pairs_tbl[[blocking_var]])
  } else {
    block_vals_raw <- rep("All pairs", length(weight_vals))
  }

  keep <- !is.na(weight_vals)
  weight_vals    <- weight_vals[keep]
  block_vals_raw <- block_vals_raw[keep]

  if (length(weight_vals) == 0) {
    stop("(*)> starling::murmuration_plot() \u2014 No non-missing weights found in '", weight_col, "'.")
  }

  block_vals_raw[is.na(block_vals_raw)] <- "(missing)"

  # ------------------------------------------------------------------
  # Per-stratum stats, computed once up front, used to decide both which
  # strata are shown individually and what order they appear in
  # ------------------------------------------------------------------
  is_above     <- weight_vals >= threshold
  total_counts <- table(block_vals_raw)
  above_counts <- tapply(is_above, block_vals_raw, sum)
  above_counts <- above_counts[names(total_counts)]  # same order/names as total_counts
  n_levels     <- length(total_counts)

  # ------------------------------------------------------------------
  # Select which strata get their own column if there are too many to show,
  # and order all shown strata left-to-right by where matches concentrate
  # ------------------------------------------------------------------
  # A stratum containing even one above-threshold pair is exactly the kind of
  # thing this plot exists to surface -- it is never pooled into "Other" just
  # because its total volume is small, the way pure frequency-based selection
  # would risk doing. Remaining slots (if any) go to the highest-volume
  # zero-match strata, for context on what a "typical" non-match stratum
  # looks like; only the least-relevant tail is pooled.
  if (use_blocks && n_levels > max_levels) {

    has_match    <- names(total_counts)[above_counts > 0]
    no_match     <- names(total_counts)[above_counts == 0]
    no_match_ord <- no_match[order(total_counts[no_match], decreasing = TRUE)]

    if (length(has_match) >= max_levels - 1) {
      # Even the matching strata alone don't all fit: keep the most
      # productive of them (some matching strata may still end up pooled)
      keep_levels <- has_match[order(above_counts[has_match], decreasing = TRUE)][seq_len(max_levels - 1)]
    } else {
      keep_levels <- c(has_match, no_match_ord[seq_len(max_levels - 1 - length(has_match))])
    }

    n_pooled <- sum(total_counts) - sum(total_counts[keep_levels])
    warning(
      "(*)> starling::murmuration_plot() \u2014 ", n_levels, " blocking strata found in '",
      blocking_var, "'; showing the ", max_levels - 1, " most relevant (any stratum with ",
      "an above-threshold pair is kept first, then the highest-volume of the rest) plus ",
      "\"Other\" (", format(n_pooled, big.mark = ","), " pairs pooled). Increase ",
      "max_levels to show more, or pass a coarser blocking_var.",
      call. = FALSE
    )
    block_vals_raw <- ifelse(block_vals_raw %in% keep_levels, block_vals_raw, "Other")
    level_order <- c(keep_levels[order(above_counts[keep_levels], decreasing = TRUE)], "Other")

  } else {
    level_order <- names(total_counts)[
      order(-above_counts[names(total_counts)], -total_counts[names(total_counts)])
    ]
  }

  full_df <- data.frame(
    weights = weight_vals,
    block   = factor(block_vals_raw, levels = level_order),
    stringsAsFactors = FALSE
  )
  full_df$match_grp <- ifelse(
    full_df$weights >= threshold, "Above threshold", "Below threshold"
  )
  n_total <- nrow(full_df)

  above_idx <- which(full_df$match_grp == "Above threshold")
  below_idx <- which(full_df$match_grp == "Below threshold")
  n_above   <- length(above_idx)
  n_below   <- length(below_idx)

  w_min <- min(full_df$weights)
  w_max <- max(full_df$weights)
  w_range_pad <- (w_max - w_min) * 0.05

  # Per-block pair counts, computed on the FULL data (not the plotted sample)
  block_summary <- as.data.frame(table(full_df$block))
  names(block_summary) <- c("block", "n_pairs")
  block_summary <- block_summary[block_summary$n_pairs > 0, , drop = FALSE]

  # ------------------------------------------------------------------
  # Subsample -- but never at the expense of above-threshold pairs
  # ------------------------------------------------------------------
  # Above-threshold pairs are typically the rare, diagnostically critical
  # class -- often well under 1% of all candidate pairs, sometimes far less.
  # Subsampling proportionally across the whole dataset (or even just within
  # each blocking stratum) can silently sample every single one of them away,
  # making a real -- if small -- signal look like no separation at all. Every
  # above-threshold pair is kept (unless there are so many that they alone
  # exceed n_sample, an unusual case); the below-threshold population, which
  # is normally overwhelming and whose *shape* rather than every point is
  # what matters, is subsampled within that remaining budget, stratified by
  # blocking stratum as before.
  all_above_shown <- TRUE

  if (n_total > n_sample) {

    if (n_above >= n_sample) {
      set.seed(42L)
      keep_above <- sample(above_idx, n_sample)
      keep_below <- integer(0)
      all_above_shown <- FALSE
    } else {
      keep_above  <- above_idx
      below_budget <- n_sample - n_above
      if (n_below > below_budget) {
        frac <- below_budget / n_below
        set.seed(42L)
        idx_by_block <- split(below_idx, full_df$block[below_idx])
        keep_below <- unlist(lapply(idx_by_block, function(idx) {
          if (length(idx) == 0) return(integer(0))
          k <- max(1L, round(length(idx) * frac))
          k <- min(k, length(idx))
          sample(idx, k)
        }), use.names = FALSE)
      } else {
        keep_below <- below_idx
      }
    }

    plot_df <- full_df[c(keep_above, keep_below), , drop = FALSE]

    message(
      "(*)> starling::murmuration_plot() \u2014 Plotting ",
      if (all_above_shown) paste0("all ", format(n_above, big.mark = ","), " above-threshold pairs")
      else paste0("a sample of ", format(length(keep_above), big.mark = ","), " above-threshold pairs"),
      " plus a ", if (use_blocks) "stratified" else "random", " sample of ",
      format(length(keep_below), big.mark = ","), " below-threshold pairs (",
      format(n_total, big.mark = ","), " total)."
    )
  } else {
    plot_df <- full_df
  }

  # ------------------------------------------------------------------
  # Build the plot
  # ------------------------------------------------------------------
  plot_title <- title %||% "Linkage Weight Distribution by Blocking Stratum"
  plot_subtitle <- paste0(
    "threshold = ", threshold, "   |   ",
    format(n_above, big.mark = ","), " above  /  ",
    format(n_below, big.mark = ","), " below",
    if (n_total > n_sample && all_above_shown) "   |   every above-threshold pair shown; below-threshold sampled" else ""
  )

  n_blocks_shown <- nlevels(plot_df$block)
  rotate_labels  <- use_blocks && n_blocks_shown > 6

  below_df <- plot_df[plot_df$match_grp == "Below threshold", , drop = FALSE]
  above_df <- plot_df[plot_df$match_grp == "Above threshold", , drop = FALSE]

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = block, y = weights, colour = match_grp)) +

    # Shaded zones (one dummy row of constant data -- avoids drawing one
    # identical rectangle per plotted point)
    ggplot2::geom_rect(
      data = data.frame(x = 1),
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = threshold, ymax = w_max + w_range_pad),
      fill = colours$zone_match, colour = NA, inherit.aes = FALSE
    ) +
    ggplot2::geom_rect(
      data = data.frame(x = 1),
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = w_min - w_range_pad, ymax = threshold),
      fill = colours$zone_nonmatch, colour = NA, inherit.aes = FALSE
    ) +

    # Below-threshold points -- the (usually enormous) background population.
    # height = 0 so a pair's actual score is never distorted.
    ggplot2::geom_jitter(
      data = below_df, width = jitter_width, height = 0,
      alpha = point_alpha, size = point_size, shape = 16
    ) +

    # Above-threshold points -- drawn on top, bigger, more opaque, a
    # different shape, so the rare/important class is never visually lost
    # against the background cloud even under heavy overplotting.
    ggplot2::geom_jitter(
      data = above_df, width = jitter_width, height = 0,
      alpha = min(point_alpha * 2.6, 0.95), size = point_size * 2, shape = 17
    ) +

    # Threshold line
    ggplot2::geom_hline(
      yintercept = threshold, colour = colours$threshold_line,
      linewidth = 0.9, linetype = "dashed"
    ) +

    # Threshold label, anchored to the left edge of the panel
    ggplot2::annotate(
      "text", x = -Inf, y = threshold + w_range_pad * 0.8,
      label = paste0("Threshold = ", threshold),
      colour = colours$threshold_text, size = 3.5, hjust = -0.05, fontface = "bold"
    ) +

    ggplot2::scale_colour_manual(
      values = c("Above threshold" = colours$match,
                 "Below threshold" = colours$nonmatch),
      name = NULL
    ) +

    ggplot2::labs(
      title    = plot_title,
      subtitle = if (show_subtitle) plot_subtitle else NULL,
      x        = if (use_blocks) paste0("Blocking stratum (", blocking_var, ")") else NULL,
      y        = "Linkage weight (Fellegi-Sunter log-likelihood ratio)"
    ) +

    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      legend.position     = "bottom",
      plot.title          = ggplot2::element_text(face = "bold"),
      plot.subtitle       = ggplot2::element_text(colour = "grey40", size = 10),
      axis.text.x         = if (!use_blocks) {
                               ggplot2::element_blank()
                             } else if (rotate_labels) {
                               ggplot2::element_text(angle = 45, hjust = 1)
                             } else {
                               ggplot2::element_text()
                             },
      axis.ticks.x        = if (use_blocks) ggplot2::element_line() else ggplot2::element_blank()
    )

  # Per-stratum pair counts, printed just below the shaded "non-match" zone.
  # No extra scale_y_continuous() expansion is needed here: the label's own y
  # position already extends the panel's computed data range, and ggplot2's
  # default 5% expansion around that is enough to keep it clear of the axis
  # without introducing extra dead space underneath it.
  if (show_counts && use_blocks) {
    p <- p + ggplot2::geom_text(
      data = block_summary,
      ggplot2::aes(x = block, y = w_min - w_range_pad * 1.8,
                   label = paste0("n=", format(n_pairs, big.mark = ","))),
      inherit.aes = FALSE, size = 3, colour = "grey40"
    )
  }

  # ------------------------------------------------------------------
  # Interactive toggle
  # ------------------------------------------------------------------
  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop(
        "(*)> starling::murmuration_plot() \u2014 Package 'plotly' is required for interactive output.\n",
        "Install it with: install.packages('plotly')"
      )
    }
    block_label <- blocking_var  # captured for use inside aes(text = ...) below
    tooltip_text <- function(d) paste0(
      if (use_blocks) paste0(block_label, ": ", d$block, "<br>") else "",
      "Weight: ", round(d$weights, 2), "<br>Group: ", d$match_grp
    )
    below_df$text <- tooltip_text(below_df)
    above_df$text <- tooltip_text(above_df)

    # p already carries the tooltip-free static layers built above; rebuilt
    # here with a `text` aesthetic so plotly's tooltip has something to show
    p_base <- ggplot2::ggplot(plot_df, ggplot2::aes(x = block, y = weights, colour = match_grp)) +
      ggplot2::geom_jitter(data = below_df, ggplot2::aes(text = text),
                           width = jitter_width, height = 0, alpha = point_alpha,
                           size = point_size, shape = 16) +
      ggplot2::geom_jitter(data = above_df, ggplot2::aes(text = text),
                           width = jitter_width, height = 0,
                           alpha = min(point_alpha * 2.6, 0.95),
                           size = point_size * 2, shape = 17) +
      ggplot2::geom_hline(yintercept = threshold,
                          colour = colours$threshold_line,
                          linewidth = 0.9, linetype = "dashed") +
      ggplot2::scale_colour_manual(
        values = c("Above threshold" = colours$match,
                   "Below threshold" = colours$nonmatch),
        name = NULL
      ) +
      ggplot2::labs(title    = plot_title,
                    subtitle = if (show_subtitle) plot_subtitle else NULL,
                    x = if (use_blocks) paste0("Blocking stratum (", blocking_var, ")") else NULL,
                    y = "Linkage weight (Fellegi-Sunter log-likelihood ratio)") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())

    return(plotly::ggplotly(p_base, tooltip = "text"))
  }

  p
}

# Null coalescing operator (avoids dependency on rlang::`%||%` here)
`%||%` <- function(x, y) if (is.null(x)) y else x
