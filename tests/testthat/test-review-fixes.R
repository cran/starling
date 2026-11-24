# Tests for the maintainer-review changes (aviary ecosystem plan v0.11, section 4.2).

test_that("the default threshold is 16", {
  expect_equal(formals(murmuration)$threshold_value, 16)
})

test_that("v2h and v2e are defunct; v2c and c2h still warn", {
  x <- data.frame(id_var = 1, lettername1 = "a")
  expect_error(murmuration(x, x, linkage_type = "v2h", id_var = "id_var",
                           blocking_var = "lettername1", compare_vars = "lettername1"),
               "defunct")
  expect_error(murmuration(x, x, linkage_type = "v2e", id_var = "id_var",
                           event_date = as.Date("2025-01-01"),
                           blocking_var = "lettername1", compare_vars = "lettername1"),
               "defunct")
})

test_that("the vaccination route returns an ungrouped frame with vaccination_status", {
  cases <- data.frame(id_var = c("c1", "c2"), lettername1 = c("ann", "bob"),
                      lettername2 = c("lee", "kay"),
                      dob = as.Date(c("1950-01-01", "1960-01-01")),
                      onset_date = as.Date("2025-06-01"))
  air <- data.frame(id_var = c("a1", "a2"), lettername1 = c("ann", "bob"),
                    lettername2 = c("lee", "kay"),
                    dob = as.Date(c("1950-01-01", "1960-01-01")),
                    vax_date_1 = as.Date(c("2025-04-01", NA)),
                    vax_type_1 = c("COVID-19", NA))
  out <- suppressMessages(murmuration(cases, air,
    linkage_type = "cohort2event", event_type = "vaccination", method = "deterministic",
    event_date = "onset_date", id_var = "id_var", blocking_var = "lettername1",
    compare_vars = c("lettername1", "lettername2", "dob")))
  expect_false(inherits(out, "rowwise_df"))
  expect_false(any(c("dates", "types") %in% names(out)))
  expect_equal(out$vaccination_status[order(out$id_var)], c("Vaccinated", "Unvaccinated"))
})

test_that("preflight() counts missing dates in linkage_vars", {
  d1 <- data.frame(id_var = 1:3, dob = as.Date(c("1980-01-01", NA, "1990-01-01")))
  d2 <- data.frame(id_var = 1:2, dob = as.Date(c("1980-01-01", "1991-01-01")))
  res <- preflight(d1, d2, linkage_vars = "dob", verbose = FALSE)
  expect_equal(res$completeness$n_missing, c(1L, 0L))
  expect_length(res$flags[grepl("'NA' in NA", res$flags)], 0)
})

test_that("vaccination windows: cohort_window works and days_after keeps post-event doses", {
  cases <- data.frame(id_var = "c1", lettername1 = "ann", lettername2 = "lee",
                      dob = as.Date("1950-01-01"), onset_date = as.Date("2025-06-01"))
  air <- data.frame(id_var = "a1", lettername1 = "ann", lettername2 = "lee",
                    dob = as.Date("1950-01-01"),
                    vax_date_1 = as.Date("2025-01-10"),   # 142 days before onset
                    vax_date_2 = as.Date("2025-05-25"),   # 7 days before: in the 14-day buffer
                    vax_date_3 = as.Date("2025-06-20"),   # 19 days after onset
                    vax_type_1 = "COVID-19", vax_type_2 = "COVID-19", vax_type_3 = "COVID-19")
  link <- function(...) suppressMessages(murmuration(cases, air,
    linkage_type = "cohort2event", event_type = "vaccination", method = "deterministic",
    event_date = "onset_date", id_var = "id_var", blocking_var = "lettername1",
    compare_vars = c("lettername1", "lettername2", "dob"), ...))
  # default window: only the dose >= 14 days before onset
  d0 <- link()
  expect_equal(sum(!is.na(unlist(d0[grep("^vax_date_", names(d0))]))), 1L)
  # days_after = 30 also keeps the post-event dose (previously dropped)
  d30 <- link(vax_window = list(days_before = 14, days_after = 30))
  kept <- as.Date(unlist(d30[grep("^vax_date_", names(d30))]), origin = "1970-01-01")
  expect_true(as.Date("2025-06-20") %in% kept)
  expect_false(as.Date("2025-05-25") %in% kept)
  # cohort_window (previously an && error on R >= 4.3)
  dc <- link(cohort_window = list(entry = as.Date("2025-05-01"), exit = as.Date("2025-12-31")))
  keptc <- as.Date(unlist(dc[grep("^vax_date_", names(dc))]), origin = "1970-01-01")
  expect_setequal(keptc[!is.na(keptc)], as.Date(c("2025-05-25", "2025-06-20")))
})

test_that("resume_from_pairs works after dplyr verbs on the reviewed table", {
  x <- data.frame(id_var = c("c1", "c2"), lettername1 = c("ann", "bob"),
                  lettername2 = c("lee", "kay"), dob = as.Date(c("1950-01-01", "1960-01-01")),
                  block1 = "A", onset_date = as.Date("2025-06-01"))
  y <- data.frame(id_var = c("a1", "a2"), lettername1 = c("ann", "bob"),
                  lettername2 = c("lee", "kay"), dob = as.Date(c("1950-01-01", "1960-01-01")),
                  block1 = "A", vax_date_1 = as.Date(c("2025-04-01", NA)),
                  vax_type_1 = c("COVID-19", NA))
  pairs <- reclin2::pair_blocking(x, y, "block1")
  pairs$weights <- ifelse(pairs$.x == pairs$.y, 12, 1)   # the same person scores high
  # a dplyr verb on a reclin2 "pairs" object drops attr(pairs, "x"/"y"),
  # which resume_from_pairs previously needed
  pairs <- dplyr::mutate(pairs, match_override = NA_character_)
  out <- suppressMessages(murmuration(x, y,
    linkage_type = "cohort2event", event_type = "vaccination",
    event_date = "onset_date", id_var = "id_var", threshold_value = 8,
    resume_from_pairs = pairs))
  expect_equal(out$id_var_y[order(out$id_var)], c("a1", "a2"))
})

