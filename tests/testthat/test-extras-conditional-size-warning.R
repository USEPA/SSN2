skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

conditional_size_design <- function(q) {
  list(observed = Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(1L, q)),
       variance = rep(1, q), info = data.frame(conditioned = rep(TRUE, q)))
}

test_that("random-effect size warnings continue past the strict threshold", {
  set.seed(2)
  # Stop before allocating the dense effect covariance.
  testthat::local_mocked_bindings(
    conditional_spatial_object = function(...) stop("past size warning"),
    .package = "SSN2"
  )
  for (q in c(10000L, 10001L, 20000L)) {
    design <- conditional_size_design(q)
    message <- if (q == 10000L) NA else paste0(
      "Conditional simulation requires factoring a covariance matrix for ",
      format(q, big.mark = ",", scientific = FALSE),
      " random-effect coefficients. This may result in exceedingly long computational times."
    )
    for (ctl in list(NULL, list(approximation = "vecchia"))) {
      expect_warning(
        expect_error(SSN2:::prepare_conditional_random(list(), design, ctl), "past size warning"),
        message, fixed = TRUE
      )
    }
  }
})

test_that("random-effect size warnings count only positive conditional variances", {
  set.seed(2)
  testthat::local_mocked_bindings(
    conditional_spatial_object = function(...) stop("past size warning"),
    .package = "SSN2"
  )
  for (ctl in list(NULL, list(approximation = "vecchia"))) {
    design <- conditional_size_design(10001L)
    design$variance[10001L] <- 0
    expect_warning(
      expect_error(SSN2:::prepare_conditional_random(list(), design, ctl), "past size warning"),
      NA
    )
    design$variance[10001L] <- 1
    design$info$conditioned[10001L] <- FALSE
    expect_warning(
      expect_error(SSN2:::prepare_conditional_random(list(), design, ctl), "past size warning"),
      NA
    )
  }
})

test_that("invalid random-effect variances fail before size warnings", {
  set.seed(2)
  testthat::local_mocked_bindings(
    conditional_spatial_object = function(...) stop("past size warning"),
    .package = "SSN2"
  )
  for (ctl in list(NULL, list(approximation = "vecchia"))) {
    design <- conditional_size_design(20000L)
    for (invalid in c(-1, NA_real_, NaN, Inf, -Inf)) {
      design$variance[1L] <- invalid
      expect_warning(
        expect_error(SSN2:::prepare_conditional_random(list(), design, ctl),
                     "Invalid random-effect variances.", fixed = TRUE),
        NA
      )
    }
  }
})
