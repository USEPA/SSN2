skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("point and block prediction use separate neighborhood defaults", {
  for (local in list(TRUE, list(method = "covariance"), list(parallel = FALSE))) {
    point <- SSN2:::get_local_list_prediction(local)
    block <- SSN2:::get_local_list_prediction_block(local)
    expect_equal(point$size, 200)
    expect_equal(block$size, 4000)
    expect_equal(block$size_new, 4000)
  }

  expect_equal(SSN2:::get_local_list_prediction(NULL)$size, 200)
  expect_equal(SSN2:::get_local_list_prediction(list(size = 35))$size, 35)
  expect_equal(SSN2:::get_local_list_prediction_block(list(size = 35))$size, 35)

  point <- SSN2:::get_local_list_prediction(FALSE)
  block <- SSN2:::get_local_list_prediction_block(FALSE)
  expect_identical(point$method, "all")
  expect_null(point$size)
  expect_identical(block$method, "all")
  expect_equal(block$size_new, Inf)
})

test_that("block prediction defaults to PID and accepts shared orderings", {
  for (local in list(TRUE, FALSE, list(size = 20), list(ordering = "pid"),
                    list(ordering = NULL))) {
    expect_identical(SSN2:::get_local_list_prediction_block(local)$ordering, "pid")
  }
  for (ordering in c("pid", "none", "random", "maxmin", "middleout",
                     "outsidein", "coordinate", "grts")) {
    expect_identical(SSN2:::get_local_list_prediction_block(list(ordering = ordering))$ordering, ordering)
  }
  for (ordering in list("network", "invalid", NA_character_,
                       character(), c("pid", "pid"), 1)) {
    expect_error(
      SSN2:::get_local_list_prediction_block(list(ordering = ordering)),
      'ordering must be "pid"', fixed = TRUE
    )
  }
})

test_that("block nodes follow the shared PID ordering across networks", {
  grid <- data.frame(netgeom = c(
    "SNETWORK (2 1 0 1 30 1)", "SNETWORK (1 1 0 1 10 2)",
    "SNETWORK (1 1 0 1 40 3)", "SNETWORK (2 1 0 1 20 4)",
    "SNETWORK (2 1 0 1 10 5)", "SNETWORK (10 1 0 1 10 6)"
  ))
  rows <- SSN2:::get_decorrelate_rows(grid)
  ordering <- SSN2:::get_decorrelate_order(rows, "pid")
  expect_identical(ordering, c(2L, 6L, 5L, 4L, 1L, 3L))
  expect_identical(SSN2:::get_block_nodes(grid, 3, "pid"), ordering[c(2, 4, 6)])
  expect_identical(SSN2:::get_block_nodes(grid, 6, "pid"), seq_len(6))
  expect_identical(SSN2:::get_block_nodes(grid, 10, "pid"), seq_len(6))
  expect_identical(SSN2:::get_block_nodes(grid, 3, "none"), seq_len(3))
  set.seed(2)
  expected_random <- sample.int(6, 3)
  set.seed(2)
  expect_identical(SSN2:::get_block_nodes(grid, 3, "random"), expected_random)
  for (invalid in c("network", "invalid")) {
    expect_error(SSN2:::get_block_nodes(grid, 3, invalid), 'ordering must be "pid"', fixed = TRUE)
  }
})

test_that("block nodes use coordinate-based ordering helpers", {
  grid <- mf04p$preds$pred1km[seq_len(12), , drop = FALSE]
  rows <- SSN2:::get_decorrelate_rows(grid)
  for (ordering in c("maxmin", "middleout", "outsidein", "coordinate", "grts")) {
    dependency <- if (ordering == "grts") "spsurvey" else "GPvecchia"
    if (!requireNamespace(dependency, quietly = TRUE)) {
      expect_error(SSN2:::get_block_nodes(grid, 4, ordering), paste("Install the", dependency))
    } else {
      set.seed(2)
      expected <- SSN2:::get_decorrelate_order(rows, ordering, grid)[seq_len(4)]
      set.seed(2)
      actual <- SSN2:::get_block_nodes(grid, 4, ordering)
      expect_identical(actual, expected)
      expect_length(unique(actual), 4)
    }
  }
})

test_that("block prediction supports each ordering through the public method", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "none", taildown_type = "none", euclid_type = "none",
    nugget_initial = nugget_initial("nugget", 1, known = "given")
  )
  fit$ssn.object$preds$pred1km <- fit$ssn.object$preds$pred1km[seq_len(12), , drop = FALSE]
  reference <- predict(fit, "pred1km", block = TRUE, local = FALSE)
  for (ordering in c("pid", "none", "random", "maxmin", "middleout",
                     "outsidein", "coordinate", "grts")) {
    if (ordering == "grts" && !requireNamespace("spsurvey", quietly = TRUE)) next
    if (ordering %in% c("maxmin", "middleout", "outsidein", "coordinate") &&
        !requireNamespace("GPvecchia", quietly = TRUE)) next
    result <- predict(fit, "pred1km", block = TRUE, se.fit = TRUE,
      local = list(method = "all", method_new = "basis", size_new = 4,
                   ordering = ordering)
    )
    expect_equal(result$fit, reference)
    expect_true(is.finite(result$se.fit))
  }
  expect_error(
    predict(fit, "pred1km", block = TRUE, local = list(ordering = "network")),
    'ordering must be "pid"', fixed = TRUE
  )
})
