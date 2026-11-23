skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("block controls resolve exact names before applying defaults", {
  withr::local_options(warnPartialMatchDollar = TRUE)
  resolve <- SSN2:::get_local_list_prediction_block
  for (method_new in c("subset", "basis")) {
    for (controls in list(list(method_new = method_new),
                         list(method_new = method_new, size_new = 3))) {
      expect_warning(actual <- resolve(controls), NA)
      expect_identical(actual, resolve(c(controls, list(method = "covariance"))))
      expect_identical(actual$method, "covariance")
      expect_identical(actual$size, 4000L)
    }
  }
  expect_identical(resolve(TRUE), resolve(list()))
  expect_identical(resolve(FALSE), resolve(list(method = "all", method_new = "basis", size_new = Inf)))
  expect_null(resolve(list(method = "all", size_new = 3))[["size"]])
  expect_identical(resolve(list(size = 5, size_new = 3))$size, 5L)
  expect_error(resolve(list(method = "subset")), "Invalid local method")
  expect_error(resolve(list(method_new = "covariance")), "local\\$method_new")
  expect_error(resolve(list(method_new = "subset", size_new = 0)), "local\\$size_new")
})

test_that("exported block controls preserve predictions and selected rows", {
  set.seed(2)
  network <- mf04p
  network$obs <- network$obs[seq_len(24), ]
  network$preds$CapeHorn <- network$preds$CapeHorn[seq_len(12), ]
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, network,
    tailup_type = "none", taildown_type = "none",
    euclid_initial = euclid_initial("exponential", 1, 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.2, known = "given"))
  withr::local_options(warnPartialMatchDollar = TRUE)
  cross_cov <- as.matrix(covmatrix(fit, "CapeHorn"))
  selected <- NULL
  local_mocked_bindings(coefficients = function(object, ...) {
    frame <- parent.frame()
    if (exists("covariance_rows", frame, inherits = FALSE)) {
      selected <<- list(
        nodes = frame$nodes, grid = frame$covariance_rows,
        obs = if (exists("cov_index", frame, inherits = FALSE)) frame$cov_index else seq_len(object$n)
      )
    }
    stats::coefficients(object, ...)
  }, .package = "SSN2")
  run <- function(controls) {
    expect_warning(result <- predict(fit, "CapeHorn", block = TRUE,
      se.fit = TRUE, local = controls), NA)
    expect_true(all(is.finite(c(result$fit, result$se.fit))))
    list(result = result, selected = selected)
  }

  for (method_new in c("subset", "basis")) {
    controls <- list(method_new = method_new, size = 6, size_new = 3,
                     ordering = "none", chunk_size = 4)
    omitted <- run(controls)
    explicit <- run(c(list(method = "covariance"), controls))
    expect_equal(omitted$result, explicit$result, tolerance = 1e-12)
    expect_identical(omitted$selected, explicit$selected)
    rows <- if (method_new == "subset") 1:3 else 1:12
    c0 <- colMeans(cross_cov[rows, , drop = FALSE])
    expected_obs <- order(c0, decreasing = TRUE)[seq_len(6)]
    expect_identical(omitted$selected$nodes, 1:3)
    expect_identical(omitted$selected$grid, rows)
    expect_identical(omitted$selected$obs, expected_obs)

    all <- run(c(list(method = "all"), controls))
    expect_identical(all$selected$nodes, 1:3)
    expect_identical(all$selected$grid, rows)
    expect_identical(all$selected$obs, seq_len(fit$n))
    full <- run(modifyList(controls, list(size = fit$n)))
    expect_equal(all$result, full$result, tolerance = 1e-10)
  }
  expect_equal(run(FALSE), run(list(method = "all", size_new = Inf)))
  expect_equal(run(TRUE), run(list()))
  expect_equal(run(TRUE)$result, run(FALSE)$result, tolerance = 1e-10)
})

test_that("block controls leave Gaussian and GLM point prediction unchanged", {
  set.seed(2)
  network <- mf04p
  network$obs <- network$obs[seq_len(24), ]
  network$obs$count <- rpois(24, 5)
  network$preds$CapeHorn <- network$preds$CapeHorn[seq_len(3), ]
  common <- list(ssn.object = network, tailup_type = "none", taildown_type = "none",
    euclid_type = "none", nugget_initial = nugget_initial("nugget", 0.2, known = "given"))
  fits <- list(
    do.call(ssn_lm, c(list(formula = Summer_mn ~ ELEV_DEM), common)),
    do.call(ssn_glm, c(list(formula = count ~ ELEV_DEM, family = "poisson"), common))
  )
  withr::local_options(warnPartialMatchDollar = TRUE)
  for (fit in fits) {
    for (method_new in c("subset", "basis")) {
      controls <- list(method_new = method_new, size = 6, size_new = 2)
      expect_warning(omitted <- predict(fit, "CapeHorn", se.fit = TRUE,
        local = controls), NA)
      expect_warning(explicit <- predict(fit, "CapeHorn", se.fit = TRUE,
        local = c(list(method = "covariance"), controls)), NA)
      expect_equal(omitted, explicit)
    }
  }
  expect_error(predict(fits[[2]], "CapeHorn", block = TRUE,
    local = list(method_new = "subset")), "Block prediction is not supported")
})
