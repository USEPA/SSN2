skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")
skip_if_not_installed("pracma")

stream_wave_learning_fixture <- function() {
  network <- mf04p
  network$obs <- network$obs[round(seq(1, NROW(network$obs), length.out = 12)), ]
  network$preds <- list(pred1km = network$preds$pred1km[1:4, ])
  network$path <- withr::local_tempdir(.local_envir = parent.frame())
  stopifnot(file.copy(file.path(mf04p$path, "binaryID.db"), network$path))
  ssn_create_distmat(network, predpts = "pred1km", among_predpts = TRUE, overwrite = TRUE)
  network
}

test_that("stream wave grids do not evaluate effective ranges or warn", {
  set.seed(2)
  network <- stream_wave_learning_fixture()
  fits <- 0L
  testthat::local_mocked_bindings(
    get_effective_range = function(...) stop("Unexpected effective-range calculation"),
    fit_decorrelate_algorithm = function(X, y, algorithm, dots) {
      fits <<- fits + 1L
      mean(y)
    },
    predict_decorrelate_algorithm = function(fit, X, algorithm) rep(fit, NROW(X)),
    .package = "SSN2"
  )
  for (component in c("up", "down", "both")) {
    before <- fits
    expect_no_warning(fit <- ssn_decorrelate(
      Summer_mn ~ ELEV_DEM, network,
      tailup_type = if (component == "down") "none" else "wave",
      taildown_type = if (component == "up") "none" else "wave",
      additive = "afvArea", dense_grid = FALSE, local = FALSE, ordering = "none",
      training = list(training_index = 1:8, test_index = 9:12)
    ))
    expect_s3_class(fit, "ssn_decorrelate")
    grid <- tidy(fit)
    expect_equal(NROW(grid), if (component == "both") 15L else 5L)
    expect_equal(fits - before, NROW(grid) + 1L)
    expect_equal(sum(grid$tailup_type == "no transformation"), 1L)
    if (component != "down") expect_true(any(grid$tailup_type == "wave"))
    if (component != "up") expect_true(any(grid$taildown_type == "wave"))
    expect_true(all(is.finite(grid$RMSPE)))
    expect_equal(fit$test$RMSPE, min(grid$RMSPE))
    expect_no_warning(prediction <- predict(fit, "pred1km", local = FALSE))
    expect_length(prediction, NROW(network$preds$pred1km))
    expect_true(all(is.finite(prediction)))
  }
})

test_that("ranger algorithm fits and predicts with fixed stream wave covariance", {
  skip_if_not_installed("ranger")
  set.seed(2)
  network <- stream_wave_learning_fixture()
  tailup <- tailup_params("wave", .4, 2500)
  taildown <- taildown_params("wave", .3, 3500)
  testthat::local_mocked_bindings(
    get_effective_range = function(...) stop("Unexpected effective-range calculation"),
    .package = "SSN2"
  )
  expect_no_warning(fit <- ssn_decorrelate(
    Summer_mn ~ ELEV_DEM, network, additive = "afvArea",
    tailup_params = tailup, taildown_params = taildown,
    euclid_params = euclid_params("none"), nugget_params = nugget_params("nugget", .25),
    algorithm = "ranger", local = FALSE, ordering = "none",
    training = list(training_index = 1:8, test_index = 9:12), evaluate_test = TRUE,
    num.trees = 20, seed = 2, num.threads = 1
  ))
  expect_s3_class(fit$fit, "ranger")
  expect_equal(fit$fit$num.trees, 20)
  expect_equal(fit$decorrelate_data$params$tailup, tailup)
  expect_equal(fit$decorrelate_data$params$taildown, taildown)
  expect_true(is.finite(fit$test$RMSPE))
  expect_no_warning(prediction <- predict(fit, "pred1km", local = FALSE))
  expect_length(prediction, NROW(network$preds$pred1km))
  expect_true(all(is.finite(prediction)))
  transformed <- ssn_decorrelate_newdata(fit$decorrelate_data, "pred1km", local = FALSE)
  learner <- predict(fit$fit, data = transformed$tX_newdata)$predictions
  expect_equal(as.vector(prediction), as.vector(ssn_recorrelate_newdata(transformed, learner)),
    tolerance = 1e-10)
})
