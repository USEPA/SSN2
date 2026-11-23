skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

parallel_decorrelation_fixture <- function() {
  network <- mf04p
  network$obs <- network$obs[rev(round(seq(1, NROW(network$obs), length.out = 12))), ]
  network$preds <- list(pred1km = network$preds$pred1km[c(5, 1, 8, 3, 6), ])
  network$obs$x <- seq(-1, 1, length.out = 12)
  network$obs$y <- seq_len(12)
  network$obs$o <- seq(.1, .5, length.out = 12)
  network$obs$g <- rep(c("a", "b", "c"), 4)
  network$obs$part <- rep(c("one", "two"), 6)
  network$preds$pred1km$x <- c(-.5, .2, .8, 0, .4)
  network$preds$pred1km$o <- seq(.1, .5, length.out = 5)
  network$preds$pred1km$g <- c("a", "new", "b", "c", "a")
  network$preds$pred1km$part <- c("one", "two", "one", "two", "one")
  network
}

parallel_decorrelation_fit <- function(network, local, ordering = "pid") {
  ssn_decorrelate_data(y ~ x + offset(o), network,
    tailup_params = tailup_params("exponential", .4, 10000),
    taildown_params = taildown_params("spherical", .3, 10000),
    euclid_params = euclid_params("exponential", .2, 12000, rotate = .3, scale = .6),
    nugget_params = nugget_params("nugget", .1), additive = "afvArea",
    randcov_params = randcov_params("(1 | g)" = .2, "(0 + x | g)" = .1),
    partition_factor = ~ part, ordering = ordering, local = local)
}

test_that("parallel decorrelation controls validate without changing serial defaults", {
  resolve <- SSN2:::get_decorrelate_local
  expect_false(resolve(FALSE, 20)$parallel)
  expect_false(resolve(TRUE, 20)$parallel)
  expect_identical(resolve(list(method = "all", parallel = TRUE, ncores = 2), 20),
                   list(method = "all", parallel = TRUE, ncores = 2L))
  for (value in list(NA, 1, "TRUE", c(TRUE, FALSE))) {
    expect_error(resolve(list(parallel = value), 20), "local\\$parallel")
  }
  for (value in list(0, NA, Inf, 1.5, "2", c(1, 2), .Machine$integer.max + 1)) {
    expect_error(resolve(list(ncores = value), 20), "local\\$ncores")
  }
  expect_error(resolve(list(parallel = TRUE, parallel = FALSE), 20), "named list")
  expect_error(resolve(list(ncore = 2), 20), "named list")
})

test_that("parallel observed and prediction transforms match serial stream algebra", {
  set.seed(2)
  network <- parallel_decorrelation_fixture()
  for (method in c("all", "covariance")) {
    serial <- parallel_decorrelation_fit(network, list(method = method, size = 4))
    parallel <- parallel_decorrelation_fit(network,
      list(method = method, size = 4, parallel = TRUE, ncores = 2))
    expect_equal(parallel$tX, serial$tX, tolerance = 1e-12)
    expect_equal(parallel$ty, serial$ty, tolerance = 1e-12)
    expect_identical(parallel$rows, serial$rows)
    prediction <- ssn_decorrelate_newdata(serial, "pred1km")
    prediction_parallel <- ssn_decorrelate_newdata(parallel, "pred1km")
    for (name in c("tX_newdata", "yscale", "yoffset", "offset", "rows")) {
      expect_equal(prediction_parallel[[name]], prediction[[name]], tolerance = 1e-12)
    }
    expect_equal(ssn_recorrelate_newdata(prediction_parallel, rep(.3, 5)),
                 ssn_recorrelate_newdata(prediction, rep(.3, 5)), tolerance = 1e-12)
    single_worker <- ssn_decorrelate_newdata(serial, "pred1km",
      local = list(method = method, size = 4, parallel = TRUE, ncores = 1))
    expect_equal(single_worker$tX_newdata, prediction$tX_newdata, tolerance = 1e-12)
  }
})

test_that("parallel decorrelation supports indexed distance files and missing responses", {
  set.seed(2)
  network <- parallel_decorrelation_fixture()
  root <- withr::local_tempdir()
  stopifnot(file.copy(network$path, root, recursive = TRUE, copy.mode = FALSE))
  network$path <- file.path(root, basename(network$path))
  ssn_create_bigdist(network, predpts = "pred1km", among_predpts = TRUE,
                     overwrite = TRUE, no_cores = 1, verbose = FALSE)
  network$obs$y[3] <- NA_real_
  serial <- parallel_decorrelation_fit(network, list(size = 3))
  parallel <- parallel_decorrelation_fit(network, list(size = 3, parallel = TRUE, ncores = 2))
  expect_identical(SSN2:::get_decorrelate_observed_backend(parallel$covariance_fit), "bigdata")
  expect_equal(parallel$tX, serial$tX, tolerance = 1e-12)
  expect_equal(parallel$ty, serial$ty, tolerance = 1e-12)
  for (name in c("pred1km", ".missing")) {
    a <- ssn_decorrelate_newdata(serial, name)
    b <- ssn_decorrelate_newdata(parallel, name)
    expect_equal(a$tX_newdata, b$tX_newdata, tolerance = 1e-12)
    expect_equal(a$yscale, b$yscale, tolerance = 1e-12)
    expect_equal(a$yoffset, b$yoffset, tolerance = 1e-12)
    expect_identical(a$rows, b$rows)
  }
})

test_that("parallel prediction chunks remain bounded and preserve empty output shapes", {
  network <- parallel_decorrelation_fixture()
  fit <- parallel_decorrelation_fit(network, list(size = 3))
  seen <- integer()
  workers <- NULL
  testthat::local_mocked_bindings(
    run_pred_dispatch = function(fn, data_list, local_list, ...) {
      seen <<- lengths(data_list)
      workers <<- local_list$ncores
      lapply(data_list, function(rows) list(tX_newdata = matrix(0, length(rows), 2),
        yscale = rep(1, length(rows)), yoffset = rep(0, length(rows))))
    }, .package = "SSN2")
  settings <- SSN2:::get_decorrelate_local(list(size = 3, parallel = TRUE, ncores = 2), 12)
  value <- SSN2:::get_decorrelate_covariance_prediction(fit, "pred1km",
    network$preds$pred1km[rep(1, 205), ], matrix(1, 205, 2), settings)
  expect_true(all(seen <= 100L))
  expect_identical(sum(seen), 205L)
  expect_identical(dim(value$tX_newdata), c(205L, 2L))
  expect_lte(workers, length(seen))
  settings$ncores <- 4L
  SSN2:::get_decorrelate_covariance_prediction(fit, "pred1km",
    network$preds$pred1km, matrix(1, 5, 2), settings)
  expect_identical(workers, 3L)
  empty <- SSN2:::get_decorrelate_covariance_prediction(fit, "pred1km",
    network$preds$pred1km[FALSE, ], matrix(numeric(), 0, 2), settings)
  expect_identical(dim(empty$tX_newdata), c(0L, 2L))
  expect_length(empty$yscale, 0)
})

test_that("parallel decorrelation closes its cluster after a worker error", {
  network <- parallel_decorrelation_fixture()
  serial <- parallel_decorrelation_fit(network, list(size = 3))
  cl <- NULL
  original <- SSN2:::make_ssn_cluster
  testthat::local_mocked_bindings(make_ssn_cluster = function(ncores) {
    cl <<- original(ncores)
    cl
  }, .package = "SSN2")
  serial$covariance_fit$ssn.object$path <- tempfile("missing-distances-")
  expect_error(ssn_decorrelate_newdata(serial, "pred1km",
    local = list(size = 3, parallel = TRUE, ncores = 2)))
  expect_false(is.null(cl))
  expect_error(parallel::clusterCall(cl, identity, TRUE))
})

test_that("parallel transforms preserve seeded ordering and coincident Euclidean sites", {
  network <- parallel_decorrelation_fixture()
  set.seed(2)
  serial <- parallel_decorrelation_fit(network, list(size = 3), ordering = "random")
  set.seed(2)
  parallel <- parallel_decorrelation_fit(network,
    list(size = 3, parallel = TRUE, ncores = 2), ordering = "random")
  expect_equal(parallel$tX, serial$tX, tolerance = 1e-12)
  expect_equal(parallel$ty, serial$ty, tolerance = 1e-12)
  sf::st_geometry(network$obs)[2] <- sf::st_geometry(network$obs)[1]
  fit <- function(local) ssn_decorrelate_data(y ~ x, network,
    euclid_params = euclid_params("exponential", .8, 12000),
    nugget_params = nugget_params("nugget", .2), ordering = "pid", local = local)
  serial <- fit(FALSE)
  parallel <- fit(list(method = "all", parallel = TRUE, ncores = 2))
  expect_equal(parallel$tX, serial$tX, tolerance = 1e-12)
  expect_equal(parallel$ty, serial$ty, tolerance = 1e-12)
})
