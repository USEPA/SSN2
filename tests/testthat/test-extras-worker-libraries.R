skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("workers inherit runtime library paths before loading SSN2", {
  set.seed(2)
  library_path <- withr::local_tempdir()
  withr::local_libpaths(c(library_path, .libPaths()))
  parent_paths <- .libPaths()
  parent_package <- getNamespaceInfo(asNamespace("SSN2"), "path")
  cl <- SSN2:::make_ssn_cluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  workers <- parallel::clusterCall(cl, function() list(
    libraries = .libPaths(),
    package = getNamespaceInfo(asNamespace("SSN2"), "path"),
    pid = Sys.getpid()
  ))
  for (worker in workers) {
    expect_identical(worker$libraries, parent_paths)
    expect_identical(worker$package, parent_package)
    expect_false(worker$pid == Sys.getpid())
  }
  expect_identical(.libPaths(), parent_paths)
})

test_that("block prediction and low-rank simulation load the parent's namespace", {
  set.seed(2)
  parent_package <- getNamespaceInfo(asNamespace("SSN2"), "path")
  if (!(requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("SSN2"))) {
    package_lib <- dirname(parent_package)
    withr::local_libpaths(setdiff(.libPaths(), package_lib))
    skip_if(package_lib %in% .libPaths(), "the package is in a mandatory base library")
  }
  parent_paths <- .libPaths()
  worker_info <- function() list(
    package = getNamespaceInfo(asNamespace("SSN2"), "path"),
    libraries = .libPaths(), pid = Sys.getpid()
  )
  environment(worker_info) <- baseenv()
  workers <- list()
  par_lapply <- parallel::parLapply
  local_mocked_bindings(parLapply = function(cl, X, fun, ...) {
    workers[[length(workers) + 1L]] <<- parallel::clusterCall(cl, worker_info)
    par_lapply(cl, X, fun, ...)
  }, .package = "parallel")

  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_initial = tailup_initial("exponential", de = 2, range = 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", nugget = 0.5, known = "given"),
    additive = "afvArea", ddf = "asymptotic")
  fit$ssn.object$preds$pred1km <- fit$ssn.object$preds$pred1km[1:12, , drop = FALSE]
  local <- list(method = "all", method_new = "basis", size_new = 4,
    ordering = "pid", chunk_size = 3, parallel = FALSE, ncores = 2)
  serial <- predict(fit, "pred1km", block = TRUE, se.fit = TRUE, local = local)
  local$parallel <- TRUE
  expect_identical(predict(fit, "pred1km", block = TRUE, se.fit = TRUE, local = local), serial)
  block_dispatches <- length(workers)
  expect_gt(block_dispatches, 0L)

  ssn <- mf04p
  ssn$obs <- ssn$obs[1:12, , drop = FALSE]
  args <- list(ssn.object = ssn, tailup_params = tailup_params("exponential", de = 2, range = 10000),
    taildown_params = taildown_params("none"), euclid_params = euclid_params("none"),
    nugget_params = nugget_params("nugget", nugget = 0.5), additive = "afvArea", samples = 3,
    local = list(approximation = "low-rank", size_base = 4, size_new = 4,
      reorder_base = "none", kmeans_new = FALSE, parallel = FALSE, ncores = 2))
  set.seed(2)
  serial <- do.call(ssn_rnorm, args)
  args$local$parallel <- TRUE
  set.seed(2)
  simulated <- do.call(ssn_rnorm, args)
  expect_identical(dim(simulated), c(12L, 3L))
  expect_true(all(is.finite(simulated)))
  expect_identical(simulated[1:4, ], serial[1:4, ])
  expect_gt(length(workers), block_dispatches)
  for (dispatch in workers) {
    expect_length(dispatch, 2L)
    for (worker in dispatch) {
      expect_identical(worker$package, parent_package)
      expect_identical(worker$libraries, parent_paths)
      expect_false(worker$pid == Sys.getpid())
    }
  }
  expect_identical(.libPaths(), parent_paths)
})

test_that("block and low-rank simulation workers are stopped after errors", {
  set.seed(2)
  stopped <- list()
  stop_cluster <- parallel::stopCluster
  local_mocked_bindings(stopCluster = function(cl) {
    stop_cluster(cl)
    stopped[[length(stopped) + 1L]] <<- cl
  }, .package = "parallel")
  fail <- function(...) stop("worker failure")
  environment(fail) <- baseenv()
  expect_error(get_block_chunk_apply(list(1, 2), fail, parallel = TRUE, ncores = 2),
    "worker failure")
  expect_length(stopped, 1L)

  ssn <- mf04p
  ssn$obs <- ssn$obs[1:12, , drop = FALSE]
  local_mocked_bindings(get_ssn_simulate_new_from_base = fail)
  expect_error(ssn_rnorm(ssn,
    tailup_params = tailup_params("exponential", de = 2, range = 10000),
    taildown_params = taildown_params("none"), euclid_params = euclid_params("none"),
    nugget_params = nugget_params("nugget", nugget = 0.5), additive = "afvArea", samples = 3,
    local = list(approximation = "low-rank", size_base = 4, size_new = 4,
      reorder_base = "none", kmeans_new = FALSE, parallel = TRUE, ncores = 2)),
    "worker failure")
  expect_length(stopped, 2L)
  for (cl in stopped) {
    for (worker in cl) expect_error(summary(worker$con), "invalid connection")
  }
})

test_that("installed workers load the parent's namespace outside library paths", {
  skip_if(requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("SSN2"),
          "requires an installed SSN2 namespace")
  set.seed(2)
  parent_package <- getNamespaceInfo(asNamespace("SSN2"), "path")
  package_lib <- dirname(parent_package)
  withr::local_libpaths(setdiff(.libPaths(), package_lib))
  skip_if(package_lib %in% .libPaths(), "the package is in a mandatory base library")
  parent_paths <- .libPaths()
  result <- SSN2:::run_pred_dispatch(function(x) x * 2, list(1, 2, 3),
    list(parallel = TRUE, ncores = 2))
  expect_equal(unlist(result), c(2, 4, 6))
  cl <- SSN2:::make_ssn_cluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  workers <- parallel::clusterCall(cl, function() list(
    libraries = .libPaths(),
    package = getNamespaceInfo(asNamespace("SSN2"), "path"),
    pid = Sys.getpid()
  ))
  for (worker in workers) {
    expect_identical(worker$libraries, parent_paths)
    expect_identical(worker$package, parent_package)
    expect_false(worker$pid == Sys.getpid())
  }
  expect_identical(.libPaths(), parent_paths)
})
