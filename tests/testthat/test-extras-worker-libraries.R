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
