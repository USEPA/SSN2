skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("unconditional local activation and Vecchia defaults respect their boundaries", {
  resolve <- function(local, n = 45L) SSN2:::get_ssn_simulate_local(local, n, NULL)
  for (n in c(1L, 5000L)) {
    expect_message(exact <- resolve(NULL, n), NA)
    expect_true(exact$exact)
  }
  expect_message(automatic <- resolve(NULL, 5001L), "exceeds 5,000")
  expect_identical(automatic, resolve(TRUE, 5001L))
  for (n in c(5000L, 5001L)) {
    expect_message(exact <- resolve(FALSE, n), NA)
    expect_true(exact$exact)
  }
  for (local in list(TRUE, list())) {
    settings <- resolve(local)
    expect_false(settings$exact)
    expect_identical(settings$method, "covariance")
    expect_equal(settings$size, 30L)
    expect_identical(settings$ordering, "pid")
  }
  for (local in list(list(method = "covariance"), list(size = 8),
                    list(method = "covariance", size = 8, ordering = "none"))) {
    settings <- resolve(local)
    explicit <- resolve(modifyList(list(method = "covariance", size = 30L, ordering = "pid"), local))
    expect_identical(settings[names(explicit)], explicit)
  }
  for (local in list(1, "vecchia", NA)) {
    expect_error(resolve(local), "local must be")
  }
  expect_error(resolve(list(approximation = "bogus")), "Unknown local option")
  expect_error(resolve(list(method = "distance")), "local\\$method must be")
  for (size in list(0, -1, 1.5, Inf, NA_real_, "30", c(1, 2))) {
    expect_error(resolve(list(size = size)), "local\\$size must be a positive integer")
  }
  for (size in list(8, 8L)) {
    expect_identical(resolve(list(size = size))$size, size)
  }
})

test_that("unconditional Vecchia uses the shared ordering validation", {
  for (ordering in list(NULL, "pid", "none", "random", "maxmin", "middleout",
                        "outsidein", "coordinate", "grts")) {
    settings <- SSN2:::get_local_vecchia_settings(list(ordering = ordering))
    expect_identical(settings$ordering, SSN2:::get_decorrelate_ordering(ordering))
  }
  expect_identical(
    SSN2:::get_local_vecchia_settings(list(), default_ordering = "none")$ordering,
    "none"
  )
  for (ordering in list("bogus", NA_character_, character(), c("pid", "none"))) {
    expected <- tryCatch(SSN2:::get_decorrelate_ordering(ordering), error = identity)
    expect_error(SSN2:::get_local_vecchia_settings(list(ordering = ordering)),
                 conditionMessage(expected), fixed = TRUE)
    expect_error(SSN2:::get_ssn_simulate_local(list(ordering = ordering), 45L, NULL),
                 conditionMessage(expected), fixed = TRUE)
  }
})

test_that("unconditional removed controls error and all dispatch stays quiet", {
  controls <- list(method_base = "all", size_base = 10, reorder_base = "none",
                   size_new = 5, kmeans_new = FALSE, parallel = TRUE, ncores = 2)
  for (control in names(controls)) {
    expect_error(SSN2:::get_ssn_simulate_local(controls[control], 45L, NULL),
                 "Unknown local option")
  }
  expect_error(SSN2:::get_ssn_simulate_local(list(parallel = FALSE), 45L, NULL),
               "Unknown local option")
  for (n in c(2000L, 2001L)) {
    expect_warning(settings <- SSN2:::get_ssn_simulate_local(list(method = "all"), n, NULL), NA)
    expect_true(settings$exact)
  }
  expect_warning(SSN2:::get_ssn_simulate_local(TRUE, 2001L, NULL), NA)
})

test_that("automatic and full-neighbor simulation dispatch to the sequential engine", {
  calls <- list()
  local_mocked_bindings(
    get_ssn_simulate_vecchia = function(covariance_fit, samples, conditioning) {
      calls[[length(calls) + 1L]] <<- conditioning
      matrix(0, NROW(covariance_fit$ssn.object$obs), samples)
    },
    get_dist_object = function(...) stop("unexpected dense dispatch"),
    .package = "SSN2"
  )
  args <- list(
    ssn.object = list(obs = data.frame(pid = seq_len(5001L))),
    tailup_params = tailup_params("none", de = 0, range = 1),
    taildown_params = taildown_params("none", de = 0, range = 1),
    euclid_params = euclid_params("exponential", de = 1, range = 1),
    nugget_params = nugget_params("nugget", nugget = 1), mean = 2, samples = 2
  )
  expect_message(automatic <- do.call(ssn_rnorm, args), "exceeds 5,000")
  expect_message(explicit <- do.call(ssn_rnorm, c(args, list(local = TRUE))), NA)
  expect_identical(automatic, explicit)
  expect_identical(dim(automatic), c(5001L, 2L))
  expect_true(all(automatic == 2))
  expect_identical(calls[[1]], calls[[2]])
  expect_false(calls[[1]]$exact)
  expect_equal(calls[[1]]$size, 30L)
  expect_identical(calls[[1]]$ordering, "pid")

  args$ssn.object$obs <- args$ssn.object$obs[1:3, , drop = FALSE]
  args$samples <- 1
  full <- do.call(ssn_rnorm, c(args, list(local = list(method = "covariance", size = 2))))
  expect_identical(full, rep(2, 3))
  expect_length(calls, 3L)
  expect_identical(calls[[3]]$method, "covariance")
  expect_equal(calls[[3]]$size, 2)
})

test_that("all unconditional families have seeded default and explicit Vecchia equivalence", {
  skip_if_not_installed("statmod")
  set.seed(2)
  copy_lsn_to_temp()
  mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"), overwrite = TRUE)
  ssn_create_distmat(mf04p, overwrite = TRUE)
  args <- list(
    ssn.object = mf04p,
    tailup_params = tailup_params("exponential", de = 0.2, range = 400),
    taildown_params = taildown_params("exponential", de = 0.1, range = 300),
    euclid_params = euclid_params("none", de = 0, range = 1),
    nugget_params = nugget_params("nugget", nugget = 0.1),
    additive = "afvArea", samples = 3
  )
  families <- c("gaussian", "poisson", "nbinomial", "binomial", "beta", "Gamma",
                "inverse.gaussian")
  for (family in families) {
    set.seed(2)
    explicit <- do.call(ssn_simulate, c(list(family = family), args,
      list(local = list(method = "covariance", size = 30,
                        ordering = "pid"))))
    expect_identical(dim(explicit), c(NROW(mf04p$obs), 3L), info = family)
    expect_true(all(is.finite(explicit)), info = family)
    for (local in list(TRUE, list(), list(method = "covariance", size = 30))) {
      set.seed(2)
      default <- do.call(ssn_simulate, c(list(family = family), args, list(local = local)))
      expect_identical(default, explicit, info = family)
    }
  }
})

test_that("all predecessors and dense simulation agree for fixed ordering and draws", {
  set.seed(2)
  copy_lsn_to_temp()
  mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"), overwrite = TRUE)
  ssn_create_distmat(mf04p, overwrite = TRUE)
  args <- list(
    ssn.object = mf04p,
    tailup_params = tailup_params("exponential", de = 0.6, range = 400),
    taildown_params = taildown_params("exponential", de = 0.4, range = 300),
    euclid_params = euclid_params("none", de = 0, range = 1),
    nugget_params = nugget_params("nugget", nugget = 0.1),
    additive = "afvArea", samples = 3, mean = 2
  )
  set.seed(2)
  dense <- do.call(ssn_rnorm, c(args, list(local = FALSE)))
  set.seed(2)
  sequential <- do.call(ssn_rnorm,
    c(args, list(local = list(method = "covariance", size = NROW(mf04p$obs) - 1L,
                             ordering = "none"))))
  expect_equal(sequential, dense, tolerance = 1e-10)
  set.seed(2)
  all <- do.call(ssn_rnorm, c(args, list(local = list(method = "all"))))
  expect_identical(all, dense)
})

test_that("co-located Euclidean covariance retains the nugget floor", {
  set.seed(2)
  copy_lsn_to_temp()
  mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"), overwrite = TRUE)
  mf04p$obs <- mf04p$obs[1:3, ]
  sf::st_geometry(mf04p$obs) <- rep(sf::st_geometry(mf04p$obs)[1], 3L)
  args <- list(
    ssn.object = mf04p,
    tailup_params = tailup_params("none", de = 0, range = 1),
    taildown_params = taildown_params("none", de = 0, range = 1),
    euclid_params = euclid_params("exponential", de = 1, range = 400),
    samples = 3
  )
  for (nugget in c(0, 0.1)) {
    args$nugget_params <- nugget_params("nugget", nugget = nugget)
    params <- list(tailup = args$tailup_params, taildown = args$taildown_params,
                   euclid = args$euclid_params, nugget = args$nugget_params, randcov = NULL)
    fit <- SSN2:::get_ssn_simulate_covariance_fit(mf04p, params, NULL, FALSE, NULL)
    covariance <- SSN2:::get_decorrelate_observed_covariance(fit, mf04p$obs)
    expect_equal(unname(covariance), matrix(1, 3, 3) + diag(max(nugget, 1e-4), 3),
                 tolerance = 1e-12)
  }
})
