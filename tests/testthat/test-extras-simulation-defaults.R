skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("unconditional local activation and Vecchia defaults respect their boundaries", {
  resolve <- function(local, n = 45L) SSN2:::get_ssn_simulate_local(local, n, NULL)
  for (n in c(1L, 5000L)) {
    expect_message(exact <- resolve(NULL, n), NA)
    expect_identical(exact$method, "exact")
  }
  expect_message(automatic <- resolve(NULL, 5001L), "exceeds 5,000")
  expect_identical(automatic, resolve(TRUE, 5001L))
  for (n in c(5000L, 5001L)) {
    expect_message(exact <- resolve(FALSE, n), NA)
    expect_identical(exact$method, "exact")
  }
  for (local in list(TRUE, list(), list(approximation = "vecchia"))) {
    settings <- resolve(local)
    expect_identical(settings$approximation, "vecchia")
    expect_identical(settings$method, "covariance")
    expect_equal(settings$size, 30L)
    expect_identical(settings$ordering, "pid")
    expect_identical(settings$parallel, FALSE)
  }
  for (local in list(list(method = "covariance"), list(size = 8),
                    list(method = "covariance", size = 8, ordering = "none"))) {
    settings <- resolve(local)
    explicit <- resolve(c(list(approximation = "vecchia"), local))
    expect_identical(settings[names(explicit)], explicit)
  }
  for (local in list(1, "vecchia", NA)) {
    expect_error(resolve(local), "local must be")
  }
  expect_error(resolve(list(approximation = "bogus")), "local\\$approximation must be")
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

test_that("unconditional mixed controls and sequential warnings match the chosen approximation", {
  controls <- list(method_base = "all", size_base = 10, reorder_base = "none",
                   size_new = 5, kmeans_new = FALSE)
  for (control in names(controls)) {
    for (approximation in list(list(), list(approximation = "vecchia"))) {
      expect_error(
        SSN2:::get_ssn_simulate_local(c(approximation, controls[control]), 45L, NULL),
        'Base and block settings require local$approximation = "low-rank".', fixed = TRUE
      )
    }
  }
  expect_warning(
    settings <- SSN2:::get_ssn_simulate_local(list(parallel = TRUE), 45L, NULL),
    'local$parallel is not used when local$approximation = "vecchia"', fixed = TRUE
  )
  expect_identical(settings$approximation, "vecchia")
  expect_identical(settings$parallel, FALSE)
  expect_warning(SSN2:::get_ssn_simulate_local(list(parallel = FALSE), 45L, NULL), NA)
  expect_warning(
    boundary <- SSN2:::get_ssn_simulate_local(list(method = "all"), 2000L, NULL), NA
  )
  expect_warning(
    large <- SSN2:::get_ssn_simulate_local(list(method = "all"), 2001L, NULL),
    'local$method = "all" can be slow for more than a few thousand locations.', fixed = TRUE
  )
  expect_identical(boundary$method, "all")
  expect_identical(large$method, "all")
  expect_identical(large$approximation, "vecchia")
  expect_equal(boundary$size, 30L)
  expect_equal(large$size, 30L)
  expect_warning(SSN2:::get_ssn_simulate_local(TRUE, 2001L, NULL), NA)

  lowrank <- SSN2:::get_ssn_simulate_local(list(approximation = "low-rank"), 45L, NULL)
  expect_identical(lowrank$approximation, "low-rank")
  expect_identical(lowrank$method_base, "all")
  expect_equal(lowrank$size_base, 5000L)
  expect_equal(lowrank$size_new, 1000L)
  expect_identical(lowrank$reorder_base, "grts")
})

test_that("automatic and all-predecessor simulation dispatch to the sequential engine", {
  calls <- list()
  local_mocked_bindings(
    get_ssn_simulate_vecchia = function(covariance_fit, samples, conditioning) {
      calls[[length(calls) + 1L]] <<- conditioning
      matrix(0, NROW(covariance_fit$ssn.object$obs), samples)
    },
    get_ssn_simulate_lowrank = function(...) stop("unexpected low-rank dispatch"),
    get_dist_object = function(...) stop("unexpected dense dispatch"),
    .package = "SSN2"
  )
  args <- list(
    ssn.object = list(obs = data.frame(pid = seq_len(5001L))),
    tailup_params = tailup_params("none", de = 0, range = 1),
    taildown_params = taildown_params("none", de = 0, range = 1),
    euclid_params = euclid_params("none", de = 0, range = 1),
    nugget_params = nugget_params("nugget", nugget = 1), mean = 2, samples = 2
  )
  expect_message(automatic <- do.call(ssn_rnorm, args), "exceeds 5,000")
  expect_message(explicit <- do.call(ssn_rnorm, c(args, list(local = TRUE))), NA)
  expect_identical(automatic, explicit)
  expect_identical(dim(automatic), c(5001L, 2L))
  expect_true(all(automatic == 2))
  expect_identical(calls[[1]], calls[[2]])
  expect_identical(calls[[1]]$approximation, "vecchia")
  expect_equal(calls[[1]]$size, 30L)
  expect_identical(calls[[1]]$ordering, "pid")

  args$ssn.object$obs <- args$ssn.object$obs[1:3, , drop = FALSE]
  args$samples <- 1
  all <- do.call(ssn_rnorm, c(args, list(local = list(method = "all"))))
  expect_identical(all, rep(2, 3))
  expect_length(calls, 3L)
  expect_identical(calls[[3]]$method, "all")
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
      list(local = list(approximation = "vecchia", method = "covariance", size = 30,
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

test_that("explicit low-rank ignores Vecchia controls and preserves seeded draws", {
  set.seed(2)
  copy_lsn_to_temp()
  mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"), overwrite = TRUE)
  ssn_create_distmat(mf04p, overwrite = TRUE)
  local <- list(approximation = "low-rank", size_base = 10, size_new = 8,
                reorder_base = "none", kmeans_new = FALSE)
  mixed <- c(local, list(method = "bogus", size = -1, ordering = "bogus"))
  settings <- SSN2:::get_ssn_simulate_local(local, NROW(mf04p$obs), mf04p$obs)
  expect_identical(settings$index$base, seq_len(10L))
  expect_equal(sort(unlist(settings$index$new, use.names = FALSE)),
               seq.int(11L, NROW(mf04p$obs)))
  args <- list(
    ssn.object = mf04p,
    tailup_params = tailup_params("exponential", de = 0.2, range = 400),
    taildown_params = taildown_params("exponential", de = 0.1, range = 300),
    euclid_params = euclid_params("none", de = 0, range = 1),
    nugget_params = nugget_params("nugget", nugget = 0.1),
    additive = "afvArea", samples = 3
  )
  for (family in c("gaussian", "poisson")) {
    set.seed(2)
    explicit <- do.call(ssn_simulate, c(list(family = family), args, list(local = local)))
    set.seed(2)
    expect_warning(ignored <- do.call(ssn_simulate,
      c(list(family = family), args, list(local = mixed))), NA)
    expect_identical(ignored, explicit)
    expect_true(all(is.finite(explicit)))
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
    c(args, list(local = list(method = "all", ordering = "none"))))
  expect_equal(sequential, dense, tolerance = 1e-10)
  set.seed(2)
  lowrank_all <- do.call(ssn_rnorm,
    c(args, list(local = list(approximation = "low-rank", method_base = "all"))))
  expect_equal(lowrank_all, dense, tolerance = 1e-10)
})

test_that("unconditional low-rank rejects singular base and block covariances", {
  set.seed(2)
  covariance <- diag(3)
  fit <- list(ssn.object = list(obs = data.frame(index = 1:3)))
  local_mocked_bindings(
    get_decorrelate_observed_covariance = function(covariance_fit, data) {
      covariance[data$index, data$index, drop = FALSE]
    },
    get_decorrelate_observed_cross_covariance = function(covariance_fit, data, newdata) {
      covariance[data$index, newdata$index, drop = FALSE]
    },
    .package = "SSN2"
  )
  local <- list(method_base = "base", index = list(base = 1L, new = list(2:3)),
                parallel = FALSE)
  expect_true(all(is.finite(SSN2:::get_ssn_simulate_lowrank(fit, 2L, local))))

  covariance[1:2, 1:2] <- 1
  singular_base <- local
  singular_base$index <- list(base = 1:2, new = list(3L))
  expect_error(SSN2:::get_ssn_simulate_lowrank(fit, 2L, singular_base),
               "leading minor of order 2 is not positive")

  covariance <- diag(3)
  covariance[2:3, 2:3] <- 1
  expect_error(SSN2:::get_ssn_simulate_lowrank(fit, 2L, local),
               "leading minor of order 2 is not positive")
})

test_that("co-located Euclidean low-rank draws retain the covariance nugget floor", {
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
    for (size_base in c(2L, 1L)) {
      local <- list(approximation = "low-rank", size_base = size_base, size_new = 2L,
                    reorder_base = "none", kmeans_new = FALSE)
      set.seed(2)
      expect_warning(draws <- do.call(ssn_rnorm, c(args, list(local = local))), NA)
      expect_identical(dim(draws), c(3L, 3L))
      expect_true(all(is.finite(draws)))
    }
  }
})
