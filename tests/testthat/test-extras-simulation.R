skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# ssn simulate local
test_that("local = list(approximation = 'vecchia', method = 'all') reproduces the exact dense path within Monte Carlo error", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, overwrite = TRUE)
  ssn_create_distmat(mf04p, overwrite = TRUE)

  tailup <- tailup_params("exponential", de = 0.6, range = 400)
  taildown <- taildown_params("exponential", de = 0.4, range = 300)
  none_eu <- euclid_params("none", de = 0, range = 1)
  nugget <- nugget_params("nugget", nugget = 0.1)

  set.seed(2)
  draws_exact <- ssn_rnorm(mf04p,
    tailup_params = tailup, taildown_params = taildown,
    euclid_params = none_eu, nugget_params = nugget, additive = "afvArea", samples = 20000
  )
  set.seed(2)
  draws_all <- ssn_rnorm(mf04p,
    tailup_params = tailup, taildown_params = taildown,
    euclid_params = none_eu, nugget_params = nugget, additive = "afvArea", samples = 20000,
    local = list(approximation = "vecchia", method = "all")
  )
  expect_equal(dim(draws_exact), dim(draws_all))
  emp_cov_exact <- cov(t(draws_exact))
  emp_cov_all <- cov(t(draws_all))
  rel_diff <- max(abs(emp_cov_exact - emp_cov_all)) / max(abs(emp_cov_exact))
  expect_true(rel_diff < 0.1)

  # exact (non-Monte-Carlo) agreement: the covariance-request helper the
  # sequential engine is built on reproduces the dense get_cov_matrix()
  # analytically, not merely empirically
  params_object <- list(tailup = tailup, taildown = taildown, euclid = none_eu, nugget = nugget, randcov = NULL)
  covariance_fit <- SSN2:::get_ssn_simulate_covariance_fit(mf04p, params_object, "afvArea", FALSE, NULL)
  full_cov_via_helper <- SSN2:::get_decorrelate_observed_covariance(covariance_fit, mf04p$obs)

  initial_object <- SSN2:::get_initial_object("exponential", "exponential", "none", "nugget", NULL, NULL, NULL, NULL)
  dist_object <- SSN2:::get_dist_object(mf04p, initial_object, "afvArea", FALSE)
  de_scale <- sum(tailup[["de"]], taildown[["de"]], none_eu[["de"]])
  full_cov_dense <- as.matrix(SSN2:::get_cov_matrix(params_object, dist_object, NULL, NULL, FALSE, de_scale, diagtol = 0))
  expect_equal(unname(as.matrix(full_cov_via_helper)), unname(full_cov_dense), tolerance = 1e-10)
})

test_that("anisotropic simulation matches the fitted covariance for both engines", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, overwrite = TRUE)
  ssn_create_distmat(mf04p, overwrite = TRUE)

  euclid_anis <- euclid_params("exponential", de = 1, range = 500, rotate = 1.0, scale = 0.4)
  nugget <- nugget_params("nugget", nugget = 0.05)
  none_tu <- tailup_params("none", de = 0, range = 1)
  none_td <- taildown_params("none", de = 0, range = 1)

  params_object <- list(tailup = none_tu, taildown = none_td, euclid = euclid_anis, nugget = nugget, randcov = NULL)
  covariance_fit <- SSN2:::get_ssn_simulate_covariance_fit(mf04p, params_object, NULL, TRUE, NULL)
  correct_cov <- SSN2:::get_decorrelate_observed_covariance(covariance_fit, mf04p$obs)

  set.seed(2)
  draws_all_anis <- ssn_rnorm(mf04p,
    tailup_params = none_tu, taildown_params = none_td,
    euclid_params = euclid_anis, nugget_params = nugget, samples = 30000,
    local = list(approximation = "vecchia", method = "all")
  )
  emp_cov_anis <- cov(t(draws_all_anis))
  rel_diff_local <- max(abs(emp_cov_anis - correct_cov)) / max(abs(correct_cov))
  expect_true(rel_diff_local < 0.1)

  set.seed(2)
  draws_dense_anis <- ssn_rnorm(mf04p,
    tailup_params = none_tu, taildown_params = none_td,
    euclid_params = euclid_anis, nugget_params = nugget, samples = 30000,
    local = FALSE
  )
  emp_cov_dense_anis <- cov(t(draws_dense_anis))
  rel_diff_dense <- max(abs(emp_cov_dense_anis - correct_cov)) / max(abs(correct_cov))
  expect_true(rel_diff_dense < 0.1)

  set.seed(2)
  expected <- t(chol(correct_cov)) %*% matrix(rnorm(NROW(mf04p$obs) * 3), NROW(mf04p$obs), 3)
  expect_equal(unname(draws_dense_anis[, 1:3]), unname(expected), tolerance = 1e-12)
})

test_that("conditional variance shrinks monotonically as neighborhood size grows (exact, no Monte Carlo)", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, overwrite = TRUE)
  ssn_create_distmat(mf04p, overwrite = TRUE)

  tailup <- tailup_params("exponential", de = 0.6, range = 400)
  taildown <- taildown_params("exponential", de = 0.4, range = 300)
  none_eu <- euclid_params("none", de = 0, range = 1)
  nugget <- nugget_params("nugget", nugget = 0.1)
  params_object <- list(tailup = tailup, taildown = taildown, euclid = none_eu, nugget = nugget, randcov = NULL)
  covariance_fit <- SSN2:::get_ssn_simulate_covariance_fit(mf04p, params_object, "afvArea", FALSE, NULL)

  rows <- SSN2:::get_decorrelate_rows(mf04p$obs)
  permutation <- SSN2:::get_decorrelate_order(rows, "pid")
  position <- 30L
  current <- permutation[[position]]
  candidates <- permutation[seq_len(position - 1L)]
  variance <- SSN2:::get_decorrelate_marginal_variance(covariance_fit, mf04p$obs[current, , drop = FALSE])
  cross_covariance <- SSN2:::get_decorrelate_observed_cross_covariance(
    covariance_fit, mf04p$obs[current, , drop = FALSE], mf04p$obs[candidates, , drop = FALSE]
  )

  cond_vars <- vapply(c(1L, 3L, 10L, 29L), function(sz) {
    keep <- SSN2:::get_decorrelate_covariance_neighbors(cross_covariance, sz)
    neighbor_index <- candidates[keep]
    covariance_neighbors <- SSN2:::get_decorrelate_observed_covariance(covariance_fit, mf04p$obs[neighbor_index, , drop = FALSE])
    chol_pool <- chol(covariance_neighbors)
    c_vec <- cross_covariance[keep]
    w <- backsolve(chol_pool, forwardsolve(t(chol_pool), c_vec))
    max(variance - sum(w * c_vec), 0)
  }, numeric(1))

  # more neighbors can only explain equal or more variance -> cond_var is
  # non-increasing as size grows
  expect_true(all(diff(cond_vars) <= 1e-8))
  # and the size = 29 (all candidates) case must be strictly less than
  # size = 1 for this fixture (genuinely informative neighbors exist)
  expect_true(cond_vars[[4]] < cond_vars[[1]])
})

test_that("local unconditional simulation supports random effects and partition factors", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, overwrite = TRUE)
  ssn_create_distmat(mf04p, overwrite = TRUE)
  mf04p$obs$netID <- as.factor(mf04p$obs$netID)

  tailup <- tailup_params("exponential", de = 0.5, range = 400)
  taildown <- taildown_params("exponential", de = 0.3, range = 300)
  none_eu <- euclid_params("none", de = 0, range = 1)
  nugget <- nugget_params("nugget", nugget = 0.1)
  randcov <- spmodel::randcov_params(netID = 0.4)

  set.seed(2)
  draws_rc <- ssn_rnorm(mf04p,
    tailup_params = tailup, taildown_params = taildown,
    euclid_params = none_eu, nugget_params = nugget,
    randcov_params = randcov, additive = "afvArea", samples = 8000,
    local = list(approximation = "vecchia", method = "covariance", size = 20)
  )
  expect_true(all(is.finite(draws_rc)))
  analytic_marginal_var <- tailup[["de"]] + taildown[["de"]] + nugget[["nugget"]] + as.numeric(randcov)
  emp_marginal_var <- mean(apply(draws_rc, 1, var))
  expect_equal(emp_marginal_var, analytic_marginal_var, tolerance = 0.05)

  set.seed(2)
  draws_pf <- ssn_rnorm(mf04p,
    tailup_params = tailup, taildown_params = taildown,
    euclid_params = none_eu, nugget_params = nugget,
    partition_factor = ~netID, additive = "afvArea", samples = 500,
    local = list(approximation = "vecchia", method = "covariance", size = 15)
  )
  expect_true(all(is.finite(draws_pf)))
  expect_equal(dim(draws_pf), c(NROW(mf04p$obs), 500))
})

test_that("local unconditional simulation: dispatch, dims, RNG, rejection sweep, duplicate locations, family forwarding", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, overwrite = TRUE)
  ssn_create_distmat(mf04p, overwrite = TRUE)

  tailup <- tailup_params("exponential", de = 0.6, range = 400)
  taildown <- taildown_params("exponential", de = 0.4, range = 300)
  none_eu <- euclid_params("none", de = 0, range = 1)
  nugget <- nugget_params("nugget", nugget = 0.1)
  n <- NROW(mf04p$obs)

  # samples = 1 shape
  d_one <- ssn_rnorm(mf04p,
    tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
    nugget_params = nugget, additive = "afvArea", samples = 1,
    local = list(approximation = "vecchia", method = "covariance", size = 10)
  )
  expect_true(is.numeric(d_one) && !is.matrix(d_one))
  expect_equal(length(d_one), n)

  # RNG reproducibility
  set.seed(2)
  d1 <- ssn_rnorm(mf04p,
    tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
    nugget_params = nugget, additive = "afvArea", samples = 5,
    local = list(approximation = "vecchia", method = "covariance", size = 10)
  )
  set.seed(2)
  d2 <- ssn_rnorm(mf04p,
    tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
    nugget_params = nugget, additive = "afvArea", samples = 5,
    local = list(approximation = "vecchia", method = "covariance", size = 10)
  )
  expect_identical(d1, d2)

  # default path unaffected by local's existence (small n, no auto-trigger)
  set.seed(2)
  d_false <- ssn_rnorm(mf04p,
    tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
    nugget_params = nugget, additive = "afvArea", samples = 3, local = FALSE
  )
  set.seed(2)
  d_default <- ssn_rnorm(mf04p,
    tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
    nugget_params = nugget, additive = "afvArea", samples = 3
  )
  expect_identical(d_false, d_default)

  # rejection sweep
  expect_error(
    ssn_rnorm(mf04p, tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
      nugget_params = nugget, additive = "afvArea", samples = 5, local = list(approximation = "vecchia", method = "bogus")),
    "local\\$method must be"
  )
  expect_error(
    ssn_rnorm(mf04p, tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
      nugget_params = nugget, additive = "afvArea", samples = 5,
      local = list(approximation = "vecchia", size = -1)),
    "local\\$size must be"
  )
  expect_error(
    ssn_rnorm(mf04p, tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
      nugget_params = nugget, additive = "afvArea", samples = 5, local = list(method = "covariance")),
    "local now defaults to a low-rank approximation"
  )

  # a partial vecchia list initializes defaults for the other list elements
  expect_identical(
    SSN2:::get_ssn_simulate_local(list(approximation = "vecchia", method = "covariance"), n, mf04p$obs),
    SSN2:::get_ssn_simulate_local(list(approximation = "vecchia"), n, mf04p$obs)
  )

  # duplicate locations do not crash
  mf04p_dup <- mf04p
  mf04p_dup$obs <- rbind(mf04p_dup$obs, mf04p_dup$obs[1, ])
  draws_dup <- ssn_rnorm(mf04p_dup,
    tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
    nugget_params = nugget, additive = "afvArea", samples = 100,
    local = list(approximation = "vecchia", method = "covariance", size = 10)
  )
  expect_true(all(is.finite(draws_dup)))
  expect_equal(dim(draws_dup), c(n + 1, 100))

  # ssn_rpois()/ssn_rbinom() forward local automatically (existing
  # match.call() re-dispatch, no new code needed in either function)
  set.seed(2)
  pois_local <- ssn_rpois(mf04p,
    tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
    nugget_params = nugget, additive = "afvArea", samples = 10,
    local = list(approximation = "vecchia", method = "covariance", size = 10)
  )
  expect_equal(dim(pois_local), c(n, 10))
  expect_true(all(pois_local == round(pois_local)) && all(is.finite(pois_local)))

  set.seed(2)
  binom_local <- ssn_rbinom(mf04p,
    tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
    nugget_params = nugget, additive = "afvArea", samples = 10, size = 5,
    local = list(approximation = "vecchia", method = "covariance", size = 10)
  )
  expect_true(all(binom_local >= 0 & binom_local <= 5))
})

test_that("local unconditional simulation: different orderings all produce valid approximations", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, overwrite = TRUE)
  ssn_create_distmat(mf04p, overwrite = TRUE)

  tailup <- tailup_params("exponential", de = 0.6, range = 400)
  taildown <- taildown_params("exponential", de = 0.4, range = 300)
  none_eu <- euclid_params("none", de = 0, range = 1)
  nugget <- nugget_params("nugget", nugget = 0.1)
  params_object <- list(tailup = tailup, taildown = taildown, euclid = none_eu, nugget = nugget, randcov = NULL)
  covariance_fit <- SSN2:::get_ssn_simulate_covariance_fit(mf04p, params_object, "afvArea", FALSE, NULL)
  true_cov <- SSN2:::get_decorrelate_observed_covariance(covariance_fit, mf04p$obs)

  for (ord in c("pid", "none", "random", "maxmin", "middleout", "outsidein", "coordinate", "grts")) {
    set.seed(2)
    d <- ssn_rnorm(mf04p,
      tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
      nugget_params = nugget, additive = "afvArea", samples = 8000,
      local = list(approximation = "vecchia", method = "covariance", size = 8, ordering = ord)
    )
    emp <- cov(t(d))
    rel <- max(abs(emp - true_cov)) / max(abs(true_cov))
    expect_true(rel < 0.15)
  }

  expect_error(
    ssn_rnorm(mf04p, tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
      nugget_params = nugget, additive = "afvArea", samples = 5,
      local = list(approximation = "vecchia", method = "covariance", size = 5, ordering = "bogus")),
    "local\\$ordering must be"
  )
})

test_that("local low-rank simulation (approximation = 'low-rank', the new default) approximates the exact distribution", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, overwrite = TRUE)
  ssn_create_distmat(mf04p, overwrite = TRUE)

  tailup <- tailup_params("exponential", de = 0.6, range = 400)
  taildown <- taildown_params("exponential", de = 0.4, range = 300)
  none_eu <- euclid_params("none", de = 0, range = 1)
  nugget <- nugget_params("nugget", nugget = 0.1)
  params_object <- list(tailup = tailup, taildown = taildown, euclid = none_eu, nugget = nugget, randcov = NULL)
  covariance_fit <- SSN2:::get_ssn_simulate_covariance_fit(mf04p, params_object, "afvArea", FALSE, NULL)
  true_cov <- SSN2:::get_decorrelate_observed_covariance(covariance_fit, mf04p$obs)

  # method_base = "all" disables subsetting entirely -- should reproduce the
  # exact distribution just like local$approximation = "vecchia"/method = "all"
  set.seed(2)
  d_all <- ssn_rnorm(mf04p, tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
    nugget_params = nugget, additive = "afvArea", samples = 8000,
    local = list(approximation = "low-rank", method_base = "all"))
  rel_all <- max(abs(cov(t(d_all)) - true_cov)) / max(abs(true_cov))
  expect_true(rel_all < 0.15)

  # genuine base+block subsetting still approximates the true covariance
  set.seed(2)
  d_block <- ssn_rnorm(mf04p, tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
    nugget_params = nugget, additive = "afvArea", samples = 8000,
    local = list(approximation = "low-rank", method_base = "base", size_base = 20, size_new = 10))
  expect_true(all(is.finite(d_block)))
  rel_block <- max(abs(cov(t(d_block)) - true_cov)) / max(abs(true_cov))
  expect_true(rel_block < 0.35)

  # local = TRUE now resolves to low-rank, matching spmodel's own default
  expect_identical(SSN2:::get_ssn_simulate_local(TRUE, 45, mf04p$obs)$approximation, "low-rank")

  # explicit vecchia still works and is unaffected by the new default
  set.seed(2)
  d_vecchia <- ssn_rnorm(mf04p, tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
    nugget_params = nugget, additive = "afvArea", samples = 8000,
    local = list(approximation = "vecchia", method = "covariance", size = 8))
  rel_vecchia <- max(abs(cov(t(d_vecchia)) - true_cov)) / max(abs(true_cov))
  expect_true(rel_vecchia < 0.15)

  # a pre-approximation local shape (method/size with no approximation) now
  # errors clearly instead of silently being treated as low-rank
  expect_error(
    ssn_rnorm(mf04p, tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
      nugget_params = nugget, additive = "afvArea", samples = 1,
      local = list(method = "covariance", size = 8)),
    "local now defaults to a low-rank approximation"
  )

  # rejection sweep for low-rank's own settings
  expect_error(
    ssn_rnorm(mf04p, tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
      nugget_params = nugget, additive = "afvArea", samples = 1,
      local = list(approximation = "low-rank", method_base = "bogus")),
    "method_base must be"
  )
  expect_error(
    ssn_rnorm(mf04p, tailup_params = tailup, taildown_params = taildown, euclid_params = none_eu,
      nugget_params = nugget, additive = "afvArea", samples = 1,
      local = list(approximation = "low-rank", reorder_base = "bogus")),
    "reorder_base must be"
  )
})

