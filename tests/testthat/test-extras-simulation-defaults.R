skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

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
    covariance <- matrix(1, 3, 3) + diag(max(nugget, 1e-4), 3)
    set.seed(2)
    expected <- t(chol(covariance)) %*% matrix(rnorm(9), 3)
    set.seed(2)
    expect_equal(unname(do.call(ssn_rnorm, args)), expected, tolerance = 1e-12)
  }
})
