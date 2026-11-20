skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("anisotropic simulation matches the fitted covariance on the exact path", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, overwrite = TRUE)
  ssn_create_distmat(mf04p, overwrite = TRUE)

  euclid_anis <- euclid_params("exponential", de = 1, range = 500, rotate = 1.0, scale = 0.4)
  nugget <- nugget_params("nugget", nugget = 0.05)
  none_tu <- tailup_params("none", de = 0, range = 1)
  none_td <- taildown_params("none", de = 0, range = 1)

  covariance_fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, euclid_type = "exponential",
    euclid_initial = euclid_initial("exponential", 1, 500, rotate = 1, scale = .4, known = "given"),
    anisotropy = TRUE, nugget_initial = nugget_initial("nugget", .05, known = "given"), ddf = "asymptotic")
  correct_cov <- covmatrix(covariance_fit)
  set.seed(2)
  draws_dense_anis <- ssn_rnorm(mf04p,
    tailup_params = none_tu, taildown_params = none_td,
    euclid_params = euclid_anis, nugget_params = nugget, samples = 30000
  )
  emp_cov_dense_anis <- cov(t(draws_dense_anis))
  rel_diff_dense <- max(abs(emp_cov_dense_anis - correct_cov)) / max(abs(correct_cov))
  expect_true(rel_diff_dense < 0.1)

  set.seed(2)
  expected <- t(chol(correct_cov)) %*% matrix(rnorm(NROW(mf04p$obs) * 3), NROW(mf04p$obs), 3)
  expect_equal(unname(draws_dense_anis[, 1:3]), unname(expected), tolerance = 1e-12)
})

