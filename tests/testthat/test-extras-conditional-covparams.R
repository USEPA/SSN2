skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# conditional covparams
test_that("conditional() simulate_covparams = TRUE works for ssn_lm()", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$preds$CapeHornSmall <- mf04p$preds$CapeHorn[1:8, ]
  ssn_create_distmat(mf04p, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE)

  fit <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    additive = "afvArea"
  )

  # basic dispatch / dimensions, output = "newdata"
  set.seed(2)
  draws_cp <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 40, simulate_covparams = TRUE)
  expect_true(is.matrix(draws_cp))
  expect_equal(dim(draws_cp), c(8, 40))
  expect_true(all(is.finite(draws_cp)))

  # output = "beta"
  set.seed(2)
  betas_cp <- conditional(fit, "CapeHornSmall", output = "beta", samples = 40, simulate_covparams = TRUE)
  expect_equal(dim(betas_cp), c(fit$p, 40))
  expect_true(all(is.finite(betas_cp)))

  R <- 3000
  set.seed(2)
  betas_cp_big <- conditional(fit, "CapeHornSmall", output = "beta", samples = R, simulate_covparams = TRUE)
  set.seed(2)
  betas_fixed_big <- conditional(fit, "CapeHornSmall", output = "beta", samples = R, simulate_covparams = FALSE)
  expect_true(all(apply(betas_cp_big, 1, sd) > apply(betas_fixed_big, 1, sd)))

  # default (simulate_covparams = FALSE) path is byte-identical whether or
  # not the argument is supplied explicitly -- unaffected by this feature's
  # existence
  set.seed(2)
  d1 <- conditional(fit, "CapeHornSmall", samples = 15)
  set.seed(2)
  d2 <- conditional(fit, "CapeHornSmall", samples = 15, simulate_covparams = FALSE)
  expect_identical(d1, d2)

  # rejection sweep
  expect_error(conditional(fit, "CapeHornSmall", simulate_covparams = "yes"), "simulate_covparams must be")
  expect_error(conditional(fit, "CapeHornSmall", output = "cov", samples = 5), "only be \"cov\", \"ssn\", \"tailup\", \"taildown\", \"euclid\", \"nugget\", or \"randcov\"")
  expect_error(conditional(fit, "CapeHornSmall", output = "bogus", samples = 5), "output must be")

  fit_fake_nonpd <- fit
  fit_fake_nonpd$ddf <- "satterthwaite"
  fit_fake_nonpd$vcov$cov <- NULL
  expect_error(
    conditional(fit_fake_nonpd, "CapeHornSmall", samples = 5, simulate_covparams = TRUE),
    "not available for this fit"
  )

  # performance warning above n = 500 observed rows (faking n avoids fitting
  # an actually-large model just for this check)
  fit_fake_big <- fit
  fit_fake_big$n <- 600
  expect_warning(
    conditional(fit_fake_big, "CapeHornSmall", samples = 2, simulate_covparams = TRUE),
    "exceedingly long"
  )
  expect_no_warning(conditional(fit, "CapeHornSmall", samples = 2, simulate_covparams = TRUE))
})

test_that("conditional() output = 'cov'/'ssn'/'tailup'/'nugget'/'randcov' match vcov()'s row names", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$obs$netID <- as.factor(mf04p$obs$netID)
  mf04p$preds$CapeHorn$netID <- as.factor(mf04p$preds$CapeHorn$netID)
  mf04p$preds$CapeHornSmall <- mf04p$preds$CapeHorn[1:8, ]
  ssn_create_distmat(mf04p, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE)

  fit <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    additive = "afvArea",
    random = ~ as.factor(netID)
  )

  vc_full <- vcov(fit, type = "cov")
  expect_false(is.null(vc_full))

  cov_out <- conditional(fit, "CapeHornSmall", output = "cov", samples = 25, simulate_covparams = TRUE)
  expect_equal(rownames(cov_out), rownames(vc_full))
  expect_equal(dim(cov_out), c(nrow(vc_full), 25))
  expect_true(all(is.finite(cov_out)))

  vc_ssn <- vcov(fit, type = "ssn")
  ssn_out <- conditional(fit, "CapeHornSmall", output = "ssn", samples = 25, simulate_covparams = TRUE)
  expect_equal(rownames(ssn_out), rownames(vc_ssn))
  expect_equal(dim(ssn_out), c(nrow(vc_ssn), 25))

  vc_tailup <- vcov(fit, type = "tailup")
  tailup_out <- conditional(fit, "CapeHornSmall", output = "tailup", samples = 25, simulate_covparams = TRUE)
  expect_equal(rownames(tailup_out), rownames(vc_tailup))
  expect_equal(dim(tailup_out), c(nrow(vc_tailup), 25))

  vc_randcov <- vcov(fit, type = "randcov")
  randcov_out <- conditional(fit, "CapeHornSmall", output = "randcov", samples = 25, simulate_covparams = TRUE)
  expect_equal(rownames(randcov_out), rownames(vc_randcov))
  expect_equal(nrow(randcov_out), nrow(vc_randcov))
  expect_equal(ncol(randcov_out), 25)

  # "cov" stacks ssn rows then randcov rows, matching vcov(type = "cov")
  expect_equal(rownames(vc_full), c(rownames(vc_ssn), rownames(vc_randcov)))
  # "tailup" is a subset of "ssn" (taildown_type/euclid_type are "none" here)
  expect_equal(rownames(vc_tailup), c("tailup_de", "tailup_range"))
})

test_that("known covariance parameters are never perturbed by simulate_covparams = TRUE", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$preds$CapeHornSmall <- mf04p$preds$CapeHorn[1:8, ]
  ssn_create_distmat(mf04p, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE)

  fit_known_range <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_initial = tailup_initial("exponential", range = 10000, known = "range"),
    additive = "afvArea"
  )

  # a genuinely fixed value is structurally excluded from the free-parameter
  # vector itself (not merely "happens to stay constant") -- confirmed
  # directly via the same context conditional() builds internally
  context <- SSN2:::get_satterthwaite_context(fit_known_range, "numeric")
  expect_false("tailup_range" %in% context$cov_names_free_orig)

  cov_out <- conditional(fit_known_range, "CapeHornSmall", output = "cov", samples = 20, simulate_covparams = TRUE)
  expect_false("tailup_range" %in% rownames(cov_out))
})

test_that("simulate_theta_draw_ssn()'s clamp-to-boundary fallback is always valid", {
  # direct unit test, mirroring spmodel's own simulate_theta_draw() test: a
  # huge covariance forces essentially every reject-and-redraw attempt to
  # fail, exercising the clamp-to-boundary fallback
  theta_hat_free <- c(tailup_de = 1, tailup_range = 500)
  huge_vcov <- diag(1e6, 2)
  dimnames(huge_vcov) <- list(names(theta_hat_free), names(theta_hat_free))
  huge_lowchol <- t(chol(huge_vcov))

  set.seed(2)
  draws <- lapply(1:30, function(i) {
    SSN2:::simulate_theta_draw_ssn(theta_hat_free, huge_lowchol, names(theta_hat_free), "exponential", max_attempts = 5)
  })
  de_vals <- vapply(draws, function(x) x$theta[["tailup_de"]], numeric(1))
  range_vals <- vapply(draws, function(x) x$theta[["tailup_range"]], numeric(1))
  exhausted_flags <- vapply(draws, function(x) x$exhausted, logical(1))

  expect_true(all(is.finite(c(de_vals, range_vals))))
  expect_true(all(de_vals >= 0) && all(range_vals >= 0))
  # a clamped range of exactly 0 must also force de to 0 (avoids
  # division-by-zero in exp(-d/range)-style formulas downstream)
  expect_true(all(de_vals[range_vals == 0] == 0))
  # a 1e6-variance draw around small fitted values should exhaust essentially
  # every attempt
  expect_true(any(exhausted_flags))
})

test_that("exhausted covariance-parameter draws are reported via an aggregated warning()", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$preds$CapeHornSmall <- mf04p$preds$CapeHorn[1:8, ]
  ssn_create_distmat(mf04p, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE)

  fit <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    additive = "afvArea"
  )

  sw <- SSN2:::get_satterthwaite_cached(fit, method = "numeric")
  expect_false(is.null(sw$vcov_theta))
  # Force rejection while retaining a positive nugget after clamping.
  sw$context$cov_val_free[c("tailup_de", "nugget")] <- c(-1e6, 1)
  sw$vcov_theta[] <- diag(1e-6, nrow(sw$vcov_theta))
  context <- SSN2:::get_conditional_context(fit, "CapeHornSmall")

  set.seed(2)
  expect_warning(
    SSN2:::draw_conditional_covparams(fit, context, sw, samples = 8, output = "newdata", max_attempts = 1),
    "exhausted 1"
  )
})

test_that("a doubly-degenerate clamped covariance-parameter draw errors clearly instead of crashing in chol()", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$preds$CapeHornSmall <- mf04p$preds$CapeHorn[1:8, ]
  ssn_create_distmat(mf04p, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE)

  fit <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    additive = "afvArea"
  )

  sw <- SSN2:::get_satterthwaite_cached(fit, method = "numeric")
  sw$context$cov_val_free[c("tailup_de", "nugget")] <- -1e6
  sw$vcov_theta[] <- diag(1e-6, nrow(sw$vcov_theta))
  context <- SSN2:::get_conditional_context(fit, "CapeHornSmall")

  expect_error(
    SSN2:::draw_conditional_covparams(fit, context, sw, samples = 2, output = "newdata", max_attempts = 1),
    "degenerate"
  )
})

test_that("simulate_covparams is ssn_lm()-only; passing it to conditional.ssn_glm() is silently ignored (matches spmodel parity)", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$preds$CapeHornSmall <- mf04p$preds$CapeHorn[1:8, ]
  ssn_create_distmat(mf04p, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE)

  fit_glm <- ssn_glm(
    formula = C16 ~ ELEV_DEM,
    ssn.object = mf04p,
    family = "poisson",
    tailup_type = "exponential",
    additive = "afvArea"
  )

  set.seed(2)
  draws1 <- conditional(fit_glm, "CapeHornSmall", type = "link", samples = 5)
  set.seed(2)
  draws2 <- conditional(fit_glm, "CapeHornSmall", type = "link", samples = 5, simulate_covparams = TRUE)
  expect_identical(draws1, draws2)
})

