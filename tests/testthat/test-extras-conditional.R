skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# conditional
test_that("conditional() Gaussian exact conditional simulation works", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$obs$netID <- as.factor(mf04p$obs$netID)
  mf04p$preds$CapeHorn$netID <- as.factor(mf04p$preds$CapeHorn$netID)

  fit <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    additive = "afvArea",
    random = ~ as.factor(netID)
  )

  # basic dispatch / dimensions
  set.seed(2)
  draws <- conditional(fit, "CapeHorn", samples = 2000)
  n_pred <- NROW(mf04p$preds$CapeHorn)
  expect_equal(dim(draws), c(n_pred, 2000))
  expect_true(is.matrix(draws))

  # mean/marginal-variance recovery vs predict(), within Monte Carlo error
  pr <- predict(fit, "CapeHorn", se.fit = TRUE)
  expect_equal(unname(rowMeans(draws)), unname(pr$fit), tolerance = 0.05)
  emp_sd <- apply(draws, 1, sd)
  expect_true(mean(emp_sd / pr$se.fit) > 0.9 && mean(emp_sd / pr$se.fit) < 1.1)

  # exact (non-Monte-Carlo) decomposition: diag(Sigma_cond) + beta-uncertainty
  # term reproduces predict()'s se.fit^2 to machine precision
  context <- SSN2:::get_conditional_context(fit, "CapeHorn")
  cond <- SSN2:::get_conditional_cov(context)
  H <- context$x0 - crossprod(context$SqrtSigInv_C0, context$SqrtSigInv_X)
  beta_unc <- diag(H %*% tcrossprod(context$cov_betahat, H))
  total_analytic <- diag(cond$Sigma_cond) + beta_unc
  expect_equal(unname(total_analytic), unname(pr$se.fit^2), tolerance = 1e-8)

  # independent solve()-based reference for conditional mean/covariance,
  # conditioning on the full observed dataset (not a subset)
  idx_pred <- 1:6
  Sigma11 <- covmatrix(fit)
  C0_full <- covmatrix(fit, "CapeHorn", cov_type = "obs.pred")
  C0 <- C0_full[, idx_pred]
  Sigma22_full <- covmatrix(fit, "CapeHorn", cov_type = "pred.pred")
  Sigma22 <- Sigma22_full[idx_pred, idx_pred]
  Sigma11_inv <- solve(Sigma11)
  reference_cond_cov <- Sigma22 - t(C0) %*% Sigma11_inv %*% C0
  reference_mu <- as.numeric(context$x0[idx_pred, , drop = FALSE] %*% context$betahat +
    t(C0) %*% Sigma11_inv %*% (context$y - context$Xmat %*% context$betahat))
  expect_equal(unname(as.matrix(reference_cond_cov)), unname(cond$Sigma_cond[idx_pred, idx_pred]), tolerance = 1e-8)
  expect_equal(reference_mu, unname(pr$fit[idx_pred]), tolerance = 1e-8)

  # joint off-diagonal covariance (not just marginal variances) matches the
  # analytic total (spatial-conditional + composition-sampled beta
  # uncertainty) within Monte Carlo error
  set.seed(2)
  draws_small <- conditional(fit, "CapeHorn", samples = 20000)
  emp_cov <- cov(t(draws_small[idx_pred, ]))
  beta_unc_cov <- H %*% tcrossprod(context$cov_betahat, H)
  analytic_total_cov <- (cond$Sigma_cond + beta_unc_cov)[idx_pred, idx_pred]
  rel_frob <- norm(emp_cov - analytic_total_cov, "F") / norm(analytic_total_cov, "F")
  expect_true(rel_frob < 0.05)

  # output = "beta": simulated fixed effects recover betahat/vcov(fit)
  set.seed(2)
  betas <- conditional(fit, "CapeHorn", output = "beta", samples = 8000)
  expect_equal(dim(betas), c(fit$p, 8000))
  expect_equal(unname(rowMeans(betas)), unname(coef(fit)), tolerance = 0.05)
  expect_equal(unname(cov(t(betas))), unname(as.matrix(vcov(fit))), tolerance = 0.05)

  # RNG reproducibility
  set.seed(2)
  d1 <- conditional(fit, "CapeHorn", samples = 10)
  set.seed(2)
  d2 <- conditional(fit, "CapeHorn", samples = 10)
  expect_identical(d1, d2)
  d3 <- conditional(fit, "CapeHorn", samples = 10)
  expect_false(identical(d1, d3))

  # rejection: no newdata, invalid newdata name, invalid output
  expect_error(conditional(fit), "requires newdata")
  expect_error(conditional(fit, "doesnotexist"), "not a valid prediction set name")
  expect_error(conditional(fit, "CapeHorn", output = "bogus"), "output must be")

  fake_other <- structure(list(), class = "not_a_supported_class")
  expect_error(conditional(fake_other, "CapeHorn"), "no applicable method")

  preds_new_level <- mf04p$preds$CapeHorn
  levels(preds_new_level$netID) <- c(levels(preds_new_level$netID), "999")
  preds_new_level$netID[1:3] <- factor("999", levels = levels(preds_new_level$netID))
  fit_new_level <- fit
  fit_new_level$ssn.object$preds$CapeHorn <- preds_new_level
  randcov_param <- coef(fit_new_level, type = "randcov")
  Sigma22_new <- covmatrix(fit_new_level, "CapeHorn", cov_type = "pred.pred")
  expect_equal(as.numeric(Sigma22_new[1, 2] - Sigma22_new[1, 4]), as.numeric(randcov_param), tolerance = 1e-6)
  draws_new_level <- conditional(fit_new_level, "CapeHorn", samples = 20)
  expect_equal(dim(draws_new_level), c(n_pred, 20))
})

test_that("conditional() output supports multiple elements, \"object\", and \"all\"", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  ssn_create_distmat(mf04p, predpts = "CapeHorn", overwrite = TRUE, among_predpts = TRUE)

  fit <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    additive = "afvArea",
    ddf = "satterthwaite"
  )

  # output = "object" returns the observed data, replicated across samples
  y_obs <- model.response(model.frame(fit))
  obj <- conditional(fit, "CapeHorn", output = "object", samples = 3)
  expect_equal(dim(obj), c(fit$n, 3))
  expect_true(all(apply(obj, 2, function(col) identical(unname(col), unname(y_obs)))))

  # a vector output returns a named list; each element matches the
  # single-output call bit-for-bit (identical RNG stream)
  set.seed(2)
  d_newdata <- conditional(fit, "CapeHorn", output = "newdata", samples = 5)
  set.seed(2)
  d_beta <- conditional(fit, "CapeHorn", output = "beta", samples = 5)
  set.seed(2)
  d_multi <- conditional(fit, "CapeHorn", output = c("newdata", "beta"), samples = 5)
  expect_type(d_multi, "list")
  expect_identical(names(d_multi), c("newdata", "beta"))
  expect_identical(d_multi$newdata, d_newdata)
  expect_identical(d_multi$beta, d_beta)

  # requested order is respected in the returned list's names
  d_reordered <- conditional(fit, "CapeHorn", output = c("object", "beta"), samples = 2)
  expect_identical(names(d_reordered), c("object", "beta"))

  # "all" expands to newdata/beta/object (not the covparam outputs, even
  # when simulate_covparams = TRUE), matching spmodel's own convention
  d_all <- conditional(fit, "CapeHorn", output = "all", samples = 2)
  expect_identical(names(d_all), c("newdata", "beta", "object"))
  d_all_cp <- conditional(fit, "CapeHorn", output = "all", samples = 2, simulate_covparams = TRUE)
  expect_identical(names(d_all_cp), c("newdata", "beta", "object"))

  # simulate_covparams = TRUE combined with a covparam output and "object"
  d_cp_multi <- conditional(fit, "CapeHorn", output = c("cov", "object"), samples = 4, simulate_covparams = TRUE)
  expect_identical(names(d_cp_multi), c("cov", "object"))
  expect_equal(dim(d_cp_multi$object), c(fit$n, 4))

  # a vector containing an invalid element still errors
  expect_error(conditional(fit, "CapeHorn", output = c("newdata", "bogus")), "output must be")
})

# conditional glm
test_that("conditional() GLM exact conditional simulation works", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)

  fit <- ssn_glm(
    formula = C16 ~ ELEV_DEM,
    ssn.object = mf04p,
    family = "poisson",
    tailup_type = "exponential",
    additive = "afvArea"
  )
  n_pred <- NROW(mf04p$preds$CapeHorn)

  # basic dispatch / dimensions, all three types
  set.seed(2)
  d_link <- conditional(fit, "CapeHorn", type = "link", samples = 500)
  set.seed(2)
  d_resp <- conditional(fit, "CapeHorn", type = "response", samples = 500)
  set.seed(2)
  d_new <- conditional(fit, "CapeHorn", type = "new", samples = 500)
  expect_equal(dim(d_link), c(n_pred, 500))
  expect_equal(dim(d_resp), c(n_pred, 500))
  expect_equal(dim(d_new), c(n_pred, 500))

  # type = "response" is the exact deterministic invlink() of type = "link"
  # (no added family-sampling noise on top of the latent draw itself)
  expect_equal(exp(d_link), d_resp)

  # type = "new" draws are integer-valued (Poisson)
  expect_true(all(d_new == round(d_new)))

  # variance ordering: response-scale variance from link-only draws is
  # strictly less than "new" draws' variance (which adds Poisson sampling
  # noise on top)
  sd_resp <- apply(d_resp, 1, sd)
  sd_new <- apply(d_new, 1, sd)
  expect_true(all(sd_new >= sd_resp - 1e-8))

  pr <- predict(fit, "CapeHorn", type = "link", se.fit = TRUE, var_correct = TRUE)
  context <- SSN2:::get_conditional_context_glm(fit, "CapeHorn")
  cond <- SSN2:::get_conditional_cov(context)
  H <- context$x0 - crossprod(context$SqrtSigInv_C0, context$SqrtSigInv_X)
  precision <- solve(covmatrix(fit))
  B <- context$cov_betahat_uncorrected %*% t(context$Xmat) %*% precision
  Q <- precision - precision %*% context$Xmat %*% B
  D <- -diag(exp(context$w))
  W <- H %*% B + t(context$C0) %*% precision
  var_adj <- W %*% solve(Q - D, t(W))
  beta_unc <- diag(H %*% tcrossprod(context$cov_betahat_uncorrected, H))
  total_analytic <- diag(cond$Sigma_cond + var_adj) + beta_unc
  expect_equal(unname(total_analytic), unname(pr$se.fit^2), tolerance = 1e-8)

  # mean recovery vs predict() (link scale)
  pr_link <- predict(fit, "CapeHorn", type = "link")
  expect_equal(unname(rowMeans(d_link)), unname(pr_link), tolerance = 0.1)

  # joint (not just marginal) link-scale covariance matches the analytic
  # total (Sigma_cond + var_adj + beta uncertainty) within Monte Carlo error
  idx <- 1:6
  beta_unc_cov <- H %*% tcrossprod(context$cov_betahat_uncorrected, H)
  analytic_total_cov <- (cond$Sigma_cond + var_adj + beta_unc_cov)[idx, idx]
  set.seed(2)
  d_small <- conditional(fit, "CapeHorn", type = "link", samples = 15000)
  emp_cov <- cov(t(d_small[idx, ]))
  rel_frob <- norm(emp_cov - analytic_total_cov, "F") / norm(analytic_total_cov, "F")
  expect_true(rel_frob < 0.05)

  # Composition sampling uses corrected fixed-effect uncertainty.
  set.seed(2)
  betas <- conditional(fit, "CapeHorn", output = "beta", samples = 6000)
  expect_equal(dim(betas), c(fit$p, 6000))
  expect_equal(unname(rowMeans(betas)), unname(coef(fit)), tolerance = 0.05)
  expect_equal(unname(diag(cov(t(betas)))), unname(diag(as.matrix(vcov(fit)))), tolerance = 0.1)

  # offset applied exactly once, on the link scale
  mf04p_off <- mf04p
  mf04p_off$obs$off_val <- 2
  mf04p_off$preds$CapeHorn$off_val <- 2
  fit_off <- ssn_glm(
    formula = C16 ~ ELEV_DEM + offset(off_val),
    ssn.object = mf04p_off,
    family = "poisson",
    tailup_type = "exponential",
    additive = "afvArea"
  )
  fit_off0 <- fit_off
  fit_off0$ssn.object$preds$CapeHorn$off_val <- 0
  set.seed(2)
  d_off2 <- conditional(fit_off, "CapeHorn", type = "link", samples = 1)
  set.seed(2)
  d_off0 <- conditional(fit_off0, "CapeHorn", type = "link", samples = 1)
  expect_equal(as.vector(d_off2) - as.vector(d_off0), rep(2, n_pred))

  # RNG reproducibility
  set.seed(2)
  a <- conditional(fit, "CapeHorn", type = "new", samples = 15)
  set.seed(2)
  b <- conditional(fit, "CapeHorn", type = "new", samples = 15)
  expect_identical(a, b)

  # rejection sweep
  expect_error(conditional(fit, "CapeHorn", type = "bogus"))
  expect_error(conditional(fit), "requires newdata")
  expect_error(conditional(fit, "CapeHorn", output = "bogus"), "output must be")

  # output = "object" returns fitted(type = "link"), replicated across
  # samples (not the observed response y, which is latent for GLM); "all"
  # and multi-element output work the same as for ssn_lm()
  w_fitted <- fitted(fit, type = "link")
  obj <- conditional(fit, "CapeHorn", output = "object", samples = 3)
  expect_equal(dim(obj), c(fit$n, 3))
  expect_true(all(apply(obj, 2, function(col) identical(unname(col), unname(w_fitted)))))

  d_all <- conditional(fit, "CapeHorn", output = "all", samples = 2)
  expect_identical(names(d_all), c("newdata", "beta", "object"))

  d_multi <- conditional(fit, "CapeHorn", output = c("beta", "newdata"), samples = 2)
  expect_identical(names(d_multi), c("beta", "newdata"))
})

test_that("conditional() for ssn_glm() binomial family respects newdata_size", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)

  set.seed(2)
  mf04p$obs$trials <- 10
  mf04p$obs$succ <- rbinom(NROW(mf04p$obs), mf04p$obs$trials, plogis(0.3 * scale(mf04p$obs$ELEV_DEM)))
  fit_bin <- ssn_glm(
    formula = cbind(succ, trials - succ) ~ ELEV_DEM,
    ssn.object = mf04p,
    family = "binomial",
    tailup_type = "exponential",
    additive = "afvArea"
  )
  n_pred <- NROW(mf04p$preds$CapeHorn)
  newdata_size <- rep(20, n_pred)

  set.seed(2)
  d_new <- conditional(fit_bin, "CapeHorn", type = "new", samples = 4000, newdata_size = newdata_size)
  set.seed(2)
  d_resp <- conditional(fit_bin, "CapeHorn", type = "response", samples = 4000, newdata_size = newdata_size)

  # "new" draws respect the trial-count bounds
  expect_true(all(d_new >= 0 & d_new <= 20))

  # law of total expectation: E[new] ~= E[response] (both averaged over the
  # same latent uncertainty), within Monte Carlo error
  mean_new <- rowMeans(d_new)
  mean_resp <- rowMeans(d_resp)
  expect_equal(unname(mean_new), unname(mean_resp), tolerance = 0.1)
})

test_that("conditional() for ssn_glm() has no method for unimplemented families is not applicable; family dispatch reuses ssn_r*() formulas", {
  # Gamma: law of total expectation and Var(new) >= Var(response) hold with
  # well-conditioned synthetic data (finite, moderate dispersion)
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)

  set.seed(2)
  mf04p$obs$posval <- rgamma(NROW(mf04p$obs), shape = 5, rate = 5 / 10)
  fit_gam <- ssn_glm(
    formula = posval ~ ELEV_DEM,
    ssn.object = mf04p,
    family = "Gamma",
    tailup_type = "exponential",
    additive = "afvArea"
  )

  set.seed(2)
  d_new <- conditional(fit_gam, "CapeHorn", type = "new", samples = 6000)
  set.seed(2)
  d_resp <- conditional(fit_gam, "CapeHorn", type = "response", samples = 6000)
  mean_new <- rowMeans(d_new)
  mean_resp <- rowMeans(d_resp)
  expect_equal(unname(mean_new), unname(mean_resp), tolerance = 0.1)
  expect_true(all(d_new > 0))
})

# conditional local
test_that("conditional() local = TRUE/list(...) works for ssn_lm()", {
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

  set.seed(2)
  d_exact <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 3000)
  set.seed(2)
  d_local_all <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 3000, local = list(approximation = "vecchia", method ="all"))
  expect_equal(dim(d_local_all), dim(d_exact))
  expect_equal(d_local_all, d_exact, tolerance = 1e-10)

  # output = "beta" never touches the spatial engine -- identical under a
  # shared RNG stream regardless of local
  set.seed(2)
  b_exact <- conditional(fit, "CapeHornSmall", output = "beta", samples = 50)
  set.seed(2)
  b_local <- conditional(fit, "CapeHornSmall", output = "beta", samples = 50, local = TRUE)
  expect_identical(b_exact, b_local)

  # "object" and multi-element output also work under local
  obj_local <- conditional(fit, "CapeHornSmall", output = "object", samples = 3, local = TRUE)
  y_obs <- model.response(model.frame(fit))
  expect_equal(dim(obj_local), c(fit$n, 3))
  expect_true(all(apply(obj_local, 2, function(col) identical(unname(col), unname(y_obs)))))

  d_local_multi <- conditional(fit, "CapeHornSmall", output = c("newdata", "object"), samples = 4, local = list(approximation = "vecchia", method ="all"))
  expect_identical(names(d_local_multi), c("newdata", "object"))

  # truncated (bounded-neighbor) approximation: relative covariance error
  # shrinks as size grows toward the full pool, and is (near-)exact once
  # size covers the whole obs+newdata pool
  set.seed(2)
  d_ref <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 8000)
  cov_ref <- cov(t(d_ref))
  rel_err <- function(sz) {
    set.seed(2)
    d <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 8000, local = list(approximation = "vecchia", method ="covariance", size = sz))
    max(abs(cov(t(d)) - cov_ref)) / max(abs(cov_ref))
  }
  errs <- vapply(c(3, 10, 45, 52), rel_err, numeric(1))
  expect_true(all(diff(errs) < 1e-6)) # non-increasing (size=52 covers the whole pool: n_obs=45 + n_new-1=7)
  expect_true(errs[[4]] < 1e-8)

  # Different factorizations need not reproduce identical individual draws.
  for (ordering in c("pid", "none", "random", "maxmin", "middleout", "outsidein", "coordinate", "grts")) {
    set.seed(2)
    d <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 200, local = list(approximation = "vecchia", method ="covariance", size = 5, ordering = ordering))
    expect_true(all(is.finite(d)))
    expect_equal(dim(d), c(8, 200))
  }

  # RNG reproducibility under the local (truncated) engine specifically
  vecchia_local <- list(approximation = "vecchia", method = "covariance", size = 5)
  set.seed(2)
  d1 <- conditional(fit, "CapeHornSmall", samples = 10, local = vecchia_local)
  set.seed(2)
  d2 <- conditional(fit, "CapeHornSmall", samples = 10, local = vecchia_local)
  expect_identical(d1, d2)
  d3 <- conditional(fit, "CapeHornSmall", samples = 10, local = vecchia_local)
  expect_false(identical(d1, d3))

  # a partial vecchia list initializes defaults for the other list elements
  expect_mapequal(
    SSN2:::get_conditional_local(list(approximation = "vecchia", size = 10), fit, "CapeHornSmall"),
    list(approximation = "vecchia", size = 10L, method = "covariance", ordering = "pid")
  )
  expect_identical(
    SSN2:::get_conditional_local(list(approximation = "vecchia", method = "covariance"), fit, "CapeHornSmall"),
    SSN2:::get_conditional_local(list(approximation = "vecchia"), fit, "CapeHornSmall")
  )

  # rejection sweep: local argument vocabulary, and local + simulate_covparams
  expect_error(conditional(fit, "CapeHornSmall", local = list(approximation = "vecchia", method = "bogus")), "method must be")
  expect_error(conditional(fit, "CapeHornSmall", local = list(approximation = "vecchia", size = -1)), "size must be")
  expect_error(conditional(fit, "CapeHornSmall", local = list(approximation = "vecchia", method = "covariance", size = 5, ordering = "bogus")), "ordering must be")
  expect_error(conditional(fit, "CapeHornSmall", local = list(method = "covariance")), "local now defaults to a low-rank approximation")

  # local + simulate_covparams = TRUE: matches spmodel -- a message and
  # Local simulation disables covariance-parameter sampling without an error.
  expect_message(
    d_local_covparams <- conditional(fit, "CapeHornSmall", local = TRUE, simulate_covparams = TRUE, samples = 5),
    "setting simulate_covparams = FALSE"
  )
  expect_true(all(is.finite(d_local_covparams)))
  # once forced to FALSE, a covparam-only output still correctly errors
  expect_error(
    suppressMessages(conditional(fit, "CapeHornSmall", local = TRUE, output = "cov", simulate_covparams = TRUE)),
    "output can only be"
  )

  # duplicate prediction location (a near-singular neighbor set) does not
  # crash and produces finite draws
  preds_dup <- mf04p$preds$CapeHornSmall
  preds_dup[2, ] <- preds_dup[1, ]
  mf04p_dup <- mf04p
  mf04p_dup$preds$CapeHornSmall <- preds_dup
  fit_dup <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p_dup,
    tailup_type = "exponential",
    additive = "afvArea"
  )
  set.seed(2)
  d_dup <- conditional(fit_dup, "CapeHornSmall", output = "newdata", samples = 50, local = list(approximation = "vecchia", method ="covariance", size = 6))
  expect_equal(dim(d_dup), c(8, 50))
  expect_true(all(is.finite(d_dup)))
})

test_that("local is automatically enabled when the observed or prediction sample size exceeds 5,000, matching spmodel", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$preds$CapeHornSmall <- mf04p$preds$CapeHorn[1:8, ]
  ssn_create_distmat(mf04p, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE)

  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea"
  )

  # faking n avoids fitting an actually-large model just for this check,
  # matching the analogous simulate_covparams performance-warning test
  fit_fake_big <- fit
  fit_fake_big$n <- 6000
  expect_message(
    conditioning <- SSN2:::get_conditional_local(NULL, fit_fake_big, "CapeHornSmall"),
    "Because the observed data size or the number of prediction locations exceeds 5,000"
  )
  expect_identical(conditioning$approximation, "low-rank")

  # small sample sizes stay exact when local is omitted (no auto-switch, no message)
  expect_no_message(
    exact_conditioning <- SSN2:::get_conditional_local(NULL, fit, "CapeHornSmall")
  )
  expect_identical(exact_conditioning$method, "exact")

  # explicit local = FALSE always stays exact, even for a large sample size --
  # only an *omitted* local auto-switches
  expect_no_message(
    forced_exact <- SSN2:::get_conditional_local(FALSE, fit_fake_big, "CapeHornSmall")
  )
  expect_identical(forced_exact$method, "exact")
})

test_that("local low-rank conditional simulation (approximation = 'low-rank', the new default) approximates the exact distribution", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$preds$CapeHornSmall <- mf04p$preds$CapeHorn[1:8, ]
  ssn_create_distmat(mf04p, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE)

  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea",
    tailup_initial = tailup_initial("exponential", de = 0.4, range = 300, known = "given"),
    nugget_initial = nugget_initial("nugget", nugget = 0.1, known = "given")
  )
  set.seed(2)
  d_exact <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 8000)
  cov_exact <- cov(t(d_exact))

  # no subsetting at all (method_base/method_new = "all") should match exact
  set.seed(2)
  d_all <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 8000,
    local = list(approximation = "low-rank", method_base = "all", method_new = "all"))
  expect_true(all(is.finite(d_all)))
  rel_all <- max(abs(cov(t(d_all)) - cov_exact)) / max(abs(cov_exact))
  expect_true(rel_all < 0.15)

  # base-only subsetting (one undivided newdata block) still approximates well
  set.seed(2)
  d_base <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 8000,
    local = list(approximation = "low-rank", method_base = "base", size_base = 20, method_new = "all"))
  expect_true(all(is.finite(d_base)))
  rel_base <- max(abs(cov(t(d_base)) - cov_exact)) / max(abs(cov_exact))
  expect_true(rel_base < 0.4)

  # genuine multi-block newdata splitting: blocks are conditionally
  # independent given the base by construction (a documented, expected
  # approximation trade-off -- see spmodel's own "blocks are assumed
  # conditionally independent" caveat), so this only checks that draws stay
  # finite with accurate means, not tight covariance agreement
  set.seed(2)
  d_multiblock <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 8000,
    local = list(approximation = "low-rank", method_base = "all", method_new = "base", size_new = 3))
  expect_true(all(is.finite(d_multiblock)))
  rel_mean <- max(abs(rowMeans(d_multiblock) - rowMeans(d_exact))) / max(abs(rowMeans(d_exact)))
  expect_true(rel_mean < 0.05)

  # local = TRUE now resolves to low-rank, matching spmodel's own default
  expect_identical(SSN2:::get_conditional_local(TRUE, fit, "CapeHornSmall")$approximation, "low-rank")

  # output = "beta"/"object" still work with low-rank
  res <- conditional(fit, "CapeHornSmall", output = c("newdata", "beta", "object"), samples = 50,
    local = list(approximation = "low-rank", method_base = "base", size_base = 20, method_new = "all"))
  expect_true(all(is.finite(res$newdata)) && all(is.finite(res$beta)) && all(is.finite(res$object)))

  # parallel = TRUE matches parallel = FALSE for a genuine multi-block config
  lowrank_multiblock <- list(approximation = "low-rank", method_base = "all",
    method_new = "base", size_new = 3, reorder_new = "none", kmeans_new = FALSE)
  set.seed(2)
  d_serial <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 8000, local = lowrank_multiblock)
  set.seed(2)
  d_par <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 8000,
    local = c(lowrank_multiblock, list(parallel = TRUE, ncores = 2)))
  expect_true(all(is.finite(d_par)))
  rel_mean_par <- max(abs(rowMeans(d_par) - rowMeans(d_serial))) / max(abs(rowMeans(d_serial)))
  expect_true(rel_mean_par < 0.05)

  # rejection sweep for low-rank's own settings
  expect_error(
    conditional(fit, "CapeHornSmall", local = list(approximation = "low-rank", method_new = "bogus")),
    "method_new must be"
  )
  expect_error(
    conditional(fit, "CapeHornSmall", local = list(approximation = "low-rank", reorder_new = "bogus")),
    "reorder_new must be"
  )
})

test_that("conditional() local = TRUE/list(...) handles random effects and partition factors", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$obs$netID <- as.factor(mf04p$obs$netID)
  mf04p$preds$CapeHorn$netID <- factor(mf04p$preds$CapeHorn$netID, levels = levels(mf04p$obs$netID))
  mf04p$preds$CapeHornSmall <- mf04p$preds$CapeHorn[1:8, ]
  ssn_create_distmat(mf04p, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE)

  fit_rand <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    additive = "afvArea",
    random = ~ as.factor(netID)
  )
  set.seed(2)
  d_exact <- conditional(fit_rand, "CapeHornSmall", output = "newdata", samples = 3000)
  set.seed(2)
  d_local <- conditional(fit_rand, "CapeHornSmall", output = "newdata", samples = 3000, local = list(approximation = "vecchia", method ="all"))
  expect_equal(d_local, d_exact, tolerance = 1e-10)

  fit_pf <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    additive = "afvArea",
    partition_factor = ~ as.factor(netID)
  )
  set.seed(2)
  d_pf_exact <- conditional(fit_pf, "CapeHornSmall", output = "newdata", samples = 3000)
  set.seed(2)
  d_pf_local <- conditional(fit_pf, "CapeHornSmall", output = "newdata", samples = 3000, local = list(approximation = "vecchia", method ="all"))
  expect_equal(d_pf_local, d_pf_exact, tolerance = 1e-10)
})

test_that("local conditional simulation accepts new partition levels", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  ssn <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  ssn$preds$CapeHornSmall <- ssn$preds$CapeHorn[1:8, ]
  ssn_create_distmat(
    ssn, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE
  )
  ssn$obs$new_partition <- factor(rep(c("a", "b"), length.out = NROW(ssn$obs)))
  ssn$preds$CapeHornSmall$new_partition <- factor(
    rep(c("c", "d"), length.out = NROW(ssn$preds$CapeHornSmall))
  )
  fit <- ssn_lm(
    Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "none", taildown_type = "none", euclid_type = "none",
    nugget_initial = nugget_initial("nugget", nugget = 1, known = "given"),
    partition_factor = ~new_partition
  )

  set.seed(2)
  draws <- conditional(
    fit, "CapeHornSmall", samples = 4,
    local = list(approximation = "vecchia", method ="covariance", size = 5)
  )
  expect_equal(dim(draws), c(8, 4))
  expect_true(all(is.finite(draws)))
})

test_that("conditional() local = TRUE/list(...) works for ssn_glm()", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$preds$CapeHornSmall <- mf04p$preds$CapeHorn[1:8, ]
  ssn_create_distmat(mf04p, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE)

  fit <- ssn_glm(
    formula = C16 ~ ELEV_DEM,
    ssn.object = mf04p,
    family = "poisson",
    tailup_type = "exponential",
    additive = "afvArea"
  )

  # All-neighbor conditioning recovers the exact joint draws.
  for (type in c("link", "response", "new")) {
    set.seed(2)
    g_exact <- conditional(fit, "CapeHornSmall", type = type, samples = 3000)
    set.seed(2)
    g_local <- conditional(fit, "CapeHornSmall", type = type, samples = 3000, local = list(approximation = "vecchia", method ="all"))
    expect_equal(g_local, g_exact, tolerance = 1e-10)
  }

  context_exact <- SSN2:::get_conditional_context_glm(fit, "CapeHornSmall")
  context_local <- SSN2:::get_conditional_context_glm(fit, "CapeHornSmall", local = TRUE)
  expect_equal(SSN2:::get_conditional_glm_joint(context_exact),
    SSN2:::get_conditional_glm_joint(context_local), tolerance = 1e-8)

  # bounded/truncated approximation: finite, correctly-typed draws across
  # every response family output, varying newdata_size, at a small size
  set.seed(2)
  g_resp <- conditional(fit, "CapeHornSmall", type = "response", samples = 500, local = list(approximation = "vecchia", method ="covariance", size = 8))
  expect_true(all(g_resp > 0) && all(is.finite(g_resp)))
  set.seed(2)
  g_new <- conditional(fit, "CapeHornSmall", type = "new", samples = 500, local = list(approximation = "vecchia", method ="covariance", size = 8))
  expect_true(all(g_new == round(g_new)))
  sd_new <- apply(g_new, 1, sd)
  sd_resp <- apply(g_resp, 1, sd)
  expect_true(all(sd_new >= sd_resp - 1e-8))

  # output = "beta"
  set.seed(2)
  betas <- conditional(fit, "CapeHornSmall", output = "beta", samples = 50, local = list(approximation = "vecchia", method ="covariance", size = 5))
  expect_equal(dim(betas), c(fit$p, 50))
  expect_true(all(is.finite(betas)))

  # binomial newdata_size threading through the local engine
  set.seed(2)
  mf04p$obs$trials <- 10
  mf04p$obs$succ <- rbinom(NROW(mf04p$obs), mf04p$obs$trials, plogis(0.3 * scale(mf04p$obs$ELEV_DEM)))
  fit_binom <- suppressWarnings(ssn_glm(
    formula = cbind(succ, trials - succ) ~ ELEV_DEM,
    ssn.object = mf04p,
    family = "binomial",
    tailup_type = "exponential",
    additive = "afvArea"
  ))
  set.seed(2)
  d_binom <- conditional(fit_binom, "CapeHornSmall", type = "new", samples = 300, newdata_size = rep(20, 8), local = list(approximation = "vecchia", method ="covariance", size = 6))
  expect_true(all(d_binom >= 0 & d_binom <= 20))
})

test_that("local low-rank GLM simulation shares latent and coefficient draws", {
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
  mf04p$preds$CapeHornSmall <- mf04p$preds$CapeHorn[1:8, ]
  ssn_create_distmat(mf04p, predpts = "CapeHornSmall", overwrite = TRUE, among_predpts = TRUE)

  fit <- ssn_glm(C16 ~ ELEV_DEM, mf04p, family = "poisson",
    euclid_initial = euclid_initial("exponential", 0.15, 12000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.12, known = "given")
  )
  set.seed(2)
  d_exact <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 4000)
  cov_exact <- cov(t(d_exact))

  set.seed(2)
  d_all <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 4000,
    local = list(approximation = "low-rank", method_base = "all", method_new = "all"))
  expect_true(all(is.finite(d_all)))
  rel_all <- max(abs(cov(t(d_all)) - cov_exact)) / max(abs(cov_exact))
  expect_true(rel_all < 0.2)

  set.seed(2)
  d_multiblock <- conditional(fit, "CapeHornSmall", output = "newdata", samples = 4000,
    local = list(approximation = "low-rank", method_base = "all", method_new = "base", size_new = 3))
  expect_true(all(is.finite(d_multiblock)))
  rel_mean <- max(abs(rowMeans(d_multiblock) - rowMeans(d_exact))) / max(abs(rowMeans(d_exact)))
  expect_true(rel_mean < 0.1)
})

# conditional local fit
conditional_local_fit_fixture <- function() {
  copy_lsn_to_temp()
  network <- ssn_import(file.path(tempdir(), "MiddleFork04.ssn"),
    predpts = "CapeHorn", overwrite = TRUE)
  network$preds$small <- network$preds$CapeHorn[c(8, 2, 6, 1, 7, 3), ]
  ssn_create_distmat(network, predpts = "small", overwrite = TRUE, among_predpts = TRUE)
  ssn_create_bigdist(network, predpts = "small", overwrite = TRUE,
    among_predpts = TRUE, verbose = FALSE)
  network$obs <- network$obs[rev(seq_len(NROW(network$obs))), ]
  network$obs$off <- seq(-0.2, 0.2, length.out = NROW(network$obs))
  network$preds$small$off <- seq(0.1, 0.3, length.out = 6)
  network$obs$group <- factor(network$obs$netID)
  network$preds$small$group <- factor(network$preds$small$netID,
    levels = levels(network$obs$group))
  network
}

conditional_local_fit_args <- function(network) {
  list(ssn.object = network, additive = "afvArea",
    tailup_initial = tailup_initial("exponential", 0.2, 10000, known = "given"),
    taildown_initial = taildown_initial("exponential", 0.15, 12000, known = "given"),
    euclid_initial = euclid_initial("exponential", 0.1, 11000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.4, known = "given"),
    local = list(index = rep(1:3, length.out = NROW(network$obs)),
      var_adjust = "none", parallel = FALSE))
}

test_that("local Gaussian fits support exact and Vecchia conditional uncertainty", {
  network <- conditional_local_fit_fixture()
  network$obs$Summer_mn[c(2, 9)] <- NA_real_
  args <- conditional_local_fit_args(network)
  fit <- do.call(ssn_lm, c(list(formula = Summer_mn ~ ELEV_DEM + offset(off),
    partition_factor = ~group), args))
  expect_false(is.null(fit$local_index))
  X <- model.matrix(fit)
  x0 <- model.matrix(~ELEV_DEM, network$preds$small)
  y <- model.response(model.frame(fit)) - model.offset(model.frame(fit))
  S <- covmatrix(fit)
  C <- covmatrix(fit, "small", cov_type = "obs.pred")
  S0 <- covmatrix(fit, "small", cov_type = "pred.pred")
  set.seed(2)
  beta <- as.vector(coef(fit)) + t(chol(vcov(fit))) %*% matrix(rnorm(fit$p * 40), fit$p)
  expected <- x0 %*% beta + t(C) %*% solve(S, as.vector(y) - X %*% beta) +
    t(chol(S0 - t(C) %*% solve(S, C))) %*% matrix(rnorm(6 * 40), 6) +
    network$preds$small$off
  set.seed(2)
  exact <- conditional(fit, "small", samples = 40)
  expect_equal(exact, expected, ignore_attr = TRUE, tolerance = 1e-9)
  set.seed(2)
  sequential <- conditional(fit, "small", samples = 40,
    local = list(approximation = "vecchia", method ="all", ordering = "none"))
  expect_equal(sequential, exact, tolerance = 1e-9)
  set.seed(2)
  missing_exact <- conditional(fit, ".missing", samples = 10)
  local({
    local_mocked_bindings(get_block_pred_backend = function(...) "dense")
    set.seed(2)
    expect_equal(conditional(fit, ".missing", samples = 10), missing_exact, tolerance = 1e-9)
  })
  set.seed(2)
  expect_equal(conditional(fit, ".missing", samples = 10,
    local = list(approximation = "vecchia", method ="all", ordering = "none")), missing_exact, tolerance = 1e-9)
  expect_equal(rownames(missing_exact), as.character(fit$missing_index))

  for (ordering in c("none", "pid", "random")) {
    set.seed(2)
    draws <- conditional(fit, "small", samples = 10,
      local = list(approximation = "vecchia", method ="covariance", size = 5, ordering = ordering))
    expect_equal(dim(draws), c(6, 10))
    expect_true(all(is.finite(draws)))
  }
  for (simulation_local in list(FALSE, TRUE, list(approximation = "vecchia", method = "all"))) {
    set.seed(2)
    expect_equal(conditional(fit, "small", output = "beta", samples = 40,
      local = simulation_local), beta, ignore_attr = TRUE)
    expect_error(conditional(fit, "small", local = simulation_local,
      simulate_covparams = TRUE), "not supported for models fitted with 'local'")
  }

  dirs <- file.path(network$path, "distance", c("obs", "small"))
  files <- unlist(lapply(dirs, list.files, pattern = "\\.RData$", full.names = TRUE))
  hidden <- paste0(files, ".hidden")
  expect_true(all(file.rename(files, hidden)))
  on.exit(file.rename(hidden, files), add = TRUE)
  set.seed(2)
  expect_equal(conditional(fit, "small", samples = 40), exact, tolerance = 1e-9)
  set.seed(2)
  expect_equal(conditional(fit, "small", samples = 40,
    local = list(approximation = "vecchia", method ="all", ordering = "none")), exact, tolerance = 1e-9)
  set.seed(2)
  expect_equal(conditional(fit, ".missing", samples = 10), missing_exact, tolerance = 1e-9)
})

test_that("local GLM fits support exact and Vecchia draws with fixed-effect uncertainty", {
  network <- conditional_local_fit_fixture()
  network$obs$C16[c(3, 10)] <- NA_real_
  fit <- do.call(ssn_glm, c(list(formula = C16 ~ ELEV_DEM + offset(off),
    family = "poisson"), conditional_local_fit_args(network)))
  expect_false(is.null(fit$local_index))
  for (type in c("link", "response", "new")) {
    set.seed(2)
    exact <- conditional(fit, "small", type = type, samples = 40)
    set.seed(2)
    sequential <- conditional(fit, "small", type = type, samples = 40,
      local = list(approximation = "vecchia", method ="all", ordering = "none"))
    expect_equal(sequential, exact, tolerance = 1e-8)
    expect_true(all(is.finite(exact)))
    draws <- conditional(fit, "small", type = type, samples = 10,
      local = list(approximation = "vecchia", method ="covariance", size = 5, ordering = "random"))
    expect_equal(dim(draws), c(6, 10))
    expect_true(all(is.finite(draws)))
    if (type == "new") expect_true(all(draws == round(draws) & draws >= 0))
  }
  context <- SSN2:::get_conditional_context_glm(fit, "small")
  set.seed(2)
  beta <- SSN2:::draw_conditional_glm_joint(context, 40)$beta
  for (simulation_local in list(FALSE, TRUE)) {
    set.seed(2)
    expect_equal(conditional(fit, "small", output = "beta", samples = 40,
      local = simulation_local), beta, ignore_attr = TRUE)
  }
})

