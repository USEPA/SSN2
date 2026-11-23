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
  expect_error(conditional(fit), "No missing data to predict")
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
  V <- solve(precision - D)
  M <- V %*% precision %*% context$Xmat
  W <- t(context$C0) %*% precision
  H <- context$x0 + W %*% (M - context$Xmat)
  var_adj <- W %*% V %*% t(W)
  beta_unc_cov <- H %*% tcrossprod(context$cov_betahat, H)
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
  expect_error(conditional(fit), "No missing data to predict")
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

test_that("local Gaussian fits support exact conditional uncertainty", {
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
  missing_exact <- conditional(fit, ".missing", samples = 10)
  local({
    local_mocked_bindings(get_block_pred_backend = function(...) "dense")
    set.seed(2)
    expect_equal(conditional(fit, ".missing", samples = 10), missing_exact, tolerance = 1e-9)
  })
  expect_equal(rownames(missing_exact), as.character(fit$missing_index))
  set.seed(2)
  expect_equal(conditional(fit, "small", output = "beta", samples = 40), beta, ignore_attr = TRUE)
  expect_error(conditional(fit, "small", simulate_covparams = TRUE),
    "not supported for models fitted with 'local'")

  dirs <- file.path(network$path, "distance", c("obs", "small"))
  files <- unlist(lapply(dirs, list.files, pattern = "\\.RData$", full.names = TRUE))
  hidden <- paste0(files, ".hidden")
  expect_true(all(file.rename(files, hidden)))
  on.exit(file.rename(hidden, files), add = TRUE)
  set.seed(2)
  expect_equal(conditional(fit, "small", samples = 40), exact, tolerance = 1e-9)
  set.seed(2)
  expect_equal(conditional(fit, ".missing", samples = 10), missing_exact, tolerance = 1e-9)
})

test_that("local GLM fits support exact draws with fixed-effect uncertainty", {
  network <- conditional_local_fit_fixture()
  network$obs$C16[c(3, 10)] <- NA_real_
  fit <- do.call(ssn_glm, c(list(formula = C16 ~ ELEV_DEM + offset(off),
    family = "poisson"), conditional_local_fit_args(network)))
  expect_false(is.null(fit$local_index))
  for (type in c("link", "response", "new")) {
    set.seed(2)
    exact <- conditional(fit, "small", type = type, samples = 40)
    expect_equal(dim(exact), c(6, 40))
    expect_true(all(is.finite(exact)))
    if (type == "new") expect_true(all(exact == round(exact) & exact >= 0))
  }
  context <- SSN2:::get_conditional_context_glm(fit, "small")
  set.seed(2)
  beta <- SSN2:::draw_conditional_glm_joint(SSN2:::get_conditional_glm_joint(context$object), 40)$beta
  set.seed(2)
  expect_equal(conditional(fit, "small", output = "beta", samples = 40), beta, ignore_attr = TRUE)
})

