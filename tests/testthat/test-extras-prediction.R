skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# predict offset forwarding
test_that("point prediction applies newdata's own offset exactly once, for every interval type (Bug A)", {
  s <- mf04p
  s$obs$my_offset <- stats::runif(nrow(s$obs), 0, 5)
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM + offset(my_offset), s, tailup_type = "exponential", additive = "afvArea")

  predict_with_offset <- function(offset_val, interval) {
    preds <- fit$ssn.object$preds$CapeHorn
    preds$my_offset <- offset_val
    fit$ssn.object$preds$CapeHorn <- preds
    predict(fit, "CapeHorn", interval = interval, se.fit = TRUE)
  }

  for (interval in c("none", "prediction", "confidence")) {
    p0 <- predict_with_offset(0, interval)
    p100 <- predict_with_offset(100, interval)
    fit0 <- unname(if (interval == "none") p0$fit[1] else p0$fit[1, "fit"])
    fit100 <- unname(if (interval == "none") p100$fit[1] else p100$fit[1, "fit"])
    expect_equal(fit100 - fit0, 100, tolerance = 1e-8, label = paste0("interval=", interval))
    # se.fit must be unaffected by a change in the offset's *value* (offset
    # is a known constant, not a random quantity, so it doesn't contribute
    # to prediction variance)
    expect_equal(p0$se.fit, p100$se.fit, tolerance = 1e-10, label = paste0("interval=", interval, " se.fit"))
  }
})

test_that("block prediction applies newdata's own (block-mean) offset exactly once, for every interval type, and no longer corrupts training residuals with newdata's offset (Bug A)", {
  s <- mf04p
  s$obs$my_offset <- stats::runif(nrow(s$obs), 0, 5)
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM + offset(my_offset), s, tailup_type = "exponential", additive = "afvArea")

  predict_with_offset <- function(offset_val, interval) {
    preds <- fit$ssn.object$preds$CapeHorn
    preds$my_offset <- offset_val
    fit$ssn.object$preds$CapeHorn <- preds
    predict(fit, "CapeHorn", block = TRUE, interval = interval, se.fit = TRUE)
  }

  for (interval in c("none", "prediction", "confidence")) {
    p0 <- predict_with_offset(0, interval)
    p100 <- predict_with_offset(100, interval)
    fit0 <- unname(if (interval == "none") p0$fit[1] else p0$fit[1, "fit"])
    fit100 <- unname(if (interval == "none") p100$fit[1] else p100$fit[1, "fit"])
    expect_equal(fit100 - fit0, 100, tolerance = 1e-8, label = paste0("interval=", interval))
    expect_equal(p0$se.fit, p100$se.fit, tolerance = 1e-10, label = paste0("interval=", interval, " se.fit"))
  }
})

test_that("block prediction for a model without an offset is unaffected (no regression from removing the buggy y <- y - offset step)", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  block_pred <- predict(fit, "CapeHorn", block = TRUE, se.fit = TRUE)
  expect_true(is.finite(block_pred$fit))
  expect_true(is.finite(block_pred$se.fit) && block_pred$se.fit > 0)
})

test_that("ssn_glm() point prediction applies newdata's own offset exactly once, for every interval type", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  s$obs$my_offset <- stats::runif(nrow(s$obs), 0, 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM + offset(my_offset), s,
    family = "poisson", tailup_type = "exponential", additive = "afvArea"
  )

  predict_with_offset <- function(offset_val, interval) {
    preds <- fitg$ssn.object$preds$CapeHorn
    preds$my_offset <- offset_val
    fitg$ssn.object$preds$CapeHorn <- preds
    predict(fitg, "CapeHorn", type = "link", interval = interval, se.fit = TRUE)
  }

  for (interval in c("none", "prediction", "confidence")) {
    p0 <- predict_with_offset(0, interval)
    p100 <- predict_with_offset(100, interval)
    fit0 <- unname(if (interval == "none") p0$fit[1] else p0$fit[1, "fit"])
    fit100 <- unname(if (interval == "none") p100$fit[1] else p100$fit[1, "fit"])
    expect_equal(fit100 - fit0, 100, tolerance = 1e-8, label = paste0("interval=", interval))
    # offset is a known constant, not a random quantity, so it must not
    # affect se.fit
    expect_equal(p0$se.fit, p100$se.fit, tolerance = 1e-10, label = paste0("interval=", interval, " se.fit"))
  }

  # type = "response" applies the inverse link to the now-correctly-shifted
  # linear predictor: for the Poisson log-link, a +100 shift on the link
  # scale multiplies the response-scale prediction by exp(100)
  preds <- fitg$ssn.object$preds$CapeHorn
  preds$my_offset <- 0
  fitg$ssn.object$preds$CapeHorn <- preds
  p0_response <- predict(fitg, "CapeHorn", type = "response", interval = "none")
  preds$my_offset <- 100
  fitg$ssn.object$preds$CapeHorn <- preds
  p100_response <- predict(fitg, "CapeHorn", type = "response", interval = "none")
  expect_equal(unname(p100_response / p0_response), rep(exp(100), length(p0_response)), tolerance = 1e-6)
})

test_that("ssn_glm() point prediction with an offset stays finite through the local/var_correct/dispersion paths", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  s$obs$my_offset <- stats::runif(nrow(s$obs), 0, 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM + offset(my_offset), s,
    family = "poisson", tailup_type = "exponential", additive = "afvArea"
  )
  preds <- fitg$ssn.object$preds$CapeHorn
  preds$my_offset <- stats::runif(nrow(preds), 0, 5)
  fitg$ssn.object$preds$CapeHorn <- preds

  p_local <- predict(fitg, "CapeHorn", type = "link", interval = "prediction", se.fit = TRUE, local = list(method = "covariance", size = 10))
  expect_true(all(is.finite(p_local$fit)))
  expect_true(all(is.finite(p_local$se.fit)))

  p_novc <- predict(fitg, "CapeHorn", type = "link", interval = "prediction", se.fit = TRUE, var_correct = FALSE)
  expect_true(all(is.finite(p_novc$fit)))
  expect_true(all(is.finite(p_novc$se.fit)))

  p_disp <- predict(fitg, "CapeHorn", type = "link", interval = "prediction", se.fit = TRUE, dispersion = 1)
  expect_true(all(is.finite(p_disp$fit)))
  expect_true(all(is.finite(p_disp$se.fit)))
})

test_that("predict.ssn_lm()'s recursive newdata = 'all' dispatch forwards scale/type/terms, not just se.fit/interval/level/local (Bug B)", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")

  direct <- predict(fit, "CapeHorn", se.fit = TRUE, scale = 3, interval = "prediction")
  via_all <- predict(fit, "all", se.fit = TRUE, scale = 3, interval = "prediction")
  expect_equal(via_all$CapeHorn, direct)

  fit_terms <- ssn_lm(Summer_mn ~ ELEV_DEM + SLOPE, mf04p, tailup_type = "exponential", additive = "afvArea")
  direct_terms <- predict(fit_terms, "CapeHorn", type = "terms")
  via_all_terms <- predict(fit_terms, "all", type = "terms")
  expect_equal(via_all_terms$CapeHorn, direct_terms)
})

test_that("predict.ssn_glm()'s recursive newdata = 'all' dispatch forwards type/terms/var_correct, not just se.fit/interval/level/local (Bug B)", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fitg <- ssn_glm(count_response ~ 1, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")

  direct <- predict(fitg, "CapeHorn", type = "response")
  via_all <- predict(fitg, "all", type = "response")
  expect_equal(via_all$CapeHorn, direct)

  direct_nc <- predict(fitg, "CapeHorn", se.fit = TRUE, var_correct = FALSE)
  via_all_nc <- predict(fitg, "all", se.fit = TRUE, var_correct = FALSE)
  expect_equal(via_all_nc$CapeHorn, direct_nc)
})

test_that("dispersion is already correctly applied through predict.ssn_glm()'s recursive dispatch (regression lock, not a new fix)", {
  s <- mf04p
  s$obs$gamma_response <- rpois(nrow(s$obs), lambda = 5) + 1
  fitg <- ssn_glm(gamma_response ~ 1, s, family = "Gamma", tailup_type = "exponential", additive = "afvArea")

  direct <- predict(fitg, "CapeHorn", dispersion = 2, type = "response")
  via_all <- predict(fitg, "all", dispersion = 2, type = "response")
  expect_equal(via_all$CapeHorn, direct)
})

test_that("newdata_size with multiple newdata sets errors instead of silently misapplying (Bug B)", {
  s <- mf04p
  s$obs$count_response <- rbinom(nrow(s$obs), size = 20, prob = 0.3)
  fitg <- ssn_glm(cbind(count_response, 20 - count_response) ~ 1, s,
    family = "binomial", tailup_type = "exponential", additive = "afvArea"
  )

  expect_error(
    predict(fitg, "all", newdata_size = rep(20, nrow(fitg$ssn.object$preds$CapeHorn))),
    "newdata_size cannot be used when predicting for multiple newdata sets"
  )
  # still works for a single, explicit dataset
  expect_no_error(predict(fitg, "CapeHorn", newdata_size = rep(20, nrow(fitg$ssn.object$preds$CapeHorn)), type = "response"))
})

# predict delta se
test_that("get_delta_se() matches a numerical (finite-difference) inverse-link derivative, all families", {
  eps <- 1e-6
  link_fits <- c(-1.5, -0.3, 0, 0.4, 1.2, 2.1)
  se_link <- 0.25

  for (family in c("poisson", "nbinomial", "Gamma", "inverse.gaussian", "binomial", "beta")) {
    size <- if (family == "binomial") 8 else NULL
    g_numeric <- (invlink(link_fits + eps, family, size) - invlink(link_fits - eps, family, size)) / (2 * eps)
    reference_se <- se_link * abs(g_numeric)
    actual_se <- get_delta_se(link_fits, se_link, family, if (family == "binomial") size else 1)
    expect_equal(actual_se, reference_se, tolerance = 1e-4, label = family)
  }
})

test_that("predict.ssn_glm() delta = TRUE matches the numerical-derivative reference, interval = 'none' and 'prediction', with and without offset", {
  eps <- 1e-6
  for (use_offset in c(FALSE, TRUE)) {
    s <- mf04p
    s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
    if (use_offset) {
      s$obs$my_offset <- stats::runif(nrow(s$obs), 0, 2)
      fitg <- ssn_glm(count_response ~ ELEV_DEM + offset(my_offset), s,
        family = "poisson", tailup_type = "exponential", additive = "afvArea"
      )
      fitg$ssn.object$preds$CapeHorn$my_offset <- stats::runif(nrow(fitg$ssn.object$preds$CapeHorn), 0, 2)
    } else {
      fitg <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")
    }

    p_link_none <- predict(fitg, "CapeHorn", type = "link", se.fit = TRUE, interval = "none")
    p_delta_none <- predict(fitg, "CapeHorn", type = "response", se.fit = TRUE, interval = "none", delta = TRUE)
    g_num <- (invlink(p_link_none$fit + eps, fitg$family, NULL) - invlink(p_link_none$fit - eps, fitg$family, NULL)) / (2 * eps)
    reference <- p_link_none$se.fit * abs(g_num)
    expect_equal(unname(p_delta_none$se.fit), unname(reference), tolerance = 1e-4, label = paste("none, offset =", use_offset))

    p_link_pi <- predict(fitg, "CapeHorn", type = "link", se.fit = TRUE, interval = "prediction")
    p_delta_pi <- predict(fitg, "CapeHorn", type = "response", se.fit = TRUE, interval = "prediction", delta = TRUE)
    link_fit_pi <- p_link_pi$fit[, "fit"]
    g_num_pi <- (invlink(link_fit_pi + eps, fitg$family, NULL) - invlink(link_fit_pi - eps, fitg$family, NULL)) / (2 * eps)
    reference_pi <- p_link_pi$se.fit * abs(g_num_pi)
    expect_equal(unname(p_delta_pi$se.fit), unname(reference_pi), tolerance = 1e-4, label = paste("prediction, offset =", use_offset))

    # delta must not change lwr/upr (already correctly invlink()-transformed)
    expect_equal(p_delta_pi$fit[, "lwr"], invlink(p_link_pi$fit[, "lwr"], fitg$family, NULL), tolerance = 1e-8)
    expect_equal(p_delta_pi$fit[, "upr"], invlink(p_link_pi$fit[, "upr"], fitg$family, NULL), tolerance = 1e-8)
  }
})

test_that("delta = FALSE (default) leaves se.fit on the link scale, unaffected by the interval = 'none' reordering", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")

  p_link <- predict(fitg, "CapeHorn", type = "link", se.fit = TRUE, interval = "none")
  p_response_nodelta <- predict(fitg, "CapeHorn", type = "response", se.fit = TRUE, interval = "none")
  expect_equal(p_response_nodelta$se.fit, p_link$se.fit)
  expect_equal(p_response_nodelta$fit, invlink(p_link$fit, fitg$family, NULL))
})

test_that("interval = 'confidence' never applies delta, se.fit always stays link scale", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")

  p_delta <- predict(fitg, "CapeHorn", type = "response", se.fit = TRUE, interval = "confidence", delta = TRUE)
  p_nodelta <- predict(fitg, "CapeHorn", type = "response", se.fit = TRUE, interval = "confidence", delta = FALSE)
  expect_equal(p_delta$se.fit, p_nodelta$se.fit)
})

test_that("delta must be TRUE or FALSE", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")
  expect_error(predict(fitg, "CapeHorn", delta = "yes"), "delta")
})

test_that("recursive newdata = 'all' dispatch forwards delta correctly", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")

  direct <- predict(fitg, "CapeHorn", type = "response", se.fit = TRUE, delta = TRUE)
  via_all <- predict(fitg, "all", type = "response", se.fit = TRUE, delta = TRUE)
  expect_equal(via_all$CapeHorn, direct)
})

# predict weight
test_that("ssn_lm() weight matrix reconstructs the ordinary point prediction exactly (no local restriction)", {
  s <- mf04p
  s$obs$my_offset <- stats::runif(nrow(s$obs), 0, 5)
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM + offset(my_offset), s, tailup_type = "exponential", additive = "afvArea")
  preds <- fit$ssn.object$preds$CapeHorn
  preds$my_offset <- stats::runif(nrow(preds), 0, 5)
  fit$ssn.object$preds$CapeHorn <- preds

  wt <- predict(fit, "CapeHorn", type = "weight")

  y_train <- model.response(model.frame(fit))
  offset_train <- model.offset(model.frame(fit))
  newdata_offset <- preds$my_offset

  reconstructed <- as.numeric(wt %*% (y_train - offset_train)) + newdata_offset
  direct <- unname(predict(fit, "CapeHorn", type = "response"))
  expect_equal(reconstructed, direct, tolerance = 1e-8)
})

test_that("ssn_lm() weight matrix has n_pred x n_obs dimensions and the documented row/column identities", {
  s <- mf04p
  s$obs$Summer_mn[c(2, 5)] <- NA
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s, tailup_type = "exponential", additive = "afvArea")

  wt_named <- predict(fit, "CapeHorn", type = "weight")
  expect_equal(dim(wt_named), c(NROW(fit$ssn.object$preds$CapeHorn), fit$n))
  expect_equal(colnames(wt_named), as.character(fit$observed_index))
  expect_null(rownames(wt_named))

  wt_missing <- predict(fit, ".missing", type = "weight")
  expect_equal(dim(wt_missing), c(length(fit$missing_index), fit$n))
  expect_equal(rownames(wt_missing), as.character(fit$missing_index))
  expect_equal(colnames(wt_missing), as.character(fit$observed_index))
})

test_that("ssn_lm() weight matrix under a local neighborhood matches an independent matrix-solve reference and zeros outside it", {
  s <- mf04p
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s, tailup_type = "exponential", additive = "afvArea")

  wt_local <- predict(fit, "CapeHorn", type = "weight", local = list(method = "covariance", size = 10))
  expect_equal(dim(wt_local), c(NROW(fit$ssn.object$preds$CapeHorn), fit$n))

  nz <- which(wt_local[1, ] != 0)
  expect_equal(length(nz), 10)
  expect_true(all(wt_local[1, -nz] == 0))

  # independent reference: rebuild row 1's weight vector from base matrix ops
  # against the exact neighborhood SSN2 chose, with no code shared with
  # get_pred()'s implementation
  cov_full <- covmatrix(fit)
  cov_local <- cov_full[nz, nz, drop = FALSE]
  chol_local <- t(chol(cov_local))
  Xmat_full <- model.matrix(fit)
  Xmat_local <- Xmat_full[nz, , drop = FALSE]
  cov_betahat <- vcov(fit)
  cov_vec_full <- covmatrix(fit, "CapeHorn")
  c0 <- as.numeric(cov_vec_full[1, nz])

  SqrtSigInv_X <- forwardsolve(chol_local, Xmat_local)
  SqrtSigInv_c0 <- forwardsolve(chol_local, c0)
  Xt_SigInv <- t(backsolve(t(chol_local), SqrtSigInv_X))
  betahat_wt <- cov_betahat %*% Xt_SigInv
  residuals_weight <- -1 * Xmat_local %*% betahat_wt
  diag(residuals_weight) <- diag(residuals_weight) + 1
  x0 <- model.matrix(delete.response(terms(fit)), fit$ssn.object$preds$CapeHorn[1, , drop = FALSE])
  x0 <- x0[, colnames(Xmat_full), drop = FALSE]
  reference_row1 <- as.numeric(x0 %*% betahat_wt + crossprod(SqrtSigInv_c0, forwardsolve(chol_local, residuals_weight)))

  expect_equal(unname(wt_local[1, nz]), reference_row1, tolerance = 1e-10)
})

test_that("ssn_lm() type = \"weight\" silently resets se.fit/interval and ignores terms, matching spmodel (T27)", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  reference <- predict(fit, "CapeHorn", type = "weight")

  expect_equal(predict(fit, "CapeHorn", type = "weight", se.fit = TRUE), reference)
  expect_equal(predict(fit, "CapeHorn", type = "weight", interval = "prediction"), reference)
  expect_equal(predict(fit, "CapeHorn", type = "weight", interval = "confidence"), reference)
  fit_terms <- ssn_lm(Summer_mn ~ ELEV_DEM + SLOPE, mf04p, tailup_type = "exponential", additive = "afvArea")
  reference_terms <- predict(fit_terms, "CapeHorn", type = "weight")
  expect_equal(predict(fit_terms, "CapeHorn", type = "weight", terms = "ELEV_DEM"), reference_terms)
  # block Kriging has no single per-observation weight vector, so this combination still errors
  expect_error(predict(fit, "CapeHorn", type = "weight", block = TRUE), "block")
})

test_that("ssn_lm() recursive newdata = 'all' dispatch forwards type = \"weight\" (Bug B regression, T08-03)", {
  s <- mf04p
  s$obs$my_offset <- stats::runif(nrow(s$obs), 0, 5)
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM + offset(my_offset), s, tailup_type = "exponential", additive = "afvArea")
  fit$ssn.object$preds$CapeHorn$my_offset <- stats::runif(nrow(fit$ssn.object$preds$CapeHorn), 0, 5)
  fit$ssn.object$preds$pred1km$my_offset <- stats::runif(nrow(fit$ssn.object$preds$pred1km), 0, 5)

  direct <- predict(fit, "CapeHorn", type = "weight")
  via_all <- predict(fit, "all", type = "weight")
  expect_equal(via_all$CapeHorn, direct)
})

test_that("ssn_glm() weight matrix reconstructs the link-scale point prediction exactly (no local restriction)", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  s$obs$my_offset <- stats::runif(nrow(s$obs), 0, 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM + offset(my_offset), s,
    family = "poisson", tailup_type = "exponential", additive = "afvArea"
  )
  preds <- fitg$ssn.object$preds$CapeHorn
  preds$my_offset <- stats::runif(nrow(preds), 0, 5)
  fitg$ssn.object$preds$CapeHorn <- preds

  wt <- predict(fitg, "CapeHorn", type = "weight")

  w_train <- fitted(fitg, type = "link")
  offset_train <- model.offset(model.frame(fitg))
  newdata_offset <- preds$my_offset

  reconstructed <- as.numeric(wt %*% (w_train - offset_train)) + newdata_offset
  direct_link <- unname(predict(fitg, "CapeHorn", type = "link"))
  expect_equal(reconstructed, direct_link, tolerance = 1e-8)

  direct_response <- unname(predict(fitg, "CapeHorn", type = "response"))
  expect_equal(invlink(reconstructed, fitg$family, NULL), direct_response, tolerance = 1e-6)
})

test_that("ssn_glm() weight matrix has n_pred x n_obs dimensions and the documented row/column identities", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")

  wt <- predict(fitg, "CapeHorn", type = "weight")
  expect_equal(dim(wt), c(NROW(fitg$ssn.object$preds$CapeHorn), fitg$n))
  expect_equal(colnames(wt), as.character(fitg$observed_index))
  expect_null(rownames(wt))
})

test_that("ssn_glm() type = \"weight\" silently resets se.fit/interval and ignores terms, matching spmodel (T27)", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")
  reference <- predict(fitg, "CapeHorn", type = "weight")

  expect_equal(predict(fitg, "CapeHorn", type = "weight", se.fit = TRUE), reference)
  expect_equal(predict(fitg, "CapeHorn", type = "weight", interval = "prediction"), reference)
  expect_equal(predict(fitg, "CapeHorn", type = "weight", interval = "confidence"), reference)
  fitg_terms <- ssn_glm(count_response ~ ELEV_DEM + SLOPE, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")
  reference_terms <- predict(fitg_terms, "CapeHorn", type = "weight")
  expect_equal(predict(fitg_terms, "CapeHorn", type = "weight", terms = "ELEV_DEM"), reference_terms)
})

test_that("ssn_glm() recursive newdata = 'all' dispatch forwards type = \"weight\"", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")

  direct <- predict(fitg, "CapeHorn", type = "weight")
  via_all <- predict(fitg, "all", type = "weight")
  expect_equal(via_all$CapeHorn, direct)
})

# augment fitted mean ci
test_that("augment.ssn_lm() interval = 'confidence' (no newdata) matches an independent fitted-mean CI reference, and se_fit = TRUE no longer crashes", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")

  Xmat <- model.matrix(fit)
  Xmat_list <- split(Xmat, seq_len(NROW(Xmat)))
  vars_reference <- as.numeric(vapply(Xmat_list, function(x) crossprod(x, vcov(fit) %*% x), numeric(1)))
  se_reference <- sqrt(vars_reference)
  tstar <- qnorm(1 - (1 - 0.95) / 2)
  lwr_reference <- unname(fitted(fit)) - tstar * se_reference
  upr_reference <- unname(fitted(fit)) + tstar * se_reference

  aug <- augment(fit, interval = "confidence")
  expect_equal(unname(aug$.lower), lwr_reference, tolerance = 1e-10)
  expect_equal(unname(aug$.upper), upr_reference, tolerance = 1e-10)

  # previously crashed unconditionally; now succeeds and matches the reference
  aug_se <- augment(fit, se_fit = TRUE)
  expect_equal(unname(aug_se$.se.fit), se_reference, tolerance = 1e-10)
})

test_that("augment.ssn_glm() interval = 'confidence' (no newdata) matches an independent fitted-mean CI reference on both scales, and se.fit never changes with type.predict", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")

  Xmat <- model.matrix(fitg)
  Xmat_list <- split(Xmat, seq_len(NROW(Xmat)))
  vars_reference <- as.numeric(vapply(Xmat_list, function(x) crossprod(x, vcov(fitg) %*% x), numeric(1)))
  se_reference <- sqrt(vars_reference)
  tstar <- qnorm(1 - (1 - 0.95) / 2)
  fitted_link <- unname(fitted(fitg, type = "link"))
  lwr_link_reference <- fitted_link - tstar * se_reference
  upr_link_reference <- fitted_link + tstar * se_reference

  aug_link <- augment(fitg, interval = "confidence", type.predict = "link")
  expect_equal(unname(aug_link$.lower), lwr_link_reference, tolerance = 1e-10)
  expect_equal(unname(aug_link$.upper), upr_link_reference, tolerance = 1e-10)

  aug_response <- augment(fitg, interval = "confidence", type.predict = "response")
  expect_equal(unname(aug_response$.lower), invlink(lwr_link_reference, fitg$family, fitg$size), tolerance = 1e-10)
  expect_equal(unname(aug_response$.upper), invlink(upr_link_reference, fitg$family, fitg$size), tolerance = 1e-10)

  # se.fit stays link-scale regardless of type.predict (matches predict()'s
  # own confidence-branch convention; delta is never applied here)
  aug_se_link <- augment(fitg, se_fit = TRUE, type.predict = "link")
  aug_se_response <- augment(fitg, se_fit = TRUE, type.predict = "response")
  expect_equal(unname(aug_se_link$.se.fit), se_reference, tolerance = 1e-10)
  expect_equal(unname(aug_se_response$.se.fit), se_reference, tolerance = 1e-10)
})

test_that("augment() warns and degrades interval = 'prediction' without newdata, matching spmodel (T27); still supports it with newdata", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  expect_warning(aug_degraded <- augment(fit, interval = "prediction"), "prediction")
  expect_false(any(c(".lower", ".upper") %in% names(aug_degraded)))
  expect_equal(aug_degraded, suppressWarnings(augment(fit, interval = "none")))
  aug_newdata <- augment(fit, newdata = "CapeHorn", interval = "prediction")
  expect_true(all(c(".lower", ".upper") %in% names(aug_newdata)))

  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fitg <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")
  expect_warning(aug_degraded_g <- augment(fitg, interval = "prediction"), "prediction")
  expect_false(any(c(".lower", ".upper") %in% names(aug_degraded_g)))
  aug_newdata_g <- augment(fitg, newdata = "CapeHorn", interval = "prediction")
  expect_true(all(c(".lower", ".upper") %in% names(aug_newdata_g)))
})

test_that("augment() default behavior (interval = 'none', se_fit = FALSE) is unaffected", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  aug <- augment(fit)
  expect_false(any(c(".lower", ".upper", ".se.fit") %in% names(aug)))
  expect_equal(unname(aug$.fitted), unname(fitted(fit)))
})

test_that("predict()'s point-level path is unaffected by local$chunk_size (T25)", {
  # regression test: point-level predict() now builds its covariance in
  # chunk_size-bounded row-chunks over newdata (get_point_pred_cov_vector_list(),
  # mirroring block prediction's chunking) instead of one unbounded call. A
  # tiny chunk_size must give identical results to the default, for both
  # ssn_lm() and ssn_glm(), and for both local$method values. ssn_glm()'s
  # se.fit = TRUE/interval = "prediction" path additionally exercises
  # get_wts_varw()'s "predvar_adjust_all" branch, which needs the full
  # observed-by-prediction covariance matrix even though it is now built from
  # chunks -- this caught a real regression (c0 accidentally resolving to the
  # unrelated cov_vector() function after chunking removed the local variable
  # that used to shadow it) before this test was added.
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea"
  )

  ref_all <- predict(fit, "CapeHorn", se.fit = TRUE, local = list(method = "all"))
  tiny_all <- predict(fit, "CapeHorn", se.fit = TRUE, local = list(method = "all", chunk_size = 2))
  expect_equal(tiny_all, ref_all)

  ref_cov <- predict(fit, "CapeHorn", se.fit = TRUE, local = list(method = "covariance", size = 50))
  tiny_cov <- predict(fit, "CapeHorn", se.fit = TRUE, local = list(method = "covariance", size = 50, chunk_size = 2))
  expect_equal(tiny_cov, ref_cov)

  fit_glm <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, family = "Gamma",
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea"
  )

  ref_glm_pred <- predict(fit_glm, "CapeHorn", se.fit = TRUE, interval = "prediction", local = list(method = "all"))
  tiny_glm_pred <- predict(fit_glm, "CapeHorn", se.fit = TRUE, interval = "prediction", local = list(method = "all", chunk_size = 2))
  expect_equal(tiny_glm_pred, ref_glm_pred)

  ref_glm_none <- predict(fit_glm, "CapeHorn", se.fit = TRUE, local = list(method = "all"))
  tiny_glm_none <- predict(fit_glm, "CapeHorn", se.fit = TRUE, local = list(method = "all", chunk_size = 2))
  expect_equal(tiny_glm_none, ref_glm_none)
})

test_that("point prediction dispatch keeps worker settings separate from dispatch settings", {
  withr::local_options(warnPartialMatchArgs = TRUE)
  worker <- function(x, local) x + local$increment
  for (parallel in c(FALSE, TRUE)) {
    expect_warning(result <- run_pred_dispatch(
      worker, list(1, 2, 3), local_list = list(parallel = parallel, ncores = 2),
      local = list(increment = 3)
    ), NA)
    expect_equal(unlist(result), c(4, 5, 6))
  }
})

test_that("point prediction's parallel dispatch matches serial output and tears down its cluster", {
  # showConnections() always includes some baseline rows in this environment
  # (unrelated to parallel clusters), so compare counts before/after rather
  # than asserting an absolute zero
  baseline_connections <- NROW(showConnections())

  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")

  serial <- predict(fit, "CapeHorn", se.fit = TRUE, local = list(method = "all", parallel = FALSE))
  parallel_res <- predict(fit, "CapeHorn", se.fit = TRUE, local = list(method = "all", parallel = TRUE, ncores = 2))
  expect_equal(parallel_res, serial)
  expect_equal(NROW(showConnections()), baseline_connections)

  fit_glm <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, family = "Gamma",
    tailup_type = "exponential", additive = "afvArea"
  )
  serial_glm <- predict(fit_glm, "CapeHorn", se.fit = TRUE, local = list(method = "all", parallel = FALSE))
  parallel_glm <- predict(fit_glm, "CapeHorn", se.fit = TRUE, local = list(method = "all", parallel = TRUE, ncores = 2))
  expect_equal(parallel_glm, serial_glm)
  expect_equal(NROW(showConnections()), baseline_connections)
})

test_that("point prediction's parallel dispatch tears down its cluster even when a worker errors", {
  baseline_connections <- NROW(showConnections())

  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  broken_get_pred <- function(...) stop("induced worker error")
  testthat::local_mocked_bindings(get_pred = broken_get_pred, .package = "SSN2")
  expect_error(
    predict(fit, "CapeHorn", local = list(method = "all", parallel = TRUE, ncores = 2)),
    "induced worker error"
  )
  expect_equal(NROW(showConnections()), baseline_connections)
})

test_that("GLM point prediction's parallel dispatch tears down its cluster even when a worker errors", {
  baseline_connections <- NROW(showConnections())

  fit_glm <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, family = "Gamma", tailup_type = "exponential", additive = "afvArea")
  broken_get_pred_glm <- function(...) stop("induced worker error")
  testthat::local_mocked_bindings(get_pred_glm = broken_get_pred_glm, .package = "SSN2")
  expect_error(
    predict(fit_glm, "CapeHorn", local = list(method = "all", parallel = TRUE, ncores = 2)),
    "induced worker error"
  )
  expect_equal(NROW(showConnections()), baseline_connections)
})

test_that("terms and confidence-interval predictions never request either covariance provider, for both model classes", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  fit_glm <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, family = "Gamma", tailup_type = "exponential", additive = "afvArea")

  broken_cov_vector_list <- function(...) stop("prediction covariance was requested")
  broken_covmatrix <- function(...) stop("observed covariance was requested")
  testthat::local_mocked_bindings(
    get_point_pred_cov_vector_list = broken_cov_vector_list,
    covmatrix = broken_covmatrix,
    .package = "SSN2"
  )

  expect_no_error(predict(fit, "CapeHorn", type = "terms"))
  expect_no_error(predict(fit, "CapeHorn", interval = "confidence"))
  expect_no_error(predict(fit_glm, "CapeHorn", type = "terms"))
  expect_no_error(predict(fit_glm, "CapeHorn", interval = "confidence"))

  # sanity check: an ordinary prediction on the same fits does hit one of
  # the mocks, confirming they are wired up correctly above
  expect_error(predict(fit, "CapeHorn"), "prediction covariance was requested")
})

test_that("predict_terms() standard errors go through the diagonal-only helper and match the full-matrix formula", {
  # Compare captured diagonal-helper inputs with the full-matrix variance formula.
  real_diag <- get_diag_XVXt
  capture_diag_calls <- function(fit, ...) {
    captured <- list()
    checking_diag <- function(X, V) {
      captured[[length(captured) + 1]] <<- list(X = X, V = V)
      real_diag(X, V)
    }
    testthat::local_mocked_bindings(get_diag_XVXt = checking_diag, .package = "SSN2")
    predict(fit, "CapeHorn", ...)
    captured
  }

  check_captured <- function(captured) {
    expect_gt(length(captured), 0)
    for (call_args in captured) {
      reference <- diag(call_args$X %*% tcrossprod(call_args$V, call_args$X))
      expect_equal(get_diag_XVXt(call_args$X, call_args$V), reference, tolerance = 1e-8, ignore_attr = TRUE)
    }
  }

  # one-column and multi-column (factor) terms together, so vcov() is
  # genuinely non-diagonal across the coefficients any one term's slice uses
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM + as.factor(netID), mf04p, tailup_type = "exponential", additive = "afvArea")
  fit_glm <- ssn_glm(Summer_mn ~ ELEV_DEM + as.factor(netID), mf04p, family = "Gamma", tailup_type = "exponential", additive = "afvArea")
  fit_no_intercept <- ssn_lm(Summer_mn ~ ELEV_DEM + as.factor(netID) - 1, mf04p, tailup_type = "exponential", additive = "afvArea")

  check_captured(capture_diag_calls(fit, type = "terms", se.fit = TRUE))
  check_captured(capture_diag_calls(fit, type = "terms", interval = "confidence"))
  check_captured(capture_diag_calls(fit_glm, type = "terms", se.fit = TRUE))
  check_captured(capture_diag_calls(fit_no_intercept, type = "terms", se.fit = TRUE))
})

test_that("an NA in a fixed-effect predictor errors instead of silently returning NA-contaminated predictions", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  broken <- fit
  broken$ssn.object$preds$CapeHorn$ELEV_DEM[1] <- NA

  for (call in list(
    quote(predict(broken, "CapeHorn")),
    quote(predict(broken, "CapeHorn", type = "terms")),
    quote(predict(broken, "CapeHorn", interval = "confidence")),
    quote(predict(broken, "CapeHorn", interval = "prediction"))
  )) {
    expect_error(eval(call), "Cannot have NA values in predictors")
  }

  fit_block <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  broken_block <- fit_block
  broken_block$ssn.object$preds$CapeHorn$ELEV_DEM[1] <- NA
  expect_error(predict(broken_block, "CapeHorn", block = TRUE), "Cannot have NA values in predictors")

  fit_glm <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, family = "Gamma", tailup_type = "exponential", additive = "afvArea")
  broken_glm <- fit_glm
  broken_glm$ssn.object$preds$CapeHorn$ELEV_DEM[1] <- NA
  expect_error(predict(broken_glm, "CapeHorn"), "Cannot have NA values in predictors")

  # NA in an unused column doesn't trigger it -- only columns the formula
  # actually references matter
  s <- mf04p
  s$obs$unused_col <- 1
  s$preds$CapeHorn$unused_col <- NA
  expect_true(anyNA(s$preds$CapeHorn$unused_col))
  fit_extra_col <- ssn_lm(Summer_mn ~ ELEV_DEM, s, tailup_type = "exponential", additive = "afvArea")
  expect_no_error(predict(fit_extra_col, "CapeHorn"))

  # a complete-predictor .missing prediction set (missing response only)
  # still works, distinguishing missing predictors from missing responses
  s2 <- mf04p
  s2$obs$Summer_mn[1] <- NA
  fit_missing <- ssn_lm(Summer_mn ~ ELEV_DEM, s2, tailup_type = "exponential", additive = "afvArea")
  expect_no_error(predict(fit_missing, ".missing"))
})


test_that("run_pred_dispatch() falls back to loadNamespace() when SSN2 is not a development install", {
  # complements the random-slope parallel-dispatch regression test (see
  # test-extras-random-effects.R), which only exercises the pkgload::load_all()
  # branch this session's dev-source setup actually takes. Forcing
  # pkgload::is_dev_package() to FALSE exercises the dev_path resolution's
  # own fallback logic (getSrcFilename()/normalizePath()/file.exists() are
  # skipped, dev_path stays NULL) and confirms dispatch still completes via
  # loadNamespace("SSN2") on the worker rather than erroring. This only
  # proves the fallback branch itself is reachable and safe, not that it
  # reflects current source -- a worker's loadNamespace("SSN2") loads
  # whatever SSN2 build happens to be installed, which is this session's
  # separately-installed package, not necessarily up to date with the
  # source under test.
  testthat::local_mocked_bindings(is_dev_package = function(...) FALSE, .package = "pkgload")
  result <- run_pred_dispatch(
    function(x) x * 2,
    list(1, 2, 3),
    list(parallel = TRUE, ncores = 2)
  )
  expect_equal(unlist(result), c(2, 4, 6))
})
