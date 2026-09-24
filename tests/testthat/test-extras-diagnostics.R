skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# convergence warnings
expect_no_warning_msg <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(force(expr), warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(msgs, 0)
}

expect_warning_msg <- function(expr, pattern) {
  msgs <- character(0)
  withCallingHandlers(force(expr), warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_true(any(grepl(pattern, msgs)), label = paste("no warning matched", pattern))
}

test_that("warn_optim_convergence() fires only for a nonzero, non-NA convergence code", {
  expect_no_warning_msg(warn_optim_convergence(0))
  expect_no_warning_msg(warn_optim_convergence(NA))
  expect_warning_msg(warn_optim_convergence(1), "did not converge")
  expect_warning_msg(warn_optim_convergence(10), "convergence code 10")
})

test_that("ssn_lm(): a forced non-convergent fit warns; a normal fit does not", {
  expect_warning_msg(
    ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea", control = list(maxit = 1)),
    "did not converge"
  )
  expect_no_warning_msg(
    ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  )
})

test_that("ssn_lm(): a fully known fit (optim() never called) never warns", {
  expect_no_warning_msg(
    ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
      tailup_type = "exponential", additive = "afvArea",
      tailup_initial = tailup_initial("exponential", de = 1, range = 100, known = "given"),
      nugget_initial = nugget_initial("nugget", nugget = 0.05, known = "given")
    )
  )
})

test_that("warn_spcov_boundary() fires only when the total (tailup_de + taildown_de + euclid_de + nugget) is near the diagtol floor", {
  near_zero <- list(tailup = c(de = 1e-20), taildown = c(de = 0), euclid = c(de = 0), nugget = c(nugget = 1e-6))
  well_identified <- list(tailup = c(de = 1.4), taildown = c(de = 0), euclid = c(de = 0), nugget = c(nugget = 0.05))
  expect_warning_msg(warn_spcov_boundary(near_zero, diagtol = 1e-4), "numerical boundary of zero")
  expect_no_warning_msg(warn_spcov_boundary(well_identified, diagtol = 1e-4))
  expect_no_warning_msg(warn_spcov_boundary(near_zero, diagtol = 0)) # diagtol <= 0 means no floor concept, per ssn_lm()'s dense path
})

test_that("ssn_glm(): warn_spcov_boundary() fires only for estmethod = 'ml' with near-zero total spatial + nugget variance", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)

  expect_warning_msg(
    ssn_glm(count_response ~ ELEV_DEM, s,
      family = "poisson", tailup_type = "none", taildown_type = "none", euclid_type = "none",
      nugget_type = "nugget", nugget_initial = nugget_initial("nugget", nugget = 1e-10, known = "given"),
      estmethod = "ml"
    ),
    "numerical boundary of zero"
  )
  expect_no_warning_msg(
    ssn_glm(count_response ~ ELEV_DEM, s,
      family = "poisson", tailup_type = "none", taildown_type = "none", euclid_type = "none",
      nugget_type = "nugget", nugget_initial = nugget_initial("nugget", nugget = 1e-10, known = "given"),
      estmethod = "reml"
    )
  )
})

test_that("warn_fitted_saturation() fires only for binomial with >= 99% saturated fitted probabilities", {
  expect_warning_msg(warn_fitted_saturation(rep(1e-10, 100), "binomial"), "Perfect separation")
  expect_no_warning_msg(warn_fitted_saturation(c(rep(1e-10, 98), 0.5, 0.5), "binomial")) # 98%, below threshold
  expect_no_warning_msg(warn_fitted_saturation(c(0.5, 0.4, 0.6), "binomial"))
  expect_no_warning_msg(warn_fitted_saturation(rep(1e-10, 100), "poisson")) # non-binomial, never fires
})

test_that("ssn_glm(): a well-behaved binomial fit does not spuriously warn about separation", {
  s <- mf04p
  s$obs$y01 <- rbinom(nrow(s$obs), 1, 0.5)
  expect_no_warning_msg(
    ssn_glm(y01 ~ ELEV_DEM, s, family = "binomial", tailup_type = "exponential", additive = "afvArea")
  )
})

test_that("binomial saturation warnings use probabilities with varying trial counts", {
  s <- mf04p
  s$obs$trials <- rep(c(10, 20, 30), length.out = NROW(s$obs))
  s$obs$successes <- rep(c(3, 8, 15), length.out = NROW(s$obs))
  expect_no_warning_msg(
    fit <- ssn_glm(cbind(successes, trials - successes) ~ 1, s,
      family = "binomial", tailup_type = "none",
      nugget_initial = nugget_initial("nugget", nugget = 0.1, known = "given")
    )
  )
  probability <- plogis(fitted(fit, type = "link"))
  expect_true(all(fitted(fit, type = "response") > 1))
  expect_true(all(probability > 0.1 & probability < 0.9))
  expect_equal(fitted(fit, type = "response"), probability * s$obs$trials,
    ignore_attr = TRUE)
})

# deviance pseudoR2 auroc audit
test_that("ssn_glm() Poisson deviance/pseudoR2 respond to a real offset and match an independent reference", {
  n <- nrow(mf04p$obs)
  mf04p_off <- mf04p
  set.seed(2)
  mf04p_off$obs$y_pois <- rpois(n, lambda = 5)
  mf04p_off$obs$off_zero <- rep(0, n)
  mf04p_off$obs$off_var <- log(seq(0.5, 3, length.out = n))

  fit_zero <- ssn_glm(y_pois ~ ELEV_DEM + offset(off_zero), mf04p_off,
    family = "poisson", tailup_type = "exponential", additive = "afvArea"
  )
  fit_var <- ssn_glm(y_pois ~ ELEV_DEM + offset(off_var), mf04p_off,
    family = "poisson", tailup_type = "exponential", additive = "afvArea"
  )

  expect_false(isTRUE(all.equal(deviance(fit_zero), deviance(fit_var))))
  expect_false(isTRUE(all.equal(pseudoR2(fit_zero), pseudoR2(fit_var))))

  y <- mf04p_off$obs$y_pois
  mu <- fitted(fit_var, type = "response")
  reference_dev_i <- 2 * (ifelse(y == 0, 0, y * log(y / mu)) - (y - mu))
  reference_dev <- sum(pmax(reference_dev_i, 0))
  expect_equal(reference_dev, deviance(fit_var), tolerance = 1e-8)
})

test_that("ssn_lm() Gaussian deviance responds to a real offset", {
  n <- nrow(mf04p$obs)
  mf04p_off <- mf04p
  set.seed(2)
  mf04p_off$obs$y_gauss <- 10 + 0.5 * mf04p_off$obs$ELEV_DEM / 100 + rnorm(n, sd = 2)
  mf04p_off$obs$off_zero <- rep(0, n)
  mf04p_off$obs$off_var <- log(seq(0.5, 3, length.out = n))

  fit_zero <- ssn_lm(y_gauss ~ ELEV_DEM + offset(off_zero), mf04p_off, tailup_type = "exponential", additive = "afvArea")
  fit_var <- ssn_lm(y_gauss ~ ELEV_DEM + offset(off_var), mf04p_off, tailup_type = "exponential", additive = "afvArea")

  expect_false(isTRUE(all.equal(deviance(fit_zero), deviance(fit_var))))
})

test_that("ssn_glm() beta deviance is finite, responds to an offset, and matches an independent reference", {
  n <- nrow(mf04p$obs)
  mf04p_off <- mf04p
  set.seed(2)
  mf04p_off$obs$y_beta <- pmin(pmax(0.3 + 0.001 * mf04p_off$obs$ELEV_DEM / 10 + rnorm(n, sd = 0.05), 0.01), 0.99)
  mf04p_off$obs$off_zero <- rep(0, n)
  mf04p_off$obs$off_var <- log(seq(0.8, 1.5, length.out = n))

  fit_zero <- ssn_glm(y_beta ~ ELEV_DEM + offset(off_zero), mf04p_off, family = "beta", tailup_type = "exponential", additive = "afvArea")
  fit_var <- ssn_glm(y_beta ~ ELEV_DEM + offset(off_var), mf04p_off, family = "beta", tailup_type = "exponential", additive = "afvArea")

  expect_true(is.finite(deviance(fit_zero)))
  expect_true(is.finite(deviance(fit_var)))
  expect_false(isTRUE(all.equal(deviance(fit_zero), deviance(fit_var))))

  y <- mf04p_off$obs$y_beta
  mu <- fitted(fit_var, type = "response")
  dispersion <- fit_var$coefficients$params_object$dispersion
  constant <- lgamma(mu * dispersion) + lgamma((1 - mu) * dispersion) - lgamma(y * dispersion) - lgamma((1 - y) * dispersion)
  reference_dev_i <- 2 * (constant + (y - mu) * dispersion * log(y) + ((1 - y) - (1 - mu)) * dispersion * log(1 - y))
  reference_dev <- sum(pmax(reference_dev_i, 0))
  expect_equal(reference_dev, deviance(fit_var), tolerance = 1e-6)
})

test_that("AUROC() enforces its family/size restrictions and matches an independent ROC-AUC reference", {
  n <- nrow(mf04p$obs)
  mf04p_auc <- mf04p
  set.seed(2)
  mf04p_auc$obs$y01 <- rbinom(n, 1, 0.5)
  fit_bin <- ssn_glm(y01 ~ ELEV_DEM, mf04p_auc, family = "binomial", tailup_type = "exponential", additive = "afvArea")

  a <- AUROC(fit_bin)
  expect_true(is.numeric(a) && length(a) == 1 && a >= 0 && a <= 1)

  mu <- fitted(fit_bin, type = "response")
  y <- as.vector(fit_bin$y)
  thresh <- sort(unique(mu))
  tpr <- sapply(thresh, function(t) sum(mu >= t & y == 1) / sum(y == 1))
  fpr <- sapply(thresh, function(t) sum(mu >= t & y == 0) / sum(y == 0))
  ord <- order(fpr, tpr)
  fpr_o <- fpr[ord]
  tpr_o <- tpr[ord]
  manual_auc <- sum(diff(fpr_o) * (head(tpr_o, -1) + tail(tpr_o, -1)) / 2)
  expect_equal(manual_auc, a, tolerance = 0.01)

  # invalid trigger: non-binomial family
  fit_gamma <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, family = "Gamma", tailup_type = "exponential", additive = "afvArea")
  expect_error(AUROC(fit_gamma), "only available when family")

  # invalid trigger: aggregated (size > 1) binomial
  mf04p_agg <- mf04p
  mf04p_agg$obs$successes <- rbinom(n, 5, 0.5)
  mf04p_agg$obs$failures <- 5 - mf04p_agg$obs$successes
  fit_agg <- ssn_glm(cbind(successes, failures) ~ ELEV_DEM, mf04p_agg, family = "binomial", tailup_type = "exponential", additive = "afvArea")
  expect_error(AUROC(fit_agg), "only available for binary models")
})

# glm residual scaling
test_that("beta leverage weights equal expected link-scale information", {
  w <- c(-2, 0, 2)
  dispersion <- 10
  mu <- plogis(w)
  expected <- dispersion^2 * (mu * (1 - mu))^2 *
    (trigamma(mu * dispersion) + trigamma((1 - mu) * dispersion))
  expect_equal(get_V(w, "beta", rep(1, length(w)), dispersion), expected)

  integrated <- vapply(seq_along(w), function(i) {
    integrate(function(y) {
      -vapply(y, function(value) {
        get_D("beta", w[[i]], value, 1, dispersion)[1, 1]
      }, numeric(1)) * dbeta(y, mu[[i]] * dispersion, (1 - mu[[i]]) * dispersion)
    }, 0, 1, rel.tol = 1e-9)$value
  }, numeric(1))
  expect_equal(integrated, expected, tolerance = 1e-7)
  expect_equal(
    get_var_y(w, "beta", rep(1, length(w)), dispersion),
    mu * (1 - mu) / (1 + dispersion)
  )
})


test_that("get_dispersion_factor() matches the documented family semantics", {
  w <- c(-1, 0, 1, 2)
  dispersion <- 3
  expect_equal(get_dispersion_factor(w, "poisson", NULL, dispersion), rep(1, 4))
  expect_equal(get_dispersion_factor(w, "binomial", 5, dispersion), rep(1, 4))
  expect_equal(get_dispersion_factor(w, "nbinomial", NULL, dispersion), rep(1, 4))
  expect_equal(get_dispersion_factor(w, "beta", NULL, dispersion), rep(1, 4))
  expect_equal(get_dispersion_factor(w, "Gamma", NULL, dispersion), rep(1 / dispersion, 4))
  expect_equal(get_dispersion_factor(w, "inverse.gaussian", NULL, dispersion), 1 / (exp(w) * dispersion))
})

test_that("Gamma hatvalues/standardized residuals/cooks.distance match an independent from-scratch reference", {
  set.seed(2)
  s <- mf04p
  s$obs$gamma_response <- rgamma(nrow(s$obs), shape = 2, rate = 2 / (5 + 0.1 * s$obs$ELEV_DEM %% 10))
  fit <- ssn_glm(gamma_response ~ ELEV_DEM, s, family = "Gamma", tailup_type = "exponential", additive = "afvArea")

  X <- model.matrix(fit)
  dispersion <- as.numeric(coef(fit, type = "dispersion"))

  V_reference <- rep(dispersion, nrow(X))
  SqrtVInv_X <- sqrt(V_reference) * X
  cov_vhat <- solve(crossprod(SqrtVInv_X, SqrtVInv_X))
  hv_reference <- unname(diag(SqrtVInv_X %*% tcrossprod(cov_vhat, SqrtVInv_X)))
  expect_equal(unname(hatvalues(fit)), hv_reference, tolerance = 1e-8)

  dev_resid <- residuals(fit, type = "deviance")
  a_phi_reference <- 1 / dispersion
  std_reference <- dev_resid / sqrt(a_phi_reference * (1 - hv_reference))
  expect_equal(unname(residuals(fit, type = "standardized")), unname(std_reference), tolerance = 1e-8)
  expect_identical(residuals(fit, type = "standardized"), rstandard(fit))

  cooks_reference <- std_reference^2 * hv_reference / (fit$p * (1 - hv_reference))
  expect_equal(unname(cooks.distance(fit)), unname(cooks_reference), tolerance = 1e-8)

  # Pearson residuals must be numerically unchanged by this fix (independent
  # reference built from mu^2 / dispersion directly, not by calling get_var_y())
  mu <- exp(fitted(fit, type = "link"))
  var_y_reference <- mu^2 / dispersion
  pearson_reference <- residuals(fit, type = "response") / sqrt(var_y_reference)
  expect_equal(unname(residuals(fit, type = "pearson")), unname(pearson_reference), tolerance = 1e-10)
})

test_that("inverse.gaussian hatvalues/standardized residuals/pearson match an independent from-scratch reference", {
  set.seed(2)
  s <- mf04p
  s$obs$ig_response <- rgamma(nrow(s$obs), shape = 3, rate = 3 / (4 + 0.05 * s$obs$ELEV_DEM %% 10))
  fit <- ssn_glm(ig_response ~ ELEV_DEM, s, family = "inverse.gaussian", tailup_type = "exponential", additive = "afvArea")

  X <- model.matrix(fit)
  dispersion <- as.numeric(coef(fit, type = "dispersion"))
  mu <- exp(fitted(fit, type = "link"))

  V_reference <- rep(dispersion + 0.5, nrow(X))
  SqrtVInv_X <- sqrt(V_reference) * X
  cov_vhat <- solve(crossprod(SqrtVInv_X, SqrtVInv_X))
  hv_reference <- unname(diag(SqrtVInv_X %*% tcrossprod(cov_vhat, SqrtVInv_X)))
  expect_equal(unname(hatvalues(fit)), hv_reference, tolerance = 1e-8)

  dev_resid <- residuals(fit, type = "deviance")
  a_phi_reference <- 1 / (mu * dispersion)
  std_reference <- dev_resid / sqrt(a_phi_reference * (1 - hv_reference))
  expect_equal(unname(residuals(fit, type = "standardized")), unname(std_reference), tolerance = 1e-8)

  var_y_reference <- mu^2 / dispersion
  pearson_reference <- residuals(fit, type = "response") / sqrt(var_y_reference)
  expect_equal(unname(residuals(fit, type = "pearson")), unname(pearson_reference), tolerance = 1e-10)
})

test_that("standardized residuals for poisson/binomial/nbinomial/beta are unaffected (a_phi = 1 regression lock)", {
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fit_pois <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")
  expect_equal(
    unname(residuals(fit_pois, type = "standardized")),
    unname(residuals(fit_pois, type = "deviance") / sqrt(1 - hatvalues(fit_pois))),
    tolerance = 1e-10
  )
})

test_that("get_hatvalues_glm() evaluates the leverage weight at the offset-inclusive predictor (any family)", {
  set.seed(2)
  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  s$obs$my_offset <- stats::runif(nrow(s$obs), 0, 3)
  fit <- ssn_glm(count_response ~ ELEV_DEM + offset(my_offset), s,
    family = "poisson", tailup_type = "exponential", additive = "afvArea"
  )

  X <- model.matrix(fit)
  w_link <- fitted(fit, type = "link") # already offset-inclusive
  mu <- exp(w_link)
  V_reference <- mu
  SqrtVInv_X <- sqrt(V_reference) * X
  cov_vhat <- solve(crossprod(SqrtVInv_X, SqrtVInv_X))
  hv_reference <- unname(diag(SqrtVInv_X %*% tcrossprod(cov_vhat, SqrtVInv_X)))
  expect_equal(unname(hatvalues(fit)), hv_reference, tolerance = 1e-8)

  # confirm it actually matters: the offset-free weight gives a different answer
  w_pre_offset <- w_link - s$obs$my_offset
  V_wrong <- exp(w_pre_offset)
  expect_false(isTRUE(all.equal(V_reference, V_wrong, tolerance = 1e-4)))
})

test_that("local (big-data) Gamma fitting with an offset succeeds and matches non-local hatvalues semantics", {
  ssn_create_bigdist(mf04p, overwrite = TRUE)
  set.seed(2)
  s <- mf04p
  s$obs$gamma_response <- rgamma(nrow(s$obs), shape = 2, rate = 2 / (5 + 0.1 * s$obs$ELEV_DEM %% 10))
  s$obs$my_offset <- stats::runif(nrow(s$obs), 0, 1)

  fit_local <- ssn_glm(gamma_response ~ ELEV_DEM + offset(my_offset), s,
    family = "Gamma", tailup_type = "exponential", additive = "afvArea",
    local = list(method = "kmeans", parallel = FALSE)
  )
  expect_s3_class(fit_local, "ssn_glm")
  expect_true(all(is.finite(hatvalues(fit_local))))
  expect_true(all(is.finite(residuals(fit_local, type = "standardized"))))
  expect_true(all(is.finite(cooks.distance(fit_local))))
})


# parity fixes
test_that("Gaussian Cook's distances agree with deletion in whitened coordinates", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, euclid_type = "exponential",
    euclid_initial = euclid_initial("exponential", 0.4, 12000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.2, known = "given"))
  eig <- eigen(as.matrix(covmatrix(fit)), symmetric = TRUE)
  whitening <- eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors)
  X <- whitening %*% model.matrix(fit)
  y <- whitening %*% model.response(model.frame(fit))
  beta <- qr.solve(X, y)
  expected <- vapply(seq_len(nrow(X)), function(i) {
    deleted <- qr.solve(X[-i, , drop = FALSE], y[-i])
    sum((X %*% (beta - deleted))^2) / ncol(X)
  }, numeric(1))
  expect_equal(unname(cooks.distance(fit)), expected, tolerance = 1e-8)
  expect_equal(unname(influence(fit)$.cooksd), expected, tolerance = 1e-8)
  expect_equal(unname(augment(fit)$.cooksd), expected, tolerance = 1e-8)
})

test_that("partitioned spatial fitted components reconstruct the residual in original order", {
  ssn_create_bigdist(mf04p, overwrite = TRUE, no_cores = 1, verbose = FALSE)
  network <- mf04p
  network$obs <- network$obs[c(45:23, 1:22), ]
  network$obs$partition <- factor(rep(1:2, length.out = nrow(network$obs)))
  network$obs$x <- as.numeric(scale(network$obs$ELEV_DEM))
  network$obs$off <- seq(-0.1, 0.1, length.out = nrow(network$obs))
  set.seed(2)
  network$obs$count <- rpois(nrow(network$obs), exp(0.5 + 0.3 * network$obs$x + network$obs$off))
  group <- rep(1:3, length.out = nrow(network$obs))
  for (glm in c(FALSE, TRUE)) for (streams in c(FALSE, TRUE)) for (local in list(NULL, list(index = group, var_adjust = "none"))) {
    common <- list(ssn.object = network, partition_factor = ~ partition, local = local,
      euclid_initial = euclid_initial("exponential", 0.4, 12000, known = "given"),
      nugget_initial = nugget_initial("nugget", 0.5, known = "given"))
    if (streams) {
      common$tailup_initial <- tailup_initial("exponential", 0.2, 12000, known = "given")
      common$taildown_initial <- taildown_initial("exponential", 0.3, 12000, known = "given")
      common$additive <- "afvArea"
    }
    fit <- if (glm) do.call(ssn_glm, c(list(count ~ x + offset(off), family = "poisson"), common)) else
      do.call(ssn_lm, c(list(Summer_mn ~ x + offset(off)), common))
    response <- if (glm) fitted(fit, "link") else model.response(model.frame(fit))
    residual <- as.numeric(response - model.offset(model.frame(fit)) - model.matrix(fit) %*% coef(fit))
    sigma <- as.matrix(covmatrix(fit))
    if (!is.null(local)) sigma <- sigma * outer(group, group, "==")
    inverse_residual <- solve(sigma, residual)
    h <- as.matrix(dist(sf::st_coordinates(network$obs)))
    mask <- outer(network$obs$partition, network$obs$partition, "==")
    if (!is.null(local)) mask <- mask * outer(group, group, "==")
    euclid_expected <- as.numeric((0.4 * exp(-h / 12000) * mask) %*% inverse_residual)
    expect_equal(unname(fitted(fit, "euclid")), euclid_expected, tolerance = 1e-7)
    expect_equal(unname(fitted(fit, "nugget")), as.numeric(0.5 * inverse_residual), tolerance = 1e-7)
    components <- fitted(fit, "euclid") + fitted(fit, "nugget")
    if (streams) components <- components + fitted(fit, "tailup") + fitted(fit, "taildown")
    expect_equal(unname(components), residual, tolerance = 1e-7)
  }
})

test_that("decorrelation resolves isotropic defaults and validates missing fixed parameters", {
  network <- mf04p
  network$preds$CapeHorn <- network$preds$CapeHorn[1:4, ]
  # euclid_params() always resolves rotate/scale to concrete values (0/1
  # when omitted) at construction time -- unlike euclid_initial(), there is
  # no separate "implicit vs explicit" state to compare downstream, so this
  # just exercises every type/local combination for a well-formed transform
  for (type in c("exponential", "circular", "matern", "cauchy", "pexponential", "none")) {
    args <- list(type, de = 0.4, range = 12000)
    if (type %in% c("matern", "cauchy", "pexponential")) args$extra <- 1
    euclid <- do.call(euclid_params, args)
    for (local in list(FALSE, list(method = "covariance", size = 12))) {
      transformed <- ssn_decorrelate_data(Summer_mn ~ ELEV_DEM, network, ordering = "none", local = local,
        euclid_params = euclid, nugget_params = nugget_params("nugget", 0.6))
      expect_true(all(is.finite(transformed$tX)))
      expect_true(all(is.finite(transformed$ty)))
      expect_true(all(is.finite(ssn_decorrelate_newdata(transformed, "CapeHorn")$tX_newdata)))
    }
  }
  # unlike euclid_initial(), a euclid_params() object cannot be constructed
  # with a missing required parameter in the first place -- the "must be
  # supplied and known" check on ssn_decorrelate_data() is unreachable via
  # euclid_params() and the validation now happens earlier, at construction
  expect_error(euclid_params("exponential", range = 12000), "de")
  expect_error(euclid_params("matern", de = 0.4, range = 12000), "extra must be specified")

  # anisotropy is no longer a ssn_decorrelate_data() argument -- a nonzero
  # rotate or a scale other than one is enough to trigger it automatically
  euclid <- euclid_params("exponential", 0.4, 12000, rotate = 0.5, scale = 0.4)
  transformed <- ssn_decorrelate_data(
    Summer_mn ~ ELEV_DEM, network, euclid_params = euclid,
    nugget_params = nugget_params("nugget", 0.6), ordering = "none", local = FALSE
  )
  expect_true(all(is.finite(transformed$tX)))
})

test_that("GLM augmentation forwards prediction variance correction on both scales", {
  network <- mf04p
  network$preds$CapeHorn <- network$preds$CapeHorn[1:4, ]
  fit <- ssn_glm(C16 ~ ELEV_DEM, network, family = "poisson",
    euclid_initial = euclid_initial("exponential", 0.15, 12000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.12, known = "given"))
  corrected <- augment(fit, newdata = "CapeHorn", se_fit = TRUE)$.se.fit
  uncorrected <- augment(fit, newdata = "CapeHorn", se_fit = TRUE, var_correct = FALSE)$.se.fit
  expect_gt(max(abs(corrected - uncorrected)), 1e-5)
  for (correction in c(FALSE, TRUE)) for (type in c("link", "response")) for (interval in c("none", "confidence", "prediction")) {
    expected <- predict(fit, "CapeHorn", type = type, se.fit = TRUE, interval = interval, var_correct = correction)
    actual <- augment(fit, newdata = "CapeHorn", type.predict = type, se_fit = TRUE,
      interval = interval, var_correct = correction)
    expect_equal(actual$.se.fit, unname(expected$se.fit))
    if (interval != "none") {
      expect_equal(actual$.fitted, unname(expected$fit[, "fit"]))
      expect_equal(actual$.lower, unname(expected$fit[, "lwr"]))
      expect_equal(actual$.upper, unname(expected$fit[, "upr"]))
      without_se <- augment(fit, newdata = "CapeHorn", type.predict = type,
        interval = interval, var_correct = correction)
      expect_equal(without_se$.lower, actual$.lower)
      expect_equal(without_se$.upper, actual$.upper)
    } else {
      expect_equal(actual$.fitted, unname(expected$fit))
    }
  }
})

