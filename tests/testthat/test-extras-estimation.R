skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# estimation helpers
test_that("orig2optim_ssn_components / optim2orig_ssn_components round-trip natural-scale values and carry the known mask through", {
  initial_object <- list(
    tailup_initial = list(initial = c(de = 2.5, range = 8000), is_known = c(de = FALSE, range = FALSE)),
    taildown_initial = list(initial = c(de = 1.2, range = 5000), is_known = c(de = FALSE, range = TRUE)),
    euclid_initial = list(
      initial = c(de = 0.8, range = 3000, rotate = 1.1, scale = 0.4),
      is_known = c(de = FALSE, range = FALSE, rotate = FALSE, scale = FALSE)
    ),
    nugget_initial = list(initial = c(nugget = 0.3), is_known = c(nugget = FALSE))
  )

  optim_val <- orig2optim_ssn_components(initial_object)
  recovered <- optim2orig_ssn_components(optim_val$value)

  expect_equal(unname(recovered[["tailup_de"]]), 2.5)
  expect_equal(unname(recovered[["tailup_range"]]), 8000)
  expect_equal(unname(recovered[["taildown_de"]]), 1.2)
  expect_equal(unname(recovered[["taildown_range"]]), 5000)
  expect_equal(unname(recovered[["euclid_de"]]), 0.8)
  expect_equal(unname(recovered[["euclid_range"]]), 3000)
  expect_equal(unname(recovered[["euclid_rotate"]]), 1.1)
  expect_equal(unname(recovered[["euclid_scale"]]), 0.4)
  expect_equal(unname(recovered[["nugget"]]), 0.3)

  # partly known: only taildown de/range's "range" is marked known, everything else free
  expect_equal(unname(optim_val$is_known), c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
})

test_that("orig2optim_ssn_components preserves a fully-known mask", {
  initial_object <- list(
    tailup_initial = list(initial = c(de = 1, range = 1), is_known = c(de = TRUE, range = TRUE)),
    taildown_initial = list(initial = c(de = 1, range = 1), is_known = c(de = TRUE, range = TRUE)),
    euclid_initial = list(
      initial = c(de = 1, range = 1, rotate = 0, scale = 1),
      is_known = c(de = TRUE, range = TRUE, rotate = TRUE, scale = TRUE)
    ),
    nugget_initial = list(initial = c(nugget = 1), is_known = c(nugget = TRUE))
  )
  optim_val <- orig2optim_ssn_components(initial_object)
  expect_true(all(optim_val$is_known))
})

test_that("orig2optim_randcov_components / optim2orig_randcov_components round-trip, report n_est, and are NULL with no random effects", {
  none_val <- orig2optim_randcov_components(list(randcov_initial = NULL))
  expect_null(none_val$value)
  expect_null(none_val$is_known)
  expect_equal(none_val$n_est, 0)
  expect_null(optim2orig_randcov_components(NULL))

  initial_object <- list(randcov_initial = list(
    initial = c(group1 = 4, group2 = 9),
    is_known = c(group1 = FALSE, group2 = TRUE)
  ))
  optim_val <- orig2optim_randcov_components(initial_object)
  expect_equal(optim_val$n_est, 1) # only group1 is free
  expect_equal(unname(optim_val$is_known), c(FALSE, TRUE))

  recovered <- optim2orig_randcov_components(optim_val$value)
  expect_equal(unname(recovered), c(4, 9))
  expect_equal(names(recovered), c("group1", "group2"))
})

test_that("clamp_optim_scale bounds only free (not known) parameters to [-50, 50]", {
  value <- c(a = 100, b = -100, c = 100, d = 10)
  is_known <- c(a = FALSE, b = FALSE, c = TRUE, d = FALSE)
  clamped <- clamp_optim_scale(value, is_known)
  expect_equal(unname(clamped), c(50, -50, 100, 10))
})

test_that("loglik_is_known_base returns the four shared is_known vectors unchanged", {
  initial_object <- list(
    tailup_initial = list(is_known = c(de = TRUE, range = FALSE)),
    taildown_initial = list(is_known = c(de = FALSE, range = FALSE)),
    euclid_initial = list(is_known = c(de = FALSE, range = TRUE, rotate = FALSE, scale = FALSE)),
    nugget_initial = list(is_known = c(nugget = FALSE))
  )
  base <- loglik_is_known_base(initial_object)
  expect_named(base, c("tailup", "taildown", "euclid", "nugget"))
  expect_equal(base$tailup, initial_object$tailup_initial$is_known)
  expect_equal(base$taildown, initial_object$taildown_initial$is_known)
  expect_equal(base$euclid, initial_object$euclid_initial$is_known)
  expect_equal(base$nugget, initial_object$nugget_initial$is_known)
})

test_that("known_optim_output_stub / trim_optim_output produce the expected optim()-shaped lists", {
  stub <- known_optim_output_stub(12.3)
  expect_equal(stub$value, 12.3)
  expect_true(is.na(stub$convergence))
  expect_named(stub, c("method", "control", "value", "counts", "convergence", "message", "hessian"))

  optim_output <- list(value = 4.5, counts = c(10, 2), convergence = 0, message = NULL, hessian = matrix(1))
  dotlist_no_hessian <- list(method = "Nelder-Mead", control = list(reltol = 1e-4), hessian = FALSE)
  trimmed <- trim_optim_output(optim_output, dotlist_no_hessian)
  expect_equal(trimmed$method, "Nelder-Mead")
  expect_equal(trimmed$value, 4.5)
  expect_false(trimmed$hessian)

  dotlist_hessian <- list(method = "Nelder-Mead", control = list(reltol = 1e-4), hessian = TRUE)
  trimmed_hess <- trim_optim_output(optim_output, dotlist_hessian)
  expect_equal(trimmed_hess$hessian, optim_output$hessian)
})

test_that("resolve_anis_rotation resolves to whichever of {rotate, pi - rotate} has the lower -2 log likelihood (Gaussian)", {
  initial_object_val <- get_initial_object(
    tailup_type = "none", taildown_type = "none", euclid_type = "exponential", nugget_type = "nugget",
    tailup_initial = NULL, taildown_initial = NULL, euclid_initial = NULL, nugget_initial = NULL
  )
  data_object <- get_data_object(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p, additive = "afvArea", anisotropy = TRUE,
    initial_object = initial_object_val, random = NULL, randcov_initial = NULL,
    partition_factor = NULL, local = NULL
  )
  cov_grid_vector <- c(
    tailup_de = 0.1, tailup_range = 1, taildown_de = 0.1, taildown_range = 1,
    euclid_de = 5, euclid_range = 3000, rotate = 0.4, scale = 0.3, nugget = 1
  )
  params_object <- get_params_object_grid(cov_grid_vector, initial_object_val)

  minustwoll_direct <- get_minustwologlik(gloglik_products(params_object, data_object, "reml"), data_object, "reml")
  params_object_q2 <- params_object
  params_object_q2$euclid[["rotate"]] <- pi - params_object_q2$euclid[["rotate"]]
  minustwoll_q2 <- get_minustwologlik(gloglik_products(params_object_q2, data_object, "reml"), data_object, "reml")
  # confirms the two rotation candidates actually differ on this fixture, so
  # the test is not vacuously true
  expect_true(abs(minustwoll_direct - minustwoll_q2) > 1e-6)

  resolved <- resolve_anis_rotation(params_object, data_object, "reml", is_glm = FALSE)
  expected_rotate <- if (minustwoll_q2 < minustwoll_direct) {
    params_object_q2$euclid[["rotate"]]
  } else {
    params_object$euclid[["rotate"]]
  }
  expect_equal(unname(resolved$euclid[["rotate"]]), unname(expected_rotate))

  resolved_minustwoll <- get_minustwologlik(gloglik_products(resolved, data_object, "reml"), data_object, "reml")
  expect_true(resolved_minustwoll <= minustwoll_direct + 1e-9)
  expect_true(resolved_minustwoll <= minustwoll_q2 + 1e-9)
})

test_that("resolve_anis_rotation dispatches to the Laplace-approximation (GLM) products/loss when is_glm = TRUE", {
  s <- mf04p
  s$obs$y_pois <- round(s$obs$Summer_mn)

  initial_object_glm_val <- get_initial_object_glm(
    tailup_type = "none", taildown_type = "none", euclid_type = "exponential", nugget_type = "nugget",
    tailup_initial = NULL, taildown_initial = NULL, euclid_initial = NULL, nugget_initial = NULL,
    family = "poisson", dispersion_initial = NULL
  )
  data_object_glm <- get_data_object_glm(
    formula = y_pois ~ ELEV_DEM, ssn.object = s, family = "poisson", additive = "afvArea", anisotropy = TRUE,
    initial_object = initial_object_glm_val, random = NULL, randcov_initial = NULL,
    partition_factor = NULL, local = NULL
  )
  cov_grid_vector <- c(
    tailup_de = 0.1, tailup_range = 1, taildown_de = 0.1, taildown_range = 1,
    euclid_de = 2, euclid_range = 3000, rotate = 0.4, scale = 0.3, nugget = 0.5, dispersion = 1
  )
  params_object <- get_params_object_grid_glm(cov_grid_vector, initial_object_glm_val)

  minustwoll_direct <- get_minustwolaploglik(
    laploglik_products(params_object, data_object_glm, "reml"), data_object_glm, "reml"
  )
  params_object_q2 <- params_object
  params_object_q2$euclid[["rotate"]] <- pi - params_object_q2$euclid[["rotate"]]
  minustwoll_q2 <- get_minustwolaploglik(
    laploglik_products(params_object_q2, data_object_glm, "reml"), data_object_glm, "reml"
  )

  resolved <- resolve_anis_rotation(params_object, data_object_glm, "reml", is_glm = TRUE)
  expected_rotate <- if (minustwoll_q2 < minustwoll_direct) {
    params_object_q2$euclid[["rotate"]]
  } else {
    params_object$euclid[["rotate"]]
  }
  expect_equal(unname(resolved$euclid[["rotate"]]), unname(expected_rotate))
})

# laploglik
test_that("beta family Hessian (get_D) matches independent finite differences", {
  w <- c(-0.8, 0.2, 1.1)
  y <- c(0.2, 0.4, 0.7)
  size <- rep(5, 3)
  dispersion <- 3
  h <- 1e-5
  numeric_D <- (get_d("beta", w + h, y, size, dispersion) - get_d("beta", w - h, y, size, dispersion)) / (2 * h)
  D <- diag(get_D("beta", w, y, size, dispersion))
  expect_equal(D, numeric_D, tolerance = 1e-6)
})

test_that("beta family Hessian matches finite differences across random interior cases", {
  set.seed(2)
  w <- runif(25, -2, 2)
  y <- runif(25, 0.05, 0.95)
  dispersion <- runif(25, 0.5, 8)
  h <- 1e-5
  for (i in seq_along(w)) {
    numeric_D <- (get_d("beta", w[i] + h, y[i], 1, dispersion[i]) - get_d("beta", w[i] - h, y[i], 1, dispersion[i])) / (2 * h)
    D <- diag(get_D("beta", w[i], y[i], 1, dispersion[i]))
    expect_equal(D, numeric_D, tolerance = 1e-6)
  }
})

test_that("other five glm family Hessians are unaffected by the beta fix", {
  w <- c(-0.8, 0.2, 1.1)
  size <- rep(5, 3)
  dispersion <- 3
  h <- 1e-5
  for (family in c("poisson", "nbinomial", "binomial", "Gamma", "inverse.gaussian")) {
    y <- if (family == "binomial") c(1, 2, 3) else c(1, 2, 4)
    numeric_D <- (get_d(family, w + h, y, size, dispersion) - get_d(family, w - h, y, size, dispersion)) / (2 * h)
    D <- diag(get_D(family, w, y, size, dispersion))
    expect_equal(D, numeric_D, tolerance = 1e-6)
  }
})

test_that("beta ssn_glm model fits and predicts with the corrected Hessian", {
  mf04p_beta <- mf04p
  mf04p_beta$obs$y_beta <- (mf04p_beta$obs$Summer_mn - 8) / 10 # rescaled into (0, 1)

  ssn_mod <- ssn_glm(y_beta ~ ELEV_DEM, mf04p_beta,
    family = "beta", tailup_type = "exponential", nugget_type = "nugget",
    additive = "afvArea"
  )
  expect_s3_class(ssn_mod, "ssn_glm")
  expect_true(is.finite(as.numeric(logLik(ssn_mod))))
  expect_vector(predict(ssn_mod, "pred1km"))
})

# floor components
test_that("component BLUPs exclude stabilization and prediction covariance includes it", {
  s <- mf04p
  s$preds$CapeHorn <- s$preds$CapeHorn[c(1, 7, 15), ]
  xy <- sf::st_coordinates(s$obs)
  xy0 <- sf::st_coordinates(s$preds$CapeHorn)
  K <- 2 * exp(-as.matrix(dist(xy)) / 12000)
  C <- 2 * exp(-sqrt(outer(xy0[, 1], xy[, 1], "-")^2 +
                       outer(xy0[, 2], xy[, 2], "-")^2) / 12000)
  for (nugget in c(0, 1e-9, 0.2)) {
    fit <- ssn_lm(Summer_mn ~ 1, s, euclid_type = "exponential",
      euclid_initial = euclid_initial("exponential", de = 2, range = 12000, known = "given"),
      nugget_initial = nugget_initial("nugget", nugget, known = "given"), ddf = "asymptotic")
    floor <- max(nugget, 2e-4)
    V <- K + diag(floor, nrow(K))
    W <- solve(V)
    X <- matrix(1, nrow(K), 1)
    X0 <- matrix(1, nrow(C), 1)
    B <- solve(t(X) %*% W %*% X)
    beta <- B %*% t(X) %*% W %*% s$obs$Summer_mn
    z <- W %*% (s$obs$Summer_mn - X %*% beta)
    expect_equal(covmatrix(fit), V, tolerance = 1e-10)
    expect_equal(unname(coef(fit, "nugget")[[1]]), nugget)
    expect_equal(unname(fitted(fit, "euclid")), as.vector(K %*% z), tolerance = 1e-8)
    expect_equal(unname(fitted(fit, "nugget")), as.vector(nugget * z), tolerance = 1e-8)
    for (size in c(nrow(K), 12)) {
      local <- if (size == nrow(K)) FALSE else list(method = "covariance", size = size)
      variance <- vapply(seq_len(nrow(C)), function(i) {
        index <- if (size == nrow(K)) seq_len(nrow(K)) else
          order(-abs(C[i, ]), -seq_len(nrow(K)))[seq_len(size)]
        Cs <- C[i, index, drop = FALSE]
        Wsub <- solve(V[index, index, drop = FALSE])
        H <- X0[i, , drop = FALSE] - Cs %*% Wsub %*% X[index, , drop = FALSE]
        as.numeric(2 + floor - Cs %*% Wsub %*% t(Cs) + H %*% B %*% t(H))
      }, numeric(1))
      actual <- predict(fit, "CapeHorn", se.fit = TRUE, local = local)
      expect_equal(unname(actual$se.fit^2), variance, tolerance = 1e-8)
    }
  }
})

test_that("the default Gaussian wrapper forwards local simulation controls", {
  args <- list(ssn.object = mf04p, tailup_params = tailup_params("none"),
    taildown_params = taildown_params("none"),
    euclid_params = euclid_params("exponential", de = 1, range = 12000),
    nugget_params = nugget_params("nugget", nugget = 0.2), samples = 2,
    local = list(approximation = "vecchia", method = "covariance", size = 10))
  set.seed(2)
  default <- do.call(ssn_simulate, args)
  set.seed(2)
  direct <- do.call(ssn_rnorm, args)
  expect_identical(default, direct)
})

# support
test_that("known independent covariance preserves linear model estimates and nobs", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "none", taildown_type = "none", euclid_type = "none",
    nugget_initial = nugget_initial("nugget", nugget = 1, known = "given")
  )
  reference <- lm(Summer_mn ~ ELEV_DEM, data = mf04p$obs)
  expect_equal(coef(fit), coef(reference), tolerance = 1e-8)
  expect_equal(nobs(fit), nobs(reference))
  expect_error(satterthwaite(fit, method = "closed"), 'method = "numeric"', fixed = TRUE)
})

test_that("beta curvature agrees with finite differences of the log density", {
  w <- c(-1.2, 0.3, 1.1)
  y <- c(0.2, 0.6, 0.7)
  phi <- 7
  log_density <- function(eta) {
    mu <- plogis(eta)
    dbeta(y, mu * phi, (1 - mu) * phi, log = TRUE)
  }
  h <- 1e-4
  expected <- (log_density(w + h) - 2 * log_density(w) + log_density(w - h)) / h^2
  expect_equal(as.numeric(diag(SSN2:::get_D("beta", w, y, NULL, phi))), expected, tolerance = 1e-5)
})

# glm offset
test_that("get_w_and_H's latent mode satisfies family score stationarity under zero and varying offsets, across families", {
  X <- cbind(1, c(-1, -0.5, 0, 0.5, 1))
  SigInv <- diag(2, 5)
  SigInv_X <- SigInv %*% X
  precision_beta <- crossprod(X, SigInv_X)
  V <- solve(precision_beta)
  Ptheta <- SigInv - SigInv_X %*% V %*% t(SigInv_X)

  family_y <- list(
    poisson = c(1, 3, 2, 4, 6),
    nbinomial = c(1, 3, 2, 4, 6),
    binomial = c(1, 3, 2, 4, 5),
    Gamma = c(0.8, 2.1, 1.4, 3.2, 4.5),
    inverse.gaussian = c(0.8, 2.1, 1.4, 3.2, 4.5),
    beta = c(0.2, 0.55, 0.4, 0.7, 0.85)
  )
  size <- rep(6, 5)

  for (family in names(family_y)) {
    y <- family_y[[family]]
    for (offset in list(rep(0, 5), c(0, 0.5, -0.3, 0.9, -0.1))) {
      d <- list(
        family = family, X_list = list(X), y_list = list(matrix(y)),
        size = size, offset = matrix(offset), diagtol = 1e-4
      )
      mode <- get_w_and_H(d, 3, list(SigInv), SigInv_X, V, precision_beta, "reml")$w
      score <- get_d(family, drop(mode) + offset, y, size, 3) - Ptheta %*% mode
      expect_equal(max(abs(score)), 0, tolerance = 1e-6)
    }
  }
})

test_that("a constant offset shifts only the fitted Poisson intercept", {
  set.seed(2)
  s <- mf04p
  s$obs$y_pois <- round(s$obs$Summer_mn)
  s$obs$const_offset <- 0.4

  tu <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
  td <- taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range"))
  ng <- nugget_initial("nugget", nugget = 0.5, known = "nugget")

  fit_no_offset <- ssn_glm(y_pois ~ ELEV_DEM, s,
    family = "poisson", tailup_initial = tu, taildown_initial = td, nugget_initial = ng,
    additive = "afvArea"
  )
  fit_offset <- ssn_glm(y_pois ~ ELEV_DEM + offset(const_offset), s,
    family = "poisson", tailup_initial = tu, taildown_initial = td, nugget_initial = ng,
    additive = "afvArea"
  )

  coef_no_offset <- coef(fit_no_offset)
  coef_offset <- coef(fit_offset)
  # a constant offset on the log-link only shifts the intercept by -offset;
  # the slope on ELEV_DEM should be essentially unchanged
  expect_equal(unname(coef_offset["(Intercept)"]), unname(coef_no_offset["(Intercept)"] - 0.4), tolerance = 1e-4)
  expect_equal(unname(coef_offset["ELEV_DEM"]), unname(coef_no_offset["ELEV_DEM"]), tolerance = 1e-4)
})

test_that("varying offset combined with a random intercept and shuffled rows converges", {
  set.seed(2)
  s <- mf04p
  s$obs <- s$obs[sample(nrow(s$obs)), ]
  s$obs$y_pois <- round(s$obs$Summer_mn)
  s$obs$audit_offset <- seq_len(nrow(s$obs)) / 20
  s$obs$audit_group <- factor(rep(1:4, length.out = nrow(s$obs)))

  fit <- ssn_glm(y_pois ~ ELEV_DEM + offset(audit_offset), s,
    family = "poisson", tailup_type = "exponential", nugget_type = "nugget",
    additive = "afvArea", random = ~audit_group
  )
  expect_s3_class(fit, "ssn_glm")
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("varying offset and varying binomial totals fit together without error", {
  set.seed(2)
  s <- mf04p
  s$obs$success <- rep(1:3, length.out = nrow(s$obs))
  s$obs$failure <- seq_len(nrow(s$obs)) + 4
  s$obs$audit_offset <- seq_len(nrow(s$obs)) / 30

  fit <- ssn_glm(cbind(success, failure) ~ ELEV_DEM + offset(audit_offset), s,
    family = "binomial", tailup_type = "exponential", nugget_type = "nugget",
    additive = "afvArea"
  )
  expect_s3_class(fit, "ssn_glm")
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("get_vcov_glm() warns about the nugget (not ie), matching the LM instability warning's terminology", {
  uncorrected <- diag(c(1, 1))
  corrected_unstable <- diag(c(-0.5, 1))
  expect_warning(
    result <- get_vcov_glm(corrected_unstable, uncorrected),
    "Consider fixing nugget \\(via nugget_initial\\)"
  )
  expect_equal(result$fixed$corrected, corrected_unstable)
  expect_equal(result$fixed$uncorrected, uncorrected)

  corrected_stable <- diag(c(0.5, 1))
  expect_no_warning(get_vcov_glm(corrected_stable, uncorrected))
})

test_that("ssn_lm() takes the closed-form IID path (no spatial dependence, no random effects) and matches lm() and the general optimizer path exactly", {
  # tailup/taildown/euclid all default to "none", no random effect -- this is
  # a plain OLS-with-normal-errors model in disguise
  fit_reml <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, estmethod = "reml")
  fit_ml <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, estmethod = "ml")
  expect_true(is.na(fit_reml$optim$method)) # confirms the closed-form path ran, not optim()
  expect_true(is.na(fit_ml$optim$method))

  lmod <- lm(Summer_mn ~ ELEV_DEM, data = mf04p$obs)
  n <- nrow(mf04p$obs)
  p <- 2
  sse <- sum(residuals(lmod)^2)

  expect_equal(unname(coef(fit_reml)), unname(coef(lmod)))
  expect_equal(unname(coef(fit_ml)), unname(coef(lmod)))
  expect_equal(fit_reml$coefficients$params_object$nugget[["nugget"]], sse / (n - p))
  expect_equal(fit_ml$coefficients$params_object$nugget[["nugget"]], sse / n)

  # cross-check against the general (numerically-optimized) path by forcing
  # euclid_type = "none" explicitly instead of relying on the default, and
  # comparing the closed-form path's result to what a direct log-likelihood
  # evaluation at the same nugget value gives via the general machinery
  known_reml <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, estmethod = "reml",
    nugget_initial = nugget_initial("nugget", nugget = sse / (n - p), known = "nugget")
  )
  expect_equal(as.numeric(logLik(fit_reml)), as.numeric(logLik(known_reml)), tolerance = 1e-8)

  known_ml <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, estmethod = "ml",
    nugget_initial = nugget_initial("nugget", nugget = sse / n, known = "nugget")
  )
  expect_equal(as.numeric(logLik(fit_ml)), as.numeric(logLik(known_ml)), tolerance = 1e-8)
  # AIC's parameter count correctly differs by one (the estimated nugget)
  # between the two, even though the log-likelihoods above match exactly
  expect_equal(AIC(fit_ml), AIC(known_ml) + 2)

  # the fixed/known-nugget case is unaffected (still goes through
  # use_gloglik_known(), not the new IID shortcut)
  fixed_fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    nugget_initial = nugget_initial("nugget", nugget = 2, known = "nugget")
  )
  expect_equal(fixed_fit$coefficients$params_object$nugget[["nugget"]], 2)

  # a genuinely spatial fit does not take the IID shortcut
  spatial_fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  expect_false(is.na(spatial_fit$optim$method))
})

test_that("ssn_lm() dispatches the IID fit's model statistics through the closed-form shortcut, avoiding the general dense-matrix path", {
  # target the general model-statistics functions themselves rather than a
  # lower-level helper like get_cov_matrix_list(): that helper is also used,
  # unrelated to this dispatch, by numerical Satterthwaite inference's own
  # gradient computation, which would otherwise trigger this mock regardless
  # of which model-statistics path actually ran
  broken_general <- function(...) stop("general model-statistics path was requested")
  testthat::local_mocked_bindings(
    get_model_stats = broken_general,
    get_model_stats_bigdata = broken_general,
    .package = "SSN2"
  )

  # exact (non-bigdata) dispatch
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, estmethod = "reml")
  expect_true(is.numeric(coef(fit)))

  # grouped/local dispatch (SSN2's local/grouped estimation shares the same
  # backend as its big-data path): IID eligibility does not depend on
  # grouping, since a truly diagonal covariance has no cross-group terms to
  # approximate away in the first place
  group_index <- rep(c(1, 2, 3), length.out = nrow(mf04p$obs))
  fit_local <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, estmethod = "reml", local = list(index = group_index))
  expect_true(is.numeric(coef(fit_local)))

  # a genuinely spatial fit still requests the general path -- the mock
  # above must actually trigger here to prove the eligibility check is not
  # simply always bypassing the general dispatch
  expect_error(
    ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea"),
    "general model-statistics path was requested"
  )

  # a random effect makes the fit ineligible even with no spatial dependence
  expect_error(
    ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, random = ~ as.factor(netID)),
    "general model-statistics path was requested"
  )
})

test_that("the IID model-statistics shortcut matches the general (known-nugget) path exactly, with reordered rows", {
  s <- mf04p
  s$obs <- s$obs[sample(nrow(s$obs)), ]

  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s, estmethod = "reml")
  nugget_val <- fit$coefficients$params_object$nugget[["nugget"]]

  # a fixed/known nugget in this all-"none" setting is IID-eligible too, so
  # without forcing it through the general dispatch here, known_fit would
  # silently take the same shortcut as fit and this comparison would check
  # the shortcut against itself rather than against the general statistics
  # implementation
  force_general <- function(cov_est_object, data_object, estmethod) {
    get_model_stats(cov_est_object, data_object, estmethod)
  }
  testthat::local_mocked_bindings(get_model_stats_iid = force_general, .package = "SSN2")
  known_fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s, estmethod = "reml",
    nugget_initial = nugget_initial("nugget", nugget = nugget_val, known = "nugget")
  )

  expect_equal(unname(coef(fit)), unname(coef(known_fit)))
  expect_equal(fit$vcov$fixed, known_fit$vcov$fixed, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(fit$fitted$response, known_fit$fitted$response, tolerance = 1e-10)
  expect_equal(fit$hatvalues, known_fit$hatvalues, tolerance = 1e-10)
  expect_equal(fit$residuals$response, known_fit$residuals$response, tolerance = 1e-10)
  expect_equal(fit$residuals$standardized, known_fit$residuals$standardized, tolerance = 1e-8)
  expect_equal(fit$cooks_distance, known_fit$cooks_distance, tolerance = 1e-8)
  expect_equal(fit$pseudoR2, known_fit$pseudoR2, tolerance = 1e-10)
  # pid-based row identities survive the shuffled input order
  expect_equal(names(fit$fitted$response), names(known_fit$fitted$response))

  lmod <- lm(Summer_mn ~ ELEV_DEM, data = s$obs)
  expect_equal(unname(fit$fitted$response), unname(fitted(lmod)), tolerance = 1e-8)
})

test_that("the IID model-statistics shortcut correctly incorporates an offset", {
  set.seed(2)
  s <- mf04p
  s$obs$off <- rnorm(nrow(s$obs))

  fit <- ssn_lm(Summer_mn ~ ELEV_DEM + offset(off), s, estmethod = "reml")
  nugget_val <- fit$coefficients$params_object$nugget[["nugget"]]

  # see the previous test: force the reference fit through the general
  # dispatch, since a fixed nugget here is IID-eligible too
  force_general <- function(cov_est_object, data_object, estmethod) {
    get_model_stats(cov_est_object, data_object, estmethod)
  }
  testthat::local_mocked_bindings(get_model_stats_iid = force_general, .package = "SSN2")
  known_fit <- ssn_lm(Summer_mn ~ ELEV_DEM + offset(off), s, estmethod = "reml",
    nugget_initial = nugget_initial("nugget", nugget = nugget_val, known = "nugget")
  )

  expect_equal(fit$fitted$response, known_fit$fitted$response, tolerance = 1e-10)
  expect_equal(fit$residuals$response, known_fit$residuals$response, tolerance = 1e-10)

  lmod <- lm(Summer_mn ~ ELEV_DEM + offset(off), data = s$obs)
  expect_equal(unname(fit$fitted$response), unname(fitted(lmod)), tolerance = 1e-8)
})

test_that("the IID model-statistics shortcut matches the forced-general path across ML/REML, exact/grouped-local, and known/estimated nugget", {
  # Cross estimation method, grouping, and nugget status against the general path.
  set.seed(2)
  s <- mf04p
  s$obs <- s$obs[sample(nrow(s$obs)), ]
  s$obs$off <- rnorm(nrow(s$obs))
  s$obs$Summer_mn[3] <- NA
  group_index <- rep(c(1, 2, 3), length.out = nrow(s$obs))

  fit_shortcut <- function(estmethod, grouped, nugget_val = NULL) {
    local_arg <- if (grouped) list(index = group_index) else NULL
    extra_args <- if (is.null(nugget_val)) {
      list()
    } else {
      list(nugget_initial = nugget_initial("nugget", nugget = nugget_val, known = "nugget"))
    }
    do.call(ssn_lm, c(
      list(Summer_mn ~ ELEV_DEM + offset(off), s, estmethod = estmethod, local = local_arg),
      extra_args
    ))
  }

  # Scope the mock to one fit so later iterations still exercise the shortcut.
  fit_forced_general <- function(estmethod, grouped, nugget_val = NULL) {
    if (grouped) {
      force_general <- function(cov_est_object, data_object, estmethod) {
        get_model_stats_bigdata(cov_est_object, data_object, estmethod)
      }
      testthat::local_mocked_bindings(get_model_stats_bigdata_iid = force_general, .package = "SSN2")
    } else {
      force_general <- function(cov_est_object, data_object, estmethod) {
        get_model_stats(cov_est_object, data_object, estmethod)
      }
      testthat::local_mocked_bindings(get_model_stats_iid = force_general, .package = "SSN2")
    }
    fit_shortcut(estmethod, grouped, nugget_val)
  }

  for (estmethod in c("ml", "reml")) {
    for (grouped in c(FALSE, TRUE)) {
      label <- sprintf("estmethod=%s grouped=%s", estmethod, grouped)
      shortcut_fit <- fit_shortcut(estmethod, grouped)
      nugget_val <- shortcut_fit$coefficients$params_object$nugget[["nugget"]]

      est_ref <- fit_forced_general(estmethod, grouped)
      known_ref <- fit_forced_general(estmethod, grouped, nugget_val)

      expect_equal(unname(coef(shortcut_fit)), unname(coef(est_ref)), tolerance = 1e-7, label = label)
      expect_equal(unname(coef(shortcut_fit)), unname(coef(known_ref)), tolerance = 1e-7, label = label)
      expect_equal(shortcut_fit$fitted$response, est_ref$fitted$response, tolerance = 1e-7, label = label)
      expect_equal(as.numeric(logLik(shortcut_fit)), as.numeric(logLik(known_ref)), tolerance = 1e-7, label = label)
    }
  }
})

# General exact/local Gaussian and GLM statistics
test_that("get_model_stats()/get_model_stats_bigdata() and their GLM counterparts all dispatch through the shared general-statistics cores", {
  gaussian_core_calls <- 0
  real_gaussian_core <- get_model_stats_core
  testthat::local_mocked_bindings(
    get_model_stats_core = function(...) {
      gaussian_core_calls <<- gaussian_core_calls + 1
      real_gaussian_core(...)
    },
    .package = "SSN2"
  )
  ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  expect_equal(gaussian_core_calls, 1)

  group_index <- rep(c(1, 2, 3), length.out = nrow(mf04p$obs))
  gaussian_core_calls <- 0
  ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea", local = list(index = group_index)
  )
  expect_equal(gaussian_core_calls, 1)

  glm_core_calls <- 0
  real_glm_core <- get_model_stats_glm_core
  testthat::local_mocked_bindings(
    get_model_stats_glm_core = function(...) {
      glm_core_calls <<- glm_core_calls + 1
      real_glm_core(...)
    },
    .package = "SSN2"
  )
  ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, family = "Gamma", tailup_type = "exponential", additive = "afvArea")
  expect_equal(glm_core_calls, 1)

  glm_core_calls <- 0
  ssn_glm(Summer_mn ~ ELEV_DEM, mf04p,
    family = "Gamma", tailup_type = "exponential", additive = "afvArea", local = list(index = group_index)
  )
  expect_equal(glm_core_calls, 1)
})

test_that("the general Gaussian statistics core matches an independent fixed-covariance GLS reference", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_initial = tailup_initial("exponential", de = 2, range = 8000, known = "given"),
    nugget_initial = nugget_initial("nugget", nugget = 1, known = "given"),
    additive = "afvArea", estmethod = "reml"
  )

  # independent GLS computation via plain solve()/crossprod() against the
  # fitted model's own (previously validated) observed covariance matrix --
  # a completely different algorithmic route than get_eigenprods()'s
  # eigendecomposition-based whitening, so this exercises the statistical
  # numerical results independently of the shared statistics calculation
  V <- as.matrix(covmatrix(fit))
  X <- model.matrix(fit)
  y <- model.response(model.frame(fit))
  Vinv <- solve(V)
  cov_betahat_ref <- solve(crossprod(X, Vinv %*% X))
  betahat_ref <- as.numeric(cov_betahat_ref %*% crossprod(X, Vinv %*% y))

  expect_equal(unname(coef(fit)), betahat_ref, tolerance = 1e-8)
  expect_equal(unname(as.matrix(fit$vcov$fixed)), unname(cov_betahat_ref), tolerance = 1e-8)

  fitted_ref <- as.numeric(X %*% betahat_ref)
  expect_equal(unname(fit$fitted$response), fitted_ref, tolerance = 1e-8)
})

test_that("a single local group reproduces the dense/exact general Gaussian and GLM fit exactly", {
  # covariance parameters are fixed (known), not numerically estimated: an
  # earlier version of this test estimated them, which left the comparison
  # exposed to optim()'s Nelder-Mead landing at very slightly different (but
  # both converged) endpoints for the dense-vs-single-group likelihood
  # evaluation -- a pre-existing property of the (untouched-by-this-change)
  # optimizer path, not of the statistics assembly this test targets. Fixing
  # the covariance parameters removes the optimizer from the comparison
  # entirely, matching the established safe pattern already used for the IID
  # case in test-extras-local.R ("... one local group must exactly reproduce
  # the dense (exact) fit").
  n <- nrow(mf04p$obs)
  one_group <- rep(1, n)
  tu <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
  ng <- nugget_initial("nugget", nugget = 0.5, known = "nugget")

  fit_dense <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_initial = tu, nugget_initial = ng, additive = "afvArea"
  )
  fit_local_one <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_initial = tu, nugget_initial = ng, additive = "afvArea",
    local = list(index = one_group)
  )
  expect_equal(coef(fit_local_one), coef(fit_dense), tolerance = 1e-8)
  expect_equal(fit_local_one$vcov$fixed, fit_dense$vcov$fixed, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(unname(fit_local_one$fitted$response), unname(fit_dense$fitted$response), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(fit_local_one)), as.numeric(logLik(fit_dense)), tolerance = 1e-6)

  fit_glm_dense <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p,
    family = "Gamma", tailup_initial = tu, nugget_initial = ng, additive = "afvArea"
  )
  fit_glm_local_one <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p,
    family = "Gamma", tailup_initial = tu, nugget_initial = ng, additive = "afvArea",
    local = list(index = one_group)
  )
  expect_equal(coef(fit_glm_local_one), coef(fit_glm_dense), tolerance = 1e-6)
  expect_equal(unname(fit_glm_local_one$fitted$response), unname(fit_glm_dense$fitted$response), tolerance = 1e-5)
})

test_that("grouped general Gaussian and GLM statistics agree between serial and parallel dispatch", {
  group_index <- rep(c(1, 2, 3), length.out = nrow(mf04p$obs))

  fit_serial <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    local = list(index = group_index)
  )
  fit_parallel <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    local = list(index = group_index, parallel = TRUE, ncores = 2)
  )
  expect_equal(coef(fit_serial), coef(fit_parallel), tolerance = 1e-10)
  expect_equal(fit_serial$fitted$response, fit_parallel$fitted$response, tolerance = 1e-10)
  expect_equal(fit_serial$vcov$fixed, fit_parallel$vcov$fixed, tolerance = 1e-10, ignore_attr = TRUE)

  fit_glm_serial <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p,
    family = "Gamma", tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    local = list(index = group_index)
  )
  fit_glm_parallel <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p,
    family = "Gamma", tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    local = list(index = group_index, parallel = TRUE, ncores = 2)
  )
  expect_equal(coef(fit_glm_serial), coef(fit_glm_parallel), tolerance = 1e-8)
  expect_equal(fit_glm_serial$fitted$response, fit_glm_parallel$fitted$response, tolerance = 1e-8)
})

test_that("general statistics preserve random effects and partition factors after consolidation", {
  ssn <- mf04p
  ssn$obs$audit_group <- factor(rep(c("a", "b", "c"), length.out = nrow(ssn$obs)))
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    random = ~audit_group, partition_factor = ~audit_group
  )
  expect_true(all(is.finite(coef(fit))))
  expect_false(is.null(fit$fitted$randcov))
  reconstructed <- fitted(fit) + residuals(fit)
  observed <- setNames(ssn$obs$Summer_mn, ssn_get_netgeom(ssn$obs, "pid")$pid)
  expect_equal(unname(reconstructed), unname(observed[names(reconstructed)]), tolerance = 1e-6)
})

