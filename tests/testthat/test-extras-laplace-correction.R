skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"))

laplace_correction_fixture <- local({
  network <- NULL
  fits <- list()
  function(family = "poisson", estmethod = "reml", local = FALSE,
           iid = FALSE, effects = FALSE) {
    key <- paste(family, estmethod, local, iid, effects)
    if (!is.null(fits[[key]])) return(fits[[key]])
    if (is.null(network)) {
      s <- mf04p
      root <- tempfile("laplace-correction-")
      dir.create(root)
      stopifnot(file.copy(s$path, root, recursive = TRUE, copy.mode = FALSE))
      s$path <- file.path(root, basename(s$path))
      rows <- round(seq(1, NROW(s$obs), length.out = 24))
      s$obs <- s$obs[rows[c(seq(2, 24, 2), seq(1, 23, 2))], ]
      s$preds <- list(correction = s$preds$pred1km[c(5, 1, 8, 3), ])
      set.seed(2)
      s$obs$x <- rnorm(24)
      s$obs$off <- seq(-.2, .3, length.out = 24)
      s$obs$trials <- rep(4:7, 6)
      s$obs$g <- factor(rep(letters[1:4], 6))
      s$obs$part <- factor(rep(c("one", "two"), 12))
      s$obs$poisson <- c(1, 3, 2, 4, 2, 5, 4, 2, 5, 3, 6, 4,
                         2, 5, 1, 4, 3, 2, 5, 4, 2, 3, 4, 1)
      s$obs$nbinomial <- s$obs$poisson
      s$obs$binomial <- rep(c(1, 3, 2, 4, 2, 1), 4)
      s$preds$correction$x <- c(-.5, .2, .8, 0)
      s$preds$correction$off <- c(.1, -.3, .2, .4)
      s$preds$correction$g <- factor(c("a", "b", "c", "d"))
      s$preds$correction$part <- factor(c("one", "two", "one", "two"))
      ssn_create_distmat(s, predpts = "correction", among_predpts = TRUE, overwrite = TRUE)
      network <<- s
    }
    s <- network
    s$obs$y <- s$obs[[family]]
    args <- list(formula = if (family == "binomial") cbind(y, trials - y) ~ x + offset(off) else y ~ x + offset(off),
      ssn.object = s, family = family, estmethod = estmethod, additive = "afvArea",
      tailup_initial = if (iid) tailup_initial("none") else tailup_initial("exponential", .1, 10000, known = "given"),
      taildown_initial = if (iid) taildown_initial("none") else taildown_initial("exponential", .1, 10000, known = "given"),
      euclid_initial = if (iid) euclid_initial("none") else euclid_initial("exponential", .2, 12000, known = "given"),
      nugget_initial = nugget_initial("nugget", .25, known = "given"),
      local = if (local) list(index = rep(c(3, 1, 2), 8), var_adjust = "none") else FALSE)
    if (family == "nbinomial") {
      args$family <- NULL
      args$dispersion_initial <- dispersion_initial("nbinomial", 3, known = "given")
    }
    if (effects) {
      args$random <- ~g
      args$randcov_initial <- randcov_initial(g = .1, known = "given")
      args$partition_factor <- ~part
    }
    fit <- do.call(ssn_glm, args)
    fits[[key]] <<- fit
    fit
  }
})

test_that("third derivatives agree with finite differences for all six families", {
  for (family in c("poisson", "nbinomial", "binomial", "Gamma", "inverse.gaussian", "beta")) {
    grid <- expand.grid(w = c(-4, -1, 0, .7, 3),
      y = if (family == "beta") c(.02, .3, .8, .98) else c(.2, 1, 5),
      size = c(1, 8), dispersion = c(.3, 2, 15))
    reference <- (diag(get_D(family, grid$w + 1e-5, grid$y, grid$size, grid$dispersion)) -
      diag(get_D(family, grid$w - 1e-5, grid$y, grid$size, grid$dispersion))) / 2e-5
    actual <- get_D3(family, grid$w, grid$y, grid$size, grid$dispersion)
    expect_type(actual, "double")
    expect_lt(max(abs(actual - reference) / pmax(1, abs(reference))), 1e-5)
  }
})

test_that("Laplace helpers match the pinned spmodel reference", {
  # Pinned spmodel 165926b values keep parity checks self-contained.
  third <- rbind(
    poisson = c(-.496585303791, -1.22140275816, -3.00416602395),
    nbinomial = c(-.348959898938, -.433169032641, .00121424663596),
    binomial = c(-.149157576881, .123347402193, .750259293455),
    Gamma = c(6.04125812241, 4.91238451847, 3.99445300438),
    inverse.gaussian = c(2.27575110552, 1.54014019061, .870664243209),
    beta = c(-.345866404855, .237080978947, .600789608996))
  shift <- rbind(
    poisson = c(-.0225437021564, -.19170484888, -.257164567377),
    nbinomial = c(-.0200397427652, -.0562253175141, -.00535351606192),
    binomial = c(-.00847965411152, .0235265913343, .0623184272399),
    Gamma = c(.305002648676, .698229049595, .37511823376),
    inverse.gaussian = c(.114327544682, .216620001378, .0859279889993),
    beta = c(-.015596336865, .0347693992472, .0526569826185))
  log_mean <- c(.522045776761, 1.49182469764, 4.48168907034)
  V <- matrix(c(.3, .04, -.02, .04, .5, .06, -.02, .06, .4), 3)
  for (family in rownames(third)) {
    y <- if (family == "beta") c(.2, .5, .8) else c(1, 2, 4)
    expect_equal(get_D3(family, c(-.7, .2, 1.1), y, c(2, 5, 8), 3),
      unname(third[family, ]), tolerance = 1e-10)
    expect_equal(get_laplace_correction(V, family, c(-.7, .2, 1.1), y, c(2, 5, 8), 3)$delta,
      unname(shift[family, ]), tolerance = 1e-10)
    expected <- switch(family,
      binomial = c(.669558521636, 2.73174468291, 5.78425739275),
      beta = c(.334779260818, .546348936582, .723032174094), log_mean)
    expect_equal(invlink_mean(c(-.7, .2, 1.1), c(.1, .4, .8), family, c(2, 5, 8)),
      expected, tolerance = 1e-10)
  }
})

test_that("stored mean shifts use projected mode curvature in observation order", {
  for (family in c("poisson", "nbinomial", "binomial")) {
    for (method in c("ml", "reml")) for (local in c(FALSE, TRUE)) {
      fit <- laplace_correction_fixture(family, method, local)
      X <- model.matrix(fit)
      Sigma <- as.matrix(covmatrix(fit))
      if (local) {
        expect_false(identical(order(fit$local_index), seq_len(fit$n)))
        Sigma[outer(fit$local_index, fit$local_index, "!=")] <- 0
      }
      Si <- solve(Sigma)
      B <- solve(crossprod(X, Si %*% X), crossprod(X, Si))
      P <- Si - Si %*% X %*% B
      eta <- fitted(fit, type = "link")
      phi <- as.vector(coef(fit, type = "dispersion"))
      D <- as.matrix(get_D(family, eta, fit$y, fit$size, phi))
      Omega <- solve(P - D)
      third <- get_D3(family, eta, fit$y, fit$size, phi)
      expected <- .5 * as.numeric(Omega %*% (third * diag(Omega)))
      expect_equal(unname(fit$laplace$delta), expected, tolerance = 1e-4)
      expect_identical(names(fit$laplace$delta), names(eta))
      expect_equal(unname(fit$laplace$betahat), as.numeric(B %*%
        (eta - model.offset(model.frame(fit)) + fit$laplace$delta)), tolerance = 1e-9)
      center <- get_latent_center(fit)
      expect_true(center$active)
      expect_equal(center$eta, eta + fit$laplace$delta)
      expect_equal(center$betahat, fit$laplace$betahat)
      reordered <- fit
      reordered$laplace$delta <- rev(reordered$laplace$delta)
      expect_identical(get_latent_center(reordered), center)
      expect_identical(get_latent_center(fit, FALSE)$eta, eta)
      expect_identical(get_latent_center(fit, FALSE)$betahat, coef(fit))
      if (method == "ml" && !local) {
        fixed <- solve(Si - D)
        expect_gt(max(abs(expected - .5 * as.numeric(fixed %*% (third * diag(fixed))))), 1e-4)
      }
    }
  }
})

test_that("inverse-link means use lognormal moments and the logit approximation", {
  grid <- expand.grid(m = c(-4, 0, 3), s = c(0, .3, 1, 2, 3, 10))
  sizes <- rep(c(1, 4, 9), length.out = nrow(grid))
  for (family in c("poisson", "nbinomial", "Gamma", "inverse.gaussian", "binomial", "beta")) {
    value <- invlink_mean(grid$m, grid$s^2, family, sizes)
    if (family %in% c("binomial", "beta")) {
      expected <- plogis(grid$m / sqrt(1 + pi * grid$s^2 / 8))
      reference <- mapply(function(m, s) integrate(function(z) plogis(m + s * z) * dnorm(z),
        -Inf, Inf, rel.tol = 1e-10)$value, grid$m, grid$s)
      expect_lt(max(abs(expected - reference)), .015)
      expect_equal(value, expected * if (family == "binomial") sizes else 1)
    } else {
      expect_identical(value, exp(grid$m + grid$s^2 / 2))
    }
    expect_equal(invlink_mean(grid$m, 0, family, sizes), invlink(grid$m, family, sizes))
  }
})

test_that("prediction changes centers but preserves mode curvature and legacy behavior", {
  for (family in c("poisson", "binomial")) for (iid in c(FALSE, TRUE)) {
    fit <- laplace_correction_fixture(family, iid = iid, effects = !iid)
    original <- fit
    old <- fit
    old$laplace <- NULL
    sizes <- c(2, 5, 8, 3)
    for (local in list(FALSE, list(method = "covariance", size = 10))) {
      args <- list(object = fit, newdata = "correction", local = local, newdata_size = sizes)
      for (correct_var in c(FALSE, TRUE)) {
        args$var_correct <- correct_var
        on <- do.call(predict, c(args, list(se.fit = TRUE)))
        off <- do.call(predict, c(args, list(se.fit = TRUE, fit_correct = FALSE)))
        expect_identical(on$se.fit, off$se.fit)
        expect_gt(max(abs(on$fit - off$fit)), 0)
        old_args <- args
        old_args$object <- old
        expect_identical(do.call(predict, c(old_args, list(se.fit = TRUE))), off)
        response <- do.call(predict, c(args, list(type = "response", se.fit = TRUE, delta = TRUE)))
        expect_equal(response$fit, invlink_mean(on$fit, on$se.fit^2, family, sizes))
        expect_equal(response$se.fit, get_delta_se(on$fit, on$se.fit, family, sizes))
        expect_identical(do.call(predict, c(args, list(type = "response"))), response$fit)
        expect_identical(do.call(predict, c(old_args, list(type = "response"))),
          do.call(predict, c(args, list(type = "response", fit_correct = FALSE))))
        link_interval <- do.call(predict, c(args, list(interval = "prediction", se.fit = TRUE)))
        response_interval <- do.call(predict, c(args, list(type = "response", interval = "prediction")))
        expect_equal(response_interval[, "fit"], response$fit, ignore_attr = TRUE)
        for (bound in c("lwr", "upr")) {
          expect_equal(response_interval[, bound], invlink(link_interval$fit[, bound], family, sizes))
        }
      }
      expect_identical(do.call(predict, c(args, list(type = "terms", se.fit = TRUE))),
        do.call(predict, c(args, list(type = "terms", se.fit = TRUE, fit_correct = FALSE))))
      weights <- do.call(predict, c(args, list(type = "weight")))
      expect_identical(weights, do.call(predict, c(args, list(type = "weight", fit_correct = FALSE))))
      if (identical(local, FALSE)) {
        expect_equal(as.numeric(on$fit - off$fit), as.numeric(weights %*% fit$laplace$delta), tolerance = 1e-9)
      }
      confidence <- do.call(predict, c(args, list(interval = "confidence", se.fit = TRUE)))
      response_confidence <- do.call(predict, c(args, list(interval = "confidence", type = "response", se.fit = TRUE)))
      Xnew <- get_newdata_model_matrix(fit, fit$ssn.object$preds$correction)$newdata_model
      expected <- as.numeric(Xnew %*% fit$laplace$betahat) + fit$ssn.object$preds$correction$off
      expect_equal(as.numeric(confidence$fit[, "fit"]), expected)
      expect_equal(response_confidence$fit[, "fit"], invlink(confidence$fit[, "fit"], family, sizes))
      expect_identical(response_confidence$se.fit, confidence$se.fit)
    }
    expect_equal(augment(fit, newdata = "correction", type.predict = "response", newdata_size = sizes)$.fitted,
      as.numeric(predict(fit, "correction", type = "response", newdata_size = sizes)))
    expect_identical(fit, original)
    for (bad in list(NA, 1, c(TRUE, FALSE), NULL, "TRUE")) {
      expect_error(predict(fit, "correction", fit_correct = bad), "fit_correct")
      expect_error(loocv(fit, fit_correct = bad), "fit_correct")
      expect_error(kcv(fit, fit_correct = bad), "fit_correct")
    }
  }
})

test_that("nonlocal CV shifts propagate through held-out kriging weights", {
  for (family in c("poisson", "binomial")) {
    fit <- laplace_correction_fixture(family)
    X <- model.matrix(fit)
    S <- as.matrix(covmatrix(fit))
    offset <- model.offset(model.frame(fit))
    w <- fitted(fit, "link") - offset
    old <- fit
    old$laplace <- NULL
    for (folds in list(seq_len(fit$n), rep(1:3, length.out = fit$n))) {
      expected <- shift <- numeric(fit$n)
      for (F in split(seq_len(fit$n), folds)) {
        Q <- solve(S[-F, -F, drop = FALSE])
        Xm <- X[-F, , drop = FALSE]
        B <- solve(crossprod(Xm, Q %*% Xm), crossprod(Xm, Q))
        L <- X[F, , drop = FALSE] %*% B + S[F, -F, drop = FALSE] %*% Q %*% (diag(nrow(Xm)) - Xm %*% B)
        expected[F] <- as.numeric(L %*% (w[-F] + fit$laplace$delta[-F])) + offset[F]
        shift[F] <- as.numeric(L %*% fit$laplace$delta[-F])
      }
      args <- list(object = fit, folds_index = folds, cv_predict = TRUE, se.fit = TRUE, local = FALSE)
      on <- do.call(kcv, args)
      off <- do.call(kcv, c(args, list(fit_correct = FALSE)))
      expect_equal(as.numeric(on$cv_predict), expected, tolerance = 1e-9)
      expect_equal(as.numeric(on$cv_predict - off$cv_predict), shift, tolerance = 1e-9)
      expect_identical(on$se.fit, off$se.fit)
      args$object <- old
      expect_identical(do.call(kcv, args), off)
      if (length(unique(folds)) == fit$n) for (correction in c(FALSE, TRUE)) {
        loo <- loocv(fit, cv_predict = TRUE, se.fit = TRUE, fit_correct = correction)
        expect_equal(loo, if (correction) on else off, tolerance = 1e-9)
        expect_identical(kcv(fit, k = fit$n, cv_predict = TRUE, se.fit = TRUE, fit_correct = correction), loo)
      }
    }
  }
})

test_that("local and nonlocal response CV preserve trial sizes and variances not requested", {
  for (family in c("poisson", "binomial")) {
    fit <- laplace_correction_fixture(family, local = TRUE, effects = TRUE)
    old <- fit
    old$laplace <- NULL
    for (local in list(FALSE, list(method = "covariance", size = 8))) for (fun in list(loocv, kcv)) {
      args <- list(object = fit, cv_predict = TRUE, local = local)
      if (identical(fun, kcv)) args$folds_index <- rep(1:3, length.out = fit$n)
      link <- do.call(fun, c(args, list(se.fit = TRUE)))
      off <- do.call(fun, c(args, list(se.fit = TRUE, fit_correct = FALSE)))
      expect_identical(link$se.fit, off$se.fit)
      response <- do.call(fun, c(args, list(type = "response", se.fit = TRUE, delta = TRUE)))
      expected <- invlink_mean(link$cv_predict, link$se.fit^2, family, fit$size)
      expect_equal(as.numeric(response$cv_predict), as.numeric(expected))
      expect_equal(response$se.fit, get_delta_se(link$cv_predict, link$se.fit, family, fit$size))
      no_se <- do.call(fun, c(args, list(type = "response")))
      expect_identical(no_se$cv_predict, response$cv_predict)
      expect_false("se.fit" %in% names(no_se))
      expect_identical(no_se$stats, response$stats)
      expect_equal(response$stats$bias, mean(fit$y - expected))
      expect_equal(response$stats$MSPE, mean((fit$y - expected)^2))
      expect_equal(response$stats$RMSPE, sqrt(mean((fit$y - expected)^2)))
      args$object <- old
      expect_identical(do.call(fun, c(args, list(se.fit = TRUE))), off)
    }
  }
})
