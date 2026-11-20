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

test_that("GLM prediction uses mode-based link and response values", {
  for (family in c("poisson", "binomial")) for (iid in c(FALSE, TRUE)) {
    fit <- laplace_correction_fixture(family, iid = iid, effects = !iid)
    original <- fit
    sizes <- c(2, 5, 8, 3)
    for (local in list(FALSE, list(method = "covariance", size = 10))) {
      args <- list(object = fit, newdata = "correction", local = local, newdata_size = sizes)
      for (correct_var in c(FALSE, TRUE)) {
        args$var_correct <- correct_var
        link <- do.call(predict, c(args, list(se.fit = TRUE)))
        response <- do.call(predict, c(args, list(type = "response", se.fit = TRUE, delta = TRUE)))
        expect_equal(response$fit, invlink(link$fit, family, sizes))
        expect_equal(response$se.fit, get_delta_se(link$fit, link$se.fit, family, sizes))
        expect_identical(do.call(predict, c(args, list(type = "response"))), response$fit)
        link_interval <- do.call(predict, c(args, list(interval = "prediction", se.fit = TRUE)))
        response_interval <- do.call(predict, c(args, list(type = "response", interval = "prediction")))
        expect_equal(response_interval[, "fit"], response$fit, ignore_attr = TRUE)
        for (bound in c("lwr", "upr")) {
          expect_equal(response_interval[, bound], invlink(link_interval$fit[, bound], family, sizes))
        }
      }
      confidence <- do.call(predict, c(args, list(interval = "confidence", se.fit = TRUE)))
      response_confidence <- do.call(predict, c(args, list(interval = "confidence", type = "response", se.fit = TRUE)))
      Xnew <- get_newdata_model_matrix(fit, fit$ssn.object$preds$correction)$newdata_model
      expected <- as.numeric(Xnew %*% coef(fit)) + fit$ssn.object$preds$correction$off
      expect_equal(as.numeric(confidence$fit[, "fit"]), expected)
      expect_equal(response_confidence$fit[, "fit"], invlink(confidence$fit[, "fit"], family, sizes))
      expect_identical(response_confidence$se.fit, confidence$se.fit)
    }
    expect_equal(augment(fit, newdata = "correction", type.predict = "response", newdata_size = sizes)$.fitted,
      as.numeric(predict(fit, "correction", type = "response", newdata_size = sizes)))
    expect_identical(fit, original)
  }
})

test_that("nonlocal GLM CV agrees with held-out kriging weights", {
  for (family in c("poisson", "binomial")) {
    fit <- laplace_correction_fixture(family)
    X <- model.matrix(fit)
    S <- as.matrix(covmatrix(fit))
    offset <- model.offset(model.frame(fit))
    w <- fitted(fit, "link") - offset
    for (folds in list(seq_len(fit$n), rep(1:3, length.out = fit$n))) {
      expected <- numeric(fit$n)
      for (F in split(seq_len(fit$n), folds)) {
        Q <- solve(S[-F, -F, drop = FALSE])
        Xm <- X[-F, , drop = FALSE]
        B <- solve(crossprod(Xm, Q %*% Xm), crossprod(Xm, Q))
        L <- X[F, , drop = FALSE] %*% B + S[F, -F, drop = FALSE] %*% Q %*% (diag(nrow(Xm)) - Xm %*% B)
        expected[F] <- as.numeric(L %*% (w[-F])) + offset[F]
      }
      args <- list(object = fit, folds_index = folds, cv_predict = TRUE, se.fit = TRUE, local = FALSE)
      cv <- do.call(kcv, args)
      expect_equal(as.numeric(cv$cv_predict), expected, tolerance = 1e-9)
      if (length(unique(folds)) == fit$n) {
        loo <- loocv(fit, cv_predict = TRUE, se.fit = TRUE)
        expect_equal(loo, cv, tolerance = 1e-9)
        expect_identical(kcv(fit, k = fit$n, cv_predict = TRUE, se.fit = TRUE), loo)
      }
    }
  }
})

test_that("local and nonlocal response CV preserve trial sizes and variances not requested", {
  for (family in c("poisson", "binomial")) {
    fit <- laplace_correction_fixture(family, local = TRUE, effects = TRUE)
    for (local in list(FALSE, list(method = "covariance", size = 8))) for (fun in list(loocv, kcv)) {
      args <- list(object = fit, cv_predict = TRUE, local = local)
      if (identical(fun, kcv)) args$folds_index <- rep(1:3, length.out = fit$n)
      link <- do.call(fun, c(args, list(se.fit = TRUE)))
      response <- do.call(fun, c(args, list(type = "response", se.fit = TRUE, delta = TRUE)))
      expected <- invlink(link$cv_predict, family, fit$size)
      expect_equal(as.numeric(response$cv_predict), as.numeric(expected))
      expect_equal(response$se.fit, get_delta_se(link$cv_predict, link$se.fit, family, fit$size))
      no_se <- do.call(fun, c(args, list(type = "response")))
      expect_identical(no_se$cv_predict, response$cv_predict)
      expect_false("se.fit" %in% names(no_se))
      expect_identical(no_se$stats, response$stats)
      expect_equal(response$stats$bias, mean(fit$y - expected))
      expect_equal(response$stats$MSPE, mean((fit$y - expected)^2))
      expect_equal(response$stats$RMSPE, sqrt(mean((fit$y - expected)^2)))
    }
  }
})
