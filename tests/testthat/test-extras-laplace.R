skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"))

laplace_ssn_data <- function(n = 36) {
  set.seed(2)
  s <- mf04p
  s$obs$.keep <- seq_len(nrow(s$obs)) <= n
  s <- ssn_subset(s, tempfile(fileext = ".ssn"), .keep)
  ssn_create_distmat(s)
  s$obs$x <- rnorm(n)
  s$obs$off <- runif(n, -.2, .2)
  s$obs$g <- factor(rep(1:6, length.out = n))
  s$obs$part <- factor(rep(1:2, length.out = n))
  eta <- .5 + .25 * s$obs$x + s$obs$off + rnorm(n, sd = .6)
  s$obs$poisson <- rpois(n, exp(eta))
  s$obs$binomial <- rbinom(n, 8, plogis(eta))
  s$obs$nbinomial <- rnbinom(n, mu = exp(eta), size = 4)
  s$obs$Gamma <- rgamma(n, shape = 4, rate = 4 / exp(eta))
  s$obs$inverse.gaussian <- exp(eta) * rep(c(.7, .9, 1.1, 1.3), length.out = n)
  s$obs$beta <- .25 + .5 * plogis(eta)
  s
}

test_that("Laplace likelihoods reproduce exact Gaussian ML and REML including constants", {
  n <- 12
  X <- cbind(1, seq(-1, 1, length.out = n))
  y <- sin(seq_len(n)) + X %*% c(.3, .7)
  offset <- seq(-.2, .2, length.out = n)^2
  logdet <- function(x) as.numeric(determinant(x, logarithm = TRUE)$modulus)
  for (blocks in c(1, 3)) {
    index <- rep(seq_len(blocks), each = n / blocks)
    rows <- split(seq_len(n), index)
    for (scale in c(.2, 2)) {
      Sigma <- scale * (.6^abs(outer(seq_len(n), seq_len(n), "-")) + diag(.3, n))
      Sigma[outer(index, index, "!=")] <- 0
      covs <- lapply(rows, function(i) Sigma[i, i, drop = FALSE])
      for (noise in c(.1, 1, 3)) {
        obj <- list(family = "gaussian", n = n, p = ncol(X),
          X_list = lapply(rows, function(i) X[i, , drop = FALSE]),
          y_list = lapply(rows, function(i) y[i, , drop = FALSE]), offset = offset,
          local_index = index, diagtol = 0, parallel = FALSE)
        V <- Sigma + diag(noise, n)
        C <- solve(crossprod(X, solve(V, X)))
        beta <- C %*% crossprod(X, solve(V, y - offset))
        residual <- y - offset - X %*% beta
        ml <- n * log(2 * pi) + logdet(V) + as.numeric(crossprod(residual, solve(V, residual)))
        for (method in c("ml", "reml")) {
          products <- with_mocked_bindings(laploglik_products(list(dispersion = noise), obj, method),
            get_cov_matrix_list = function(...) covs,
            get_w_init = function(family, y, dispersion) rep(0, length(y)),
            get_d = function(family, w, y, size, dispersion) (y - w) / dispersion,
            get_D = function(family, w, y, size, dispersion) Matrix::Diagonal(n, -1 / dispersion),
            get_l00 = function(family, w, y, size, dispersion) -2 * sum(dnorm(as.numeric(y), as.numeric(w), sqrt(dispersion), log = TRUE)))
          expected <- ml + if (method == "reml") logdet(solve(C)) - ncol(X) * log(2 * pi) else 0
          expect_equal(get_minustwolaploglik(products, obj, method), expected,
            tolerance = 1e-10, info = paste(blocks, scale, noise, method))
        }
      }
    }
  }
})

test_that("independent-error ML approaches ordinary GLM likelihoods as variance shrinks", {
  s <- laplace_ssn_data()
  set.seed(2)
  eta <- .3 + .4 * s$obs$x + s$obs$off
  s$obs$poisson <- rpois(nrow(s$obs), exp(eta))
  s$obs$binomial <- rbinom(nrow(s$obs), 8, plogis(eta))
  for (family in c("poisson", "binomial")) {
    s$obs$response <- s$obs[[family]]
    formula <- if (family == "poisson") response ~ x + offset(off) else cbind(response, 8 - response) ~ x + offset(off)
    ordinary <- glm(formula, family = family, data = s$obs)
    gaps <- numeric()
    for (v in c(.01, .001, .0001)) {
      fit <- ssn_glm(formula, s, family = as.character(family), estmethod = "ml", tailup_type = "none",
        nugget_initial = nugget_initial("nugget", v, known = "given"))
      gaps <- c(gaps, abs(-2 * as.numeric(logLik(fit) - logLik(ordinary))))
      expect_lt(tail(gaps, 1), nrow(s$obs) * v)
      expect_equal(as.numeric(AIC(fit) - AIC(ordinary)), -2 * as.numeric(logLik(fit) - logLik(ordinary)), tolerance = 1e-9)
    }
    expect_true(all(diff(gaps) < 0))
    none <- ssn_glm(formula, s, family = as.character(family), estmethod = "ml", tailup_type = "none", nugget_type = "none")
    expect_equal(as.numeric(logLik(none)), as.numeric(logLik(fit)), tolerance = 1e-10)
  }
})

test_that("Laplace determinants preserve fixed-parameter inference across relevant covariance paths", {
  s <- laplace_ssn_data(24)
  original <- get_w_and_H
  original_D <- get_D
  logdet <- function(x) as.numeric(determinant(x, logarithm = TRUE)$modulus)
  for (family in c("poisson", "binomial", "nbinomial", "Gamma", "inverse.gaussian", "beta")) {
    s$obs$response <- s$obs[[family]]
    formula <- if (family == "binomial") cbind(response, 8 - response) ~ x + offset(off) else response ~ x + offset(off)
    for (path in c("network", "blocks", "random", "anisotropic", "partition")) {
      args <- list(formula = formula, ssn.object = s, family = family, additive = "afvArea",
        tailup_initial = tailup_initial("exponential", de = .1, range = 1000, known = "given"),
        taildown_initial = taildown_initial("exponential", de = .1, range = 1000, known = "given"),
        euclid_initial = euclid_initial("exponential", de = .1, range = 1000,
          rotate = if (path == "anisotropic") .3 else 0,
          scale = if (path == "anisotropic") .7 else 1, known = "given"),
        nugget_initial = nugget_initial("nugget", .2, known = "given"),
        local = if (path == "blocks") list(index = rep(1:3, each = 8)) else FALSE)
      if (path == "random") {
        args$random <- ~g
        args$randcov_initial <- randcov_initial(g = .1, known = "given")
      }
      if (path == "anisotropic") args$anisotropy <- TRUE
      if (path == "partition") args$partition_factor <- ~part
      if (!family %in% c("poisson", "binomial")) {
        args$dispersion_initial <- dispersion_initial(as.character(family), 4, known = "given")
        args$family <- NULL
      }
      fits <- list()
      for (method in c("ml", "reml")) {
        args$estmethod <- method
        checked <- FALSE
        fits[[method]] <- with_mocked_bindings(do.call(ssn_glm, args),
          get_w_and_H = function(data_object, dispersion, SigInv_list, SigInv_X, cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = FALSE) {
            D <- NULL
            result <- with_mocked_bindings(original(data_object, dispersion, SigInv_list,
              SigInv_X, cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = TRUE),
              get_D = function(...) { D <<- original_D(...); D })
            if (!checked) {
              checked <<- TRUE
              expected <- if (method == "ml") logdet(Matrix::bdiag(SigInv_list) - D) else -logdet(result$mHInv)
              if (length(SigInv_list) > 1 && method == "reml") {
                diag(cov_betahat_Inv) <- diag(cov_betahat_Inv) + data_object$diagtol
                expected <- expected + logdet(cov_betahat) + logdet(cov_betahat_Inv)
              }
              expect_equal(result$mHldet, expected, tolerance = 1e-9, info = paste(family, path, method))
            }
            result
          })
        expect_true(checked)
      }
      expect_identical(coef(fits$ml), coef(fits$reml))
      expect_identical(fitted(fits$ml, type = "link"), fitted(fits$reml, type = "link"))
      expect_identical(vcov(fits$ml), vcov(fits$reml))
      expect_identical(vcov(fits$ml, var_correct = FALSE), vcov(fits$reml, var_correct = FALSE))
      set.seed(2)
      ml_draws <- conditional(fits$ml, output = c("latent", "beta"), samples = 5)
      set.seed(2)
      reml_draws <- conditional(fits$reml, output = c("latent", "beta"), samples = 5)
      expect_identical(ml_draws, reml_draws)
    }
  }
})

test_that("estimated nuggets recover independently bracketed fixed-variance optimums", {
  set.seed(2)
  s <- mf04p
  s$obs <- s$obs[seq_len(40), ]
  s$obs$x <- rnorm(40)
  s$obs$y <- rpois(40, exp(.5 + .3 * s$obs$x + rnorm(40, sd = .8)))
  for (method in c("reml", "ml")) {
    fixed <- function(theta) ssn_glm(y ~ x, s, family = "poisson", tailup_type = "none",
      estmethod = method, nugget_initial = nugget_initial("nugget", exp(theta), known = "given"))
    reference <- optimize(function(theta) -2 * as.numeric(logLik(fixed(theta))), log(c(.01, 10)), tol = 1e-7)
    for (start in c(1e-8, .3, 10)) {
      fit <- ssn_glm(y ~ x, s, family = "poisson", tailup_type = "none", estmethod = method,
        nugget_initial = nugget_initial("nugget", start))
      nugget <- coef(fit, type = "nugget")[["nugget"]]
      ref <- fixed(log(nugget))
      expect_equal(nugget, exp(reference$minimum), tolerance = 1e-4)
      expect_equal(fit$optim$value, reference$objective, tolerance = 1e-7)
      expect_equal(fit$optim$value, -2 * as.numeric(logLik(ref)), tolerance = 1e-10)
      expect_equal(coef(fit), coef(ref), tolerance = 1e-10)
      expect_equal(vcov(fit), vcov(ref), tolerance = 1e-10)
      expect_lt(fit$optim$value, -2 * as.numeric(logLik(fixed(log(1e-4)))))
      expect_false(fit$is_known$nugget[["nugget"]])
    }
  }
})

test_that("sole-nugget GLM searches reconstruct fixed fits across families and contexts", {
  s <- laplace_ssn_data()
  for (family in c("poisson", "binomial", "nbinomial", "Gamma", "inverse.gaussian", "beta")) {
    s$obs$response <- s$obs[[family]]
    formula <- if (family == "binomial") cbind(response, 8 - response) ~ x + offset(off) else response ~ x + offset(off)
    paths <- if (family == "poisson") c("iid", "network", "blocks", "random", "partition", "anisotropic") else "iid"
    for (path in paths) for (method in c("ml", "reml")) {
      args <- list(formula = formula, ssn.object = s, family = family, estmethod = method,
        tailup_type = "none", nugget_initial = nugget_initial("nugget", .3))
      if (!family %in% c("poisson", "binomial")) {
        args$dispersion_initial <- dispersion_initial(as.character(family), 4, known = "given")
        args$family <- NULL
      }
      if (path == "blocks") args$local <- list(index = rep(1:3, each = 12))
      if (path == "partition") args$partition_factor <- ~part
      if (path == "random") {
        args$random <- ~g
        args$randcov_initial <- randcov_initial(g = .1, known = "given")
      }
      if (path == "network") {
        args$tailup_initial <- tailup_initial("exponential", de = .1, range = 1000, known = "given")
        args$taildown_initial <- taildown_initial("exponential", de = .1, range = 1000, known = "given")
        args$additive <- "afvArea"
      }
      if (path == "anisotropic") {
        args$euclid_initial <- euclid_initial("exponential", de = .1, range = 1000, rotate = .3, scale = .7, known = "given")
        args$anisotropy <- TRUE
      }
      fit <- do.call(ssn_glm, args)
      nugget <- coef(fit, type = "nugget")[["nugget"]]
      fixed <- function(value) {
        args$nugget_initial <- nugget_initial("nugget", value, known = "given")
        do.call(ssn_glm, args)
      }
      ref <- fixed(nugget)
      info <- paste(family, path, method)
      expect_equal(fit$optim$value, -2 * as.numeric(logLik(ref)), tolerance = 1e-9, info = info)
      expect_equal(coef(fit), coef(ref), tolerance = 1e-9, info = info)
      expect_equal(vcov(fit), vcov(ref), tolerance = 1e-9, info = info)
      expect_equal(fit$optim$value, -2 * as.numeric(logLik(fit)), tolerance = 1e-9, info = info)
      for (candidate in c(1e-4, .3, nugget / 1.05, nugget * 1.05)) {
        expect_true(fit$optim$value <= -2 * as.numeric(logLik(fixed(candidate))) + 1e-5, info = info)
      }
    }
  }
})

test_that("fixed, other-scalar and multiple-parameter GLM fits retain optimizer behavior", {
  s <- laplace_ssn_data(24)
  for (path in c("fixed", "range", "multiple", "dispersion")) {
    args <- list(formula = poisson ~ x + offset(off), ssn.object = s, family = "poisson",
      tailup_type = "none", nugget_initial = nugget_initial("nugget", .3,
        known = if (path == "multiple") NULL else "given"))
    args$euclid_initial <- euclid_initial("exponential", de = .1, range = 1000,
      known = if (path %in% c("range", "multiple")) "de" else "given")
    if (path == "dispersion") {
      args$formula <- Gamma ~ x + offset(off)
      args$family <- NULL
      args$dispersion_initial <- dispersion_initial("Gamma", 4)
    }
    fit <- do.call(ssn_glm, args)
    ref <- with_mocked_bindings(do.call(ssn_glm, args),
      optim_with_nugget_bound = function(optim_args, optim_dotlist, lower) do.call("optim", c(optim_args, optim_dotlist)))
    for (field in c("coefficients", "optim", "is_known", "fitted", "vcov")) {
      expect_identical(fit[[field]], ref[[field]], info = paste(path, field))
    }
  }
})
