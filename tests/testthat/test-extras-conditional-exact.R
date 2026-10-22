skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

conditional_exact_fixture <- function(glm = FALSE, stream = FALSE, missing = FALSE) {
  network <- mf04p
  network$obs <- network$obs[round(seq(1, NROW(network$obs), length.out = 12)), ]
  network$preds <- list(pred1km = network$preds$pred1km[c(5, 1, 8, 3), ])
  network$obs$x <- seq(-1, 1, length.out = 12)
  network$obs$off <- seq(-.2, .2, length.out = 12)
  network$obs$y <- c(1, 3, 2, 4, 2, 5, 4, 2, 5, 3, 6, 4)
  network$obs$g <- factor(rep(c("a", "b"), 6))
  network$obs$part <- factor(rep(c("one", "two"), each = 6))
  network$preds$pred1km$x <- c(-.5, .2, .8, 0)
  network$preds$pred1km$off <- c(.1, -.3, .2, .4)
  network$preds$pred1km$g <- factor(c("a", "a", "new", "new"))
  network$preds$pred1km$part <- factor(c("one", "two", "one", "one"))
  if (stream) {
    root <- tempfile("conditional-exact-")
    dir.create(root)
    stopifnot(file.copy(network$path, root, recursive = TRUE, copy.mode = FALSE))
    network$path <- file.path(root, basename(network$path))
    ssn_create_distmat(network, predpts = "pred1km", among_predpts = TRUE, overwrite = TRUE)
  }
  if (missing) network$obs$y[c(2, 7)] <- NA
  args <- list(formula = y ~ x + offset(off), ssn.object = network,
    additive = "afvArea", random = ~g, partition_factor = ~part,
    randcov_initial = randcov_initial(g = .3, known = "given"),
    tailup_initial = if (stream) tailup_initial("exponential", .3, 10000, known = "given") else tailup_initial("none"),
    taildown_initial = if (stream) taildown_initial("spherical", .2, 10000, known = "given") else taildown_initial("none"),
    euclid_initial = euclid_initial("exponential", .4, 12000, known = "given"),
    nugget_initial = nugget_initial("nugget", .25, known = "given"))
  if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
  do.call(if (glm) ssn_glm else ssn_lm, args)
}

expect_exact_moments <- function(draws, mean, covariance) {
  n <- NCOL(draws)
  expect_lt(max(abs(rowMeans(draws) - mean) / sqrt(diag(covariance) / n)), 6)
  mcse <- sqrt((outer(diag(covariance), diag(covariance)) + covariance^2) / (n - 1))
  expect_lt(max(abs(cov(t(draws)) - covariance) / mcse), 6)
}

test_that("shared-effect exact draws recover independent total-covariance moments", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) for (stream in c(FALSE, TRUE)) {
    fit <- conditional_exact_fixture(glm, stream)
    X <- model.matrix(fit)
    pred <- fit$ssn.object$preds$pred1km
    Xp <- cbind(1, pred$x)
    Sigma <- as.matrix(covmatrix(fit))
    cross <- as.matrix(covmatrix(fit, "pred1km", cov_type = "obs.pred"))
    prior <- as.matrix(covmatrix(fit, "pred1km", cov_type = "pred.pred"))
    C <- as.matrix(vcov(fit))
    W <- t(cross) %*% solve(Sigma)
    R <- prior - W %*% cross
    beta <- as.vector(if (glm) fit$laplace$betahat else coef(fit))
    eta <- if (glm) fitted(fit, type = "link") else model.response(model.frame(fit))
    center <- if (glm) as.vector(eta + fit$laplace$delta[names(eta)]) else eta
    center <- center - model.offset(model.frame(fit))
    V <- if (glm) solve(solve(Sigma) + diag(exp(fitted(fit, type = "link")))) else matrix(0, 12, 12)
    M <- if (glm) V %*% solve(Sigma, X) else matrix(0, 12, NCOL(X))
    B <- Xp + W %*% (M - X)
    covariance <- R + W %*% V %*% t(W) + B %*% C %*% t(B)
    mean <- as.vector(Xp %*% beta + W %*% (center - X %*% beta) + pred$off)
    set.seed(2)
    draws <- conditional(fit, "pred1km", output = c("all", "random", if (glm) "latent"),
                         samples = 6000, local = FALSE)
    expect_exact_moments(rbind(draws$beta, draws$newdata), c(beta, mean),
      rbind(cbind(C, C %*% t(B)), cbind(B %*% C, covariance)))
    expect_exact_moments(matrix(colMeans(draws$newdata), 1), mean(mean), matrix(sum(covariance) / 16, 1))
    expect_equal(draws$object, matrix(eta, 12, 6000), ignore_attr = TRUE)
    if (glm) {
      expect_exact_moments(rbind(draws$beta, draws$latent),
        c(beta, center + model.offset(model.frame(fit))),
        rbind(cbind(C, C %*% t(M)), cbind(M %*% C, V + M %*% C %*% t(M))))
      expect_identical(rownames(draws$latent), rownames(X))
    }
    design <- get_conditional_random_design(fit, pred)
    Z <- as.matrix(design$observed)
    S <- Sigma - Z %*% diag(design$variance) %*% t(Z)
    H <- solve(diag(1 / design$variance) + t(Z) %*% solve(S, Z))
    residual <- if (glm) draws$latent - model.offset(model.frame(fit)) - X %*% draws$beta else
      center - X %*% draws$beta
    innovation <- draws$random - H %*% t(Z) %*% solve(S, residual)
    expect_exact_moments(innovation, rep(0, NROW(H)), H)
    expect_identical(attr(draws$random, "effect_info"), design$info)
    expect_equal(sum(!design$info$observed), 1)
    expect_equal(as.numeric(H[!design$info$observed, !design$info$observed]), .3)
    expect_equal(length(unique(design$info$partition)), 2)
  }
})

test_that("effect preparation preserves slopes, zero variances and unconditioned levels", {
  set.seed(2)
  fit <- conditional_exact_fixture()
  obs <- fit$ssn.object$obs
  pred <- fit$ssn.object$preds$pred1km
  design <- get_randcov_design(obs, c("1 | g" = .3, "x | g" = .2, "1 | part" = 0),
                              ~part, pred)
  S <- as.matrix(covmatrix(conditional_spatial_object(fit)))
  effect <- prepare_conditional_random(fit, design, covariance = S)
  active <- effect$active
  Z <- as.matrix(design$observed[, active, drop = FALSE])
  H <- solve(diag(1 / design$variance[active]) + t(Z) %*% solve(S, Z))
  expect_equal(effect$covariance, H, ignore_attr = TRUE, tolerance = 1e-10)
  set.seed(2)
  draws <- draw_conditional_random(effect, matrix(0, 12, 6000))
  expect_true(all(draws[design$variance == 0, ] == 0))
  positive <- which(design$variance > 0)
  prior <- diag(design$variance)
  prior[active, active] <- H
  expect_exact_moments(draws[positive, ], rep(0, length(positive)), prior[positive, positive])
})

test_that("exact outputs avoid prediction data and unnecessary factors", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) {
    fit <- conditional_exact_fixture(glm)
    fit$ssn.object$preds <- list()
    set.seed(2)
    value <- conditional(fit, output = c("beta", "object", "random", if (glm) "latent"),
                         samples = 1, local = FALSE)
    expect_equal(dim(value$beta), c(2, 1))
    expect_equal(dim(value$random), c(4, 1))
    if (glm) expect_equal(dim(value$latent), c(12, 1))
    expect_error(conditional(fit, local = FALSE), "No missing data to predict")
  }
  fit <- conditional_exact_fixture(TRUE)
  local_mocked_bindings(covmatrix = function(...) stop("Unexpected covariance read"),
    get_prediction_newdata = function(...) stop("Unexpected prediction read"), .package = "SSN2")
  expect_equal(dim(conditional(fit, output = "beta", samples = 2, local = FALSE)), c(2, 2))
  expect_equal(dim(conditional(fit, output = "object", samples = 2, local = FALSE)), c(12, 2))
})

test_that("missing predictions, output combinations and offsets retain row orientation", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) {
    fit <- conditional_exact_fixture(glm, missing = TRUE)
    set.seed(2)
    implicit <- conditional(fit, samples = 7, local = FALSE)
    set.seed(2)
    explicit <- conditional(fit, ".missing", samples = 7, local = FALSE)
    expect_identical(implicit, explicit)
    expect_identical(rownames(implicit), as.character(fit$missing_index))
    expect_equal(dim(implicit), c(2, 7))
  }
  fit <- conditional_exact_fixture(TRUE)
  set.seed(2)
  link <- conditional(fit, "pred1km", output = c("all", "latent", "random"), samples = 7, local = FALSE)
  set.seed(2)
  response <- conditional(fit, "pred1km", output = c("all", "latent", "random"),
                          samples = 7, local = FALSE, type = "response")
  expect_equal(response$newdata, exp(link$newdata))
  expect_identical(response$latent, link$latent)
  expect_identical(response$random, link$random)
  expect_named(link, c("newdata", "beta", "object", "latent", "random"))
})

test_that("exact prediction reuses the GLM model covariance factor", {
  set.seed(2)
  fit <- conditional_exact_fixture(TRUE)
  fit <- conditional_spatial_object(fit)
  original <- covmatrix
  observed_calls <- 0L
  local_mocked_bindings(covmatrix = function(object, ...) {
    if (!length(list(...))) observed_calls <<- observed_calls + 1L
    original(object, ...)
  }, .package = "SSN2")
  conditional(fit, "pred1km", samples = 2, local = FALSE)
  expect_equal(observed_calls, 1L)
})

test_that("covariance-parameter draws retain paired realized effects without prediction data", {
  set.seed(2)
  fit <- conditional_exact_fixture()
  fit$ssn.object$preds <- list()
  params <- fit$coefficients$params_object
  theta <- c(nugget = .25, "1 | g" = .3)
  sw <- list(context = list(cov_names_free_orig = names(theta), cov_val_free = theta,
                             params_object = params), vcov_theta = diag(.01, 2))
  count <- 0L
  local_mocked_bindings(
    get_satterthwaite_cached = function(...) sw,
    simulate_theta_draw_ssn = function(...) {
      count <<- count + 1L
      list(theta = theta * count, exhausted = FALSE, invalid_names = character())
    },
    fill_perturbed_params_object = function(theta, ...) {
      params$nugget[["nugget"]] <- theta[["nugget"]]
      params$randcov[["1 | g"]] <- theta[["1 | g"]]
      params
    }, .package = "SSN2")
  set.seed(2)
  paired <- conditional(fit, output = c("cov", "beta", "random"),
                        samples = 3, simulate_covparams = TRUE, local = FALSE)
  expect_equal(paired$cov, outer(theta, 1:3), ignore_attr = TRUE)
  set.seed(2)
  X <- model.matrix(fit)
  y <- model.response(model.frame(fit)) - model.offset(model.frame(fit))
  for (b in 1:3) {
    object <- fit
    object$coefficients$params_object$nugget[["nugget"]] <- theta[["nugget"]] * b
    object$coefficients$params_object$randcov[["1 | g"]] <- theta[["1 | g"]] * b
    P <- solve(as.matrix(covmatrix(object)))
    C <- solve(t(X) %*% P %*% X)
    beta <- as.vector(coef(fit)) + as.vector(t(chol(C)) %*% rnorm(NCOL(X)))
    design <- get_conditional_random_design(object)
    Z <- as.matrix(design$observed)
    S <- as.matrix(covmatrix(conditional_spatial_object(object)))
    H <- solve(diag(1 / design$variance) + t(Z) %*% solve(S, Z))
    effect <- H %*% t(Z) %*% solve(S, y - X %*% beta) +
      t(chol(H)) %*% rnorm(NCOL(Z))
    expect_equal(paired$beta[, b], beta, ignore_attr = TRUE, tolerance = 1e-10)
    expect_equal(paired$random[, b], as.vector(effect), ignore_attr = TRUE, tolerance = 1e-10)
  }
  count <- 0L
  set.seed(2)
  only <- conditional(fit, output = "random", samples = 3, simulate_covparams = TRUE, local = FALSE)
  expect_identical(only, paired$random)
})

test_that("joint factors use stored coefficient covariance and permit positive second derivitave entries", {
  fit <- conditional_exact_fixture(TRUE)
  C <- as.matrix(vcov(fit))
  original_vcov <- vcov
  D <- diag(c(.01, rep(-1, 11)))
  local_mocked_bindings(
    vcov = function(object, ...) 2 * original_vcov(object, ...),
    get_D = function(...) D, .package = "SSN2")
  joint <- get_conditional_glm_joint(fit)
  P <- solve(as.matrix(covmatrix(fit)))
  expect_equal(tcrossprod(joint$cov_betahat_lowchol), 2 * C, ignore_attr = TRUE)
  expect_equal(crossprod(joint$cond_prec_upchol), P - D, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(joint$wts_latent, solve(P - D, P %*% model.matrix(fit)),
               ignore_attr = TRUE, tolerance = 1e-10)
})

test_that("effect covariance guards reject invalid variances", {
  fit <- conditional_exact_fixture()
  design <- get_conditional_random_design(fit)
  design$variance[1] <- -1
  expect_error(prepare_conditional_random(fit, design), "Invalid random-effect variances", fixed = TRUE)
})
