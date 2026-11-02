skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")
skip_if_not_installed("pracma")

conditional_wave_network <- function() {
  network <- mf04p
  network$obs <- network$obs[rev(round(seq(1, NROW(network$obs), length.out = 12))), ]
  network$preds <- list(pred1km = network$preds$pred1km[c(5, 1, 8, 3), ])
  network$obs$x <- seq(-1, 1, length.out = 12)
  network$obs$off <- seq(-.2, .2, length.out = 12)
  network$obs$y <- c(1, 3, 2, 4, 2, 5, 4, 2, 5, 3, 6, 4)
  network$obs$g <- factor(rep(c("a", "b"), 6))
  network$obs$part <- factor(rep(c("one", "two"), each = 6))
  network$preds$pred1km$x <- c(-.5, .2, .8, 0)
  network$preds$pred1km$off <- c(.1, -.3, .2, .4)
  network$preds$pred1km$g <- factor(c("a", "b", "new", "a"))
  network$preds$pred1km$part <- factor(c("one", "two", "one", "two"))
  network$path <- withr::local_tempdir(.local_envir = parent.frame())
  stopifnot(file.copy(file.path(mf04p$path, "binaryID.db"), network$path))
  ssn_create_distmat(network, predpts = "pred1km", among_predpts = TRUE, overwrite = TRUE)
  network
}

conditional_wave_fit <- function(network, glm = FALSE, component = "both", effects = FALSE) {
  args <- list(formula = y ~ x + offset(off), ssn.object = network, additive = "afvArea",
    tailup_initial = if (component == "down") tailup_initial("none") else
      tailup_initial("wave", .4, 2500, known = "given"),
    taildown_initial = if (component == "up") taildown_initial("none") else
      taildown_initial("wave", .3, 3500, known = "given"),
    euclid_initial = euclid_initial("none"),
    nugget_initial = nugget_initial("nugget", .25, known = "given"))
  if (effects) {
    args$euclid_initial <- euclid_initial("exponential", .2, 12000, known = "given")
    args$random <- ~g
    args$randcov_initial <- randcov_initial(g = .15, known = "given")
    args$partition_factor <- ~part
  }
  if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
  do.call(if (glm) ssn_glm else ssn_lm, args)
}

conditional_wave_reference <- function(fit, corrected = TRUE) {
  glm <- inherits(fit, "ssn_glm")
  X <- model.matrix(fit)
  pred <- fit$ssn.object$preds$pred1km
  Xp <- cbind(1, pred$x)
  S <- covmatrix(fit)
  C <- covmatrix(fit, "pred1km", cov_type = "obs.pred")
  P <- covmatrix(fit, "pred1km", cov_type = "pred.pred")
  B <- as.matrix(vcov(fit))
  beta <- as.vector(if (glm && corrected) fit$laplace$betahat else coef(fit))
  eta <- if (glm) fitted(fit, type = "link") else model.response(model.frame(fit))
  center <- eta
  if (glm && corrected) center <- center + fit$laplace$delta[names(eta)]
  center <- as.vector(center - model.offset(model.frame(fit)))
  V <- if (glm) solve(solve(S) + diag(exp(eta))) else matrix(0, fit$n, fit$n)
  M <- if (glm) V %*% solve(S, X) else matrix(0, fit$n, NCOL(X))
  W <- t(C) %*% solve(S)
  H <- Xp + W %*% (M - X)
  mean <- as.vector(Xp %*% beta + W %*% (center - X %*% beta) + pred$off)
  variance <- P - W %*% C + W %*% V %*% t(W) + H %*% B %*% t(H)
  list(mean = c(beta, mean), covariance = rbind(cbind(B, B %*% t(H)),
      cbind(H %*% B, variance)),
    latent_mean = c(beta, center + model.offset(model.frame(fit))),
    latent_covariance = rbind(cbind(B, B %*% t(M)), cbind(M %*% B, V + M %*% B %*% t(M))))
}

expect_conditional_wave_moments <- function(draws, mean, covariance) {
  samples <- NCOL(draws)
  expect_true(all(is.finite(draws)))
  expect_lt(max(abs(rowMeans(draws) - mean) / sqrt(diag(covariance) / samples)), 6)
  mcse <- sqrt((outer(diag(covariance), diag(covariance)) + covariance^2) / (samples - 1))
  expect_lt(max(abs(cov(t(draws)) - covariance) / mcse), 6)
}

test_that("exact wave conditional draws recover joint coefficient and prediction moments", {
  set.seed(2)
  network <- conditional_wave_network()
  for (glm in c(FALSE, TRUE)) for (component in c("up", "down", "both")) {
    fit <- conditional_wave_fit(network, glm, component, effects = component == "both")
    reference <- conditional_wave_reference(fit)
    set.seed(2)
    draws <- conditional(fit, "pred1km", output = c("all", if (glm) "latent"),
      samples = 3000, local = FALSE)
    expect_equal(dim(draws$newdata), c(4L, 3000L))
    expect_conditional_wave_moments(rbind(draws$beta, draws$newdata),
      reference$mean, reference$covariance)
    h <- c(0, 0, rep(.25, 4))
    expect_conditional_wave_moments(matrix(colMeans(draws$newdata), 1),
      sum(h * reference$mean), matrix(as.numeric(t(h) %*% reference$covariance %*% h)))
    eta <- if (glm) fitted(fit, type = "link") else model.response(model.frame(fit))
    expect_equal(draws$object, matrix(eta, 12, 3000), ignore_attr = TRUE)
    if (glm) expect_conditional_wave_moments(rbind(draws$beta, draws$latent),
      reference$latent_mean, reference$latent_covariance)
  }
})

test_that("wave conditional local paths recover full-neighborhood moments", {
  set.seed(2)
  network <- conditional_wave_network()
  for (glm in c(FALSE, TRUE)) {
    fit <- conditional_wave_fit(network, glm)
    reference <- conditional_wave_reference(fit)
    for (local in list(
      list(approximation = "vecchia", method = "covariance", size = 16, ordering = "none"),
      list(approximation = "low-rank", size_base = 12, reorder_base = "none",
           size_new = 4, reorder_new = "none", kmeans_new = FALSE))) {
      set.seed(2)
      draws <- conditional(fit, "pred1km", output = c("beta", "newdata"), samples = 3000, local = local)
      expect_conditional_wave_moments(rbind(draws$beta, draws$newdata),
        reference$mean, reference$covariance)
    }
    for (local in list(list(size = 3, ordering = "none"),
      list(approximation = "low-rank", size_base = 6, size_new = 2,
           reorder_base = "none", reorder_new = "none", kmeans_new = FALSE))) {
      set.seed(2)
      draws <- conditional(fit, "pred1km", samples = 7, local = local)
      expect_equal(dim(draws), c(4L, 7L))
      expect_true(all(is.finite(draws)))
    }
  }
})

test_that("Poisson wave conditionals distinguish latent means and response draws", {
  set.seed(2)
  network <- conditional_wave_network()
  fit <- conditional_wave_fit(network, glm = TRUE)
  for (corrected in c(FALSE, TRUE)) {
    for (local in list(FALSE, list(size = 3, ordering = "none"),
      list(approximation = "low-rank", size_base = 6, size_new = 2,
           reorder_base = "none", reorder_new = "none", kmeans_new = FALSE))) {
      set.seed(2)
      link <- conditional(fit, "pred1km", type = "link", samples = 23,
        fit_correct = corrected, local = local)
      expected <- vapply(seq_len(NCOL(link)), function(j) rpois(4, exp(link[, j])), numeric(4))
      set.seed(2)
      response <- conditional(fit, "pred1km", type = "response", samples = 23,
        fit_correct = corrected, local = local)
      expect_equal(response, exp(link), tolerance = 1e-12)
      set.seed(2)
      new <- conditional(fit, "pred1km", type = "new", samples = 23,
        fit_correct = corrected, local = local)
      expect_equal(new, expected, ignore_attr = TRUE, tolerance = 1e-12)
      expect_true(all(is.finite(new) & new >= 0 & new == round(new)))
    }
    reference <- conditional_wave_reference(fit, corrected)
    set.seed(2)
    draws <- conditional(fit, "pred1km", output = c("beta", "newdata", "latent"),
      samples = 3000, local = FALSE, fit_correct = corrected)
    expect_conditional_wave_moments(rbind(draws$beta, draws$newdata),
      reference$mean, reference$covariance)
    expect_conditional_wave_moments(rbind(draws$beta, draws$latent),
      reference$latent_mean, reference$latent_covariance)
  }
})
