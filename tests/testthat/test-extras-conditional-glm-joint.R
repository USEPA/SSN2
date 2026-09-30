skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
  "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

conditional_glm_joint_fixture <- local({
  fixture <- NULL
  function() {
    if (!is.null(fixture)) return(fixture)
    withr::local_preserve_seed()
    set.seed(2)
    root <- tempfile("conditional-glm-joint-")
    dir.create(root)
    stopifnot(file.copy(
      system.file("lsndata/MiddleFork04.ssn", package = "SSN2"),
      root, recursive = TRUE, copy.mode = FALSE
    ))
    network <- ssn_import(file.path(root, "MiddleFork04.ssn"),
      predpts = "CapeHorn", overwrite = TRUE, verbose = FALSE)
    network$obs <- network$obs[seq_len(16), ]
    network$preds <- list(joint = network$preds$CapeHorn[c(4, 1, 3, 2), ])
    network$obs$x <- seq(-1, 1, length.out = 16)
    network$preds$joint$x <- c(-0.8, 0.4, 1.4, -1.3)
    network$obs$off <- seq(-0.3, 0.3, length.out = 16)
    network$preds$joint$off <- c(0.2, -0.4, 0.1, 0.3)

    xy <- rbind(sf::st_coordinates(network$obs),
      sf::st_coordinates(network$preds$joint))
    covariance <- 0.4 * exp(-as.matrix(dist(xy)) / 12000) + diag(0.25, 20)
    X <- cbind("(Intercept)" = 1, x = network$obs$x)
    Xnew <- cbind("(Intercept)" = 1, x = network$preds$joint$x)
    eta <- as.vector(X %*% c(0.3, 0.2) + network$obs$off +
      t(chol(covariance[1:16, 1:16])) %*% rnorm(16))
    network$obs$y <- rpois(16, exp(eta))
    formula_env <- new.env(parent = baseenv())
    formula_env$offset <- stats::offset
    formula <- as.formula("y ~ x + offset(off)", env = formula_env)
    fit <- ssn_glm(formula, network, family = "poisson", local = FALSE,
      euclid_initial = euclid_initial("exponential", 0.4, 12000, known = "given"),
      nugget_initial = nugget_initial("nugget", 0.25, known = "given"))

    Sigma <- covariance[1:16, 1:16]
    C <- covariance[17:20, 1:16]
    Knew <- covariance[17:20, 17:20]
    P <- solve(Sigma)
    V0 <- solve(crossprod(X, P %*% X))
    B <- V0 %*% t(X) %*% P
    Q <- P - P %*% X %*% B
    w <- as.vector(fitted(fit, type = "link"))
    w_free <- w - network$obs$off
    D <- -diag(exp(w))
    L <- solve(Q - D)
    H <- Xnew - C %*% P %*% X
    R <- Knew - C %*% P %*% t(C)
    W <- H %*% B + C %*% P
    Vbeta <- V0 + B %*% L %*% t(B)
    Veta <- R + H %*% V0 %*% t(H) + W %*% L %*% t(W)
    Vcross <- V0 %*% t(H) + B %*% L %*% t(W)
    beta <- as.vector(coef(fit))
    mean_new <- as.vector(Xnew %*% beta +
      C %*% P %*% (w_free - X %*% beta) + network$preds$joint$off)
    fixture <<- list(fit = fit, Sigma = Sigma, C = C, Knew = Knew,
      X = X, Xnew = Xnew, V0 = V0, B = B, Q = Q, D = D, L = L,
      H = H, R = R, W = W, Vbeta = Vbeta, Veta = Veta, Vcross = Vcross,
      beta = beta, w = w, w_free = w_free, mean_new = mean_new)
    fixture
  }
})

conditional_glm_joint_methods <- function() {
  list(exact = FALSE,
    allbase = list(approximation = "low-rank", method_base = "all", method_new = "all"),
    allneighbor = list(approximation = "vecchia", method = "all", ordering = "none"))
}

expect_conditional_glm_joint_moments <- function(draws, mean, covariance, label) {
  n <- ncol(draws)
  expect_true(all(is.finite(draws)), info = label)
  mean_error <- abs(rowMeans(draws) - mean) / sqrt(diag(covariance) / n)
  # Gaussian covariance Monte Carlo Standard Error (MCSE) covers diagonal and off-diagonal entries separately.
  covariance_mcse <- sqrt((outer(diag(covariance), diag(covariance)) +
    covariance^2) / (n - 1))
  covariance_error <- abs(cov(t(draws)) - covariance) / covariance_mcse
  expect_lt(max(mean_error), 6, label = paste(label, "mean / MCSE"))
  expect_lt(max(covariance_error), 6, label = paste(label, "covariance / MCSE"))
}

test_that("GLM joint factors match independent Euclidean Laplace matrices", {
  f <- conditional_glm_joint_fixture()
  expect_equal(model.matrix(f$fit), f$X, ignore_attr = TRUE)
  expect_equal(unname(covmatrix(f$fit)), unname(f$Sigma), tolerance = 1e-10)
  expect_equal(unname(covmatrix(f$fit, "joint", cov_type = "obs.pred")),
    unname(t(f$C)), tolerance = 1e-10)
  expect_equal(unname(covmatrix(f$fit, "joint", cov_type = "pred.pred")),
    unname(f$Knew), tolerance = 1e-10)
  expect_equal(unname(vcov(f$fit, var_correct = FALSE)), unname(f$V0), tolerance = 1e-8)
  # Stored vcov uses the preceding Newton Hessian with a 1e-4 stopping tolerance.
  expect_lt(max(abs(vcov(f$fit) - f$Vbeta)), 1e-5)
  prediction <- predict(f$fit, "joint", type = "link", se.fit = TRUE)
  expect_equal(as.vector(prediction$fit), f$mean_new, tolerance = 1e-8)
  expect_equal(as.vector(prediction$se.fit^2), unname(diag(f$Veta)), tolerance = 1e-7)

  for (local in c(FALSE, TRUE)) {
    context <- SSN2:::get_conditional_context_glm(f$fit, "joint", local = local)
    joint <- SSN2:::get_conditional_glm_joint(context)
    expect_equal(unname(joint$weights_beta), unname(f$B), tolerance = 1e-8)
    expect_equal(unname(joint$mH_upchol), unname(chol(f$Q - f$D)), tolerance = 1e-8)
    expect_equal(unname(joint$beta_lowchol), unname(t(chol(f$V0))), tolerance = 1e-8)
    expect_equal(unname(crossprod(joint$mH_upchol)), unname(f$Q - f$D), tolerance = 1e-8)
    if (!local) {
      cond <- SSN2:::get_conditional_cov(context)
      expect_equal(unname(cond$Sigma_cond), unname(f$R), tolerance = 1e-8)
    }
  }
})

test_that("joint beta and offset-free latent draws preserve their dependence", {
  f <- conditional_glm_joint_fixture()
  context <- SSN2:::get_conditional_context_glm(f$fit, "joint")
  set.seed(2)
  draws <- SSN2:::draw_conditional_glm_joint(context, samples = 10000)
  expect_equal(dim(draws$beta), c(2L, 10000L))
  expect_equal(dim(draws$w), c(16L, 10000L))
  BL <- f$B %*% f$L
  covariance <- rbind(cbind(f$Vbeta, BL), cbind(t(BL), f$L))
  expect_conditional_glm_joint_moments(rbind(draws$beta, draws$w),
    c(f$beta, f$w_free), covariance, "beta and offset-free w")
  set.seed(2)
  expect_identical(SSN2:::draw_conditional_glm_joint(context, samples = 10000), draws)
  set.seed(2)
  one <- SSN2:::draw_conditional_glm_joint(context, samples = 1)
  expect_equal(dim(one$beta), c(2L, 1L))
  expect_equal(dim(one$w), c(16L, 1L))
})

test_that("exact-limit GLM paths match beta, new-site and cross moments", {
  f <- conditional_glm_joint_fixture()
  covariance <- rbind(cbind(f$Vbeta, f$Vcross), cbind(t(f$Vcross), f$Veta))
  excess <- f$H %*% (f$Vbeta - f$V0) %*% t(f$H)
  expect_gt(max(diag(excess) / (diag(f$Veta) * sqrt(2 / 9999))), 12)
  for (method in names(conditional_glm_joint_methods())) {
    local <- conditional_glm_joint_methods()[[method]]
    set.seed(2)
    draws <- conditional(f$fit, "joint", output = c("beta", "newdata"),
      type = "link", samples = 10000, local = local)
    expect_named(draws, c("beta", "newdata"))
    expect_equal(dim(draws$beta), c(2L, 10000L))
    expect_equal(dim(draws$newdata), c(4L, 10000L))
    expect_conditional_glm_joint_moments(rbind(draws$beta, draws$newdata),
      c(f$beta, f$mean_new), covariance, method)
    expect_conditional_glm_joint_moments(matrix(colMeans(draws$newdata), nrow = 1),
      mean(f$mean_new), matrix(sum(f$Veta) / 16, 1, 1), paste(method, "joint mean"))
  }
})

test_that("locally fitted GLMs retain fitting-block coefficient uncertainty", {
  f <- conditional_glm_joint_fixture()
  set.seed(2)
  fit <- ssn_glm(f$fit$formula, f$fit$ssn.object, family = "poisson",
    euclid_initial = euclid_initial("exponential", 0.4, 12000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.25, known = "given"),
    local = list(index = rep(1:2, each = 8), var_adjust = "none", parallel = FALSE))
  expect_equal(unname(fit$local_index), rep(1:2, each = 8))
  expect_equal(model.matrix(fit), f$X, ignore_attr = TRUE)
  P <- matrix(0, 16, 16)
  for (index in list(1:8, 9:16)) P[index, index] <- solve(f$Sigma[index, index])
  information <- crossprod(f$X, P %*% f$X)
  V0 <- solve(information)
  B <- V0 %*% t(f$X) %*% P
  # Local fitting regularizes the latent Hessian, not var_adjust = "none" beta V0.
  regularized <- solve(information + diag(fit$diagtol, 2))
  Q <- P - P %*% f$X %*% regularized %*% t(f$X) %*% P
  w <- as.vector(fitted(fit, type = "link"))
  D <- -diag(exp(w))
  L <- solve(Q - D)
  Vbeta <- V0 + B %*% L %*% t(B)
  expect_equal(unname(vcov(fit, var_correct = FALSE)), unname(V0), tolerance = 1e-8)
  expect_lt(max(abs(vcov(fit) - Vbeta)), 1e-5)
  for (local in c(FALSE, TRUE)) {
    context <- SSN2:::get_conditional_context_glm(fit, "joint", local = local)
    joint <- SSN2:::get_conditional_glm_joint(context)
    expect_equal(unname(joint$weights_beta), unname(B), tolerance = 1e-8)
    expect_equal(unname(joint$mH_upchol), unname(chol(Q - D)), tolerance = 1e-8)
    expect_equal(unname(joint$beta_lowchol), unname(t(chol(V0))), tolerance = 1e-8)
  }
  for (method in names(conditional_glm_joint_methods())) {
    set.seed(2)
    beta <- conditional(fit, "joint", output = "beta", samples = 10000,
      local = conditional_glm_joint_methods()[[method]])
    expect_equal(dim(beta), c(2L, 10000L))
    expect_conditional_glm_joint_moments(beta, as.vector(coef(fit)), Vbeta,
      paste("local fit beta", method))
  }
})

test_that("GLM object snapshots and beta-only outputs retain their contract", {
  f <- conditional_glm_joint_fixture()
  for (local in conditional_glm_joint_methods()) {
    for (samples in c(1L, 7L)) {
      snapshot <- matrix(f$w, nrow = 16, ncol = samples)
      set.seed(2)
      object <- conditional(f$fit, "joint", output = "object", samples = samples, local = local)
      expect_equal(unname(object), snapshot)
      set.seed(2)
      all <- conditional(f$fit, "joint", output = "all", samples = samples, local = local)
      expect_named(all, c("newdata", "beta", "object"))
      expect_equal(unname(all$object), snapshot)
      set.seed(2)
      beta <- conditional(f$fit, "joint", output = "beta", samples = samples, local = local)
      expect_equal(dim(beta), c(2L, samples))
      expect_identical(rownames(beta), names(coef(f$fit)))
      expect_identical(beta, all$beta)
      set.seed(2)
      expect_identical(conditional(f$fit, "joint", output = "all", samples = samples,
        local = local), all)
    }
  }
})

test_that("GLM prediction offsets enter once and response scales share link draws", {
  f <- conditional_glm_joint_fixture()
  no_new_offset <- f$fit
  off <- no_new_offset$ssn.object$preds$joint$off
  no_new_offset$ssn.object$preds$joint$off <- 0
  for (local in conditional_glm_joint_methods()) {
    set.seed(2)
    link <- conditional(f$fit, "joint", type = "link", samples = 31, local = local)
    set.seed(2)
    without <- conditional(no_new_offset, "joint", type = "link", samples = 31, local = local)
    expect_equal(link - without, matrix(off, nrow = 4, ncol = 31), ignore_attr = TRUE,
      tolerance = 1e-10)
    set.seed(2)
    response <- conditional(f$fit, "joint", type = "response", samples = 31, local = local)
    expect_equal(response, exp(link), tolerance = 1e-12)
    set.seed(2)
    counts <- conditional(f$fit, "joint", type = "new", samples = 31, local = local)
    expect_equal(dim(counts), c(4L, 31L))
    expect_true(all(is.finite(counts) & counts >= 0 & counts == round(counts)))
    set.seed(2)
    expect_identical(conditional(f$fit, "joint", type = "new", samples = 31, local = local), counts)
  }
})

test_that("one prediction and one GLM draw keep matrix dimensions", {
  f <- conditional_glm_joint_fixture()
  fit <- f$fit
  fit$ssn.object$preds$joint <- fit$ssn.object$preds$joint[2, , drop = FALSE]
  for (local in conditional_glm_joint_methods()) {
    for (type in c("link", "response", "new")) {
      set.seed(2)
      draws <- conditional(fit, "joint", output = "all", type = type,
        samples = 1, local = local)
      expect_equal(dim(draws$newdata), c(1L, 1L))
      expect_equal(dim(draws$beta), c(2L, 1L))
      expect_equal(dim(draws$object), c(16L, 1L))
      expect_true(all(is.finite(unlist(draws))))
      expect_equal(as.vector(draws$object), f$w)
      set.seed(2)
      expect_identical(conditional(fit, "joint", type = type,
        samples = 1, local = local), draws$newdata)
    }
  }
})

test_that("invalid GLM latent precision gives an actionable error", {
  f <- conditional_glm_joint_fixture()
  # Force invalid curvature without depending on an optimizer failure.
  local_mocked_bindings(get_D = function(...) diag(1e6, 16), .package = "SSN2")
  for (local in c(FALSE, TRUE)) {
    context <- SSN2:::get_conditional_context_glm(f$fit, "joint", local = local)
    expect_error(SSN2:::get_conditional_glm_joint(context),
      "latent-process precision is not positive definite")
    expect_error(SSN2:::draw_conditional_glm_joint(context, samples = 1),
      "latent-process precision is not positive definite")
  }
})

for (family in c("Gamma", "nbinomial", "binomial", "beta", "inverse.gaussian")) {
  test_that(paste("joint conditional smoke respects", family, "response semantics"), {
    if (family == "inverse.gaussian") skip_if_not_installed("statmod")
    f <- conditional_glm_joint_fixture()
    network <- f$fit$ssn.object
    dispersion <- switch(family, beta = 20, binomial = 1, 5)
    bounded <- family %in% c("binomial", "beta")
    sample_response <- function(mu, size) {
      switch(family,
        Gamma = rgamma(length(mu), shape = dispersion, scale = mu / dispersion),
        nbinomial = rnbinom(length(mu), mu = mu, size = dispersion),
        binomial = rbinom(length(mu), size = size, prob = mu),
        beta = pmin(pmax(rbeta(length(mu), mu * dispersion,
          (1 - mu) * dispersion), 1e-4), 1 - 1e-4),
        inverse.gaussian = statmod::rinvgauss(length(mu), mean = mu,
          dispersion = 1 / (mu * dispersion)))
    }
    eta <- 0.3 + 0.2 * network$obs$x + network$obs$off
    network$obs$size <- rep(c(4, 7), each = 8)
    set.seed(2)
    network$obs$y <- sample_response(if (bounded) plogis(eta) else exp(eta), network$obs$size)
    formula <- if (family == "binomial") {
      as.formula("cbind(y, size - y) ~ x + offset(off)", env = environment(f$fit$formula))
    } else {
      f$fit$formula
    }
    initial <- do.call(dispersion_initial,
      list(family = family, dispersion = dispersion, known = "given"))
    fit <- ssn_glm(formula, network,
      dispersion_initial = initial,
      euclid_initial = euclid_initial("exponential", 0.4, 12000, known = "given"),
      nugget_initial = nugget_initial("nugget", 0.25, known = "given"), local = FALSE)
    expect_identical(fit$family, family)
    expect_equal(as.vector(coef(fit, type = "dispersion")), dispersion)
    size <- c(1, 3, 7, 11)
    for (local in conditional_glm_joint_methods()) {
      set.seed(2)
      link <- conditional(fit, "joint", type = "link", samples = 23,
        newdata_size = size, local = local)
      mu <- if (bounded) plogis(link) else exp(link)
      expected_new <- vapply(seq_len(ncol(mu)), function(j) sample_response(mu[, j], size),
        numeric(4))
      set.seed(2)
      response <- conditional(fit, "joint", type = "response", samples = 23,
        newdata_size = size, local = local)
      expect_equal(response, if (family == "binomial") mu * size else mu, tolerance = 1e-12)
      set.seed(2)
      new <- conditional(fit, "joint", type = "new", samples = 23,
        newdata_size = size, local = local)
      expect_equal(new, expected_new, ignore_attr = TRUE, tolerance = 1e-12)
      expect_equal(dim(new), c(4L, 23L))
      expect_true(all(is.finite(link) & is.finite(response) & is.finite(new)))
      if (family %in% c("binomial", "nbinomial")) {
        expect_true(all(new >= 0 & new == round(new)))
        if (family == "binomial") expect_true(all(new <= size))
      } else {
        expect_true(all(new > 0))
        if (family == "beta") expect_true(all(new < 1))
      }
    }
  })
}
