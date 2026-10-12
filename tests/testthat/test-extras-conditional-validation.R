skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

conditional_validation_fixture <- function(glm = FALSE) {
  set.seed(2)
  network <- mf04p
  network$obs <- network$obs[round(seq(1, NROW(network$obs), length.out = 12)), ]
  network$preds <- list(validation = network$preds$pred1km[c(5, 1, 8), ])
  network$obs$x <- seq(-1, 1, length.out = 12)
  network$obs$off <- seq(-.2, .2, length.out = 12)
  network$obs$y <- c(1, 3, 2, 4, 2, 5, 4, 2, 5, 3, 6, 4)
  network$obs$g <- factor(rep(c("a", "b"), 6))
  network$preds$validation$x <- c(-.5, .2, .8)
  network$preds$validation$off <- c(.1, -.3, .2)
  network$preds$validation$g <- factor(c("a", "b", "new"))
  args <- list(formula = y ~ x + offset(off), ssn.object = network,
    random = ~g, randcov_initial = randcov_initial(g = .3, known = "given"),
    tailup_initial = tailup_initial("none"), taildown_initial = taildown_initial("none"),
    euclid_initial = euclid_initial("exponential", .4, 12000, known = "given"),
    nugget_initial = nugget_initial("nugget", .25, known = "given"))
  if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
  do.call(if (glm) ssn_glm else ssn_lm, args)
}

test_that("public samples validation precedes output shortcuts and dispatch", {
  set.seed(2)
  invalid <- list(NULL, numeric(), NA_real_, NaN, Inf, -Inf, 0, -1, 1.5,
                  c(1, 2), TRUE, "2")
  for (glm in c(FALSE, TRUE)) {
    fit <- conditional_validation_fixture(glm)
    for (local in list(FALSE, list(size = 2, ordering = "none"))) {
      for (samples in invalid) {
        expect_error(conditional(fit, output = "object", samples = samples, local = local),
                     "samples must be a positive integer", fixed = TRUE)
      }
      expect_identical(dim(conditional(fit, output = "object", samples = 1L, local = local)),
                       c(12L, 1L))
    }
  }
})

test_that("conditional controls reject missing nonfinite and nonscalar values", {
  set.seed(2)
  fit <- conditional_validation_fixture()
  invalid <- list(numeric(), NA_real_, NaN, Inf, -Inf, 0, -1, 1.5, c(1, 2), TRUE, "2")
  for (name in c("size", "chunk_size", "size_random", "groups_random")) {
    for (value in invalid) {
      local <- list(size = 2, ordering = "none")
      local[name] <- list(value)
      expect_error(conditional(fit, output = "object", samples = 1, local = local),
                   paste0("local$", name, " must be a positive integer"), fixed = TRUE)
    }
  }
  for (name in c("approximation", "method", "ordering", "method_random")) {
    for (value in list(character(), NA_character_, c("none", "none"), "invalid")) {
      local <- list(size = 2, ordering = "none")
      local[name] <- list(value)
      expect_error(conditional(fit, output = "object", samples = 1, local = local), name)
    }
  }
  for (value in list(NA, 0, "vecchia", c(TRUE, FALSE))) {
    expect_error(conditional(fit, output = "object", samples = 1, local = value),
                 "local must be FALSE, TRUE, or a list", fixed = TRUE)
  }
  defaults <- get_local_list_conditional(list(size = NULL, chunk_size = NULL), fit)
  expect_equal(defaults$size, 30L)
  expect_equal(defaults$chunk_size, 1000L)
  expect_error(get_local_list_conditional(list(groups_random = 13), fit),
               "cannot exceed the number of observed sites", fixed = TRUE)
  expect_error(get_local_list_conditional(list(size_response = NULL), fit), "size_response")
  expect_true(get_local_list_conditional(list(method = "all", size = NA, chunk_size = NA), fit)$exact)
})

test_that("invalid output and covariance simulation flags fail explicitly", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) {
    fit <- conditional_validation_fixture(glm)
    for (output in list(NULL, character(), NA_character_, c("beta", NA), 1, "invalid")) {
      expect_error(conditional(fit, output = output, samples = 1, local = FALSE),
                   "output must be", fixed = TRUE)
    }
    if (!glm) {
      for (flag in list(NULL, logical(), NA, c(TRUE, FALSE), 1, "TRUE")) {
        expect_error(conditional(fit, output = "object", samples = 1,
                                 simulate_covparams = flag, local = FALSE),
                     "simulate_covparams must be TRUE or FALSE", fixed = TRUE)
      }
    }
  }
})

test_that("singleton predictions and prediction-free output combinations keep matrix shapes", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) {
    fit <- conditional_validation_fixture(glm)
    fit$ssn.object$preds$validation <- fit$ssn.object$preds$validation[3, ]
    for (local in list(FALSE, list(size = 2, ordering = "none", chunk_size = 1))) {
      output <- c("all", "random", if (glm) "latent")
      for (samples in c(1L, 2L)) {
        value <- conditional(fit, "validation", output = output, samples = samples, local = local)
        expect_setequal(names(value), c("newdata", "beta", "object", "random", if (glm) "latent"))
        expect_identical(dim(value$newdata), c(1L, samples))
        expect_identical(dim(value$beta), c(2L, samples))
        expect_identical(dim(value$object), c(12L, samples))
        expect_identical(dim(value$random), c(3L, samples))
        expect_true(all(is.finite(unlist(value[c("newdata", "beta", "object", "random")]))))
        if (glm) expect_identical(dim(value$latent), c(12L, samples))
      }
      empty <- fit
      empty$ssn.object$preds$validation <- empty$ssn.object$preds$validation[FALSE, ]
      expect_error(conditional(empty, "validation", samples = 1, local = local), "has no rows")
      empty$ssn.object$preds <- list()
      output <- c("object", "beta", "random", if (glm) "latent")
      set.seed(2)
      absent <- conditional(empty, output = output, samples = 1, local = local)
      set.seed(2)
      ignored <- conditional(empty, "not-a-prediction-set", output = output, samples = 1, local = local)
      expect_identical(absent, ignored)
      expect_identical(names(absent), output)
      expect_identical(dim(absent$random), c(2L, 1L))
      expect_identical(dim(absent$beta), c(2L, 1L))
      expect_error(conditional(empty, samples = 1, local = local), "No missing data to predict")
    }
  }
})

test_that("user return paths reject a missing requested random output", {
  set.seed(2)
  calls <- c(exact = 0L, vecchia = 0L, lowrank = 0L, covparams = 0L)
  missing_draw <- function(name) function(...) {
    calls[name] <<- calls[name] + 1L
    list(beta = matrix(0, 2, 1))
  }
  local_mocked_bindings(
    draw_conditional_exact = missing_draw("exact"),
    draw_conditional_vecchia_model = missing_draw("vecchia"),
    draw_conditional_lowrank = missing_draw("lowrank"),
    draw_conditional_covparams = missing_draw("covparams"),
    get_satterthwaite_cached = function(...) list(vcov_theta = diag(1)),
    .package = "SSN2")
  for (glm in c(FALSE, TRUE)) {
    fit <- conditional_validation_fixture(glm)
    for (local in list(FALSE, list(size = 2, ordering = "none"),
                      list(approximation = "low-rank", size_base = 3, reorder_base = "none"))) {
      for (output in list("random", c("beta", "random"))) {
        expect_error(conditional(fit, output = output, samples = 1, local = local),
                     "Requested random-effect draws were not generated", fixed = TRUE)
      }
    }
    if (!glm) for (output in list("random", c("beta", "random"))) {
      expect_error(conditional(fit, output = output, samples = 1,
                               local = FALSE, simulate_covparams = TRUE),
                   "Requested random-effect draws were not generated", fixed = TRUE)
    }
  }
  expect_identical(calls, c(exact = 4L, vecchia = 4L, lowrank = 4L, covparams = 2L))
})

test_that("local latent factors and conditional weights respect covariance scale", {
  A <- matrix(c(2, .4, .2, .4, 1.5, .3, .2, .3, 1), 3)
  X <- cbind(1, c(-1, 0, 1))
  D <- c(-.2, .1, -.3)
  V <- solve(solve(A) - diag(D))
  M <- V %*% solve(A, X)
  for (scale in c(1e-12, 1, 1e12)) {
    gaussian <- prepare_conditional_local_gaussian(scale * A, X, 1:3, D / scale)
    expect_equal(gaussian$V / scale, V, tolerance = 1e-12)
    expect_equal(gaussian$M, M, tolerance = 1e-12)
    expect_equal(gaussian$mean, as.numeric(1:3))
    weights <- prepare_conditional_local_weights(gaussian$V, 3, 1:2, "scale", gaussian$factor)
    oracle <- solve(V[1:2, 1:2], V[1:2, 3])
    expect_equal(weights$weights, as.vector(oracle), tolerance = 1e-12)
    expect_equal(weights$variance / scale, V[3, 3] - sum(oracle * V[1:2, 3]), tolerance = 1e-12)
    rho <- 1 - 1e-9
    almost <- scale * matrix(c(1, rho, rho, 1), 2)
    precision <- conditional_local_factor(solve(almost), "precision")
    for (factor in list(NULL, precision)) {
      qr_weights <- prepare_conditional_local_weights(almost, 2, 1, "cancellation", factor)
      expect_gt(qr_weights$variance, 0)
      expect_equal(qr_weights$variance / (scale * (1 - rho) * (1 + rho)), 1, tolerance = 1e-6)
    }
  }
})

test_that("local solves reject singular indefinite and illegal covariance matrices", {
  invalid <- list(matrix(1, 2, 2), matrix(c(1, 2, 2, 1), 2),
                  matrix(NA_real_, 2, 2),
                  matrix(Inf, 2, 2), matrix(1, 2, 3), matrix(numeric(), 0, 0))
  for (A in invalid) {
    expect_error(prepare_conditional_local_gaussian(A, matrix(1, NROW(A), 1),
                                                   rep(0, NROW(A)), rep(0, NROW(A))),
                 "spatial covariance", fixed = TRUE)
  }
  for (D in list(c(1, 0), c(2, 0), c(NA, 0), c(Inf, 0))) {
    expect_error(prepare_conditional_local_gaussian(diag(2), matrix(1, 2, 1), c(0, 0), D),
                 "latent precision", fixed = TRUE)
  }
})

test_that("empty conditional neighborhoods return the marginal variance", {
  for (variance in c(0, 1e-12, 1, 1e12)) {
    value <- prepare_conditional_local_weights(matrix(variance), 1, integer(), "empty")
    expect_identical(value, list(weights = numeric(), variance = variance))
  }
})

test_that("exact GLM preparation rejects non-positive-definite factors", {
  set.seed(2)
  fit <- conditional_validation_fixture(TRUE)
  n <- fit$n
  covariance <- diag(n)
  curvature <- -diag(n)
  local_mocked_bindings(covmatrix = function(...) covariance,
                        get_D = function(...) curvature, .package = "SSN2")
  covariance <- matrix(1, n, n)
  expect_error(get_conditional_glm_joint(fit), "not positive", fixed = TRUE)
  covariance <- diag(n)
  for (value in c(1, 2)) {
    curvature <- -diag(n)
    curvature[1, 1] <- value
    expect_error(get_conditional_glm_joint(fit),
                 "The conditional GLM latent precision is not positive definite", fixed = TRUE)
  }
  curvature <- -diag(n)
  joint <- get_conditional_glm_joint(fit)
  expect_equal(crossprod(joint$cond_prec_upchol), 2 * diag(n), tolerance = 1e-12)
})

test_that("reduced GLM preparation rejects non-positive-definite factors", {
  set.seed(2)
  fit <- conditional_validation_fixture(TRUE)
  index <- c(7L, 2L, 10L)
  covariance <- diag(3)
  curvature <- -diag(3)
  local_mocked_bindings(conditional_vecchia_covariance = function(...) covariance,
                        get_D = function(...) curvature, .package = "SSN2")
  covariance <- matrix(1, 3, 3)
  expect_error(get_conditional_glm_base(fit, index, covariance = list()), "not positive", fixed = TRUE)
  covariance <- diag(3)
  for (value in c(1, 2)) {
    curvature <- -diag(3)
    curvature[1, 1] <- value
    expect_error(get_conditional_glm_base(fit, index, covariance = list()),
                 "The simulation base latent precision is not positive definite", fixed = TRUE)
  }
  curvature <- -diag(3)
  base <- get_conditional_glm_base(fit, index, covariance = list())
  expect_equal(tcrossprod(base$cov_lowchol_mH), 2 * diag(3), tolerance = 1e-12)
})

test_that("exact prediction draws reject singular conditional covariance", {
  context <- list(Sigma22 = matrix(1, 2, 2), SqrtSigInv_C0 = matrix(0, 1, 2),
                  newdata_name = "validation")
  expect_error(get_conditional_cov(context), "not positive", fixed = TRUE)
})
