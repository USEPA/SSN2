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
    for (samples in invalid) {
      expect_error(conditional(fit, output = "object", samples = samples),
                   "samples must be a positive integer", fixed = TRUE)
    }
    expect_identical(dim(conditional(fit, output = "object", samples = 1L)),
                     c(12L, 1L))
  }
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
    output <- c("all", "random", if (glm) "latent")
    for (samples in c(1L, 2L)) {
      value <- conditional(fit, "validation", output = output, samples = samples)
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
    expect_error(conditional(empty, "validation", samples = 1), "has no rows")
    empty$ssn.object$preds <- list()
    output <- c("object", "beta", "random", if (glm) "latent")
    set.seed(2)
    absent <- conditional(empty, output = output, samples = 1)
    set.seed(2)
    ignored <- conditional(empty, "not-a-prediction-set", output = output, samples = 1)
    expect_identical(absent, ignored)
    expect_identical(names(absent), output)
    expect_identical(dim(absent$random), c(2L, 1L))
    expect_identical(dim(absent$beta), c(2L, 1L))
    expect_error(conditional(empty, samples = 1), "No missing data to predict")
  }
})

test_that("user return paths reject a missing requested random output", {
  set.seed(2)
  calls <- c(exact = 0L, covparams = 0L)
  missing_draw <- function(name) function(...) {
    calls[name] <<- calls[name] + 1L
    list(beta = matrix(0, 2, 1))
  }
  local_mocked_bindings(
    draw_conditional_exact = missing_draw("exact"),
    draw_conditional_covparams = missing_draw("covparams"),
    get_satterthwaite_cached = function(...) list(vcov_theta = diag(1)),
    .package = "SSN2")
  for (glm in c(FALSE, TRUE)) {
    fit <- conditional_validation_fixture(glm)
    for (output in list("random", c("beta", "random"))) {
      expect_error(conditional(fit, output = output, samples = 1),
                   "Requested random-effect draws were not generated", fixed = TRUE)
    }
    if (!glm) for (output in list("random", c("beta", "random"))) {
      expect_error(conditional(fit, output = output, samples = 1,
                               simulate_covparams = TRUE),
                   "Requested random-effect draws were not generated", fixed = TRUE)
    }
  }
  expect_identical(calls, c(exact = 4L, covparams = 2L))
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

test_that("exact prediction draws reject singular conditional covariance", {
  context <- list(Sigma22 = matrix(1, 2, 2), SqrtSigInv_C0 = matrix(0, 1, 2),
                  newdata_name = "validation")
  expect_error(get_conditional_cov(context), "not positive", fixed = TRUE)
})
